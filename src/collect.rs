use std::{sync::Arc, time::Duration, collections::HashMap};
use tokio::{sync::{mpsc, Mutex}, time::sleep};
use anyhow::Result;
use regex::Regex;
use serialport;
use std::io;
use simple_logger::SimpleLogger;
use libgtr;
use log::{info, warn, error, debug, trace};

#[derive(serde::Serialize, Debug)]
pub struct Record {
    pub timestamp: chrono::DateTime<chrono::Utc>,
    pub weight: f64,
    pub cog_1: f64,
    pub cog_2: f64,
    pub cog_3: f64,
    pub cog_4: f64,
    pub cog_5: f64,
    pub cog_6: f64,
    pub cog_7: f64,
    pub cog_8: f64,
}

#[derive(Debug, PartialEq)]
enum Phase {
    StartPump,
    WaitStableDecrease,
    StopPump,
    TareScale,
    StartFill,
    WaitStableIncrease,
    StopFill,
    StoreData,
    Done,
}

struct Context {
    weight_rx: mpsc::Receiver<f64>,
    psu: Arc<Mutex<mxpsu::MxSeries>>,
    capture: Vec<Record>,
    stable_secs: u64,
    gtr: Option<libgtr::GtrSerialReader>,
    output_file: String,
}

impl Context {
    fn new(
        weight_rx: mpsc::Receiver<f64>, 
        psu: Arc<Mutex<mxpsu::MxSeries>>, 
        stable_secs: u64,
        output_file: String
    ) -> Self {
        info!("Creating new Context with stable_secs={}, output={}", stable_secs, output_file);
        Context {
            weight_rx,
            psu,
            capture: Vec::new(),
            stable_secs,
            gtr: None,
            output_file,
        }
    }

    async fn init_gtr(&mut self, port: &str) -> Result<()> {
        info!("Initializing GTR on port: {}", port);
        self.gtr = Some(libgtr::GtrSerialReader::new(port, 115200)?);
        info!("GTR initialized successfully");
        Ok(())
    }

    async fn stop_gtr(&mut self) -> Result<()> {
        if let Some(mut gtr) = self.gtr.take() {
            info!("Stopping GTR reader...");
            tokio::task::spawn_blocking(move || {
                if let Err(e) = gtr.stop() {
                    error!("Error stopping GTR: {:?}", e);
                } else {
                    info!("GTR reader stopped successfully");
                }
            }).await?;
        }
        Ok(())
    }

    async fn collect_cog_data(&mut self, weight: f64) -> Result<()> {
        if let Some(ref mut gtr) = self.gtr {
            let mut packets_collected = 0;
            
            // Drain ALL available packets from the queue
            loop {
                match gtr.try_recv_packet() {
                    Ok(Some(packet)) => {
                        let record = Record {
                            timestamp: chrono::Utc::now(),
                            weight, // Use current weight for ALL packets in this timeframe
                            cog_1: packet.sensors[0].cog_value as f64,
                            cog_2: packet.sensors[1].cog_value as f64,
                            cog_3: packet.sensors[2].cog_value as f64,
                            cog_4: packet.sensors[3].cog_value as f64,
                            cog_5: packet.sensors[4].cog_value as f64,
                            cog_6: packet.sensors[5].cog_value as f64,
                            cog_7: packet.sensors[6].cog_value as f64,
                            cog_8: packet.sensors[7].cog_value as f64,
                        };

                        self.capture.push(record);
                        packets_collected += 1;
                        
                        debug!(
                            "Packet #{}: Weight={:.2}lb, Counter={}, COGs=[{:.2}, {:.2}, {:.2}, {:.2}, {:.2}, {:.2}, {:.2}, {:.2}]",
                            packets_collected,
                            weight,
                            packet.header.counter,
                            packet.sensors[0].cog_value as f64,
                            packet.sensors[1].cog_value as f64,
                            packet.sensors[2].cog_value as f64,
                            packet.sensors[3].cog_value as f64,
                            packet.sensors[4].cog_value as f64,
                            packet.sensors[5].cog_value as f64,
                            packet.sensors[6].cog_value as f64,
                            packet.sensors[7].cog_value as f64
                        );
                    }
                    Ok(None) => {
                        // No more packets available - break out of loop
                        break;
                    }
                    Err(libgtr::GtrError::ThreadComm(msg)) => {
                        error!("GTR communication error: {}", msg);
                        break;
                    }
                    Err(e) => {
                        error!("GTR error: {:?}", e);
                        break;
                    }
                }
            }
            
            if packets_collected > 0 {
                info!("Collected {} GTR packets at weight {:.2}lb (Total records: {})", 
                      packets_collected, weight, self.capture.len());
            } else {
                debug!("No GTR packets available at weight {:.2}lb", weight);
            }
        } else {
            warn!("GTR not initialized, skipping COG data collection");
        }
        Ok(())
    }

    async fn flush_to_csv(&self) -> Result<()> {
        info!("Saving {} records to CSV file: {}", self.capture.len(), self.output_file);
        let mut wtr = csv::Writer::from_path(&self.output_file)?;
        for rec in &self.capture {
            wtr.serialize(rec)?;
        }
        wtr.flush()?;
        info!("Data saved to {} successfully", self.output_file);
        Ok(())
    }
}

// ... (copy all the other functions: psu_set, wait_for_stability, scale_reader, etc.)

async fn psu_set(ctx: Arc<Mutex<mxpsu::MxSeries>>, ch: u8, on: bool) -> Result<()> {
    let action = if on { "ON" } else { "OFF" };
    info!("Setting PSU channel {} {}", ch, action);

    let psu = ctx.clone();
    tokio::task::spawn_blocking(move || {
        let mut guard = psu.blocking_lock();
        if on {
            guard.turn_on(ch);
        } else {
            guard.turn_off(ch);
        }
    })
    .await?;

    debug!("PSU channel {} set {} successfully", ch, action);
    Ok(())
}

async fn wait_for_stability(rx: &mut mpsc::Receiver<f64>, increasing: bool, secs: u64) {
    let mut stable_count = 0;
    let trend_str = if increasing { "increasing" } else { "decreasing" };

    info!("Waiting for {} trend to stabilize for {} seconds...", trend_str, secs);

    let mut last = match rx.recv().await {
        Some(w) => {
            info!("Initial weight: {:.2}lb", w);
            w
        }
        None => {
            error!("No weight data received");
            return;
        }
    };

    while stable_count < secs {
        sleep(Duration::from_secs(1)).await;

        let mut current_weight = last;
        while let Ok(w) = rx.try_recv() {
            current_weight = w;
        }

        let delta = current_weight - last;
        let ok = if increasing { delta >= 0.0 } else { delta <= 0.0 };

        let status = if ok { "✓" } else { "✗" };
        let reset_msg = if !ok { " (reset counter)" } else { "" };

        debug!(
            "Weight: {:.2}lb, Delta: {:+.2}lb, Stable: {}/{} {}{}",
            current_weight,
            delta,
            stable_count + if ok { 1 } else { 0 },
            secs,
            status,
            reset_msg
        );

        stable_count = if ok { stable_count + 1 } else { 0 };
        last = current_weight;
    }

    info!("Stability achieved! Final weight: {:.2}lb", last);
}

async fn scale_reader(port_name: &str, tx: mpsc::Sender<f64>) -> Result<()> {
    info!("Opening serial port {}...", port_name);
    let mut port = serialport::new(port_name, 9600)
        .timeout(Duration::from_millis(1000))
        .open()?;
    info!("Scale reader started on {}", port_name);

    let mut extra_charmap = HashMap::new();
    extra_charmap.insert(0xa0u8, 0x20u8); // space
    extra_charmap.insert(0xb0u8, 0x30u8); // 0
    extra_charmap.insert(0xb1u8, 0x31u8); // 1
    extra_charmap.insert(0xb2u8, 0x32u8); // 2
    extra_charmap.insert(0xb3u8, 0x33u8); // 3
    extra_charmap.insert(0xb4u8, 0x34u8); // 4
    extra_charmap.insert(0xb5u8, 0x35u8); // 5
    extra_charmap.insert(0xb6u8, 0x36u8); // 6
    extra_charmap.insert(0xb7u8, 0x37u8); // 7
    extra_charmap.insert(0xb8u8, 0x38u8); // 8
    extra_charmap.insert(0xb9u8, 0x39u8); // 9

    let mut buffer = String::new();
    let pattern = Regex::new(r"\d{1,3}\.\d{2}(?:lb|l)")?;
    let mut last_weight: Option<f64> = None;

    loop {
        let mut read_buffer = vec![0; 1024];
        match port.read(&mut read_buffer) {
            Ok(bytes_read) if bytes_read > 0 => {
                read_buffer.truncate(bytes_read);
                let decoded = decode_ascii_with_extra(&read_buffer, &extra_charmap);
                let filtered = filter_printable(&decoded);

                if !filtered.trim().is_empty() {
                    buffer.push_str(&filtered);
                    buffer = remove_chars(&buffer, &['\r', '\n']);

                    if let Some(captures) = pattern.find(&buffer) {
                        let weight_str = captures.as_str().replace("lb", "").replace("l", "");
                        buffer.clear();
                        if let Ok(weight) = weight_str.parse::<f64>() {
                            if last_weight.map_or(true, |last| (weight - last).abs() > 0.01) {
                                debug!("Scale: {:.2}lb -> {:.2}lb", last_weight.unwrap_or(0.0), weight);
                                last_weight = Some(weight);
                            }

                            if tx.send(weight).await.is_err() {
                                warn!("Scale reader: Channel closed, exiting");
                                break;
                            }
                        }
                    }
                }
            }
            Ok(_) => {
                trace!("No data read from scale");
            }
            Err(ref e) if e.kind() == io::ErrorKind::TimedOut => {
                trace!("Scale read timeout");
            }
            Err(e) => {
                error!("Error reading from serial port: {}", e);
                break;
            }
        }
    }

    info!("Scale reader task terminated");
    Ok(())
}

fn decode_ascii_with_extra(data: &[u8], charmap: &HashMap<u8, u8>) -> String {
    let dot_pos = data.iter().position(|&b| b == b'.');
    if dot_pos.map_or(true, |pos| pos < 3) {
        return String::new();
    }

    let mut mapped_data = data.to_vec();
    for byte in &mut mapped_data {
        if let Some(&replacement) = charmap.get(byte) {
            *byte = replacement;
        }
    }

    String::from_utf8_lossy(&mapped_data).to_string()
}

fn filter_printable(s: &str) -> String {
    s.chars()
        .filter(|c| c.is_ascii_graphic() || c.is_ascii_whitespace())
        .collect()
}

fn remove_chars(s: &str, chars: &[char]) -> String {
    s.chars().filter(|c| !chars.contains(c)).collect()
}

impl Phase {
    async fn run(self, ctx: &mut Context) -> Result<Phase> {
        info!("Entering phase: {:?}", self);

        match self {
            Phase::StartPump => {
                info!("Starting pump...");
                psu_set(ctx.psu.clone(), 2, true).await?;
                info!("Pump ON (channel 2)");
                Ok(Phase::WaitStableDecrease)
            }
            Phase::WaitStableDecrease => {
                info!("Waiting for pumping to stabilize...");
                wait_for_stability(&mut ctx.weight_rx, true, ctx.stable_secs).await;
                Ok(Phase::StopPump)
            }
            Phase::StopPump => {
                info!("Stopping pump...");
                psu_set(ctx.psu.clone(), 2, false).await?;
                info!("Pump OFF (channel 2)");
                Ok(Phase::TareScale)
            }
            Phase::TareScale => {
                info!("Taring scale");
                Ok(Phase::StartFill)
            }
            Phase::StartFill => {
                info!("Starting fill process...");
                psu_set(ctx.psu.clone(), 1, true).await?;
                psu_set(ctx.psu.clone(), 3, true).await?;
                info!("Filling ON (channels 1 & 3)");
                Ok(Phase::WaitStableIncrease)
            }
            Phase::WaitStableIncrease => {
                info!("Collecting data during fill process...");
                wait_for_stability_with_collection(ctx, false, ctx.stable_secs).await;
                Ok(Phase::StopFill)
            }
            Phase::StopFill => {
                info!("Stopping fill process...");
                psu_set(ctx.psu.clone(), 1, false).await?;
                psu_set(ctx.psu.clone(), 3, false).await?;
                info!("Filling OFF (channels 1 & 3)");
                Ok(Phase::StoreData)
            }
            Phase::StoreData => {
                info!("Storing collected data...");
                ctx.flush_to_csv().await?;
                Ok(Phase::Done)
            }
            Phase::Done => {
                info!("Process completed successfully");
                Ok(Phase::Done)
            }
        }
    }
}

async fn wait_for_stability_with_collection(
    ctx: &mut Context,
    increasing: bool,
    secs: u64,
) {
    let mut stable_count = 0;
    let trend_str = if increasing { "increasing" } else { "decreasing" };

    info!(
        "Waiting for {} trend to stabilize for {} seconds while collecting data...",
        trend_str, secs
    );

    let mut last = match ctx.weight_rx.recv().await {
        Some(w) => {
            info!("Initial weight: {:.2}lb", w);
            w
        }
        None => {
            error!("No weight data received");
            return;
        }
    };

    let mut collection_cycle = 0;

    while stable_count < secs {
        sleep(Duration::from_secs(1)).await;

        let mut current_weight = last;
        while let Ok(w) = ctx.weight_rx.try_recv() {
            current_weight = w;
        }

        collection_cycle += 1;
        debug!("Collection cycle #{} at weight {:.2}lb", collection_cycle, current_weight);
        if let Err(e) = ctx.collect_cog_data(current_weight).await {
            error!("Error collecting COG data: {}", e);
        }

        let delta = current_weight - last;
        let ok = if increasing { delta >= 0.0 } else { delta <= 0.0 };

        let status = if ok { "✓" } else { "✗" };
        let reset_msg = if !ok { " (reset counter)" } else { "" };

        debug!(
            "Weight: {:.2}lb, Delta: {:+.2}lb, Stable: {}/{} {}{}",
            current_weight,
            delta,
            stable_count + if ok { 1 } else { 0 },
            secs,
            status,
            reset_msg
        );

        stable_count = if ok { stable_count + 1 } else { 0 };
        last = current_weight;
    }

    info!("Stability achieved! Final weight: {:.2}lb", last);
    info!("Total collection cycles: {}, Total records captured: {}", 
          collection_cycle, ctx.capture.len());
}

pub async fn run(
    psu_port: String,
    scale_port: String, 
    gtr_port: String,
    stable_secs: u64,
    output_file: String
) -> Result<()> {
    SimpleLogger::new()
        .with_level(log::LevelFilter::Info)
        .env()
        .init() 
        .unwrap();

    info!("Starting data collection");
    info!("PSU: {}, Scale: {}, GTR: {}", psu_port, scale_port, gtr_port);

    let (wt_tx, wt_rx) = mpsc::channel(128);
    tokio::spawn(async move {
        if let Err(e) = scale_reader(&scale_port, wt_tx).await {
            error!("Scale reader error: {}", e);
        }
    });

    info!("Connecting to PSU...");
    let psu = Arc::new(Mutex::new(mxpsu::MxSeries::connect_serial(&psu_port, 9600)?));
    info!("PSU connected successfully");

    let mut ctx = Context::new(wt_rx, psu.clone(), stable_secs, output_file);

    info!("Initializing GTR...");
    ctx.init_gtr(&gtr_port).await?;

    let mut state = Phase::StartPump;

    info!("Starting state machine...");
    while state != Phase::Done {
        state = state.run(&mut ctx).await?;
    }

    info!("Shutting down GTR...");
    ctx.stop_gtr().await?;

    info!("Data collection completed successfully");
    Ok(())
}