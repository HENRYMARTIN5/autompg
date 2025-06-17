use anyhow::Result;
use plotters::prelude::*;
use csv::Reader;
use serde::Deserialize;
use std::f64::consts::PI;

#[derive(Debug, Deserialize)]
struct Record {
    timestamp: String,
    weight: f64,
    cog_1: f64,
    cog_2: f64,
    cog_3: f64,
    cog_4: f64,
    cog_5: f64,
    cog_6: f64,
    cog_7: f64,
    cog_8: f64,
}

pub fn run(input_file: String, output_file: String, width: u32, height: u32) -> Result<()> {
    println!("Loading data from: {}", input_file);
    
    // Read CSV data
    let mut rdr = Reader::from_path(&input_file)?;
    let mut records = Vec::new();
    
    for result in rdr.deserialize() {
        let record: Record = result?;
        records.push(record);
    }
    
    if records.is_empty() {
        println!("No data found in {}", input_file);
        return Ok(());
    }
    
    println!("Loaded {} records", records.len());
    
    // Create main plot
    create_time_series_plot(&records, &output_file, width, height)?;
    
    // Create FFT surface plots for each COG sensor with data
    create_fft_surface_plots(&records, width, height)?;
    
    // Print summary statistics
    print_summary(&records);
    
    Ok(())
}

fn create_time_series_plot(records: &[Record], output_file: &str, width: u32, height: u32) -> Result<()> {
    let root = BitMapBackend::new(output_file, (width, height)).into_drawing_area();
    root.fill(&WHITE)?;
    
    // Split the drawing area - weight on top, COG sensors below (fix: use (2, 1) for 2 rows, 1 column)
    let areas = root.split_evenly((2, 1));
    let upper = &areas[0];
    let lower = &areas[1];
    
    // Upper plot: Weight
    let mut weight_chart = ChartBuilder::on(upper)
        .caption("Weight Over Time", ("sans-serif", 30))
        .margin(5)
        .x_label_area_size(40)
        .y_label_area_size(60)
        .build_cartesian_2d(
            0f64..records.len() as f64,
            0f64..records.iter().map(|r| r.weight).fold(0.0, f64::max) + 0.5
        )?;
    
    weight_chart
        .configure_mesh()
        .x_desc("Time (samples)")
        .y_desc("Weight (lb)")
        .draw()?;
    
    weight_chart.draw_series(
        LineSeries::new(
            records.iter().enumerate().map(|(i, r)| (i as f64, r.weight)),
            BLUE.stroke_width(2),
        )
    )?
    .label("Weight")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 10, y)], BLUE));
    
    weight_chart.configure_series_labels().draw()?;
    
    // Lower plot: COG sensors with individual scaling
    create_cog_subplots(lower, records)?;
    
    root.present()?;
    println!("Time series plot saved to: {}", output_file);
    Ok(())
}

fn create_cog_subplots(area: &DrawingArea<BitMapBackend, plotters::coord::Shift>, records: &[Record]) -> Result<()> {
    let colors = [RED, GREEN, MAGENTA, CYAN, BLACK, RGBColor(255, 165, 0), RGBColor(128, 0, 128), RGBColor(255, 255, 0)];
    
    // Create a 2x4 grid for the 8 COG sensors
    let areas = area.split_evenly((2, 4));
    
    for sensor in 1..=8 {
        let cog_values: Vec<f64> = records.iter().map(|r| get_cog_value(r, sensor)).collect();
        
        // Skip sensors with no data or no variation
        if cog_values.iter().all(|&v| v == 0.0) {
            continue;
        }
        
        let min_cog = cog_values.iter().fold(f64::INFINITY, |a, &b| a.min(b));
        let max_cog = cog_values.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
        let range = max_cog - min_cog;
        
        // Add 5% padding to the range for better visualization
        let padding = range * 0.05;
        let y_min = min_cog - padding;
        let y_max = max_cog + padding;
        
        // If there's virtually no variation, create a narrow range around the mean
        let (final_y_min, final_y_max) = if range < min_cog * 0.001 { // Less than 0.1% variation
            let center = (min_cog + max_cog) / 2.0;
            let artificial_range = center * 0.001; // 0.1% of the center value
            (center - artificial_range, center + artificial_range)
        } else {
            (y_min, y_max)
        };
        
        let area_idx = sensor - 1;
        let mut chart = ChartBuilder::on(&areas[area_idx])
            .caption(&format!("COG {}", sensor), ("sans-serif", 16))
            .margin(2)
            .x_label_area_size(25)
            .y_label_area_size(40)
            .build_cartesian_2d(
                0f64..records.len() as f64,
                final_y_min..final_y_max
            )?;
        
        chart
            .configure_mesh()
            .x_desc("Samples")
            .y_desc("Value")
            .x_label_formatter(&|x| format!("{:.0}", x))
            .y_label_formatter(&|y| format!("{:.0}", y))
            .draw()?;
        
        let cog_data: Vec<(f64, f64)> = records.iter().enumerate().map(|(i, r)| {
            (i as f64, get_cog_value(r, sensor))
        }).collect();
        
        chart.draw_series(
            LineSeries::new(cog_data, colors[(sensor - 1) % colors.len()].stroke_width(1))
        )?;
    }
    
    Ok(())
}

fn create_fft_surface_plots(records: &[Record], width: u32, height: u32) -> Result<()> {
    let window_size = 256;
    let overlap = window_size / 2;
    let step_size = window_size - overlap;
    
    for sensor in 1..=8 {
        let cog_data: Vec<f64> = records.iter().map(|r| get_cog_value(r, sensor)).collect();
        
        if cog_data.iter().all(|&v| v == 0.0) {
            continue;
        }
        
        println!("Creating FFT surface plot for COG {}", sensor);
        
        // Just remove the mean, don't do complex detrending
        let mean = cog_data.iter().sum::<f64>() / cog_data.len() as f64;
        let centered_data: Vec<f64> = cog_data.iter().map(|&x| x - mean).collect();
        
        let mut fft_data = Vec::new();
        let mut time_points = Vec::new();
        
        for start in (0..centered_data.len().saturating_sub(window_size)).step_by(step_size) {
            let window = &centered_data[start..start + window_size];
            let fft_result = compute_fft(window);
            fft_data.push(fft_result);
            time_points.push(start + window_size / 2);
        }
        
        if fft_data.is_empty() {
            continue;
        }
        
        let output_file = format!("cog_{}_fft_surface.png", sensor);
        create_surface_plot(&fft_data, &time_points, &output_file, width, height, sensor)?;
    }
    
    Ok(())
}

fn compute_fft(data: &[f64]) -> Vec<f64> {
    let n = data.len();
    
    // Apply Hanning window to reduce spectral leakage
    let windowed: Vec<f64> = data.iter().enumerate().map(|(i, &x)| {
        let window = 0.5 * (1.0 - (2.0 * PI * i as f64 / (n - 1) as f64).cos());
        x * window
    }).collect();
    
    // Simple DFT implementation
    let mut magnitude = Vec::new();
    
    for k in 1..n/2 { // Skip DC (k=0) and only positive frequencies
        let mut real = 0.0;
        let mut imag = 0.0;
        
        for j in 0..n {
            let angle = -2.0 * PI * (k as f64) * (j as f64) / (n as f64);
            real += windowed[j] * angle.cos();
            imag += windowed[j] * angle.sin();
        }
        
        magnitude.push((real * real + imag * imag).sqrt());
    }
    
    magnitude
}

fn create_surface_plot(
    fft_data: &[Vec<f64>], 
    time_points: &[usize], 
    output_file: &str, 
    width: u32, 
    height: u32,
    sensor: usize
) -> Result<()> {
    let root = BitMapBackend::new(output_file, (width, height)).into_drawing_area();
    root.fill(&WHITE)?;
    
    if fft_data.is_empty() {
        return Ok(());
    }
    
    let freq_bins = fft_data[0].len();
    let time_bins = fft_data.len();
    
    // Collect all magnitude values for statistical analysis
    let all_magnitudes: Vec<f64> = fft_data.iter().flat_map(|row| row.iter().copied()).collect();
    
    // Calculate percentiles to filter outliers
    let mut sorted_magnitudes = all_magnitudes.clone();
    sorted_magnitudes.sort_by(|a, b| a.partial_cmp(b).unwrap());
    
    let p95_idx = (sorted_magnitudes.len() as f64 * 0.95) as usize;
    let p5_idx = (sorted_magnitudes.len() as f64 * 0.05) as usize;
    
    let min_magnitude = sorted_magnitudes[p5_idx];
    let max_magnitude = sorted_magnitudes[p95_idx.min(sorted_magnitudes.len() - 1)];
    
    println!("COG {} FFT: Min={:.1}, Max={:.1}, P95={:.1}, {} outliers clipped", 
             sensor, 
             sorted_magnitudes[0], 
             sorted_magnitudes[sorted_magnitudes.len()-1],
             max_magnitude,
             all_magnitudes.iter().filter(|&&x| x > max_magnitude).count());
    
    let mut chart = ChartBuilder::on(&root)
        .caption(&format!("COG {} - FFT Spectrogram (5%-95% Range: {:.1} - {:.1})", sensor, min_magnitude, max_magnitude), ("sans-serif", 30))
        .margin(10)
        .x_label_area_size(60)
        .y_label_area_size(80)
        .build_cartesian_2d(
            0f64..time_bins as f64,
            0f64..freq_bins as f64
        )?;
    
    chart
        .configure_mesh()
        .x_desc("Time Window")
        .y_desc("Frequency Bin")
        .draw()?;
    
    // Create high-resolution interpolated surface
    let interpolation_factor = 4; // 4x interpolation for smoother surface
    let interp_time_bins = time_bins * interpolation_factor;
    let interp_freq_bins = freq_bins * interpolation_factor;
    
    // Calculate the size of each interpolated cell in chart coordinates
    let time_step = time_bins as f64 / interp_time_bins as f64;
    let freq_step = freq_bins as f64 / interp_freq_bins as f64;
    
    for it in 0..interp_time_bins {
        for ift in 0..interp_freq_bins {
            // Map interpolated coordinates to original data coordinates
            let time_coord = (it as f64 + 0.5) * time_step;
            let freq_coord = (ift as f64 + 0.5) * freq_step;
            
            // Bilinear interpolation
            let magnitude = bilinear_interpolate(&fft_data, time_coord, freq_coord);
            
            // Clip outliers
            let clipped_magnitude = magnitude.min(max_magnitude).max(min_magnitude);
            
            let normalized_mag = if max_magnitude > min_magnitude {
                (clipped_magnitude - min_magnitude) / (max_magnitude - min_magnitude)
            } else {
                0.0
            };
            
            // Improved color gradient with gamma correction for better perception
            let gamma_corrected = normalized_mag.powf(0.5); // Gamma correction
            
            let color = if gamma_corrected < 0.25 {
                let t = gamma_corrected / 0.25;
                RGBColor(0, 0, (128.0 + 127.0 * t) as u8) // Dark blue to blue
            } else if gamma_corrected < 0.5 {
                let t = (gamma_corrected - 0.25) / 0.25;
                RGBColor(0, (255.0 * t) as u8, 255) // Blue to cyan
            } else if gamma_corrected < 0.75 {
                let t = (gamma_corrected - 0.5) / 0.25;
                RGBColor((255.0 * t) as u8, 255, (255.0 * (1.0 - t)) as u8) // Cyan to yellow
            } else {
                let t = (gamma_corrected - 0.75) / 0.25;
                RGBColor(255, (255.0 * (1.0 - t * 0.8)) as u8, 0) // Yellow to red
            };
            
            // Draw filled rectangle for this interpolated cell
            let x1 = it as f64 * time_step;
            let x2 = (it + 1) as f64 * time_step;
            let y1 = ift as f64 * freq_step;
            let y2 = (ift + 1) as f64 * freq_step;
            
            chart.draw_series(std::iter::once(Rectangle::new([
                (x1, y1), 
                (x2, y2)
            ], color.filled())))?;
        }
    }
    
    // Add outlier information
    let outlier_count = all_magnitudes.iter().filter(|&&x| x > max_magnitude).count();
    chart.draw_series(std::iter::once(Text::new(
        format!("{} outliers clipped", outlier_count), 
        (time_bins as f64 * 0.7, freq_bins as f64 * 0.05), 
        ("sans-serif", 16).into_font().color(&BLACK)
    )))?;
    
    root.present()?;
    println!("FFT surface plot saved to: {} (high-resolution interpolated surface)", output_file);
    Ok(())
}

fn bilinear_interpolate(data: &[Vec<f64>], x: f64, y: f64) -> f64 {
    let x_max = data.len() as f64 - 1.0;
    let y_max = data[0].len() as f64 - 1.0;
    
    // Clamp coordinates
    let x_clamped = x.max(0.0).min(x_max);
    let y_clamped = y.max(0.0).min(y_max);
    
    // Get integer and fractional parts
    let x0 = x_clamped.floor() as usize;
    let y0 = y_clamped.floor() as usize;
    let x1 = (x0 + 1).min(data.len() - 1);
    let y1 = (y0 + 1).min(data[0].len() - 1);
    
    let dx = x_clamped - x0 as f64;
    let dy = y_clamped - y0 as f64;
    
    // Get corner values
    let v00 = data[x0][y0];
    let v10 = data[x1][y0];
    let v01 = data[x0][y1];
    let v11 = data[x1][y1];
    
    // Bilinear interpolation
    let v0 = v00 * (1.0 - dx) + v10 * dx;
    let v1 = v01 * (1.0 - dx) + v11 * dx;
    
    v0 * (1.0 - dy) + v1 * dy
}

fn get_cog_value(record: &Record, sensor: usize) -> f64 {
    match sensor {
        1 => record.cog_1,
        2 => record.cog_2,
        3 => record.cog_3,
        4 => record.cog_4,
        5 => record.cog_5,
        6 => record.cog_6,
        7 => record.cog_7,
        8 => record.cog_8,
        _ => 0.0,
    }
}

fn print_summary(records: &[Record]) {
    if records.is_empty() {
        return;
    }
    
    let weights: Vec<f64> = records.iter().map(|r| r.weight).collect();
    let min_weight = weights.iter().fold(f64::INFINITY, |a, &b| a.min(b));
    let max_weight = weights.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
    let avg_weight = weights.iter().sum::<f64>() / weights.len() as f64;
    
    println!("\n=== DATA SUMMARY ===");
    println!("Total records: {}", records.len());
    println!("Weight range: {:.2} - {:.2} lb", min_weight, max_weight);
    println!("Average weight: {:.2} lb", avg_weight);
    println!("Weight change: {:.2} lb", max_weight - min_weight);
    
    // COG analysis
    println!("\n=== COG ANALYSIS ===");
    for sensor in 1..=8 {
        let cog_values: Vec<f64> = records.iter().map(|r| get_cog_value(r, sensor)).collect();
        
        let non_zero_count = cog_values.iter().filter(|&&v| v != 0.0).count();
        if non_zero_count > 0 {
            let min_cog = cog_values.iter().fold(f64::INFINITY, |a, &b| if b != 0.0 { a.min(b) } else { a });
            let max_cog = cog_values.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
            let avg_cog = cog_values.iter().filter(|&&v| v != 0.0).sum::<f64>() / non_zero_count as f64;
            let variation = ((max_cog - min_cog) / avg_cog) * 100.0;
            
            println!("COG {}: {:.0} samples, range {:.0} - {:.0}, avg {:.0}, variation {:.3}%", 
                     sensor, non_zero_count, min_cog, max_cog, avg_cog, variation);
        } else {
            println!("COG {}: No data", sensor);
        }
    }
    
    println!("\n=== FFT SURFACE PLOTS ===");
    println!("Generated FFT spectrograms for active COG sensors");
    println!("Files: cog_N_fft_surface.png (where N = sensor number)");
    println!("Features: Mean-centered data, linear scale, individual COG subplots");
}