use clap::{Parser, Subcommand};
use anyhow::Result;

mod collect;
mod viz;

#[derive(Parser)]
#[command(name = "autompg")]
#[command(about = "Automatic MPG tank measurement and analysis")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Collect weight and COG data from hardware
    Collect {
        /// PSU serial port
        #[arg(long, default_value = "/dev/ttyACM0")]
        psu_port: String,
        
        /// Scale serial port
        #[arg(long, default_value = "/dev/ttyUSB2")]
        scale_port: String,
        
        /// GTR serial port
        #[arg(long, default_value = "/dev/ttyUSB0")]
        gtr_port: String,
        
        /// Stability time in seconds
        #[arg(long, default_value = "5")]
        stable_secs: u64,
        
        /// Output CSV file
        #[arg(long, default_value = "tank_run.csv")]
        output: String,
    },
    /// Visualize collected data
    Viz {
        /// Input CSV file
        #[arg(long, default_value = "tank_run.csv")]
        input: String,
        
        /// Output image file
        #[arg(long, default_value = "tank_analysis.png")]
        output: String,
        
        /// Image width
        #[arg(long, default_value = "1200")]
        width: u32,
        
        /// Image height
        #[arg(long, default_value = "800")]
        height: u32,
    },
}

#[tokio::main]
async fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Collect { 
            psu_port, 
            scale_port, 
            gtr_port, 
            stable_secs, 
            output 
        } => {
            collect::run(psu_port, scale_port, gtr_port, stable_secs, output).await
        }
        Commands::Viz { 
            input, 
            output, 
            width, 
            height 
        } => {
            viz::run(input, output, width, height)
        }
    }
}