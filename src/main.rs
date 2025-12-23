use anyhow::Result;
use clap::{Parser, ValueEnum};

// avoid vec!, use ndarray.
// use ndarray::{Array1, Array2, Axis};

#[derive(Parser)]
#[command(author, version, about)]
struct Args {
    /// advective scheme
    #[arg(
        //required = false,
        short = 's',
        long = "scheme",
        default_value_t = Scheme::Central,
        value_enum,
    )]
    scheme: Scheme,

    /// Mass of the pollutant
    #[arg(
        short,
        long,
        default_value_t = 5000.
    )]
    mass: f64,

    /// Width of the channel
    #[arg(
        short,
        long,
        default_value_t = 60.
    )]
    width: f64,

    /// Height of the channel
    #[arg(
        short('H'),
        long("height"),
        default_value_t = 1.
    )]
    height: f64,

    /// Length of the channel
    #[arg(
        short,
        long,
        default_value_t = 3500.
    )]
    length: f64,

    /// Target time
    #[arg(
        short,
        long,
        default_value_t = 0.2
    )]
    t_target: f64,

    /// Dispersion coefficient
    #[arg(
        short,
        long,
        default_value_t = 150000.
    )]
    k_x: f64,

    /// Mean flow speed
    #[arg(
        short,
        long,
        default_value_t = 2400.
    )]
    u: f64,
}

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
enum Scheme {
    Central,
    Backward
}

fn calc_dispersion (
    k_x: f64,
    alpha: f64,
    u: f64,
    dx: f64,
    dt: f64) -> f64 {

    return k_x + u*(-dx*(alpha-0.5) + u*dt/2.0);
}

fn check_stability (dx: f64, dt: f64, u: f64, k_x: f64) {

}

fn run(args: Args) -> Result<()> {
    //
    Ok(())
}


fn main() {
    let args = Args::parse();
    if let Err(e) = run(args) {
        eprintln!("Error: {e}.");
        std::process::exit(1);
    }
}
