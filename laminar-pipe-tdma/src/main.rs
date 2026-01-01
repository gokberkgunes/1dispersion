use clap::{Parser, ValueEnum};

use ndarray::{Array};

/// Solves laminar pipe flow with 1D grid usign Thomas algorithm (TDMA),
/// d^2u/dy^2 = 1/mu * dp/dx.
#[derive(Parser)]
#[command(author, version, about)]
struct Args {
    #[arg(
        //required = false,
        short = 's',
        long = "scheme",
        default_value_t = BoundaryCondition::Dirichlet,
        value_enum,
    )]
    scheme: BoundaryCondition,

}

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
enum BoundaryCondition {
    Neumann,
    Dirichlet
}

const DP_DX: f64 = -1.0;
const H: f64 = 0.1;
const MU: f64 = 1e-3;
const N: [usize; 4] =[10, 1000, 3000, 5000];

fn tdma(bc: BoundaryCondition, n: usize) -> Result<String, String> {
    let dy: f64 = H/(n as f64);
    let a = 1.0;
    let b = 2.0;
    let c = 1.0;
    let d = -DP_DX*dy.powi(2)/MU;
    let mut e = Array::zeros(n);
    let mut f = Array::zeros(n);
    let mut u = Array::zeros(n);

    if bc == BoundaryCondition::Neumann {
    }
    // Report CFL and Diffusion numbers
    Ok(format!(""))
}

fn run(args: Args) -> Result<(), String> {

    Ok(())
}


fn main() {
    let args = Args::parse();
    if let Err(e) = run(args) {
        eprintln!("Error: {e}.");
        std::process::exit(1);
    }
}
