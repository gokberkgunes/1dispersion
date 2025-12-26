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
    /// Spill Location
    #[arg(
        short('p'),
        long,
        default_value_t = 500.
    )]
    spill_location: f64,
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

// Von Neumann Analysis.
// If this is invalid, solution may diverge locally. We must have amplification
// factor sigma <= 1.
// fn check_stability (dx: f64, dt: f64, u: f64, k_x: f64) {
//
// }

// TODO: Avoid calculating dx and dt, take them from user or enable program
// for suggestions only. Calculating dx and dt should not be a task of this
// program at this point.
fn run(args: Args) -> Result<()> {
    let alpha: f64;
    let dx: f64;
    let dt: f64;
    let k_x: f64;
    let safety_factor = 0.8;

    if args.scheme == Scheme::Central {
        alpha = 0.5;
        // Cell Peclet number should be less than 2. Below uses this condition.
        //dx = safety_factor * args.k_x/((1.0-alpha) * args.u);
    } else {
        alpha = 1.0;
        // For backward difference method, spatial discritization is unlimited.
        // This is because we do not use downstream information to calculate 
        // advection anymore, i.e., upwinding.
        //dx = 2.0 * args.k_x * args.u;
    }
    // Trial Values
    dx = 50.0;
    dt = 0.005;
    // Diffusion number should be less than 0.5; however, we use 0.25 to avoid
    // oscillations.
    //dt = safety_factor * 0.25*dx.powi(2) / args.k_x;

    k_x = calc_dispersion(args.k_x, alpha, args.u, dx, dt);

    // number of points
    let n = (args.length/dx).round() as usize;
    // number of timestep
    let num_timestep =(args.t_target/dt) as usize;
    //let grid: Vec<f64> = (0..=n).map(|i| (i as f64)*dx).collect();
    // injection point of the total mass
    let injection_idx = (args.spill_location/dx) as usize;

    let mut c: Vec<f64> = (0..=n).map(|_i| 0.0).collect();
    c[injection_idx] = args.mass / (dx * args.width * args.height);

    let mut c_new: Vec<f64> = (0..=n).map(|i| (i as f64)*dx).collect();

    for _ in 0..=num_timestep {
        for j in 1..n {
            c_new[j] = c[j] - args.u*dt*(c[j+1] - c[j-1])/(2.*dx)
                + k_x*dt*(c[j+1]-2.0*c[j]+c[j-1])/(dx.powi(2));
        }

        // Inlet, dirichlet b.c.
        //c_new[0] = 0.0;

        // Inlet Ghost cell b.c. (central)
        c_new[0] = c[0] - args.u*dt * (c[0] - 0.0)/dx + k_x*dt*(c[1] -2.0*c[0] + 0.0)/(dx.powi(2));
        c_new[n] = c_new[n-1];

        // Swap c_new with c. Data of c_new isn't needed anyway.
        std::mem::swap(&mut c, &mut c_new);
    }

    println!("{:?}", c);

    Ok(())
}


fn main() {
    let args = Args::parse();
    if let Err(e) = run(args) {
        eprintln!("Error: {e}.");
        std::process::exit(1);
    }
}
