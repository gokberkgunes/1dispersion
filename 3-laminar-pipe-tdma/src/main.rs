use anyhow::Result;
use clap::{error::ErrorKind, CommandFactory, Parser, ValueEnum};

use ndarray::Array1;

/// Solves laminar flow between two flat plates using central differences with 1D grid with
/// Thomas algorithm (TDMA),
/// d^2u/dy^2 = 1/mu * dp/dx.
#[derive(Parser)]
#[command(author, version, about)]
struct Args {
    #[arg(
        //required = false,
        short = 'b',
        long = "boundary-condition-type",
        default_value_t = BoundaryCondition::Neumann,
        value_enum,
    )]
    boundary_condition_type: BoundaryCondition,
    /// Boundary condition value, should only be
    #[arg(
        required_if_eq("boundary_condition_type", "dirichlet"),
        short = 'v',
        long = "boundary-value"
    )]
    boundary_value: Option<f64>,

    /// Number of points between the wall and midpoint between walls
    #[arg(
        //required = false,
        short = 'n',
        long = "number-points",
        default_value_t = 10,
    )]
    n_points: usize,

    /// pressure drop gradient, dp/dx or C_p.
    #[arg(
        //required = false,
        short = 'd',
        long = "pressure-gradient",
        default_value_t = -1.0,
    )]
    dpdx: f64,

    /// Dynamic viscosity of the fluid.
    #[arg(
        //required = false,
        short = 'm',
        long = "viscosity",
        default_value_t = 1e-3,
    )]
    mu: f64,

    /// Distance between plate to centerline.
    #[arg(
        //required = false,
        short = 'H',
        long = "half-height",
        default_value_t = 0.1,
    )]
    half_height: f64,
}

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
enum BoundaryCondition {
    Neumann,
    Dirichlet,
}

fn tdma_prepare_flat_plate(
    bc_type: BoundaryCondition,
    bc_value: f64,
    n: usize,
    dpdx: f64,
    h: f64,
    mu: f64,
) -> Result<String> {
    let dy: f64 = h / ((n - 1) as f64);
    let a_weight = -1.;
    let b_weight = 2.;
    let c_weight = -1.;
    let d_weight = -dpdx * dy.powi(2) / mu;
    let mut a = Array1::<f64>::from_elem(n, a_weight); // first a doesn't exist.
    let mut b = Array1::<f64>::from_elem(n, b_weight);
    let mut c = Array1::<f64>::from_elem(n, c_weight); // last c doesn't exist.
    let mut d = Array1::<f64>::from_elem(n, d_weight);

    // Boundary condition for wall
    // Accept that bottom is wall, u_1 = 0.
    b[0] = 1.;
    c[0] = 0.;
    d[0] = 0.;

    // We modify last equation according to the boundary condition.
    // For Neumann we have (u_N - u_(N-1)) / dy = B.C.
    if bc_type == BoundaryCondition::Neumann {
        // Second order with ghost cell.
        // u_(N-1) = u_(N+1). Then, [u_(N+1) - u_(N-1)]/(2dy) = B.C., where B.C. should be 0 at
        // pipe centerline.
        //
        // Then, A_N * u_(N-1) + B_N u_N + C_N u_(N+1) = d can be rewritten as
        // (A_N + C_N) u_(N-1) + B_N u_N = d - C_N * B.C. * 2 * dy
        a[n-1] = a_weight + c_weight;
        b[n-1] = b_weight;
        d[n-1] = d_weight - c_weight * bc_value * 2. * dy;

        // Below is first order backward.
        // u_N - u_(N-1) = B.C.
        // a = -1, b = 1, d = B.C. * dy
        //d_prime[n-1] = (bc_value*dy - (-1.)*d_prime[n-2])/(1. - (-1.)*c_prime[n-2])
    } else {
        d[n-1] = bc_value;
    }

    let u = tdma(&a, &b, &c, &d);
    println!("{:?}", u);

    // Report CFL and Diffusion numbers
    Ok(format!("Solved."))
}

/// TDMA has a Russian variation called Progorka method that envisions that
/// solutions can be reduced to a form of u_i = E_i u_(i+1) + F_i.
///
/// However, TDMA is much more straightforward to understand, and reduce to
/// the same method.
/// This method does row-reduction starting from first row to the last row.
/// Then do reverse sweep: find velocities from end one by one.
/// This algorithm does not require to calculate A and B of weights matrix:
/// A reduces to 0, B reduces to 1.
// Solves: A[i]*x[i-1] + B[i]*x[i] + C[i]*x[i+1] = D[i]
fn tdma(
    a: &Array1<f64>,
    b: &Array1<f64>,
    c: &Array1<f64>,
    d: &Array1<f64>
) -> Array1<f64>
{
    let n = a.len();
    let mut x = Array1::<f64>::zeros(n);
    let mut c_prime = Array1::<f64>::zeros(n);
    let mut d_prime = Array1::<f64>::zeros(n);
    let mut inv_denom;

    // a[0] doesn't exist for first row.
    // Normalize c[0] and d[0]. b[0] reduces to 1 anyway.
    c_prime[0] = c[0]/b[0];
    d_prime[0] = d[0]/b[0];

    for i in 1..n {
        inv_denom = 1./(b[i] - a[i] * c_prime[i-1]);
        c_prime[i] = c[i] * inv_denom; // c[n-1] is unneeded
        d_prime[i] = (d[i] - a[i] * d_prime[i-1]) * inv_denom;
    }

    x[n-1] = d_prime[n-1];
    for i in (0..n-1).rev() {
        x[i] = d_prime[i] - c_prime[i] * x[i+1];
    }

    return x;
}

fn run(args: Args) -> Result<()> {
    let n = args.n_points;
    let bc_type = args.boundary_condition_type;
    let dpdx = args.dpdx;
    let h = args.half_height;
    let mu = args.mu;

    let bc_value = match args.boundary_condition_type {
        BoundaryCondition::Dirichlet => args.boundary_value.unwrap(),
        BoundaryCondition::Neumann => 0.0,
    };
    if args.boundary_condition_type == BoundaryCondition::Neumann && args.boundary_value.is_some() {
        Args::command()
            .error(
                ErrorKind::ArgumentConflict,
                "The argument '--boundary-condition-value' cannot be used with 'Neumann' type.",
            )
            .exit();
    }

    tdma_prepare_flat_plate(bc_type, bc_value, n, dpdx, h, mu)?;
    Ok(())
}

fn main() {
    let args = Args::parse();
    if let Err(e) = run(args) {
        eprintln!("Error: {e}.");
        std::process::exit(1);
    }
}
