use std::fs::File;
use std::io::{BufWriter, Write};

/// Calculates friction factor in a pipe using Reynolds number, roughness of
/// pipe and pipe diameter.
fn swamee_jean(roughness: f64, u: f64, diameter: f64, rho: f64, mu: f64) -> f64 {
    let reynolds_no = u * diameter * rho / mu;
    return 0.25
        / (((roughness / (3.7 * diameter) + 5.74 / reynolds_no.powf(0.9)).log10()).powi(2));
}

// Perhaps prefer some library? Newton-Rhapson, Brent, Secant, or Bisection
// Apply bisection method to get steady state velocity.
// Iterate to solve dV/dt = H - (1 + f * L/D) * V^2/(2g)) * g/L, for dV/dt = 0.
// Equation reduces to 0 = H - (1 + f * L/D) * V^2/(2g)
// Iterations are needed since f and V are unknown and dependent on each other.
fn bisection<F>(mut a: f64, mut b: f64, f: F, tol: f64, max_iter: usize) -> Result<f64, String>
where
    F: Fn(f64) -> f64,
{
    let mut mid = 0.0;
    let mut fmid;
    let mut fa = f(a);

    if fa * f(b) > 0. {
        return Err("Root may not be found since f(a) * f(b) > 0".into());
    }

    for _ in 0..max_iter {
        mid = (a + b) / 2.0;
        fmid = f(mid);

        // If our function is close to 0, return the root.
        if fmid.abs() < tol {
            return Ok(mid);
        }

        if fa.signum() == fmid.signum() {
            a = mid;
            fa = fmid; // avoids double f(a) calculation
        } else {
            b = mid;
        }
    }

    Err(format!(
        "Reached maximum iterations, current root guess: {}",
        mid
    ))
}

fn solve<F>(f: F, dt: f64, total_time: f64, u_steady: f64) -> std::io::Result<()>
where
    F: Fn(f64) -> f64,
{
    let mut time = 0.0;
    let mut vel = 0.0;
    let mut variation = (u_steady - vel) / u_steady;
    let file = File::create(format!("{:.2}.out", u_steady))?;
    let mut writer = BufWriter::new(file); // use BufWriter. buffered should be faster

    writeln!(writer, "{vel:.8} {time:.2} {variation:.4}")?;
    while time <= total_time {
        vel += f(vel);
        variation = (u_steady - vel) / u_steady;
        time += dt;
        writeln!(writer, "{vel:.8} {time:.2} {variation:.4}")?;
    }

    Ok(())
}

fn main() {
    let g = 9.81;
    let h = 8.0;
    let rho = 1000.0;
    let diam = 0.3;
    let dt = 0.02;
    let t = 60.0;
    let ks = 0.0001;
    let mu = 0.001;
    let mut v = vec![0., 0., 0., 0.];
    let length = vec![50., 100., 200., 500.];

    for i in 0..length.len() {
        let equation = |u: f64| {
            h - (1. + swamee_jean(ks, u, diam, rho, mu) * length[i] / diam) * u.powi(2) / (2. * g)
        };
        match bisection(0.01, (2.0 * g * h).sqrt(), equation, 1e-5, 100) {
            Ok(root) => v[i] = root,
            Err(e) => println!("ERROR: Bisection {e}."),
        }
    }

    for i in 0..length.len() {
        let equation = |u: f64| {
            dt * (h
                - (1. + swamee_jean(ks, u, diam, rho, mu) * length[i] / diam) * u.powi(2)
                    / (2. * g))
                * g
                / length[i]
        };
        if let Err(e) = solve(equation, dt, t, v[i]) {
            eprintln!("Error while solving: {}", e);
        }
    }
}
