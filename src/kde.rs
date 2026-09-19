//! Binned kernel density estimate and its derivatives (R `drvkde_1d`, `symconv_1d`).
//!
//! The density and its **second** derivative are what the peak test needs: a peak is a point
//! where the curvature is significantly negative. Both come from one convolution of the binned
//! counts with a Gaussian-derivative kernel, done through an FFT exactly as R does it.
use std::sync::Arc;

use rustfft::{FftPlanner, num_complex::Complex};

use crate::binning;

/// Symmetric convolution through the FFT (R `symconv_1d`).
///
/// `rr` is the half-kernel (index 0 is the centre), `ss` the binned counts. `skewflag` is
/// `(-1)^drv`: the kernel's mirrored half changes sign for odd derivatives.
pub fn symconv(rr: &[f64], ss: &[f64], skewflag: f64) -> Vec<f64> {
    let l = rr.len() - 1;
    let m = ss.len();
    let p = (m + l).next_power_of_two();

    let mut rp = vec![Complex::new(0.0, 0.0); p];
    for (i, &v) in rr.iter().enumerate() {
        rp[i] = Complex::new(v, 0.0);
    }
    if l > 0 {
        // R: rp[(P - L + 1):P] <- skewflag * rr[(L + 1):2] — the mirrored half, reversed.
        for k in 1..=l {
            rp[p - k] = Complex::new(skewflag * rr[k], 0.0);
        }
    }
    let mut sp = vec![Complex::new(0.0, 0.0); p];
    for (i, &v) in ss.iter().enumerate() {
        sp[i] = Complex::new(v, 0.0);
    }

    let mut planner = FftPlanner::new();
    let fwd: Arc<dyn rustfft::Fft<f64>> = planner.plan_fft_forward(p);
    let inv: Arc<dyn rustfft::Fft<f64>> = planner.plan_fft_inverse(p);
    fwd.process(&mut rp);
    fwd.process(&mut sp);
    let mut prod: Vec<Complex<f64>> = rp.iter().zip(&sp).map(|(a, b)| a * b).collect();
    inv.process(&mut prod);
    // R's fft(inverse = TRUE) is unnormalised, and the port divides by P.
    prod[..m].iter().map(|c| c.re / p as f64).collect()
}

pub struct Drvkde {
    pub grid: Vec<f64>,
    pub est: Vec<f64>,
}

/// Standard normal density.
fn dnorm(x: f64) -> f64 {
    (-0.5 * x * x).exp() / (2.0 * std::f64::consts::PI).sqrt()
}

/// Binned KDE of derivative order `drv` (R `drvkde_1d`, `binned = TRUE`).
pub fn drvkde(gcounts: &[f64], drv: usize, h: f64, a: f64, b: f64) -> Drvkde {
    let m = gcounts.len();
    let tau = 4.0 + drv as f64;
    let n: f64 = gcounts.iter().sum();
    let gpoints = binning::grid(a, b, m);

    // Kernel support, in grid steps: R clamps into [1, M].
    let l = ((tau * h * (m as f64 - 1.0) / (b - a)).floor() as i64).clamp(1, m as i64) as usize;
    let fac = (b - a) / (h * (m as f64 - 1.0));
    let arg: Vec<f64> = (0..=l).map(|i| i as f64 * fac).collect();

    // Hermite polynomial of order `drv`, by the same recurrence R uses.
    let hm: Vec<f64> = match drv {
        0 => vec![1.0; arg.len()],
        1 => arg.clone(),
        _ => {
            let mut old0 = vec![1.0; arg.len()];
            let mut old1 = arg.clone();
            let mut new = old1.clone();
            for ihm in 2..=drv {
                for k in 0..arg.len() {
                    new[k] = arg[k] * old1[k] - (ihm as f64 - 1.0) * old0[k];
                }
                old0 = std::mem::replace(&mut old1, new.clone());
            }
            new
        }
    };

    let sign = if drv % 2 == 0 { 1.0 } else { -1.0 };
    let kappam: Vec<f64> = arg
        .iter()
        .zip(&hm)
        .map(|(a, h_)| h_ * (dnorm(*a) / h.powi(drv as i32 + 1)) * sign / n)
        .collect();

    let est = symconv(&kappam, gcounts, sign);
    Drvkde { grid: gpoints, est }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn convolution_of_a_delta_returns_the_kernel() {
        let mut ss = vec![0.0; 16];
        ss[0] = 1.0;
        let rr = [1.0, 0.5, 0.25];
        let got = symconv(&rr, &ss, 1.0);
        for (i, want) in rr.iter().enumerate() {
            assert!((got[i] - want).abs() < 1e-12, "tap {i}");
        }
    }

    #[test]
    fn a_density_integrates_to_about_one() {
        // 4000 points from a standard normal, binned, then estimated: the grid sum times the
        // step should be close to 1.
        let x: Vec<f64> = (0..4000)
            .map(|i| {
                // deterministic pseudo-normal: Box-Muller on a low-discrepancy pair
                let u1 = ((i as f64 + 0.5) / 4000.0).clamp(1e-9, 1.0 - 1e-9);
                let u2 = (((i * 7919) % 4000) as f64 + 0.5) / 4000.0;
                (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
            })
            .collect();
        let (a, b) = (-6.0, 6.0);
        let g = binning::linbin(&x, a, b, 401);
        let d = drvkde(&g, 0, 0.3, a, b);
        let step = (b - a) / 400.0;
        let mass: f64 = d.est.iter().sum::<f64>() * step;
        assert!((mass - 1.0).abs() < 0.05, "mass {mass}");
    }
}
