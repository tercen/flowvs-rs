//! R's `stats::density()` with a Gaussian kernel, ported faithfully.
//!
//! flowVS locates a peak by maximising **this** density inside each significant region, not the
//! binned estimate used for the significance test, so the two have to coexist. R's version has
//! particulars that matter: `bw.nrd0`, a grid rounded up to a power of two and at least 512, an
//! extension of four bandwidths beyond the requested range, the kernel wrapped into the second
//! half of the array, and a final linear interpolation back onto the user's grid.
use std::sync::Arc;

use rustfft::{FftPlanner, num_complex::Complex};

use crate::stats;

/// R `bw.nrd0`: `0.9 · min(sd, IQR/1.34) · n^(-1/5)`, with R's degenerate-case fallbacks.
///
/// The divisor is **1.34**, not the 1.349 that the same rule uses elsewhere in flowVS and that
/// older R releases used here. It only matters when the IQR rule binds rather than the standard
/// deviation, which is why it can hide: on a two-population mixture `sd` is usually the smaller
/// of the two, and the wrong constant costs nothing; on a single tight population it shifts the
/// bandwidth by about 0.7% and moves every peak with it.
pub fn bw_nrd0(x: &[f64]) -> f64 {
    let n = x.len();
    assert!(n >= 2, "bw.nrd0 needs at least two observations");
    let mut sorted = x.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let iqr = stats::quantile_sorted(&sorted, 0.75) - stats::quantile_sorted(&sorted, 0.25);
    let sd = stats::sd(x);
    let mut lo = sd.min(iqr / 1.34);
    if lo == 0.0 {
        lo = if sd > 0.0 {
            sd
        } else if x[0].abs() > 0.0 {
            x[0].abs()
        } else {
            1.0
        };
    }
    0.9 * lo * (n as f64).powf(-0.2)
}

/// R's `BinDist`: linear binning of `x` onto `n` points over `[lo, up]`, returning `2n` values
/// with the second half zero (the padding the convolution needs). Weights are `1/len(x)`.
fn bin_dist(x: &[f64], lo: f64, up: f64, n: usize) -> Vec<f64> {
    let mut y = vec![0.0; 2 * n];
    let xdelta = (up - lo) / (n as f64 - 1.0);
    let w = 1.0 / x.len() as f64;
    for &xi in x {
        if !xi.is_finite() {
            continue;
        }
        let xpos = (xi - lo) / xdelta;
        if !xpos.is_finite() {
            continue;
        }
        let ix = xpos.floor() as i64;
        let fx = xpos - ix as f64;
        if 0 <= ix && ix < n as i64 - 1 {
            y[ix as usize] += (1.0 - fx) * w;
            y[ix as usize + 1] += fx * w;
        } else if ix == -1 {
            y[0] += fx * w;
        } else if ix == n as i64 - 1 {
            y[ix as usize] += (1.0 - fx) * w;
        }
    }
    y
}

pub struct Density {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub bw: f64,
}

/// R's `density(x)` with every default: `n = 512`, and the range extended by `cut = 3`
/// bandwidths past the data. flowVS uses this one to place valleys, and a different grid puts a
/// valley in a different place, so the defaults are part of the algorithm.
pub fn density_default(x: &[f64]) -> Density {
    let bw = bw_nrd0(x);
    let from = x.iter().cloned().fold(f64::INFINITY, f64::min) - 3.0 * bw;
    let to = x.iter().cloned().fold(f64::NEG_INFINITY, f64::max) + 3.0 * bw;
    density(x, 512, from, to)
}

/// `density(x, n = n_user, from, to)` with the Gaussian kernel and `bw.nrd0`.
pub fn density(x: &[f64], n_user: usize, from: f64, to: f64) -> Density {
    let bw = bw_nrd0(x);
    let mut n = n_user.max(512);
    if n > 512 {
        n = n.next_power_of_two();
    }
    let ext = 4.0;
    let lo = from - ext * bw;
    let up = to + ext * bw;

    let y = bin_dist(x, lo, up, n); // length 2n
    // R: kords <- seq.int(0, (2n-1)/(n-1) * (up - lo), length.out = 2n); then the second half is
    // mirrored to negative lags.
    let span = (2.0 * n as f64 - 1.0) / (n as f64 - 1.0) * (up - lo);
    let mut kords: Vec<f64> = (0..2 * n)
        .map(|i| span * i as f64 / (2.0 * n as f64 - 1.0))
        .collect();
    for i in (n + 1)..(2 * n) {
        // R indices: kords[(n+2):(2n)] <- -kords[n:2]
        kords[i] = -kords[2 * n - i];
    }
    let kords: Vec<f64> = kords
        .iter()
        .map(|k| (-0.5 * (k / bw) * (k / bw)).exp() / (bw * (2.0 * std::f64::consts::PI).sqrt()))
        .collect();

    // conv = Re(ifft(fft(y) * Conj(fft(kords)))) / length(y)
    let p = 2 * n;
    let mut fy: Vec<Complex<f64>> = y.iter().map(|v| Complex::new(*v, 0.0)).collect();
    let mut fk: Vec<Complex<f64>> = kords.iter().map(|v| Complex::new(*v, 0.0)).collect();
    let mut planner = FftPlanner::new();
    let fwd: Arc<dyn rustfft::Fft<f64>> = planner.plan_fft_forward(p);
    let inv: Arc<dyn rustfft::Fft<f64>> = planner.plan_fft_inverse(p);
    fwd.process(&mut fy);
    fwd.process(&mut fk);
    let mut prod: Vec<Complex<f64>> = fy.iter().zip(&fk).map(|(a, b)| a * b.conj()).collect();
    inv.process(&mut prod);
    let conv: Vec<f64> = prod[..n]
        .iter()
        .map(|c| (c.re / p as f64).max(0.0))
        .collect();

    let xords: Vec<f64> = crate::binning::grid(lo, up, n);
    let xs: Vec<f64> = crate::binning::grid(from, to, n_user);
    let ys = approx(&xords, &conv, &xs);
    Density { x: xs, y: ys, bw }
}

/// R's `approx`: linear interpolation, NaN outside the range.
pub fn approx(x: &[f64], y: &[f64], xout: &[f64]) -> Vec<f64> {
    xout.iter()
        .map(|&t| {
            if t < x[0] || t > x[x.len() - 1] {
                return f64::NAN;
            }
            match x.binary_search_by(|v| v.partial_cmp(&t).unwrap()) {
                Ok(i) => y[i],
                Err(0) => y[0],
                Err(i) if i >= x.len() => y[y.len() - 1],
                Err(i) => {
                    let (x0, x1) = (x[i - 1], x[i]);
                    let (y0, y1) = (y[i - 1], y[i]);
                    y0 + (y1 - y0) * (t - x0) / (x1 - x0)
                }
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bandwidth_matches_r() {
        // R: bw.nrd0(1:10) -> 1.7192864
        let x: Vec<f64> = (1..=10).map(f64::from).collect();
        assert!(
            (bw_nrd0(&x) - 1.719286404692283).abs() < 1e-12,
            "{}",
            bw_nrd0(&x)
        );
    }

    #[test]
    fn interpolation_is_linear_between_knots() {
        let got = approx(&[0.0, 1.0, 2.0], &[0.0, 10.0, 0.0], &[0.5, 1.5]);
        assert!((got[0] - 5.0).abs() < 1e-12 && (got[1] - 5.0).abs() < 1e-12);
    }
}
