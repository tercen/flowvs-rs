//! Bandwidth selection (R `compute_bandwidth`).
//!
//! `bwNS = min(sd, IQR/1.349) · (4/(7n))^(1/9)`, then `h = bwFac · bwNS` with `bwFac = 2`.
//! The `(4/(7n))^(1/9)` exponent is the normal-scale rule for a **second** derivative, which is
//! what the peak test estimates — not the familiar `n^(-1/5)` of a density.
use crate::stats;

pub const BW_FAC: f64 = 2.0;

/// `h` for a sample already transformed to asinh space.
pub fn compute(y: &[f64], bw_fac: f64) -> f64 {
    let n = y.len();
    if n == 0 {
        return 0.0;
    }
    let mut sorted = y.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let sigma = stats::sd(y).min(stats::iqr_scaled_sorted(&sorted));
    bw_fac * sigma * (4.0 / (7.0 * n as f64)).powf(1.0 / 9.0)
}

/// Grid range for the binned estimate (R `dflt_counts`): the data range, widened by `supp · h`.
pub const SUPP: f64 = 3.7;

pub fn grid_range(y: &[f64], h: f64) -> (f64, f64) {
    let min = y.iter().cloned().fold(f64::INFINITY, f64::min);
    let max = y.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    (min - SUPP * h, max + SUPP * h)
}
