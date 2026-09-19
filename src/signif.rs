//! Which grid points are significant local modes (R `signif_feature_region_1d`, a port of
//! flowStats' `SignifFeatureRegion` for d = 1).
//!
//! Three conditions have to hold at a grid point for it to count as part of a peak:
//!
//! 1. the curvature is significantly non-zero, by a Wald test on the second derivative with a
//!    Holm-type sequential correction over the whole grid;
//! 2. the curvature is **negative**, so the point is a mode rather than a trough;
//! 3. the effective sample size there is at least 5, so the estimate rests on real data.
use crate::kde::{self, Drvkde};
use crate::stats;

pub struct Significance {
    /// One flag per grid point.
    pub curv: Vec<bool>,
    /// The standardised curvature statistic, kept for diagnostics.
    pub lambda: Vec<f64>,
}

/// `n` is the number of observations, `dest` the density (drv = 0) on the same grid.
pub fn signif_curvature(
    n: usize,
    gcounts: &[f64],
    dest: &Drvkde,
    h: f64,
    a: f64,
    b: f64,
    signif_level: f64,
) -> Significance {
    let m = gcounts.len();
    let nf = n as f64;

    // R clips the density at zero before using it as a variance scalar.
    let est0: Vec<f64> = dest.est.iter().map(|v| v.max(0.0)).collect();
    let ess: Vec<f64> = est0
        .iter()
        .map(|e| nf * e * h * (2.0 * std::f64::consts::PI).sqrt())
        .collect();
    let sig2: Vec<f64> = est0
        .iter()
        .map(|e| e / (8.0 * std::f64::consts::PI.sqrt() * nf * h))
        .collect();

    let d2 = kde::drvkde(gcounts, 2, h, a, b);
    let lambda: Vec<f64> = sig2
        .iter()
        .zip(&d2.est)
        .map(|(s2, f2)| {
            let inv = 1.0 / (s2 * 3.0 * h.powi(-4)).sqrt();
            if inv.is_finite() { inv * f2 } else { 0.0 }
        })
        .collect();

    let pvals: Vec<f64> = lambda.iter().map(|l| stats::chisq1_sf(l * l)).collect();

    // Holm-type sequential rejection, exactly as the R port writes it: sort the p-values, reject
    // while p_(i) <= level / (m + 1 - i), and mark every grid point sharing a rejected value.
    let mut ord: Vec<f64> = pvals.clone();
    ord.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
    let num_test = ord.len();
    let mut rejected: Vec<f64> = Vec::new();
    for (i, p) in ord.iter().enumerate() {
        let threshold = signif_level / (num_test as f64 + 1.0 - (i as f64 + 1.0));
        if *p <= threshold && *p > 0.0 {
            rejected.push(*p);
        }
    }

    let mut curv = vec![false; m];
    for (k, p) in pvals.iter().enumerate() {
        // A p-value of exactly zero is significant by construction (R sets it before the loop).
        if *p == 0.0 || rejected.iter().any(|r| r == p) {
            curv[k] = true;
        }
    }
    for k in 0..m {
        curv[k] = curv[k] && lambda[k] < 0.0 && ess[k] >= 5.0;
    }
    Significance { curv, lambda }
}
