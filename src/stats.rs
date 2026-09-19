//! The handful of R statistics the algorithm depends on, with R's conventions.
//!
//! These look trivial and are not: R's `quantile` defaults to **type 7**, its `sd` and `var`
//! divide by `n - 1`, and `IQR` is itself a type-7 quantile difference. Getting any of them
//! subtly wrong moves the bandwidth, which moves every peak.

/// `qnorm(0.75) - qnorm(0.25)`, the constant R's `IQR`-to-sigma scaling divides by.
pub const IQR_TO_SIGMA: f64 = 1.3489795003921634;

/// R `quantile(x, p)` with the default type 7, on a slice that is already sorted ascending.
pub fn quantile_sorted(sorted: &[f64], p: f64) -> f64 {
    assert!(!sorted.is_empty(), "quantile of an empty sample");
    let n = sorted.len();
    if n == 1 {
        return sorted[0];
    }
    let h = (n as f64 - 1.0) * p.clamp(0.0, 1.0);
    let lo = h.floor() as usize;
    let hi = h.ceil() as usize;
    if lo == hi {
        sorted[lo]
    } else {
        sorted[lo] + (h - lo as f64) * (sorted[hi] - sorted[lo])
    }
}

/// Sorts a copy and takes the quantile — convenience for one-off calls.
pub fn quantile(x: &[f64], p: f64) -> f64 {
    let mut v = x.to_vec();
    v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    quantile_sorted(&v, p)
}

/// R `var` / `sd`: denominator `n - 1`.
pub fn sd(x: &[f64]) -> f64 {
    let n = x.len();
    if n < 2 {
        return 0.0;
    }
    let mean = x.iter().sum::<f64>() / n as f64;
    let ss: f64 = x.iter().map(|v| (v - mean) * (v - mean)).sum();
    (ss / (n as f64 - 1.0)).sqrt()
}

/// R `IQR(x) / (qnorm(0.75) - qnorm(0.25))` — the robust sigma estimate flowVS uses.
pub fn iqr_scaled_sorted(sorted: &[f64]) -> f64 {
    (quantile_sorted(sorted, 0.75) - quantile_sorted(sorted, 0.25)) / IQR_TO_SIGMA
}

/// Survival function of chi-squared with one degree of freedom: `1 - pchisq(w, 1)`.
///
/// `pchisq(w, 1) = erf(sqrt(w/2))`, so the upper tail is `erfc(sqrt(w/2))`. Computing it this
/// way keeps the precision R has in the tail, where every p-value that matters lives.
pub fn chisq1_sf(w: f64) -> f64 {
    if !w.is_finite() || w <= 0.0 {
        return 1.0;
    }
    libm::erfc((w / 2.0).sqrt())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quantile_matches_r_type_7() {
        // R: quantile(c(1,2,3,4), c(0, 0.25, 0.5, 0.75, 1)) -> 1.00 1.75 2.50 3.25 4.00
        let x = [1.0, 2.0, 3.0, 4.0];
        for (p, want) in [
            (0.0, 1.0),
            (0.25, 1.75),
            (0.5, 2.5),
            (0.75, 3.25),
            (1.0, 4.0),
        ] {
            assert!((quantile(&x, p) - want).abs() < 1e-12, "p={p}");
        }
    }

    #[test]
    fn sd_uses_the_r_denominator() {
        // R: sd(c(1,2,3,4)) -> 1.290994
        assert!((sd(&[1.0, 2.0, 3.0, 4.0]) - 1.2909944487358056).abs() < 1e-12);
    }

    #[test]
    fn chisq_tail_matches_r() {
        // R: 1 - pchisq(c(0.5, 3.841459, 10), df = 1)
        for (w, want) in [
            (0.5, 0.479500122186953),
            (3.841458820694124, 0.05),
            (10.0, 0.001565402258392),
        ] {
            let got = chisq1_sf(w);
            assert!((got - want).abs() < 1e-9, "w={w}: {got} vs {want}");
        }
    }
}
