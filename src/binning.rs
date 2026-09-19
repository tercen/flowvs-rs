//! Linear binning onto the evaluation grid (R `linbin_1d`).
//!
//! Each value's weight is split between the two grid points it falls between, which is what
//! makes the binned KDE agree with the direct one to O(1/M²) rather than O(1/M).

/// `M` equally spaced grid points from `a` to `b`.
pub fn grid(a: f64, b: f64, m: usize) -> Vec<f64> {
    if m == 1 {
        return vec![a];
    }
    let step = (b - a) / (m as f64 - 1.0);
    (0..m).map(|i| a + step * i as f64).collect()
}

/// Weights of `x` distributed over the grid `[a, b]` with `m` points.
pub fn linbin(x: &[f64], a: f64, b: f64, m: usize) -> Vec<f64> {
    let mut counts = vec![0.0; m];
    if m < 2 {
        return counts;
    }
    let delta = (b - a) / (m as f64 - 1.0);
    for &xi in x {
        // R works in 1-based grid positions; keep that arithmetic and subtract at the index.
        let pos = (xi - a) / delta + 1.0;
        if !(pos >= 1.0 && pos <= m as f64) {
            continue; // outside the grid: R skips it
        }
        let lo = pos.floor();
        let hi = pos.ceil();
        if lo == hi {
            counts[lo as usize - 1] += 1.0;
        } else {
            let frac = pos - lo;
            counts[lo as usize - 1] += 1.0 - frac;
            counts[hi as usize - 1] += frac;
        }
    }
    counts
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_point_on_a_grid_node_lands_whole() {
        let c = linbin(&[0.5], 0.0, 1.0, 3);
        assert_eq!(c, vec![0.0, 1.0, 0.0]);
    }

    #[test]
    fn a_point_between_nodes_is_split_linearly() {
        let c = linbin(&[0.25], 0.0, 1.0, 3);
        assert!((c[0] - 0.5).abs() < 1e-12 && (c[1] - 0.5).abs() < 1e-12);
    }

    #[test]
    fn total_weight_is_preserved_and_outside_points_are_dropped() {
        let c = linbin(&[0.1, 0.9, 2.0, -1.0], 0.0, 1.0, 11);
        assert!((c.iter().sum::<f64>() - 2.0).abs() < 1e-12);
    }
}
