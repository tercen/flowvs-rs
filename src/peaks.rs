//! Peak regions and peaks (R `curv1_filter` boundaries, then `curv_peaks`).
//!
//! The significance test marks grid points; a **region** is a run of them, and its boundaries sit
//! half a grid step outside the first and last point of the run. The peak itself is then the
//! maximum of R's ordinary `density()` inside that region — a second, smoother estimate than the
//! one used to decide significance, which is what flowVS does and is worth not "simplifying".
use crate::optimize;
use crate::rdensity;

/// A significant region, in data coordinates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Region {
    pub left: f64,
    pub right: f64,
}

/// A peak: where it is and how high the density is there.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Peak {
    pub x: f64,
    pub y: f64,
}

/// Runs of significant grid points → regions (R `curv1_filter`'s boundary arithmetic).
pub fn regions(significant: &[bool], grid: &[f64]) -> Vec<Region> {
    let m = significant.len();
    let mut out = Vec::new();
    let mut i = 0;
    while i < m {
        if !significant[i] {
            i += 1;
            continue;
        }
        let start = i;
        while i + 1 < m && significant[i + 1] {
            i += 1;
        }
        let end = i; // inclusive
        // R: lowLims = (xGrid[k] + xGrid[k-1]) / 2 with 1-based k, and the first grid point when
        // the run starts at the very beginning; the upper limit is the midpoint past the run.
        let left = if start == 0 {
            grid[0]
        } else {
            (grid[start] + grid[start - 1]) / 2.0
        };
        let right = if end + 1 >= m {
            grid[m - 1]
        } else {
            (grid[end + 1] + grid[end]) / 2.0
        };
        out.push(Region { left, right });
        i += 1;
    }
    out
}

/// Peaks inside each region (R `curv_peaks`): maximise `density(y, n = 201, from, to)` there.
///
/// Regions that lie in the outer `border_quant` of the data range are dropped, as R does, since
/// a "peak" at the very edge is an artefact of the tail rather than a population.
pub fn curv_peaks(regions: &[Region], y: &[f64], border_quant: f64) -> Vec<Peak> {
    if regions.is_empty() || y.len() < 2 {
        return Vec::new();
    }
    let from = y.iter().cloned().fold(f64::INFINITY, f64::min);
    let to = y.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    // R: quantile(c(from, to), p) on two points is linear interpolation between them.
    let lo_cut = from + border_quant * (to - from);
    let hi_cut = from + (1.0 - border_quant) * (to - from);

    let dens = rdensity::density(y, 201, from, to);
    let f = |t: f64| -> f64 {
        let t = t.clamp(from, to); // approxfun is NA outside; a region never should be
        rdensity::approx(&dens.x, &dens.y, &[t])[0]
    };

    let mut peaks = Vec::new();
    for r in regions {
        if r.right <= lo_cut || r.left >= hi_cut {
            continue;
        }
        let (x, fx) = optimize::brent_fmax(r.left, r.right, f, optimize::default_tol());
        peaks.push(Peak { x, y: fx });
    }
    peaks
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_run_becomes_one_region_with_half_step_edges() {
        let grid: Vec<f64> = (0..11).map(|i| i as f64).collect();
        let mut sig = vec![false; 11];
        sig[4] = true;
        sig[5] = true;
        let r = regions(&sig, &grid);
        assert_eq!(r.len(), 1);
        // R: low 3.5, upp 5.5 for a run at 1-based grid positions 5..6
        assert!((r[0].left - 3.5).abs() < 1e-12 && (r[0].right - 5.5).abs() < 1e-12);
    }

    #[test]
    fn two_runs_become_two_regions() {
        let grid: Vec<f64> = (0..11).map(|i| i as f64).collect();
        let mut sig = vec![false; 11];
        sig[1] = true;
        sig[8] = true;
        assert_eq!(regions(&sig, &grid).len(), 2);
    }

    #[test]
    fn a_run_touching_the_edge_uses_the_grid_end() {
        let grid: Vec<f64> = (0..5).map(|i| i as f64).collect();
        let r = regions(&[true, true, false, true, true], &grid);
        assert_eq!(r.len(), 2);
        assert_eq!(r[0].left, 0.0);
        assert_eq!(r[1].right, 4.0);
    }
}
