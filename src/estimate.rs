//! Stages 7–12: populations, the Bartlett objective, and the cofactor search.
//!
//! This is where flowVS's idea lives. Split each sample into populations at the valleys between
//! peaks, measure each population's variance, and choose the cofactor that makes those variances
//! as equal as possible — Bartlett's statistic being the measure of "as equal as possible".
//!
//! Every constant here is from the R reference and is normative
//! (`flowvs-rust-plan.md` §2): border trim 0.1%, population trim 1%, peaks merged closer than 5%
//! of the range, samples and populations trimmed at 10%, and the search over `[e^i, e^(i+1)]`
//! for `i = -1..9` followed by a local refinement in tenths.
use crate::peaks::{self, Peak};
use crate::rdensity;
use crate::stats;
use crate::{bandwidth, binning, kde, signif};

pub const BORDER_QUANT: f64 = 0.001;
pub const POPULATION_QUANT: f64 = 0.01;
pub const PEAK_DISTANCE_THR: f64 = 0.05;
pub const MAX_BT: f64 = 1e9;

/// One population's summary, R's `peaksStats` row.
#[derive(Debug, Clone, Copy)]
pub struct Population {
    pub n: usize,
    pub mean: f64,
    pub median: f64,
    pub variance: f64,
    /// Which sample it came from (R's fifth column).
    pub sample: usize,
}

/// R's `round`: half to even. `round(2.5)` is 2, not 3, and the 10% trims depend on it.
fn r_round(x: f64) -> f64 {
    let r = x.round();
    if (x - x.trunc()).abs() == 0.5 && r % 2.0 != 0.0 {
        r - x.signum()
    } else {
        r
    }
}

fn median_sorted(sorted: &[f64]) -> f64 {
    let n = sorted.len();
    if n == 0 {
        return f64::NAN;
    }
    if n % 2 == 1 {
        sorted[n / 2]
    } else {
        (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
    }
}

fn median(x: &[f64]) -> f64 {
    let mut v = x.to_vec();
    v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    median_sorted(&v)
}

fn variance(x: &[f64]) -> f64 {
    let s = stats::sd(x);
    s * s
}

/// Peaks and the valleys between them, for one transformed sample
/// (R `find_density_peaks_valleys` with `use_curv1filter = TRUE`).
pub fn peaks_and_valleys(y: &[f64], signif_level: f64, bw_corr: f64) -> (Vec<Peak>, Vec<f64>) {
    if y.len() < 10 {
        return (Vec::new(), Vec::new());
    }
    let h = bandwidth::compute(y, bandwidth::BW_FAC) * bw_corr;
    let (a, b) = bandwidth::grid_range(y, h);
    let gcounts = binning::linbin(y, a, b, crate::GRIDSIZE);
    let d0 = kde::drvkde(&gcounts, 0, h, a, b);
    let sig = signif::signif_curvature(y.len(), &gcounts, &d0, h, a, b, signif_level);
    let regions = peaks::regions(&sig.curv, &d0.grid);
    let mut pk = peaks::curv_peaks(&regions, y, 0.01);

    let dens = rdensity::density_default(y);
    if pk.is_empty() {
        // R falls back to the highest point of the default density.
        let mut best = 0usize;
        for i in 1..dens.y.len() {
            if dens.y[i] > dens.y[best] {
                best = i;
            }
        }
        pk.push(Peak {
            x: dens.x[best],
            y: dens.y[best],
        });
    }
    pk.sort_by(|p, q| p.x.partial_cmp(&q.x).unwrap_or(std::cmp::Ordering::Equal));
    let valleys = valleys_between(&pk, &dens);
    (pk, valleys)
}

/// The lowest point of the density strictly between each pair of consecutive peaks.
fn valleys_between(pk: &[Peak], dens: &rdensity::Density) -> Vec<f64> {
    let mut valleys = Vec::new();
    for w in pk.windows(2) {
        let mut best: Option<(f64, f64)> = None;
        for (x, y) in dens.x.iter().zip(&dens.y) {
            if *x > w[0].x && *x < w[1].x && best.is_none_or(|(_, by)| *y < by) {
                best = Some((*x, *y));
            }
        }
        if let Some((x, _)) = best {
            valleys.push(x);
        }
    }
    valleys
}

/// Merge peaks closer than 5% of the data range, keeping the taller (R's `peak.rm.idx` walk).
fn merge_close_peaks(pk: &[Peak], mrange: f64) -> Vec<Peak> {
    if pk.len() < 2 {
        return pk.to_vec();
    }
    let mut remove: Vec<usize> = Vec::new();
    let (mut p1, mut p2) = (0usize, 1usize);
    while p2 < pk.len() {
        if (pk[p1].x - pk[p2].x).abs() < mrange * PEAK_DISTANCE_THR {
            if pk[p1].y < pk[p2].y {
                remove.push(p1);
                p1 = p2;
            } else {
                remove.push(p2);
            }
            p2 += 1;
        } else {
            p1 = p2;
            p2 += 1;
        }
    }
    let mut kept = pk.to_vec();
    let mut i = 0;
    kept.retain(|_| {
        let keep = !remove.contains(&i);
        i += 1;
        keep
    });
    kept
}

/// Populations of one transformed sample (R `density_peaks`).
pub fn density_peaks(y: &[f64], signif_level: f64, bw_corr: f64) -> Vec<Population> {
    // R trims the extreme values twice: the exact min/max, then the 0.1% band.
    let ymin = y.iter().cloned().fold(f64::INFINITY, f64::min);
    let ymax = y.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let mut y: Vec<f64> = y
        .iter()
        .cloned()
        .filter(|v| *v > ymin && *v < ymax)
        .collect();
    if y.len() < 10 {
        return Vec::new();
    }
    let lo = stats::quantile(&y, BORDER_QUANT);
    let hi = stats::quantile(&y, 1.0 - BORDER_QUANT);
    y.retain(|v| *v > lo && *v < hi);
    if y.len() < 10 {
        return Vec::new();
    }

    let (pk, valleys) = peaks_and_valleys(&y, signif_level, bw_corr);
    if pk.is_empty() {
        return Vec::new();
    }
    let ymin = y.iter().cloned().fold(f64::INFINITY, f64::min);
    let ymax = y.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let mrange = ymax - ymin;

    let merged = merge_close_peaks(&pk, mrange);
    let dens = rdensity::density_default(&y);
    let valleys = if merged.len() == pk.len() {
        valleys
    } else {
        valleys_between(&merged, &dens)
    };
    if merged.is_empty() {
        return Vec::new();
    }

    let mut out = Vec::new();
    let npeaks = merged.len();
    let mut left = ymin;
    for (i, p) in merged.iter().enumerate() {
        let mut right = if i + 1 == npeaks {
            ymax
        } else if valleys.len() > i {
            valleys[i]
        } else {
            // R: the density minimum between this peak and the next, else their midpoint.
            let mut best: Option<(f64, f64)> = None;
            for (x, yv) in dens.x.iter().zip(&dens.y) {
                if *x > p.x && *x < merged[i + 1].x && best.is_none_or(|(_, by)| *yv < by) {
                    best = Some((*x, *yv));
                }
            }
            best.map(|(x, _)| x)
                .unwrap_or_else(|| (p.x + merged[i + 1].x) / 2.0)
        };

        let mut sel: Vec<f64> = y
            .iter()
            .cloned()
            .filter(|v| *v >= left && *v <= right)
            .collect();
        if !sel.is_empty() {
            let mut lo_p = stats::quantile(&sel, POPULATION_QUANT);
            let mut hi_p = stats::quantile(&sel, 1.0 - POPULATION_QUANT);
            // Asymmetry rule: an outermost peak with a long tail is capped at three times its
            // short half-width, so a skirt does not become the population's variance.
            let r = (hi_p - p.x) / (p.x - lo_p + 1e-10);
            if r >= 3.0 && i + 1 == npeaks {
                right = hi_p.min(p.x + 3.0 * (p.x - lo_p));
                sel = y
                    .iter()
                    .cloned()
                    .filter(|v| *v >= left && *v <= right)
                    .collect();
                if !sel.is_empty() {
                    lo_p = stats::quantile(&sel, POPULATION_QUANT);
                    hi_p = stats::quantile(&sel, 1.0 - POPULATION_QUANT);
                }
            } else if r <= 1.0 / 3.0 && i == 0 {
                left = lo_p.max(p.x - 3.0 * (hi_p - p.x));
                sel = y
                    .iter()
                    .cloned()
                    .filter(|v| *v >= left && *v <= right)
                    .collect();
                if !sel.is_empty() {
                    lo_p = stats::quantile(&sel, POPULATION_QUANT);
                    hi_p = stats::quantile(&sel, 1.0 - POPULATION_QUANT);
                }
            }
            let population: Vec<f64> = sel
                .into_iter()
                .filter(|v| *v >= lo_p && *v <= hi_p)
                .collect();
            if population.len() > 1 {
                out.push(Population {
                    n: population.len(),
                    mean: population.iter().sum::<f64>() / population.len() as f64,
                    median: median(&population),
                    variance: variance(&population),
                    sample: 0,
                });
            }
        }
        left = right;
    }
    out
}

/// Bartlett's statistic over populations (R `bartlett_test`), or `MAX_BT` when it cannot be formed.
pub fn bartlett(pops: &[Population]) -> f64 {
    let kept: Vec<&Population> = pops.iter().filter(|p| p.n > 1).collect();
    if kept.len() <= 1 {
        return MAX_BT;
    }
    if kept
        .iter()
        .any(|p| p.variance <= 0.0 || !p.variance.is_finite())
    {
        return MAX_BT;
    }
    let sum1: f64 = kept
        .iter()
        .map(|p| (p.n as f64 - 1.0) * p.variance.ln())
        .sum();
    let sum2: f64 = kept.iter().map(|p| (p.n as f64 - 1.0) * p.variance).sum();
    let n_total: f64 = kept.iter().map(|p| p.n as f64).sum();
    let sum3: f64 = kept.iter().map(|p| 1.0 / (p.n as f64 - 1.0)).sum();
    let k = kept.len() as f64;
    let bt = ((n_total - k) * (sum2 / (n_total - k)).ln() - sum1)
        / (1.0 + (sum3 - 1.0 / (n_total - k)) / (3.0 * (k - 1.0)));
    if bt.is_finite() { bt } else { MAX_BT }
}

/// The objective at one cofactor, over every sample of a channel (R `flowvs_objective`).
pub fn objective(cofactor: f64, samples: &[Vec<f64>], signif_level: f64, bw_corr: f64) -> f64 {
    let mut pops: Vec<Population> = Vec::new();
    for (i, s) in samples.iter().enumerate() {
        let y = crate::transform(s, cofactor);
        for mut p in density_peaks(&y, signif_level, bw_corr) {
            p.sample = i;
            pops.push(p);
        }
    }
    if pops.len() <= 1 {
        return MAX_BT;
    }

    // Drop the 10% of samples whose peak count is furthest from the median count.
    let mut ids: Vec<usize> = pops.iter().map(|p| p.sample).collect();
    ids.sort_unstable();
    ids.dedup();
    let counts: Vec<f64> = ids
        .iter()
        .map(|id| pops.iter().filter(|p| p.sample == *id).count() as f64)
        .collect();
    let med = median(&counts);
    let mut order: Vec<usize> = (0..ids.len()).collect();
    // R's `order(..., decreasing = TRUE)` is stable, so ties keep ascending index order.
    order.sort_by(|a, b| {
        (counts[*b] - med)
            .abs()
            .partial_cmp(&(counts[*a] - med).abs())
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let rm = r_round(ids.len() as f64 / 10.0) as usize;
    if rm > 0 && ids.len() > rm + 1 {
        let drop: Vec<usize> = order[..rm].iter().map(|i| ids[*i]).collect();
        pops.retain(|p| !drop.contains(&p.sample));
    }
    if pops.len() <= 1 {
        return MAX_BT;
    }

    // Drop the 10% of populations whose variance is furthest from the median variance.
    let vars: Vec<f64> = pops.iter().map(|p| p.variance).collect();
    let medv = median(&vars);
    let mut vorder: Vec<usize> = (0..pops.len()).collect();
    vorder.sort_by(|a, b| {
        (vars[*b] - medv)
            .abs()
            .partial_cmp(&(vars[*a] - medv).abs())
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let rm = r_round(pops.len() as f64 / 10.0) as usize;
    if rm > 0 && pops.len() > rm + 1 {
        let drop: Vec<usize> = vorder[..rm].to_vec();
        pops = pops
            .into_iter()
            .enumerate()
            .filter(|(i, _)| !drop.contains(i))
            .map(|(_, p)| p)
            .collect();
    }
    if pops.len() <= 1 {
        return MAX_BT;
    }
    bartlett(&pops)
}

/// The cofactor search (R `optim_cofactor`): Brent inside each `[e^i, e^(i+1)]`, then a local
/// refinement in tenths around the best one.
pub fn optim_cofactor(samples: &[Vec<f64>], signif_level: f64, bw_corr: f64) -> f64 {
    let (cf_low, cf_high) = (-1i32, 10i32);
    let mut cfopt = Vec::new();
    let mut btopt = Vec::new();
    for i in cf_low..cf_high {
        let low = (i as f64).exp();
        let high = ((i + 1) as f64).exp();
        let tol = (high - low) / 10.0;
        let (x, fx) = crate::optimize::brent_fmin(
            low,
            high,
            |c| objective(c, samples, signif_level, bw_corr),
            tol,
        );
        cfopt.push(x);
        btopt.push(fx);
    }
    let mut best = 0usize;
    for i in 1..btopt.len() {
        if btopt[i] < btopt[best] {
            best = i;
        }
    }
    let centre = cfopt[best];
    let del = centre / 10.0;
    let mut cf_local = [0.0f64; 11];
    let mut bt_local = [0.0f64; 11];
    for k in 0..5 {
        cf_local[k] = centre - (5 - k) as f64 * del;
        cf_local[6 + k] = centre + (k + 1) as f64 * del;
    }
    cf_local[5] = centre;
    bt_local[5] = btopt[best];
    for k in (0..11).filter(|k| *k != 5) {
        bt_local[k] = objective(cf_local[k], samples, signif_level, bw_corr);
    }
    let mut bestk = 0usize;
    for k in 1..11 {
        if bt_local[k] < bt_local[bestk] {
            bestk = k;
        }
    }
    cf_local[bestk]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn r_round_is_half_to_even() {
        assert_eq!(r_round(2.5), 2.0);
        assert_eq!(r_round(3.5), 4.0);
        assert_eq!(r_round(2.4), 2.0);
        assert_eq!(r_round(-2.5), -2.0);
    }

    #[test]
    fn bartlett_is_max_without_two_usable_populations() {
        let p = Population {
            n: 10,
            mean: 0.0,
            median: 0.0,
            variance: 1.0,
            sample: 0,
        };
        assert_eq!(bartlett(&[p]), MAX_BT);
        let bad = Population { variance: 0.0, ..p };
        assert_eq!(bartlett(&[p, bad]), MAX_BT);
    }

    #[test]
    fn bartlett_is_zero_for_identical_variances() {
        let p = Population {
            n: 100,
            mean: 0.0,
            median: 0.0,
            variance: 2.0,
            sample: 0,
        };
        assert!(bartlett(&[p, p, p]).abs() < 1e-9);
    }

    #[test]
    fn merging_keeps_the_taller_of_two_close_peaks() {
        let pk = [
            Peak { x: 0.0, y: 1.0 },
            Peak { x: 0.01, y: 5.0 },
            Peak { x: 2.0, y: 3.0 },
        ];
        let m = merge_close_peaks(&pk, 10.0);
        assert_eq!(m.len(), 2);
        assert_eq!(m[0].y, 5.0);
    }
}

/// How the estimator is run.
///
/// `threads` is explicit on purpose. Channels are independent, so the search parallelises almost
/// perfectly across them, but peak memory grows with the number in flight: each worker holds the
/// transformed copy of one sample plus its grids. An operator has to book a fixed amount of
/// memory before it runs (see the `create-rust-operator` skill §7), and a default of "every core
/// on the machine" makes that booking depend on the machine. So: set it.
#[derive(Debug, Clone, Copy)]
pub struct Options {
    pub signif_level: f64,
    pub bw_corr: f64,
    /// Channels estimated at once. `0` means "as many as the machine has", which is convenient
    /// for a CLI and wrong for an operator.
    pub threads: usize,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            signif_level: 0.05,
            bw_corr: 1.0,
            threads: 1,
        }
    }
}

/// What the search found for one channel.
#[derive(Debug, Clone, Copy)]
pub struct ChannelEstimate {
    pub cofactor: f64,
    /// Bartlett's statistic there — `MAX_BT` means no usable populations were found, which is
    /// the "flowVS has nothing to say about this channel" case.
    pub objective: f64,
}

/// Estimate one cofactor per channel (R `est_param_flowvs`), optionally in parallel.
///
/// `channels[c][s]` is the raw values of channel `c` in sample `s`. The result is in channel
/// order and does not depend on `threads`: each channel is independent and deterministic.
pub fn estimate_cofactors(channels: &[Vec<Vec<f64>>], opts: Options) -> Vec<ChannelEstimate> {
    let one = |samples: &Vec<Vec<f64>>| -> ChannelEstimate {
        let cofactor = optim_cofactor(samples, opts.signif_level, opts.bw_corr);
        let objective = objective(cofactor, samples, opts.signif_level, opts.bw_corr);
        ChannelEstimate {
            cofactor,
            objective,
        }
    };
    if opts.threads == 1 || channels.len() < 2 {
        return channels.iter().map(one).collect();
    }
    use rayon::prelude::*;
    let n = if opts.threads == 0 {
        std::thread::available_parallelism()
            .map(|v| v.get())
            .unwrap_or(1)
    } else {
        opts.threads
    };
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(n)
        .build()
        .expect("build the estimator thread pool");
    pool.install(|| channels.par_iter().map(one).collect())
}
