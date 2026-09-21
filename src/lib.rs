//! flowVS cofactor estimation, ported from Tercen's validated pure-R implementation
//! (`tercen/asinh_operator/flowvs_standalone.R`), which is itself a port of Azad et al. (2016)
//! and the flowStats peak machinery it depends on.
//!
//! Every constant here is normative: it comes from the reference and changing it changes the
//! answer. `STATUS.md` records the agreement stage by stage.
//!
//! Stage by stage, per channel and candidate cofactor:
//!
//! 1. transform `y = asinh(x / c)`  — [`transform`]
//! 2. bandwidth                      — [`bandwidth`]
//! 3. grid and linear binning        — [`binning`]
//! 4. binned KDE and 2nd derivative  — [`kde`]
//! 5. significant local modes        — [`signif`]
//! 6. peak regions and peaks         — [`peaks`]
//! 7. populations, Bartlett objective, cofactor search — [`estimate`]
//!
//! All stages are implemented and checked against fixtures dumped from the R reference
//! (`fixtures/gen_fixtures.R`).
pub mod bandwidth;
pub mod binning;
pub mod estimate;
pub mod kde;
pub mod optimize;
pub mod peaks;
pub mod rdensity;
pub mod signif;
pub mod stats;

/// `asinh(x / cofactor)` — the transform whose variance the cofactor stabilises.
pub fn transform(x: &[f64], cofactor: f64) -> Vec<f64> {
    x.iter().map(|v| (v / cofactor).asinh()).collect()
}

/// Everything stages 1–5 produce for one sample at one cofactor.
pub struct Stages {
    pub y: Vec<f64>,
    pub h: f64,
    pub range: (f64, f64),
    pub gcounts: Vec<f64>,
    pub density: Vec<f64>,
    pub curvature: Vec<f64>,
    pub significant: Vec<bool>,
    pub grid: Vec<f64>,
}

/// Run stages 1–5 with the reference defaults (`gridsize` 401, `bwFac` 2, `signifLevel` 0.05).
pub fn stages(x: &[f64], cofactor: f64, signif_level: f64) -> Stages {
    let y = transform(x, cofactor);
    let h = bandwidth::compute(&y, bandwidth::BW_FAC);
    let (a, b) = bandwidth::grid_range(&y, h);
    let gcounts = binning::linbin(&y, a, b, GRIDSIZE);
    let d0 = kde::drvkde(&gcounts, 0, h, a, b);
    let d2 = kde::drvkde(&gcounts, 2, h, a, b);
    let sig = signif::signif_curvature(y.len(), &gcounts, &d0, h, a, b, signif_level);
    Stages {
        y,
        h,
        range: (a, b),
        gcounts,
        density: d0.est,
        curvature: d2.est,
        significant: sig.curv,
        grid: d0.grid,
    }
}

/// Grid size used throughout flowVS.
pub const GRIDSIZE: usize = 401;
