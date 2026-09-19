//! End-to-end parity for stages 7–12: populations, the Bartlett objective, and the cofactor the
//! search lands on. Fixtures come from `fixtures/gen_estimate_fixtures.R`, three synthetic
//! two-population samples with different positive means.
use flowvs::estimate;

fn read_csv(path: &str) -> (Vec<String>, Vec<Vec<f64>>) {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("read {path}: {e}"));
    let mut lines = text.lines();
    let header: Vec<String> = lines
        .next()
        .unwrap()
        .split(',')
        .map(|h| h.trim().trim_matches('"').to_string())
        .collect();
    let rows = lines
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            l.split(',')
                .map(|v| v.trim().trim_matches('"').parse().unwrap_or(f64::NAN))
                .collect()
        })
        .collect();
    (header, rows)
}

fn col(header: &[String], rows: &[Vec<f64>], name: &str) -> Vec<f64> {
    let i = header
        .iter()
        .position(|h| h == name)
        .unwrap_or_else(|| panic!("no column {name} in {header:?}"));
    rows.iter().map(|r| r[i]).collect()
}

fn samples() -> Vec<Vec<f64>> {
    let dir = concat!(env!("CARGO_MANIFEST_DIR"), "/fixtures/");
    (1..=3)
        .map(|i| {
            let (h, r) = read_csv(&format!("{dir}sample{i}_raw.csv"));
            col(&h, &r, "value")
        })
        .collect()
}

#[test]
fn populations_match_the_r_reference() {
    let dir = concat!(env!("CARGO_MANIFEST_DIR"), "/fixtures/");
    let s = samples();
    let y = flowvs::transform(&s[0], 150.0);
    let mine = estimate::density_peaks(&y, 0.05, 1.0);

    let (h, r) = read_csv(&format!("{dir}populations_sample1_cf150.csv"));
    let (rn, rmean, rmed, rvar) = (
        col(&h, &r, "n"),
        col(&h, &r, "mean"),
        col(&h, &r, "median"),
        col(&h, &r, "variance"),
    );
    assert_eq!(
        mine.len(),
        rn.len(),
        "{} populations, R found {}",
        mine.len(),
        rn.len()
    );
    for (i, p) in mine.iter().enumerate() {
        assert_eq!(p.n as f64, rn[i], "population {i}: size");
        for (got, want, what) in [
            (p.mean, rmean[i], "mean"),
            (p.median, rmed[i], "median"),
            (p.variance, rvar[i], "variance"),
        ] {
            let rel = (got - want).abs() / want.abs().max(1e-12);
            assert!(rel < 1e-9, "population {i} {what}: {got} vs R {want}");
        }
    }
    println!("populations: {} matching R exactly", mine.len());
}

#[test]
fn the_objective_matches_the_r_reference() {
    let dir = concat!(env!("CARGO_MANIFEST_DIR"), "/fixtures/");
    let s = samples();
    let (h, r) = read_csv(&format!("{dir}objective.csv"));
    let (cfs, objs) = (col(&h, &r, "cofactor"), col(&h, &r, "objective"));
    let mut worst = 0.0f64;
    for (cf, want) in cfs.iter().zip(&objs) {
        let got = estimate::objective(*cf, &s, 0.05, 1.0);
        let rel = (got - want).abs() / want.abs().max(1e-12);
        println!("cofactor {cf}: rust {got:.6} R {want:.6} ({rel:.2e})");
        worst = worst.max(rel);
    }
    // The plan asks for 1e-6 relative on the objective at a fixed cofactor (§3).
    assert!(worst < 1e-6, "worst relative difference {worst:e}");
}

#[test]
fn the_estimated_cofactor_matches_the_r_reference() {
    let dir = concat!(env!("CARGO_MANIFEST_DIR"), "/fixtures/");
    let (h, r) = read_csv(&format!("{dir}optim_cofactor.csv"));
    let want = col(&h, &r, "cofactor")[0];
    let got = estimate::optim_cofactor(&samples(), 0.05, 1.0);
    let ratio = (got / want).ln().abs();
    println!("cofactor: rust {got:.6} R {want:.6}  |ln ratio| = {ratio:.3e}");
    // §3: within 10% for a channel with two clear populations. This fixture has them, so hold it
    // to something much tighter and let the plan's bound be the fallback.
    assert!(
        ratio < 0.01f64.ln_1p(),
        "rust {got} vs R {want} (|ln ratio| = {ratio})"
    );
}

/// How long the search takes, which is the reason for the port: the R implementation needs
/// 21.5 minutes for 32 channels on 66k cells, and a review loop has to be interactive.
#[test]
fn search_timing() {
    let s = samples();
    let t = std::time::Instant::now();
    let cf = estimate::optim_cofactor(&s, 0.05, 1.0);
    let secs = t.elapsed().as_secs_f64();
    println!(
        "optim_cofactor over {} samples x {} cells: {:.2} s (cofactor {cf:.3})",
        s.len(),
        s[0].len(),
        secs
    );
}
