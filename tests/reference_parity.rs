//! Stage-by-stage parity with the R reference.
//!
//! `fixtures/gen_fixtures.R` runs `flowvs_standalone.R` on two deterministic synthetic samples —
//! a two-population mixture, which is what flowVS is designed for, and a single population,
//! where it has nothing to stabilise — and dumps every intermediate. Checking those, rather than
//! only the final cofactor, is what makes a disagreement locatable: a wrong bandwidth and a
//! wrong significance rule both move the answer, and look identical at the end.
use flowvs::{GRIDSIZE, stages};

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
    let i = header.iter().position(|h| h == name).unwrap();
    rows.iter().map(|r| r[i]).collect()
}

fn rel(a: f64, b: f64) -> f64 {
    (a - b).abs() / b.abs().max(1e-12)
}

fn check_case(tag: &str) {
    let dir = concat!(env!("CARGO_MANIFEST_DIR"), "/fixtures/");
    let (rh, rr) = read_csv(&format!("{dir}{tag}_raw.csv"));
    let raw = col(&rh, &rr, "value");
    let (sh, sr) = read_csv(&format!("{dir}{tag}_scalars.csv"));
    let scalar = |n: &str| col(&sh, &sr, n)[0];
    let (gh, gr) = read_csv(&format!("{dir}{tag}_stages.csv"));

    let st = stages(&raw, scalar("cofactor"), 0.05);

    assert_eq!(st.y.len(), scalar("n") as usize, "{tag}: sample size");
    assert!(
        rel(st.h, scalar("bandwidth")) < 1e-12,
        "{tag}: bandwidth {} vs R {}",
        st.h,
        scalar("bandwidth")
    );
    assert!(
        rel(st.range.0, scalar("range_lo")) < 1e-12 && rel(st.range.1, scalar("range_hi")) < 1e-12,
        "{tag}: grid range ({}, {}) vs R ({}, {})",
        st.range.0,
        st.range.1,
        scalar("range_lo"),
        scalar("range_hi")
    );

    let r_grid = col(&gh, &gr, "grid");
    let r_counts = col(&gh, &gr, "gcounts");
    let r_est0 = col(&gh, &gr, "est0");
    let r_est2 = col(&gh, &gr, "est2");
    let r_sig = col(&gh, &gr, "signif_curv");
    assert_eq!(r_grid.len(), GRIDSIZE);

    let worst = |a: &[f64], b: &[f64]| -> f64 {
        a.iter()
            .zip(b)
            .map(|(x, y)| (x - y).abs() / y.abs().max(1e-9))
            .fold(0.0f64, f64::max)
    };
    let w_grid = worst(&st.grid, &r_grid);
    let w_counts = worst(&st.gcounts, &r_counts);
    let w_est0 = worst(&st.density, &r_est0);
    let w_est2 = worst(&st.curvature, &r_est2);
    println!(
        "{tag}: grid {w_grid:.2e}  binning {w_counts:.2e}  density {w_est0:.2e}  curvature {w_est2:.2e}"
    );
    assert!(w_grid < 1e-12, "{tag}: grid points differ by {w_grid:e}");
    assert!(w_counts < 1e-10, "{tag}: linear binning differs by {w_counts:e}");
    assert!(w_est0 < 1e-8, "{tag}: density differs by {w_est0:e}");
    assert!(w_est2 < 1e-8, "{tag}: second derivative differs by {w_est2:e}");

    let mine: Vec<f64> = st.significant.iter().map(|b| f64::from(*b)).collect();
    let disagree: Vec<usize> = mine
        .iter()
        .zip(&r_sig)
        .enumerate()
        .filter(|(_, (a, b))| a != b)
        .map(|(i, _)| i)
        .collect();
    println!(
        "{tag}: significant points R={} rust={} disagreements={}",
        r_sig.iter().sum::<f64>(),
        mine.iter().sum::<f64>(),
        disagree.len()
    );
    assert!(
        disagree.is_empty(),
        "{tag}: significance differs at grid points {:?}",
        &disagree[..disagree.len().min(12)]
    );
}

#[test]
fn stages_match_the_r_reference_on_a_two_population_mixture() {
    check_case("mixture");
}

#[test]
fn stages_match_the_r_reference_on_a_single_population() {
    check_case("single");
}
