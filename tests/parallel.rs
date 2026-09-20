//! Channel-level parallelism: it must change the time and nothing else.
use flowvs::estimate::{self, Options};

fn synth(seed: u64, n: usize, mu_pos: f64, sd_pos: f64) -> Vec<f64> {
    // Deterministic pseudo-normals (Box-Muller on a linear congruential stream), so the test
    // needs no RNG crate and gives the same data on every machine.
    let mut s = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
    let mut next = || {
        s = s
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((s >> 11) as f64 + 0.5) / (1u64 << 53) as f64
    };
    (0..n)
        .map(|i| {
            let (u1, u2) = (next().max(1e-12), next());
            let z = (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos();
            if i % 10 < 7 {
                z * 30.0
            } else {
                mu_pos + z * sd_pos
            }
        })
        .collect()
}

fn channels(n_channels: usize, n_samples: usize, n_cells: usize) -> Vec<Vec<Vec<f64>>> {
    (0..n_channels)
        .map(|c| {
            (0..n_samples)
                .map(|s| {
                    synth(
                        (c * 100 + s) as u64 + 1,
                        n_cells,
                        1200.0 + 120.0 * c as f64,
                        400.0 + 40.0 * s as f64,
                    )
                })
                .collect()
        })
        .collect()
}

#[test]
fn parallel_and_sequential_agree_exactly() {
    let data = channels(6, 3, 1500);
    let seq = estimate::estimate_cofactors(&data, Options::default());
    let par = estimate::estimate_cofactors(
        &data,
        Options {
            threads: 4,
            ..Options::default()
        },
    );
    assert_eq!(seq.len(), par.len());
    for (i, (a, b)) in seq.iter().zip(&par).enumerate() {
        assert_eq!(
            a.cofactor.to_bits(),
            b.cofactor.to_bits(),
            "channel {i}: {} sequential vs {} parallel",
            a.cofactor,
            b.cofactor
        );
    }
    println!(
        "6 channels, identical bit for bit: {:?}",
        seq.iter().map(|e| e.cofactor).collect::<Vec<_>>()
    );
}

/// Run with `--ignored --nocapture` to see the speedup; it is a measurement, not an assertion,
/// because a loaded machine would make it flaky.
#[test]
#[ignore]
fn speedup_across_channels() {
    let data = channels(32, 4, 3000);
    for threads in [1usize, 4, 8, 16] {
        let t = std::time::Instant::now();
        let r = estimate::estimate_cofactors(
            &data,
            Options {
                threads,
                ..Options::default()
            },
        );
        println!(
            "{:>2} threads: {:>6.2} s for {} channels (first cofactor {:.2})",
            threads,
            t.elapsed().as_secs_f64(),
            r.len(),
            r[0].cofactor
        );
    }
}

/// End-to-end against the R reference's own real-data answers, without that data entering this
/// repository. Point `FLOWVS_INPUT_CSV` at a CSV with `sample_id` and channel columns:
///
/// ```text
/// FLOWVS_INPUT_CSV=~/tercen/flowvs/flowvs_input.csv \
/// FLOWVS_EXPECT_CSV=~/tercen/flowvs/flowvs_output.csv \
///   cargo test --release real_data -- --nocapture
/// ```
#[test]
fn real_data() {
    let Ok(input) = std::env::var("FLOWVS_INPUT_CSV") else {
        return;
    };
    let text = std::fs::read_to_string(&input).expect("read FLOWVS_INPUT_CSV");
    let mut lines = text.lines();
    let header: Vec<String> = lines
        .next()
        .unwrap()
        .split(',')
        .map(|h| h.trim().trim_matches('"').to_string())
        .collect();
    let i_sample = header.iter().position(|h| h == "sample_id").unwrap();
    let chan_idx: Vec<usize> = (0..header.len())
        .filter(|i| *i != i_sample && header[*i] != "cell_id")
        .collect();

    let mut sample_ids: Vec<String> = Vec::new();
    let mut per_channel: Vec<Vec<Vec<f64>>> = vec![Vec::new(); chan_idx.len()];
    for line in lines.filter(|l| !l.trim().is_empty()) {
        let f: Vec<&str> = line.split(',').collect();
        let sid = f[i_sample].trim().trim_matches('"').to_string();
        let s = match sample_ids.iter().position(|x| *x == sid) {
            Some(k) => k,
            None => {
                sample_ids.push(sid);
                for c in per_channel.iter_mut() {
                    c.push(Vec::new());
                }
                sample_ids.len() - 1
            }
        };
        for (c, i) in chan_idx.iter().enumerate() {
            if let Ok(v) = f[*i].trim().trim_matches('"').parse::<f64>() {
                per_channel[c][s].push(v);
            }
        }
    }
    println!(
        "{} channels x {} samples, {} values",
        per_channel.len(),
        sample_ids.len(),
        per_channel.iter().flatten().map(|v| v.len()).sum::<usize>()
    );

    let t = std::time::Instant::now();
    let got = estimate::estimate_cofactors(
        &per_channel,
        Options {
            threads: 0,
            ..Options::default()
        },
    );
    let secs = t.elapsed().as_secs_f64();

    for (c, i) in chan_idx.iter().enumerate() {
        println!("  {:<6} {:>14.4}", header[*i], got[c].cofactor);
    }
    println!("estimated in {secs:.2} s");

    if let Ok(expect) = std::env::var("FLOWVS_EXPECT_CSV") {
        let want = std::fs::read_to_string(&expect).expect("read FLOWVS_EXPECT_CSV");
        let mut worst: f64 = 0.0;
        for line in want.lines().skip(1).filter(|l| !l.trim().is_empty()) {
            let f: Vec<&str> = line.split(',').collect();
            let name = f[0].trim().trim_matches('"');
            let w: f64 = f[1].trim().parse().unwrap();
            let c = chan_idx
                .iter()
                .position(|i| header[*i] == name)
                .unwrap_or_else(|| panic!("channel {name} not in the input"));
            let ratio = (got[c].cofactor / w).ln().abs();
            println!(
                "  {name}: rust {:.4} vs R {:.4}  |ln ratio| {ratio:.3e}",
                got[c].cofactor, w
            );
            worst = worst.max(ratio);
        }
        // flowvs-rust-plan.md §3: within 10% for channels with real populations.
        assert!(
            worst < 1.10f64.ln(),
            "worst |ln ratio| {worst:e} exceeds the plan's 10% bound"
        );
    }
}
