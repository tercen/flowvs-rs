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

    // FLOWVS_FLOOR=2.5 turns on the negative-spread floor, to check it leaves a well-behaved
    // channel alone — the cofactor check found flowVS and σ_neg agree on bright markers.
    let floor = std::env::var("FLOWVS_FLOOR")
        .ok()
        .and_then(|v| v.parse::<f64>().ok())
        .map(|factor| estimate::Floor::SigmaNeg { factor })
        .unwrap_or(estimate::Floor::None);
    let t = std::time::Instant::now();
    let got = estimate::estimate_cofactors(
        &per_channel,
        Options {
            threads: 0,
            floor,
            ..Options::default()
        },
    );
    let secs = t.elapsed().as_secs_f64();

    for (c, i) in chan_idx.iter().enumerate() {
        println!(
            "  {:<6} {:>12.2}  flowVS {:>12.2}  sigma_neg {:>10}  {:<9} runner-up {}",
            header[*i],
            got[c].cofactor,
            got[c].flowvs_cofactor,
            got[c]
                .sigma_neg_cofactor
                .map(|v| format!("{v:.2}"))
                .unwrap_or_else(|| "-".into()),
            got[c].status.as_str(),
            got[c]
                .runner_up
                .map(|(cf, bt)| format!("{cf:.2} (bt {bt:.2})"))
                .unwrap_or_else(|| "none".into())
        );
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

/// The floor is the principled answer to a degenerate minimum: a cofactor below the spread of
/// the negative population cannot be stabilising it. Off by default, so flowVS parity holds.
#[test]
fn the_negative_spread_floor_rejects_a_degenerate_cofactor() {
    use flowvs::estimate::{Floor, Status, sigma_neg};

    // one channel, four samples, a tight negative population and a distant positive one — the
    // shape that makes the objective bimodal
    let data = channels(1, 4, 4000);
    let s0 = &data[0][0];
    let sn = sigma_neg(s0).expect("a negative population to measure");
    println!("sigma_neg = {sn:.2}, floor = {:.2}", 2.5 * sn);

    let plain = estimate::estimate_cofactors(&data, Options::default());
    let floored = estimate::estimate_cofactors(
        &data,
        Options {
            floor: Floor::SigmaNeg { factor: 2.5 },
            ..Options::default()
        },
    );
    println!(
        "flowVS {:.2} ({:?}) -> floored {:.2} ({:?})",
        plain[0].cofactor, plain[0].status, floored[0].cofactor, floored[0].status
    );
    // the floor never lowers a cofactor, and it reports what flowVS said either way
    assert!(floored[0].cofactor >= plain[0].cofactor - 1e-9);
    assert_eq!(floored[0].flowvs_cofactor, plain[0].flowvs_cofactor);
    if floored[0].status == Status::Floored {
        assert!(floored[0].cofactor > plain[0].cofactor);
        assert_eq!(
            floored[0].cofactor,
            floored[0].sigma_neg_cofactor.unwrap(),
            "a floored channel takes the floor exactly"
        );
    }
}

/// A channel with no negative population — a viability-gated CD45 — has no σ_neg, and the floor
/// must leave it alone rather than invent one.
#[test]
fn a_channel_without_negatives_has_no_floor() {
    use flowvs::estimate::sigma_neg;
    let all_positive: Vec<f64> = (0..5000).map(|i| 1000.0 + (i % 700) as f64).collect();
    assert!(sigma_neg(&all_positive).is_none());
}

/// σ_neg reads the negative half as a half-normal: median|x| = 0.6745 σ.
#[test]
fn sigma_neg_recovers_a_known_spread() {
    use flowvs::estimate::sigma_neg;
    // a clean symmetric sample with sd 200 — no positive population, which is what σ_neg
    // assumes it is looking at
    let mut st = 12345u64;
    let mut next = || {
        st = st
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((st >> 11) as f64 + 0.5) / (1u64 << 53) as f64
    };
    let x: Vec<f64> = (0..40000)
        .map(|_| {
            let (u1, u2) = (next().max(1e-12), next());
            200.0 * (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
        })
        .collect();
    let sn = sigma_neg(&x).expect("negatives");
    println!("sigma_neg {sn:.1} for a sample with sd 200");
    assert!((sn - 200.0).abs() / 200.0 < 0.15, "got {sn}");
}
