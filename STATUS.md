# flowvs-rs — status, 2026-09-21

Built overnight against the goal in `~/tercen/goals/2026-09-19-asinh-flowvs.md`, and published to
`github.com/tercen/flowvs-rs` on 2026-09-21 so that `asinh_rust_operator` can depend on it as a
normal git dependency rather than a path. Six commits.

## Where it got to

All twelve stages of `estParamFlowVS` are ported and checked against the reference
(`tercen/asinh_operator/flowvs_standalone.R`, 937 lines of pure R), stage by stage rather than
only at the end.

| what | agreement with R |
|---|---|
| bandwidth, grid range | 1e-12 relative |
| linear binning | 2.7e-12 |
| density, second derivative | 1.3e-10 |
| significant grid points | identical: 113 of 401, and 69 of 401 |
| peak regions, peak positions and heights | identical count; positions inside 1% of a grid cell |
| R's own `density()` (used to place valleys) | 5e-14 |
| populations (n, mean, median, variance) | identical to 1e-9 |
| Bartlett objective at four cofactors | 5e-15 relative |
| **estimated cofactor** | 79.872038 vs 79.872038, `|ln ratio|` 2.2e-16 |

The plan (`flowvs-rust-plan.md` §3) asked for 1e-6 on the objective and 10% on the cofactor. This
is far inside that, on data with two clear populations — which is the easy case, and the one the
tolerance was *not* written for.

### Real data

The reference's own three-channel case (`~/tercen/flowvs/flowvs_input.csv`, 12 samples, 230,013
values, with R's published answers alongside it) reproduces exactly:

| channel | rust | R | `|ln ratio|` |
|---|---:|---:|---:|
| CD4 | 6450.3437 | 6450.3437 | 6.4e-15 |
| CD8 | 4766.6984 | 4766.6984 | 2.2e-16 |
| CD3 | 6317.3435 | 6317.3435 | 4.4e-16 |

The plan allows 10% on a real channel. That data is **not** in this repository and will not be:
the test reads `FLOWVS_INPUT_CSV` and does nothing when it is unset.

```bash
FLOWVS_INPUT_CSV=~/tercen/flowvs/flowvs_input.csv \
FLOWVS_EXPECT_CSV=~/tercen/flowvs/flowvs_output.csv \
  cargo test --release real_data -- --nocapture
```

### Speed

That whole three-channel estimate takes **2.1 s**. On the synthetic 32-channel benchmark:

| threads | 32 channels |
|---:|---:|
| 1 | 11.21 s |
| 4 | 3.48 s |
| 8 | 2.17 s |
| 16 | 1.74 s |

R needs 21.5 minutes for 32 channels on 66k cells, so this is the interactive review loop the
plan asks for. Parallelism is **across channels** and the thread count is an explicit option that
**defaults to 1**: peak memory grows with the number in flight, and an operator has to book a
fixed amount before it runs. Results are bit-identical whatever the thread count, which is a test.

## Two traps worth remembering

- **R 4.6's `bw.nrd0` divides the IQR by 1.34**, not the 1.349 the same rule uses elsewhere in
  flowVS and that older R used. It only bites when the IQR rule binds rather than the standard
  deviation, so the mixture fixture passed with the wrong constant and the single-population one
  did not. Exactly the tight channels where estimation is already fragile.
- **R's `optimize` replaces a NaN with the worst possible value.** The peak search hits this
  whenever a significant region reaches past the data, because the significance grid is widened
  by 3.7 bandwidths. Clamping instead invents a plateau and moves the peak.

Also: R's `round` is half-to-even, which the 10% sample and population trims depend on.

## Beyond flowVS, and off by default

Two guards against the degenerate cofactor a dim channel produces, both added after the parity
work and both inert unless asked for:

- **`Floor::SigmaNeg { factor }`** refuses a cofactor below `factor · σ_neg`, measured from the
  negative population's median absolute deviation and taken as the median across samples. 2.5 is
  the value the cofactor check used. `Options::floor` defaults to `Floor::None`.
- **The runner-up.** The search keeps the best optimum from any other interval, and a channel whose
  runner-up is nearly as good at a very different cofactor is reported `Fragile` rather than
  silently resolved.

Every estimate carries a `Status` — `Resolved`, `Fragile`, `Floored` or `Unstable` — so a caller
can show the three that need a human instead of a table of equally confident numbers.

## What is not here

- **The rest of the plan's §5**: the per-batch mode and machine-readable diagnostics beyond
  `Status`.
- **Committed real-data fixtures.** Everything in the repository is synthetic and deterministic
  by design. The real-data check above runs from a path outside it.
- **A CLI, a Sarno operator, a wasm build.** The plan's §4 architecture beyond the library.
- **Parity against flowVS itself** (the C-backed original), as opposed to Tercen's R port. The
  plan notes the two already differ by up to 12.7% on a borderline channel.

## Suggested next steps

1. Widen the real-data check: more channels, and the public `omip69_1k_donor` files, which could
   be committed.
2. Then §5's improvements, each behind an option that defaults to flowVS behaviour.
3. ~~Decide where this crate lives~~ — `tercen/flowvs-rs`, public, AGPL-3.0, as the plan's §10 Q1
   proposed.
