# flowvs-rs — status, morning of 2026-09-20

Built overnight against the goal in `~/tercen/goals/2026-09-19-asinh-flowvs.md`. **Local git only:
no remote, nothing published.** Three commits.

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

Speed, which is the reason the port exists: 0.17 s against R's 1.1 s on the same three-sample
problem, release build. R needs 21.5 minutes for 32 channels on 66k cells.

## Two traps worth remembering

- **R 4.6's `bw.nrd0` divides the IQR by 1.34**, not the 1.349 the same rule uses elsewhere in
  flowVS and that older R used. It only bites when the IQR rule binds rather than the standard
  deviation, so the mixture fixture passed with the wrong constant and the single-population one
  did not. Exactly the tight channels where estimation is already fragile.
- **R's `optimize` replaces a NaN with the worst possible value.** The peak search hits this
  whenever a significant region reaches past the data, because the significance grid is widened
  by 3.7 bandwidths. Clamping instead invents a plateau and moves the peak.

Also: R's `round` is half-to-even, which the 10% sample and population trims depend on.

## What is not here

- **The improvements the plan proposes** (§5): the negative-spread estimator as a floor, the
  guardrails for dim and single-peak channels, the per-batch mode, machine-readable diagnostics.
  Faris chose plain flowVS parity first; these are the part that makes it better than flowVS
  rather than equal to it.
- **Real-data fixtures.** Everything here is synthetic and deterministic, by design — no patient
  data in this repository, ever. The public `omip69_1k_donor` files are the obvious next fixture.
- **A CLI, a Sarno operator, a wasm build.** The plan's §4 architecture beyond the library.
- **Parity against flowVS itself** (the C-backed original), as opposed to Tercen's R port. The
  plan notes the two already differ by up to 12.7% on a borderline channel.

## Suggested next steps

1. Add the public OMIP files as a real-data fixture and check a handful of channels end to end.
2. Then §5's improvements, each behind an option that defaults to flowVS behaviour.
3. Decide where this crate lives (`tercen/flowvs-rs` per the plan's §10 Q1) before it grows.
