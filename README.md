# flowvs-rs

`flowVS` channel-specific variance stabilisation in Rust: given the raw values of a channel across
several samples, estimate the `asinh` cofactor that makes the channel's populations share a
variance.

It is a port of Tercen's validated pure-R implementation of Azad et al. (2016), stage for stage,
and it agrees with that implementation to within floating-point noise — `|ln ratio|` of 2e-16 on
real three-channel data with published answers. What R needs 21.5 minutes for (32 channels,
66k cells) this does in 1.7 seconds, which is the difference between a batch job and a review loop.

```toml
[dependencies]
flowvs = { git = "https://github.com/tercen/flowvs-rs" }
```

```rust
use flowvs::estimate::{estimate_cofactors, Options, Status};

// channels[c][s] — the raw values of channel c in sample s
let est = estimate_cofactors(&channels, Options { threads: 8, ..Default::default() });

for e in &est {
    println!("{:.1}  ({})", e.cofactor, e.status.as_str());
}
```

Results do not depend on `threads`: channels are independent and each is deterministic. That is a
test, not a hope.

## What it estimates, and what it tells you when it cannot

flowVS splits each sample into populations at the valleys between density peaks, measures each
population's variance after `asinh(x / c)`, and searches for the `c` that minimises Bartlett's
statistic. On a channel with two clean populations that is a well-posed question. On a dim channel
it is not, and the original returns a number anyway.

So every estimate carries a `Status`:

| status | meaning |
|---|---|
| `Resolved` | one clear minimum |
| `Fragile` | another interval came close with a very different cofactor — a histogram should decide |
| `Floored` | flowVS went below the negative population's spread and the floor took over |
| `Unstable` | fewer than two usable populations; flowVS has nothing to say about this channel |

`Floored` needs `Options::floor`, which is **off by default** because the point of this crate is to
be flowVS. `Floor::SigmaNeg { factor: 2.5 }` refuses any cofactor below `2.5 · σ_neg`, where
`σ_neg` comes from the negative population's median absolute deviation. A cofactor smaller than the
noise it is meant to stabilise is a degenerate minimum: `asinh` is then a logarithm of the negative
population, whose variance can match another population's by coincidence.

The same caution applies to subsampling. The same channel can yield a cofactor of 2 or of 99
depending on how many cells per sample are drawn, with a healthy-looking Bartlett statistic at
both. Check the spread before freezing a cofactor table.

## Parity

Twelve stages, each checked against fixtures dumped from the R reference rather than only the final
number — bandwidth 1e-12, linear binning 2.7e-12, density and second derivative 1.3e-10, the
significant grid points identical, populations identical to 1e-9, the Bartlett objective 5e-15.
`fixtures/gen_fixtures.R` and `gen_estimate_fixtures.R` regenerate them; everything committed is
synthetic and seeded, and no study data enters this repository.

The real-data check reads its input from the environment and does nothing when it is unset:

```bash
FLOWVS_INPUT_CSV=… FLOWVS_EXPECT_CSV=… cargo test --release real_data -- --nocapture
```

Two traps the port paid for, both worth knowing if you compare against R yourself:

- **R 4.6's `bw.nrd0` divides the IQR by 1.34**, not the 1.349 the same rule uses elsewhere. It
  only bites when the IQR rule binds rather than the standard deviation — exactly the tight
  channels where estimation is already fragile.
- **R's `optimize` replaces a NaN with the worst possible value.** Clamping instead invents a
  plateau and moves the peak.

`STATUS.md` has the full table and what is deliberately missing.

## Licence

GPL-2.0-or-later (changed from AGPL-3.0 in 0.1.1). The R flowVS package it is ported from is
Artistic-2.0, which is GPL-compatible; GPL-2.0-or-later lets this crate link into both the
GPL-2.0-only FlowSOM operator family and GPL-3 code. `LICENSE` holds the GPL-2 text.
