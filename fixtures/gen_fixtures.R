#!/usr/bin/env Rscript
# Dump per-stage intermediates from the reference implementation
# (tercen/asinh_operator/flowvs_standalone.R) so the Rust port can be checked stage by stage
# rather than only on the final cofactor.
#
# Inputs are synthetic and deterministic: no patient data ever enters this repository.
#
#   Rscript fixtures/gen_fixtures.R
suppressWarnings(source(file.path(Sys.getenv("HOME"), "tercen/asinh_operator/flowvs_standalone.R")))

set.seed(20260920)
out <- file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), "")
if (out == "") out <- "fixtures/"

# --- case 1: a two-population mixture, the shape flowVS is designed for -----------------------
n <- 4000
neg <- rnorm(round(n * 0.7), mean = 0, sd = 30)
pos <- rnorm(n - length(neg), mean = 1800, sd = 600)
raw <- c(neg, pos)
# --- case 2: a single population, where flowVS has nothing to stabilise ------------------------
raw1 <- rnorm(2500, mean = 200, sd = 80)

dump_case <- function(tag, x, cofactor) {
  y <- asinh(x / cofactor)
  write.csv(data.frame(value = x), paste0(out, tag, "_raw.csv"), row.names = FALSE)

  h <- compute_bandwidth(y, bwFac = 2)
  dc <- dflt_counts(y, gridsize = 401, h = h, supp = 3.7)
  rng <- dc$range.x[[1]]
  gcounts <- dc$counts
  d0 <- drvkde_1d(gcounts, drv = 0, bandwidth = h, binned = TRUE,
                  range.x = dc$range.x, gridsize = 401)
  d2 <- drvkde_1d(gcounts, drv = 2, bandwidth = h, binned = TRUE,
                  range.x = dc$range.x, gridsize = 401)
  sfr <- signif_feature_region_1d(length(y), gcounts, 401, d0, h,
                                  signifLevel = 0.05, range.x = dc$range.x, curv = TRUE)

  write.csv(
    data.frame(
      grid = d0$x.grid[[1]],
      gcounts = gcounts,
      est0 = d0$est,
      est2 = d2$est,
      signif_curv = as.integer(sfr$curv)
    ),
    paste0(out, tag, "_stages.csv"), row.names = FALSE
  )
  scalars <- data.frame(
    cofactor = cofactor, n = length(y), bandwidth = h,
    range_lo = rng[1], range_hi = rng[2],
    y_min = min(y), y_max = max(y),
    sd = sd(y), iqr_scaled = IQR(y) / (qnorm(0.75) - qnorm(0.25)),
    q001 = unname(quantile(y, 0.001)), q999 = unname(quantile(y, 0.999)),
    q25 = unname(quantile(y, 0.25)), q75 = unname(quantile(y, 0.75))
  )
  write.csv(scalars, paste0(out, tag, "_scalars.csv"), row.names = FALSE)
  cat(sprintf("%s: n=%d h=%.10g range=[%.10g, %.10g] signif=%d\n",
              tag, length(y), h, rng[1], rng[2], sum(sfr$curv)))
}

dump_case("mixture", raw, 150)
dump_case("single", raw1, 150)
cat("done\n")
