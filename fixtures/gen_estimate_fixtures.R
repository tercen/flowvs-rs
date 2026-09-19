#!/usr/bin/env Rscript
# Stages 7-12 from the reference: populations, the Bartlett objective at fixed cofactors, and the
# cofactor the search lands on. Synthetic data only.
suppressWarnings(source(file.path(Sys.getenv("HOME"), "tercen/asinh_operator/flowvs_standalone.R")))
out <- "fixtures/"

set.seed(20260920)
mk <- function(n, mu_pos, sd_pos, frac_neg = 0.7) {
  neg <- rnorm(round(n * frac_neg), mean = 0, sd = 30)
  pos <- rnorm(n - length(neg), mean = mu_pos, sd = sd_pos)
  c(neg, pos)
}
samples <- list(mk(2000, 1800, 600), mk(2000, 2100, 700), mk(2000, 1500, 500))
for (i in seq_along(samples)) {
  write.csv(data.frame(value = samples[[i]]),
            sprintf("%ssample%d_raw.csv", out, i), row.names = FALSE)
}

# populations for sample 1 at cofactor 150
ps <- density_peaks(asinh(samples[[1]] / 150), bwFac = 2, signifLevel = 0.05, bwCorr = 1.0)
write.csv(as.data.frame(ps), paste0(out, "populations_sample1_cf150.csv"), row.names = FALSE)
cat("populations:", nrow(ps), "\n")

# objective at a few cofactors
cfs <- c(50, 150, 500, 1500)
obj <- sapply(cfs, function(c) flowvs_objective(c, samples, signifLevel = 0.05, bwCorr = 1.0))
write.csv(data.frame(cofactor = cfs, objective = obj),
          paste0(out, "objective.csv"), row.names = FALSE)
print(data.frame(cofactor = cfs, objective = obj))

t0 <- Sys.time()
cf <- optim_cofactor(samples, verbose = FALSE, signifLevel = 0.05, bwCorr = 1.0)
cat(sprintf("optim_cofactor = %.10g  (%.1f s)\n", cf, as.numeric(Sys.time() - t0, units = "secs")))
write.csv(data.frame(cofactor = cf), paste0(out, "optim_cofactor.csv"), row.names = FALSE)
