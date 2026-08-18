## Six approaches fitted to one simulated dataset, for the vignette figure.
##
## Scenario 8 (delta = 4, top_factor = 2.0, R = 2.3) is the central case in the
## reaching regime: the design runs well past the zero crossing, so every
## approach has something to act on, and the true ErC10/ErC50 are known.
##
## The point of the figure is the fitted CURVE, not the endpoint: it shows
## directly how flooring and pinning bend the fit in the region the endpoints
## are read from.

source("R/setup.R"); source("R/data_prep.R"); source("R/diagnostics.R")
source("R/simulate.R"); load_bayesnec(); source("R/arms.R"); source("R/metrics.R")
use_compile_cache()
options(mc.cores = 4L)

cal <- calibrate_sigma(read_sgr("c_proliferum"))
SIGMA_0 <- sigma_0_at(cal$cv_control, R_ref = 2.3, t = 7)
truth <- sim_truth(R = 2.3, delta = 4, t = 7)
design <- sim_design(truth, top_factor = 2.0)
sim <- simulate_dataset(truth, design, cv_control = cal$cv_control,
                        sigma_ratio = cal$sigma_ratio, seed = 20260818,
                        sigma_mode = "absolute", sigma_0_abs = SIGMA_0)
saveRDS(list(truth = truth, sim = sim), "analysis/vignette_example_data.rds")

ap <- prepare_sgr(sim, "raw", meta = sim_meta())
prior <- arm_prior(ap$x, ap$y)
cat("true top", signif(truth$top,4), "bot", signif(truth$bot,4),
    "ErC10", signif(true_ecx(truth,10),4), "ErC50", signif(true_ecx(truth,50),4), "\n")

## Arm D needs the crossing from arm A, so A is fitted first.
fits <- list()
fits$A <- fit_arm("A", sim, prior = prior, meta = sim_meta())
crossing <- zero_crossing(fits$A$fit)
cat("estimated zero crossing:", signif(crossing, 4), "\n")
for (a in c("B1", "B2", "B3", "C", "D")) {
  cat("fitting", a, "\n")
  fits[[a]] <- fit_arm(a, sim, prior = prior, crossing = crossing,
                       meta = sim_meta())
}

## Fitted curves on a common grid, plus endpoints, for the vignette.
xx <- seq(0, max(sim$x), length.out = 200)
curves <- do.call(rbind, lapply(names(fits), function(a) {
  p <- fitted(fits[[a]]$fit, newdata = data.frame(x = xx), re_formula = NA)
  data.frame(arm = a, x = xx, est = p[, "Estimate"],
             lo = p[, "Q2.5"], hi = p[, "Q97.5"])
}))
ends <- do.call(rbind, lapply(names(fits), function(a)
  cbind(arm = a, endpoint_table(fits[[a]]$fit, a, "sim"))))
pars <- do.call(rbind, lapply(names(fits), function(a)
  cbind(arm = a, parameter_table(fits[[a]]$fit, a, "sim"))))

write.csv(curves, "analysis/vignette_example_curves.csv", row.names = FALSE)
write.csv(ends, "analysis/vignette_example_endpoints.csv", row.names = FALSE)
write.csv(pars, "analysis/vignette_example_parameters.csv", row.names = FALSE)

options(width = 200)
cat("\n=== endpoints, this dataset ===\n"); print(ends, digits = 4, row.names = FALSE)
cat("\n=== bot ===\n")
print(pars[grepl("^bot", pars$parameter), c("arm","mean","sd","q2.5","q97.5")],
      digits = 4, row.names = FALSE)
cat("\n=== sigma ===\n")
print(pars[pars$parameter == "sigma", c("arm","mean","sd")], digits = 4, row.names = FALSE)
