## Is `true_ecx()` the same quantity that `bayesnec::ecx()` estimates?
##
## The coverage validation showed arm A covering ErC50 at 0.817 even with the
## variance model correctly specified, driven by a bias of -0.225 against an
## interval half-width of about 0.4. Before that is reported as a property of
## the fitting, it has to be shown not to be a mismatch between the estimand and
## the estimator.
##
## `check_ecx_estimator()` does not settle this: it compares `true_ecx()` to a
## grid evaluation of the SAME analytic curve, so it validates the formula
## against itself. The open question is whether `bayesnec::ecx()` targets the
## same number. `ecx_x_absolute()` sets `range_y <- c(0, max(y))` from FITTED
## values, whereas `true_ecx()` solves against the `top` PARAMETER.
##
## Design: simulate at a negligible noise level, so the posterior sits on the
## true curve and any remaining discrepancy cannot be sampling error. Then
## compare three numbers:
##
##   1. true_ecx()                  -- the study's estimand
##   2. ecx() on the fit            -- what the study actually records
##   3. analytic ECx computed from the POSTERIOR DRAWS of top/bot/beta/nec
##
## (2) vs (3) isolates `ecx()`'s numerics from the fit: if (3) matches (1) and
## (2) does not, the estimator and the estimand are different quantities and the
## bias column is an artefact. If (3) also misses, the fit itself is biased.
##
## Run: Rscript analysis/phase5_estimand_check.R

source("R/setup.R")
source("R/data_prep.R")
source("R/diagnostics.R")
source("R/simulate.R")
load_bayesnec()
source("R/arms.R")
source("R/metrics.R")

use_compile_cache()
options(mc.cores = 4L)

MCMC_Q <- list(chains = 4L, iter = 2000L, warmup = 1000L, adapt_delta = 0.95,
               max_treedepth = 12L, seed = 20260812L)

truth <- sim_truth(R = 2.3, delta = 4, t = 7)
design <- sim_design(truth, top_factor = 1.0, n_conc = 12, n_rep = 5,
                     n_control = 10)

## Analytic ECx from a single parameter draw, same algebra as true_ecx() but
## applied to the posterior rather than to the truth.
ecx_from_draw <- function(top, bot, beta, nec, p) {
  target <- top * (1 - p / 100)
  nec + log((top - bot) / (target - bot)) / exp(beta)
}

report <- function(sigma_0, label) {
  cat("\n================ ", label, " (sigma_0 = ", signif(sigma_0, 3), ") ================\n",
      sep = "")
  sim <- simulate_dataset(truth, design, cv_control = NA, sigma_ratio = 1,
                          seed = 424242, sigma_mode = "absolute",
                          sigma_0_abs = sigma_0)
  ap <- prepare_sgr(sim, "raw", meta = sim_meta())
  pr <- arm_prior(ap$x, ap$y)
  f <- fit_arm("A", sim, prior = pr, mcmc = MCMC_Q, meta = sim_meta())

  dr <- brms::as_draws_df(f$fit$fit)
  gp <- function(n) dr[[paste0("b_", n, "_Intercept")]]
  top <- gp("top"); bot <- gp("bot"); beta <- gp("beta"); nec <- gp("nec")

  cat("\nparameter recovery (posterior mean vs truth):\n")
  print(data.frame(
    parameter = c("top", "bot", "beta", "nec"),
    truth = round(c(truth$top, truth$bot, truth$beta, truth$nec), 4),
    posterior_mean = round(c(mean(top), mean(bot), mean(beta), mean(nec)), 4),
    row.names = NULL))

  out <- do.call(rbind, lapply(c(10, 50), function(p) {
    pkg <- bayesnec::ecx(f$fit, ecx_val = p, type = "absolute")
    drawn <- ecx_from_draw(top, bot, beta, nec, p)
    data.frame(
      endpoint = paste0("ErC", p),
      true_ecx = round(true_ecx(truth, p), 4),
      pkg_ecx = round(unname(pkg[1]), 4),
      pkg_lower = round(unname(pkg[2]), 4),
      pkg_upper = round(unname(pkg[3]), 4),
      posterior_analytic = round(median(drawn), 4),
      pa_lower = round(unname(stats::quantile(drawn, 0.025)), 4),
      pa_upper = round(unname(stats::quantile(drawn, 0.975)), 4))
  }))
  out$pkg_minus_true <- round(out$pkg_ecx - out$true_ecx, 4)
  out$analytic_minus_true <- round(out$posterior_analytic - out$true_ecx, 4)
  cat("\n")
  print(out, row.names = FALSE)

  # NSEC for completeness: its target is the true nec only approximately.
  ns <- bayesnec::nsec(f$fit)
  cat("\nNSEC:", round(unname(ns[1]), 4), "[", round(unname(ns[2]), 4), ",",
      round(unname(ns[3]), 4), "]   true nec:", truth$nec, "\n")
  invisible(out)
}

## Negligible noise: the posterior should sit on the true curve, so anything
## left is definitional or numerical.
tiny <- report(1e-4, "NEGLIGIBLE NOISE")
## The sweep's actual noise level, for comparison.
real <- report(sigma_0_at(0.0956, R_ref = 2.3, t = 7), "SWEEP NOISE LEVEL")

saveRDS(list(tiny = tiny, real = real), "analysis/phase5_estimand_check.rds")

cat("\n--- how to read this ---\n")
cat("At negligible noise, pkg_minus_true near zero means ecx() and true_ecx()\n")
cat("are the same quantity and the coverage shortfall is a fitting property.\n")
cat("A systematic non-zero pkg_minus_true with analytic_minus_true near zero\n")
cat("means the estimator and the estimand differ, and the bias/coverage\n")
cat("columns are measuring that mismatch rather than anything about the arms.\n")
