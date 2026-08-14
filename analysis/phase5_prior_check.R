## Does holding the prior fixed within a cell change the answer?
##
## Phase 5 fits every iteration of a cell with one reference prior rather than
## rebuilding `bnec()`'s data-dependent default from each simulated dataset.
## The reason is mechanical -- a per-dataset prior writes different constants
## into the Stan program, so every fit recompiles -- but the question of whether
## it changes the inference is scientific, and is answered here rather than
## asserted.
##
## Design: one cell, N iterations, arm A only. Each dataset is fitted twice,
## once with its own default prior and once with the cell's reference prior,
## with the same seed. The comparison is the endpoint difference against the
## Monte Carlo spread of the endpoint across iterations -- a prior effect that
## is small relative to sampling variation cannot matter for coverage or bias,
## which are the study's outputs.
##
## Run: Rscript analysis/phase5_prior_check.R   (N_CHECK iterations, default 5)

source("R/setup.R")
source("R/data_prep.R")
source("R/diagnostics.R")
source("R/simulate.R")
load_bayesnec()
source("R/arms.R")
source("R/metrics.R")

use_compile_cache()
options(mc.cores = min(4L, STUDY_CORES))

N_CHECK <- as.integer(Sys.getenv("N_CHECK", "5"))
MCMC_SIM <- list(chains = 4L, iter = 2000L, warmup = 1000L, adapt_delta = 0.95,
                 max_treedepth = 12L, seed = 20260812L)

cal <- calibrate_sigma(read_sgr("c_proliferum"))
SIGMA_0 <- sigma_0_at(cal$cv_control, R_ref = 2.3, t = 7)

## The mid cell of the core design: delta = 4, top_factor = 1.0, R = 2.3. Chosen
## because the transition is fully inside the design, so the endpoints are all
## estimable and any prior effect has somewhere to show.
truth <- sim_truth(R = 2.3, delta = 4, t = 7)
design <- sim_design(truth, top_factor = 1.0, n_conc = 12, n_rep = 5,
                     n_control = 10)

ref <- simulate_dataset(truth, design, cv_control = cal$cv_control,
                        sigma_ratio = cal$sigma_ratio, seed = 6e5 + 5,
                        sigma_mode = "absolute", sigma_0_abs = SIGMA_0)
prior_ref <- arm_prior(prepare_sgr(ref, "raw", meta = sim_meta())$x,
                       prepare_sgr(ref, "raw", meta = sim_meta())$y)

cat("reference prior:\n"); print(as.data.frame(prior_ref)[, c("prior", "nlpar")])
cat("\n")

rows <- list()
for (i in seq_len(N_CHECK)) {
  sim <- simulate_dataset(truth, design, cv_control = cal$cv_control,
                          sigma_ratio = cal$sigma_ratio, seed = 7e5 + i,
                          sigma_mode = "absolute", sigma_0_abs = SIGMA_0)
  ap <- prepare_sgr(sim, "raw", meta = sim_meta())
  pr_own <- arm_prior(ap$x, ap$y)
  mcmc <- MCMC_SIM; mcmc$seed <- MCMC_SIM$seed + i

  for (mode in c("own", "fixed")) {
    pr <- if (mode == "own") pr_own else prior_ref
    f <- try(suppressMessages(fit_arm("A", sim, prior = pr, mcmc = mcmc,
                                      meta = sim_meta())), silent = TRUE)
    if (inherits(f, "try-error")) {
      cat("iteration", i, mode, "FAILED\n"); next
    }
    et <- endpoint_table(f$fit, "A", "sim")
    rows[[length(rows) + 1]] <- data.frame(iteration = i, mode = mode,
                                           endpoint = et$endpoint,
                                           estimate = et$estimate,
                                           lower = et$lower, upper = et$upper)
  }
  cat("iteration", i, "done\n")
}

res <- do.call(rbind, rows)
saveRDS(res, "analysis/phase5_prior_check.rds")

w <- stats::reshape(res[, c("iteration", "mode", "endpoint", "estimate")],
                    idvar = c("iteration", "endpoint"), timevar = "mode",
                    direction = "wide")
w$diff <- w$estimate.fixed - w$estimate.own
w$rel <- w$diff / w$estimate.own

cat("\n--- per-iteration differences (fixed - own) ---\n")
print(w, row.names = FALSE)

cat("\n--- summary by endpoint ---\n")
summ <- do.call(rbind, lapply(split(w, w$endpoint), function(d) {
  data.frame(endpoint = d$endpoint[1],
             mean_abs_rel_diff = mean(abs(d$rel)),
             max_abs_rel_diff = max(abs(d$rel)),
             # The comparator: how much the endpoint moves between iterations
             # under the fixed prior. A prior effect well inside this is
             # invisible to bias and coverage.
             mc_sd_of_estimate = stats::sd(d$estimate.own),
             ratio_priordiff_to_mcsd = mean(abs(d$diff)) /
               stats::sd(d$estimate.own))
}))
print(summ, row.names = FALSE)
utils::write.csv(summ, "analysis/phase5_prior_check.csv", row.names = FALSE)

cat("\nInterpretation: ratio_priordiff_to_mcsd well below 1 means switching the\n",
    "prior moves an endpoint by much less than the endpoint moves between\n",
    "simulated datasets, so the fixed prior cannot materially change bias or\n",
    "coverage. Above about 0.5 it should be reported, and PRIOR_MODE =\n",
    "'per_iteration' reruns the sweep the slow, faithful way.\n")
