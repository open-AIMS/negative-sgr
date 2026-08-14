## Does arm A reach nominal coverage when the model is CORRECTLY specified?
##
## The sweep cells generate heteroscedastic data (sigma_ratio ~ 8, calibrated
## from c_proliferum) and fit it homoscedastic, because `bnec()` has no route to
## a distributional sigma (bayesnec issue #191). That is deliberate and
## realistic, but it means every arm's absolute coverage is degraded by a
## misspecification shared across arms -- the first sweep cell gives arm A a 95%
## interval covering ErC50 in 77% of iterations.
##
## Before any of that is reported as a property of the *arms*, it has to be
## shown to be a property of the *variance model*. This script refits the same
## machinery with sigma_ratio = 1, where the fitted model matches the generating
## model exactly. Arm A must then reach roughly nominal coverage. If it does
## not, the estimand or the interval extraction is wrong and the sweep's
## coverage numbers mean nothing.
##
## Run: Rscript analysis/phase5_coverage_validation.R   (N_VAL, default 60)

source("R/setup.R")
source("R/data_prep.R")
source("R/diagnostics.R")
source("R/simulate.R")
load_bayesnec()
source("R/arms.R")
source("R/metrics.R")

use_compile_cache()
options(mc.cores = 1L)

N_VAL <- as.integer(Sys.getenv("N_VAL", "60"))
WORKERS_VAL <- as.integer(Sys.getenv("WORKERS_VAL", "4"))
MCMC_SIM <- list(chains = 4L, iter = 2000L, warmup = 1000L, adapt_delta = 0.95,
                 max_treedepth = 12L, seed = 20260812L)

cal <- calibrate_sigma(read_sgr("c_proliferum"))
SIGMA_0 <- sigma_0_at(cal$cv_control, R_ref = 2.3, t = 7)

## Mid cell, transition fully inside the design so all three endpoints are
## estimable, run at both variance settings.
truth <- sim_truth(R = 2.3, delta = 4, t = 7)
design <- sim_design(truth, top_factor = 1.0, n_conc = 12, n_rep = 5,
                     n_control = 10)
te <- true_endpoints(truth)

ref <- simulate_dataset(truth, design, cv_control = cal$cv_control,
                        sigma_ratio = 1, seed = 6e5 + 99,
                        sigma_mode = "absolute", sigma_0_abs = SIGMA_0)
rp <- prepare_sgr(ref, "raw", meta = sim_meta())
prior_ref <- arm_prior(rp$x, rp$y)

one <- function(i, sigma_ratio) {
  sim <- simulate_dataset(truth, design, cv_control = cal$cv_control,
                          sigma_ratio = sigma_ratio, seed = 8e5 + i,
                          sigma_mode = "absolute", sigma_0_abs = SIGMA_0)
  mcmc <- MCMC_SIM; mcmc$seed <- MCMC_SIM$seed + i
  f <- try(suppressMessages(fit_arm("A", sim, prior = prior_ref, mcmc = mcmc,
                                    meta = sim_meta())), silent = TRUE)
  if (inherits(f, "try-error")) return(NULL)
  et <- endpoint_table(f$fit, "A", "sim")
  et$sigma_ratio <- sigma_ratio
  et$iteration <- i
  et
}

results <- list()
for (sr in c(1, cal$sigma_ratio)) {
  cat("\n=== sigma_ratio =", round(sr, 3),
      if (sr == 1) "(correctly specified)" else "(misspecified, as in the sweep)",
      "===\n")
  # Warm the cache serially first, for the reason in phase5_run.R.
  invisible(one(0L, sr))
  res <- parallel::mclapply(seq_len(N_VAL), one, sigma_ratio = sr,
                            mc.cores = WORKERS_VAL, mc.preschedule = FALSE)
  ok <- !vapply(res, is.null, logical(1)) &
    !vapply(res, inherits, logical(1), "try-error")
  cat(sum(!ok), "failures of", N_VAL, "\n")
  results[[length(results) + 1]] <- do.call(rbind, res[ok])
}

all <- do.call(rbind, results)
all <- merge(all, te, by = "endpoint", all.x = TRUE)
saveRDS(all, "analysis/phase5_coverage_validation.rds")

summ <- do.call(rbind, lapply(
  split(all, list(all$endpoint, all$sigma_ratio), drop = TRUE),
  function(d) data.frame(
    endpoint = d$endpoint[1],
    sigma_ratio = round(d$sigma_ratio[1], 3),
    n = nrow(d),
    truth = round(d$truth[1], 3),
    mean_estimate = round(mean(d$estimate), 3),
    bias = round(mean(d$estimate) - d$truth[1], 4),
    coverage = round(mean(d$lower <= d$truth & d$upper >= d$truth), 3),
    coverage_mcse = round(sqrt(0.95 * 0.05 / nrow(d)), 3),
    mean_width = round(mean(d$upper - d$lower), 3))))

cat("\n--- arm A, nominal 95% ---\n")
print(summ[order(summ$endpoint, summ$sigma_ratio), ], row.names = FALSE)
utils::write.csv(summ, "analysis/phase5_coverage_validation.csv",
                 row.names = FALSE)

cat("\nRead: at sigma_ratio = 1 the fitted model IS the generating model, so\n",
    "coverage should sit near 0.95 within about 2 MCSE. If it does, the sweep's\n",
    "low absolute coverage is the heteroscedasticity misspecification and the\n",
    "arm CONTRASTS remain interpretable. If it does not, the estimand or the\n",
    "interval extraction is at fault and the coverage column must be withdrawn.\n")
