## Phase 5 pilot -- measure before committing.
##
## The plan requires a budget to be published before the factorial is run. This
## script fits one cell a handful of times under two MCMC settings, records the
## wall-clock cost of every fit, and converts it into core-hours for candidate
## factorial sizes and iteration counts. It does NOT run the study.

source("R/setup.R")
source("R/data_prep.R")
source("R/diagnostics.R")
source("R/simulate.R")
load_bayesnec()
source("R/arms.R")
source("R/metrics.R")

use_compile_cache()
options(mc.cores = min(STUDY_CORES, 4L))

N_PILOT <- as.integer(Sys.getenv("N_PILOT", "8"))
ARMS <- c("A", "B1", "B2", "B3", "C")

## Calibrate the generating variance from the assay the truth is modelled on.
cal <- calibrate_sigma(read_sgr("c_proliferum"))
cat("sigma calibration from c_proliferum:\n"); str(cal)

truth <- sim_truth(R = 2.3, delta = 7.97, t = 7)
cat("\ntruth:\n"); str(truth)
cat("\ntrue endpoints:\n"); print(true_endpoints(truth))
cat("\nanalytic vs grid ErCx (must agree to within the grid spacing):\n")
print(check_ecx_estimator(truth))

design <- sim_design(truth, top_factor = 2, n_conc = 12, n_rep = 5,
                     n_control = 10)
cat("\ndesign: n =", nrow(design), "over", length(unique(design$x)),
    "concentrations; top =", signif(max(design$x), 4),
    "vs zero-crossing", signif(truth$zero_crossing, 4), "\n")
res <- resolves_transition(truth, design)
cat("concentrations between ErC10 and the zero crossing:", res$n_in_transition,
    if (res$adequate) "-- adequate\n" else "-- TOO FEW, the curve is a step\n")
if (!res$adequate) stop("Design does not resolve the transition; recalibrate.")

## Two settings: the Phase 3 settings, and a leaner set for the simulation.
SETTINGS <- list(
  full = list(chains = 4L, iter = 4000L, warmup = 2000L, adapt_delta = 0.99,
              max_treedepth = 12L, seed = 1L),
  lean = list(chains = 4L, iter = 2000L, warmup = 1000L, adapt_delta = 0.95,
              max_treedepth = 12L, seed = 1L)
)

rows <- list()
for (setting in names(SETTINGS)) {
  mcmc <- SETTINGS[[setting]]
  for (i in seq_len(N_PILOT)) {
    sim <- simulate_dataset(truth, design, cv_control = cal$cv_control,
                            sigma_ratio = cal$sigma_ratio, seed = 1000 + i)
    ## One prior per iteration, built from that iteration's arm-A data and
    ## handed to every arm -- the same rule as the real-data phases, for the
    ## same reason: the default prior is a function of the response vector, so
    ## letting each arm take its own default would confound the likelihood
    ## contrast with a prior contrast.
    a_prep <- prepare_sgr(sim, "raw", meta = sim_meta())
    pr <- arm_prior(a_prep$x, a_prep$y)
    ## Arm D needs this iteration's own arm-A crossing, never the truth.
    fits <- list()
    for (arm in ARMS) {
      mcmc_i <- mcmc; mcmc_i$seed <- mcmc$seed + i
      f <- try(suppressMessages(
        fit_arm(arm, sim, prior = pr, mcmc = mcmc_i, meta = sim_meta())
      ), silent = TRUE)
      if (inherits(f, "try-error")) {
        rows[[length(rows) + 1]] <- data.frame(
          setting = setting, iteration = i, arm = arm, elapsed = NA,
          divergences = NA, max_rhat = NA, ok = FALSE,
          error = as.character(f), stringsAsFactors = FALSE)
        next
      }
      fits[[arm]] <- f
      dg <- fit_diagnostics(f$fit)
      rows[[length(rows) + 1]] <- data.frame(
        setting = setting, iteration = i, arm = arm, elapsed = f$elapsed,
        divergences = dg$divergences, max_rhat = dg$max_rhat, ok = TRUE,
        error = NA_character_, stringsAsFactors = FALSE)
    }
    if (!is.null(fits$A)) {
      xc <- zero_crossing(fits$A$fit)
      f <- try(suppressMessages(
        fit_arm("D", sim, prior = pr, crossing = xc, mcmc = mcmc,
                meta = sim_meta())), silent = TRUE)
      ok <- !inherits(f, "try-error")
      rows[[length(rows) + 1]] <- data.frame(
        setting = setting, iteration = i, arm = "D",
        elapsed = if (ok) f$elapsed else NA,
        divergences = if (ok) fit_diagnostics(f$fit)$divergences else NA,
        max_rhat = if (ok) fit_diagnostics(f$fit)$max_rhat else NA,
        ok = ok, error = if (ok) NA_character_ else as.character(f),
        stringsAsFactors = FALSE)
    }
    cat("."); utils::flush.console()
  }
  cat(" ", setting, "done\n")
}

timing <- do.call(rbind, rows)
utils::write.csv(timing, "analysis/phase5_pilot_timing.csv", row.names = FALSE)

cat("\n===== PER-FIT COST =====\n")
agg <- aggregate(elapsed ~ setting + arm, timing, function(z)
  c(mean = mean(z), sd = stats::sd(z)))
print(agg)
per_iteration <- tapply(timing$elapsed, timing$setting,
                        function(z) sum(z, na.rm = TRUE)) /
  tapply(timing$iteration, timing$setting, function(z) length(unique(z)))
cat("\nseconds per simulation iteration (all 6 arms):\n")
print(round(per_iteration, 1))

cat("\n===== FAILURES =====\n")
print(aggregate(ok ~ setting + arm, timing, function(z) mean(!z)))

cat("\n===== MCSE ARITHMETIC =====\n")
for (m in c(0.01, 0.02, 0.04)) {
  cat(sprintf("MCSE %.0f%% on 95%% coverage needs %d iterations per cell\n",
              m * 100, iterations_for_coverage_mcse(m)))
}

cat("\n===== BUDGET =====\n")
## Fits run one at a time here, each using `chains` cores. In the real sweep the
## cheaper arrangement is one core per chain-free fit and many fits at once, so
## core-hours -- not wall-clock -- is the transferable number.
cores_per_fit <- min(STUDY_CORES, SETTINGS$lean$chains)
for (setting in names(per_iteration)) {
  s_per_iter <- per_iteration[[setting]]
  for (n_cells in c(9, 12, 18, 72)) {
    for (n_iter in c(50, 240, 950)) {
      ch <- s_per_iter * n_cells * n_iter * cores_per_fit / 3600
      cat(sprintf("%-5s %2d cells x %4d iter: %8.1f core-hours (%6.1f h on %d cores)\n",
                  setting, n_cells, n_iter, ch, ch / STUDY_CORES, STUDY_CORES))
    }
  }
}
cat("\nMCSE at those iteration counts: 50 -> +/-3.1%, 240 -> +/-1.4%,",
    "950 -> +/-0.7% on 95% coverage.\n")
