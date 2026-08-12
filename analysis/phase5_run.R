## Phase 5 -- the simulation sweep.
##
## Scale is set by environment variables so the same script runs the pilot-sized
## and the full-sized sweep:
##
##   N_ITER   iterations per cell (default 50)
##   DESIGN   "core" (9 cells: delta x top_factor at R = 2.3)
##            "rsep" (12 cells: core plus an R sweep that separates R from Delta)
##            "full" (72 cells: the complete factorial)
##   WORKERS  parallel iterations (default STUDY_CORES)
##
## Results are written per cell to analysis/phase5/<cell>.rds as they complete,
## so an interrupted sweep resumes rather than restarts.

source("R/setup.R")
source("R/data_prep.R")
source("R/diagnostics.R")
source("R/simulate.R")
load_bayesnec()
source("R/arms.R")
source("R/metrics.R")

use_compile_cache()

N_ITER <- as.integer(Sys.getenv("N_ITER", "50"))
DESIGN <- Sys.getenv("DESIGN", "core")
WORKERS <- as.integer(Sys.getenv("WORKERS", as.character(STUDY_CORES)))
ARMS <- c("A", "B1", "B2", "B3", "C", "D")
OUT <- "analysis/phase5"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## Each worker runs one whole iteration (all six arms). Chains run sequentially
## inside a worker so that WORKERS, not WORKERS x chains, is the core count.
options(mc.cores = 1L)
MCMC_SIM <- list(chains = 4L, iter = 2000L, warmup = 1000L, adapt_delta = 0.95,
                 max_treedepth = 12L, seed = 20260812L)

cal <- calibrate_sigma(read_sgr("c_proliferum"))

## ----------------------------------------------------------------- cells ----
## Delta levels are the discarded effect fractions measured on the real
## datasets, so a simulation cell and a real dataset sit on the same axis:
## the group-mean lower bound (3.7) and the arm-A fitted value (8.0) on
## c_proliferum, plus a mild level. The fitted asymptote is deeper than any
## group mean because the curve never plateaus in these designs, so both ends
## of that range are represented.
DELTAS <- c(2, 4, 8)
TOPS <- c(0.8, 1.0, 2.0)
RS <- c(2.3, 3.3, 17, 73)

cells <- switch(
  DESIGN,
  core = expand.grid(delta = DELTAS, top_factor = TOPS, R = 2.3,
                     sigma_ratio = cal$sigma_ratio,
                     KEEP.OUT.ATTRS = FALSE),
  rsep = rbind(
    expand.grid(delta = DELTAS, top_factor = TOPS, R = 2.3,
                sigma_ratio = cal$sigma_ratio, KEEP.OUT.ATTRS = FALSE),
    # R varied with Delta and design held fixed: the only cells that can tell
    # the two apart, because the four real datasets confound them.
    expand.grid(delta = 4, top_factor = 2.0, R = setdiff(RS, 2.3),
                sigma_ratio = cal$sigma_ratio, KEEP.OUT.ATTRS = FALSE)
  ),
  full = expand.grid(delta = DELTAS, top_factor = TOPS, R = RS,
                     sigma_ratio = c(1, cal$sigma_ratio),
                     KEEP.OUT.ATTRS = FALSE),
  stop("Unknown DESIGN: ", DESIGN)
)
cells$cell <- sprintf("d%.1f_t%.1f_R%.1f_s%.1f", cells$delta, cells$top_factor,
                      cells$R, cells$sigma_ratio)

cat("DESIGN =", DESIGN, "|", nrow(cells), "cells x", N_ITER, "iterations x",
    length(ARMS), "arms =", nrow(cells) * N_ITER * length(ARMS), "fits\n")
cat("WORKERS =", WORKERS, "\n\n")

## ------------------------------------------------------------ one iteration --
run_iteration <- function(i, truth, design, sigma_ratio) {
  sim <- simulate_dataset(truth, design, cv_control = cal$cv_control,
                          sigma_ratio = sigma_ratio, seed = 7e5 + i)
  a_prep <- prepare_sgr(sim, "raw", meta = sim_meta())
  pr <- arm_prior(a_prep$x, a_prep$y)
  mcmc <- MCMC_SIM; mcmc$seed <- MCMC_SIM$seed + i

  out <- list()
  fitA <- NULL
  for (arm in c("A", setdiff(ARMS, c("A", "D")))) {
    f <- try(suppressMessages(
      fit_arm(arm, sim, prior = pr, mcmc = mcmc, meta = sim_meta())),
      silent = TRUE)
    if (inherits(f, "try-error")) {
      out[[arm]] <- data.frame(arm = arm, endpoint = NA, estimate = NA,
                               lower = NA, upper = NA, divergences = NA,
                               max_rhat = NA, ok = FALSE,
                               error = as.character(f))
      next
    }
    if (arm == "A") fitA <- f
    dgn <- fit_diagnostics(f$fit)
    et <- endpoint_table(f$fit, arm, "sim")
    out[[arm]] <- data.frame(arm = arm, endpoint = et$endpoint,
                             estimate = et$estimate, lower = et$lower,
                             upper = et$upper, divergences = dgn$divergences,
                             max_rhat = dgn$max_rhat, ok = TRUE,
                             error = NA_character_)
  }
  ## Arm D truncates at THIS iteration's arm-A crossing. Using the true
  ## crossing would leak the truth into the design and make arm D look better
  ## than any analyst could make it.
  if (!is.null(fitA)) {
    xc <- try(zero_crossing(fitA$fit), silent = TRUE)
    f <- if (inherits(xc, "try-error") || !is.finite(xc)) xc else
      try(suppressMessages(
        fit_arm("D", sim, prior = pr, crossing = xc, mcmc = mcmc,
                meta = sim_meta())), silent = TRUE)
    if (inherits(f, "try-error")) {
      out[["D"]] <- data.frame(arm = "D", endpoint = NA, estimate = NA,
                               lower = NA, upper = NA, divergences = NA,
                               max_rhat = NA, ok = FALSE,
                               error = as.character(f))
    } else {
      dgn <- fit_diagnostics(f$fit)
      et <- endpoint_table(f$fit, "D", "sim")
      out[["D"]] <- data.frame(arm = "D", endpoint = et$endpoint,
                               estimate = et$estimate, lower = et$lower,
                               upper = et$upper, divergences = dgn$divergences,
                               max_rhat = dgn$max_rhat, ok = TRUE,
                               error = NA_character_)
    }
  }
  cbind(iteration = i, do.call(rbind, out), row.names = NULL)
}

## ------------------------------------------------------------------ sweep ----
for (k in seq_len(nrow(cells))) {
  cel <- cells[k, ]
  path <- file.path(OUT, paste0(cel$cell, ".rds"))
  if (file.exists(path)) {
    cat("[skip]", cel$cell, "\n"); next
  }
  truth <- sim_truth(R = cel$R, delta = cel$delta, t = 7)
  design <- sim_design(truth, top_factor = cel$top_factor, n_conc = 12,
                       n_rep = 5, n_control = 10)
  resolution_check <- resolves_transition(truth, design)
  if (!resolution_check$adequate) {
    cat(sprintf("[skip] %s -- only %d concentrations between ErC10 and the",
                cel$cell, resolution_check$n_in_transition),
        "zero crossing; the design cannot resolve the curve.\n")
    next
  }
  t0 <- Sys.time()
  res <- parallel::mclapply(seq_len(N_ITER), run_iteration, truth = truth,
                            design = design, sigma_ratio = cel$sigma_ratio,
                            mc.cores = WORKERS, mc.preschedule = FALSE)
  bad <- vapply(res, inherits, logical(1), "try-error")
  res <- do.call(rbind, res[!bad])
  res <- cbind(cel[rep(1, nrow(res)), ], res, row.names = NULL)
  res <- merge(res, true_endpoints(truth), by = "endpoint", all.x = TRUE)
  attr(res, "truth") <- truth
  attr(res, "worker_failures") <- sum(bad)
  saveRDS(res, path)
  cat(sprintf("[done] %s  %.1f min  %d worker failures\n", cel$cell,
              as.numeric(difftime(Sys.time(), t0, units = "mins")), sum(bad)))
}

cat("\nSweep complete. Summarise with analysis/phase5_report.R\n")
