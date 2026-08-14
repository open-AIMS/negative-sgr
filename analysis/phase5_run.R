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
PRIOR_MODE <- Sys.getenv("PRIOR_MODE", "fixed")
stopifnot(PRIOR_MODE %in% c("fixed", "per_iteration"))
OUT <- "analysis/phase5"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## Each worker runs one whole iteration (all six arms). Chains run sequentially
## inside a worker so that WORKERS, not WORKERS x chains, is the core count.
options(mc.cores = 1L)
MCMC_SIM <- list(chains = 4L, iter = 2000L, warmup = 1000L, adapt_delta = 0.95,
                 max_treedepth = 12L, seed = 20260812L)

cal <- calibrate_sigma(read_sgr("c_proliferum"))

## Residual scale is held ABSOLUTE across cells, anchored at R = 2.3.
##
## Under the "cv" mode the whole generating model is exactly scale-equivariant
## in the growth rate, so R cannot change any endpoint and the R sweep would be
## a null manipulation -- see the long note in R/simulate.R. Anchoring here means
## the nine R = 2.3 core cells are numerically identical to the pilot, while the
## R sweep varies signal-to-noise the way a real change in control growth does.
SIGMA_0 <- sigma_0_at(cal$cv_control, R_ref = 2.3, t = 7)

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
cat("WORKERS =", WORKERS, "| PRIOR_MODE =", PRIOR_MODE, "\n\n")

## ------------------------------------------------------------ one iteration --
run_iteration <- function(i, truth, design, sigma_ratio, prior_ref = NULL) {
  sim <- simulate_dataset(truth, design, cv_control = cal$cv_control,
                          sigma_ratio = sigma_ratio, seed = 7e5 + i,
                          sigma_mode = "absolute", sigma_0_abs = SIGMA_0)
  # The fraction of responses the convention actually discards, recorded per
  # iteration because it is the quantity R is expected to act through.
  f_neg <- mean(sim$sgr < 0)
  a_prep <- prepare_sgr(sim, "raw", meta = sim_meta())
  ## Prior held fixed within a cell (PRIOR_MODE = "fixed"), or rebuilt from each
  ## simulated dataset as `bnec()` would (PRIOR_MODE = "per_iteration").
  ##
  ## Fixed is the default because of compilation, and the size of that effect is
  ## easy to underestimate. brms writes prior constants into the Stan program as
  ## literals, and `bnec()`'s gaussian defaults are functions of the response --
  ## `normal(quantile(y, 0.9), 2.5 * sd(y))` for `top`, likewise for `bot`. Two
  ## simulated datasets from the SAME cell therefore produce two different Stan
  ## programs: measured, 40 iterations gave 40 distinct prior strings and so 40
  ## distinct model hashes. Every fit then recompiles from scratch at 3-5 minutes
  ## a time, and the first sweep left 1805 files and 2.6 GB in the compile cache
  ## for one cell. Holding the prior fixed collapses that to one program per arm
  ## per cell.
  ##
  ## It is also the cleaner comparison. The study already gives every arm within
  ## an iteration the same prior (see `arm_prior()`), so that an arm contrast is
  ## a likelihood contrast; fixing it across iterations extends the same
  ## principle, and removes a nuisance source of variation shared by all arms.
  ## The cost is fidelity to what a practitioner would get, since their prior
  ## does move with their data -- checked in `analysis/phase5_prior_check.R`
  ## rather than assumed.
  pr <- if (identical(PRIOR_MODE, "fixed") && !is.null(prior_ref)) prior_ref else
    arm_prior(a_prep$x, a_prep$y)
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
  ##
  ## `zero_crossing()` returns Inf when the fitted curve never reaches zero
  ## growth inside the tested range, and `truncate_at_crossing(dat, Inf)` keeps
  ## every row -- so arm D is then *identical* to arm A, as `zero_crossing()`'s
  ## own documentation says. It is reported as arm A's result rather than
  ## refitted, which is both correct and saves the most expensive fit in the
  ## study in exactly the cells where it would be redundant.
  ##
  ## This is the branch that destroyed the first sweep. The old code did
  ## `f <- xc` for the non-finite case and then tested only
  ## `inherits(f, "try-error")`, which is FALSE for a bare Inf, so it fell
  ## through to `f$fit` and every worker died with "$ operator is invalid for
  ## atomic vectors". It bit hardest in the `top_factor = 0.8` cells, where the
  ## top concentration lies *below* the true crossing by construction and Inf is
  ## the expected return, not the exception: 239 of 240 iterations were lost.
  if (!is.null(fitA)) {
    xc <- try(zero_crossing(fitA$fit), silent = TRUE)
    d_degenerate <- !inherits(xc, "try-error") && !is.finite(xc)
    f <- if (inherits(xc, "try-error")) xc else if (d_degenerate) NULL else
      try(suppressMessages(
        fit_arm("D", sim, prior = pr, crossing = xc, mcmc = mcmc,
                meta = sim_meta())), silent = TRUE)
    if (d_degenerate) {
      a_rows <- out[["A"]]
      a_rows$arm <- "D"
      a_rows$error <- "crossing not reached in design; arm D == arm A"
      out[["D"]] <- a_rows
    } else if (inherits(f, "try-error")) {
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
  cbind(iteration = i, f_neg = f_neg, do.call(rbind, out), row.names = NULL)
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
  ## The cell's reference prior, built from a dataset drawn with a seed no
  ## analysis iteration uses, so the prior is never derived from a dataset it is
  ## then used to fit.
  ref <- simulate_dataset(truth, design, cv_control = cal$cv_control,
                          sigma_ratio = cel$sigma_ratio, seed = 6e5 + k,
                          sigma_mode = "absolute", sigma_0_abs = SIGMA_0)
  ref_prep <- prepare_sgr(ref, "raw", meta = sim_meta())
  prior_ref <- arm_prior(ref_prep$x, ref_prep$y)

  ## Warm the compile cache SERIALLY before forking. Eighteen workers each
  ## meeting an uncompiled model at once means eighteen concurrent C++ compiles
  ## at roughly 1.5 GB apiece, which is how the first sweep exhausted a 31 GB
  ## box and lost its workers with no R-level error to show for it. One
  ## iteration up front compiles every distinct program the cell needs; the rest
  ## then reuse the executables.
  cat("  warming compile cache for", cel$cell, "... ")
  tw <- Sys.time()
  warm <- try(run_iteration(0L, truth, design, cel$sigma_ratio, prior_ref),
              silent = TRUE)
  cat(sprintf("%.1f min%s\n", as.numeric(difftime(Sys.time(), tw, units = "mins")),
              if (inherits(warm, "try-error")) " [FAILED]" else ""))
  if (inherits(warm, "try-error")) {
    cat("[abort]", cel$cell, "-- warm-up failed:\n", as.character(warm), "\n")
    next
  }

  t0 <- Sys.time()
  res <- parallel::mclapply(seq_len(N_ITER), run_iteration, truth = truth,
                            design = design, sigma_ratio = cel$sigma_ratio,
                            prior_ref = prior_ref,
                            mc.cores = WORKERS, mc.preschedule = FALSE)
  bad <- vapply(res, inherits, logical(1), "try-error")
  if (any(bad)) {
    # Silently dropping these is what let a 239/240 failure rate look like a
    # completed cell. Write them down.
    msgs <- unique(vapply(res[bad], as.character, character(1)))
    cat(sprintf("  %d/%d worker failures in %s; distinct messages:\n",
                sum(bad), length(bad), cel$cell))
    for (m in utils::head(msgs, 5)) cat("    ", trimws(m), "\n")
    writeLines(msgs, file.path(OUT, paste0(cel$cell, ".failures.txt")))
  }
  res <- do.call(rbind, res[!bad])
  if (is.null(res) || !nrow(res)) {
    cat("[abort]", cel$cell, "-- every iteration failed; cell not written.\n")
    next
  }
  res <- cbind(cel[rep(1, nrow(res)), ], res, row.names = NULL)
  res <- merge(res, true_endpoints(truth), by = "endpoint", all.x = TRUE)
  attr(res, "truth") <- truth
  attr(res, "worker_failures") <- sum(bad)
  saveRDS(res, path)
  cat(sprintf("[done] %s  %.1f min  %d worker failures\n", cel$cell,
              as.numeric(difftime(Sys.time(), t0, units = "mins")), sum(bad)))
}

cat("\nSweep complete. Summarise with analysis/phase5_report.R\n")
