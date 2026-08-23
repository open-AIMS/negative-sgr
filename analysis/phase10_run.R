## Phase 10 -- the simulation under bayesnec's model-averaged workflow.
##
## Phases 5 and 7 fixed one equation in advance. That is right for scoring a
## convention -- the generating model is fitted, so the only misspecification is
## the one each arm deliberately introduces -- and it is not the workflow
## bayesnec recommends. Under B3, E and F the constrained curve visibly cannot
## follow the data, so it flattens its descent and pushes `nec` right, and the
## toxicity estimates follow. A more flexible equation might follow the floored
## data adequately OVER THE REGION ErC10 and ErC50 are read from, without ever
## going below zero. If so, the reported bias and the collapsed coverage are
## artefacts of fixing `nec4param` rather than properties of the convention.
##
## Five arms: A (reference, intact), C (left-censored), B3 (floored and pinned),
## E (floored, scaled, Beta), F (floored, Gamma). Three cells, run 12 -> 8 -> 9.
##
## THE PAIRING IS THE POINT. Seeds are reused from `phase8_run.R` -- `6e5 + k`
## for the cell's reference-prior dataset and `7e5 + i` for iteration i -- so
## every simulated dataset is bit-identical to the one the single-model arms saw
## and averaged-versus-single is a PAIRED contrast rather than two independent
## estimates. `k` is the ORIGINAL cell index into the twelve-cell grid and is
## never renumbered when cells are subset; `phase8_run.R` renumbers and would
## give cell 12 a different reference prior when run alone. `p10_cells()` carries
## the index as a column for exactly this reason, and the pilot verifies the
## pairing against the stored Phase 8 output rather than trusting the arithmetic.
##
## Writes are chunked and additive, as in Phases 7 and 8: each block goes to its
## own file under a temporary name and is renamed, so a file that exists is
## complete, a resume skips it by filename, and an interrupt costs at most one
## block. Nothing already on disk is ever rewritten.
##
## Run:    CELLS=12,8,9 N_ITER_FROM=1 N_ITER_TO=100 CHUNK=20 WORKERS=18 \
##           Rscript analysis/phase10_run.R
## Resume: re-issue the identical command.
## Pilot:  PHASE10_PILOT=1 Rscript analysis/phase10_run.R      (Gate 2)

source("R/setup.R"); source("R/setup_phase10.R")
source("R/data_prep.R"); source("R/diagnostics.R"); source("R/simulate.R")
load_bayesnec_p10(); source("R/arms.R"); source("R/metrics.R")
source("R/model_average.R")
p10_assert_arms_ready()
use_compile_cache()
options(mc.cores = 1L)          # one whole iteration per worker, chains serial

PILOT     <- identical(Sys.getenv("PHASE10_PILOT", "0"), "1")
ITER_FROM <- as.integer(Sys.getenv("N_ITER_FROM", if (PILOT) "1" else "1"))
ITER_TO   <- as.integer(Sys.getenv("N_ITER_TO",   if (PILOT) "2" else "100"))
CHUNK     <- as.integer(Sys.getenv("CHUNK", if (PILOT) "2" else "20"))
## The pilot's cache warm-up (iteration 0) is serial regardless -- it has to be,
## or 18 concurrent C++ compiles exhaust the box -- and that is where the pilot's
## time goes. The two analysis iterations after it may as well run together.
## Per-arm cost comes from the recorded `elapsed` column, not from wall-clock, so
## parallelism here does not blur the timing.
WORKERS   <- as.integer(Sys.getenv("WORKERS",
                                   if (PILOT) "2" else as.character(STUDY_CORES)))
OUT       <- Sys.getenv("OUT", if (PILOT) "analysis/phase10_pilot" else
                               "analysis/phase10")
TAG       <- Sys.getenv("TAG", if (PILOT) "pilot" else "ma")
CELLS     <- Sys.getenv("CELLS", "")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

cal <- calibrate_sigma(read_sgr("c_proliferum"))
SIGMA_0 <- sigma_0_at(cal$cv_control, R_ref = 2.3, t = 7)
all_cells <- p10_cells(cal$sigma_ratio)

wanted <- if (nzchar(CELLS)) {
  as.integer(strsplit(CELLS, ",")[[1]])
} else if (PILOT) {
  ## Cell 12 alone. It is where the single-model failure is most extreme (B3 and
  ## E contained the true ErC50 in none of 500 datasets) AND where averaging has
  ## its best chance, because at a 1.9% control CV the data resolve the stacking
  ## weights sharply. The hypothesis's best case is the right place to look first.
  12L
} else {
  P10_CELL_ORDER
}
cells <- all_cells[match(wanted, all_cells$cell_index), , drop = FALSE]
stopifnot(!anyNA(cells$cell_index))

## Refuse to start if another instance is already running against this OUT.
##
## Two runners sharing an output directory is not hypothetical: it happened
## during the pilot launch, when a `head` on a not-yet-created log returned a
## non-zero exit that read as a failed launch and prompted a relaunch. Blocks are
## written to a temporary name and renamed, so the loser's work is discarded
## rather than corrupted, but both processes redirect to the same log, which
## interleaves and becomes unreadable, and the duplicated cache warm-up wastes
## the best part of an hour.
##
## A PID LOCKFILE, NOT `pgrep`. The first attempt at this matched
## `--file=analysis/phase10_run` with pgrep and refused to start *anything*: at
## startup a transient wrapper process shares that command line, and excluding
## only `Sys.getpid()` leaves it looking like a second instance. Ancestry is
## fiddly to get right and the failure mode is silent refusal. A lock naming its
## own PID has no such problem, and it self-heals -- a lock left by a crashed run
## names a PID that is either gone or belongs to something else, and both are
## detectable. The cmdline check guards against PID reuse.
lock_path <- file.path(OUT, ".phase10.lock")
if (file.exists(lock_path)) {
  other <- suppressWarnings(as.integer(readLines(lock_path, warn = FALSE)[1]))
  cmd <- if (!is.na(other)) {
    tryCatch(paste(readLines(file.path("/proc", other, "cmdline"),
                             warn = FALSE), collapse = " "),
             error = function(e) "")
  } else ""
  if (!is.na(other) && other != Sys.getpid() && grepl("phase10_run", cmd)) {
    stop("Another phase10_run.R is already running (pid ", other, ").\n",
         "Stop it with  pkill -f \"[-]-file=analysis/phase10\"  or wait; ",
         "a resume skips completed blocks by filename.")
  }
  ## Stale lock: the PID is gone, or has been reused by something unrelated.
  unlink(lock_path)
}
writeLines(as.character(Sys.getpid()), lock_path)
## Released on any exit, including an error, so a failed run does not leave a
## lock that blocks the resume it is meant to permit.
invisible(reg.finalizer(environment(), function(e) unlink(lock_path),
                        onexit = TRUE))

iters <- ITER_FROM:ITER_TO
cat("Phase 10 |", nrow(cells), "cells x", length(iters), "iterations x",
    length(P10_ARMS), "arms\n")
cat("cells", paste(cells$cell_index, collapse = ", "),
    "| iterations", ITER_FROM, "-", ITER_TO, "| WORKERS =", WORKERS,
    "| CHUNK =", CHUNK, "| OUT =", OUT, "\n")
cat("pin:", BAYESNEC_WORKTREE_P10, "@", BAYESNEC_SHA_P10, "\n\n")

## ------------------------------------------------------------ one iteration --
run_iteration <- function(i, truth, design, specs, x_grid) {
  sim <- simulate_dataset(truth, design, cv_control = cal$cv_control,
                          sigma_ratio = cal$sigma_ratio, seed = 7e5 + i,
                          sigma_mode = "absolute", sigma_0_abs = SIGMA_0)
  ## `< 0`, matching phase5_run.R and phase8_run.R. On a continuous response the
  ## strict and non-strict forms differ only on an event of probability zero, but
  ## this column is the pairing check against the stored sweep, so it has to be
  ## computed the same way.
  f_neg <- mean(sim$sgr < 0)
  mcmc <- P10_MCMC; mcmc$seed <- P10_MCMC$seed + i

  ends <- list(); wts <- list(); crv <- list(); rhs <- list()

  ## One row-builder for both stages. `stage` is "all_models" (average over the
  ## whole candidate set) or "converged" (average over the models that passed
  ## the R-hat cutoff). The converged stage is the study's primary result --
  ## dropping a non-converged component before averaging is what the workflow
  ## does -- and the all-models stage is carried alongside as the sensitivity to
  ## that choice, at the cost of an extra ecx()/nsec() and no extra sampling.
  row_for <- function(obj, arm, stage, requested, elapsed, n_kept) {
    w <- p10_weights(obj, requested)
    dgn <- fit_diagnostics(p10_pull(obj))
    et <- endpoint_table(p10_pull(obj), arm, "sim")
    list(
      ends = data.frame(arm = arm, stage = stage, endpoint = et$endpoint,
                        estimate = et$estimate, lower = et$lower,
                        upper = et$upper, divergences = dgn$divergences,
                        max_rhat = dgn$max_rhat, n_requested = length(requested),
                        n_fitted = w$n_fitted[1], n_kept = n_kept,
                        elapsed = elapsed, ok = TRUE, error = NA_character_),
      wts = cbind(iteration = i, arm = arm, stage = stage, w, row.names = NULL))
  }
  fail_row <- function(arm, stage, requested, msg, n_kept = NA_integer_) {
    data.frame(arm = arm, stage = stage, endpoint = NA_character_,
               estimate = NA_real_, lower = NA_real_, upper = NA_real_,
               divergences = NA_real_, max_rhat = NA_real_,
               n_requested = length(requested), n_fitted = NA_integer_,
               n_kept = n_kept, elapsed = NA_real_, ok = FALSE, error = msg)
  }

  for (arm in P10_ARMS) {
    spec <- specs[[arm]]
    f <- try(suppressMessages(p10_fit(arm, sim, spec, mcmc = mcmc)),
             silent = TRUE)
    if (inherits(f, "try-error")) {
      ## Recorded, never filtered. A failed fit writes one row carrying the
      ## message rather than three endpoint rows, which is why a cell file's row
      ## count is not iterations x arms x stages x 3 (Trap 11).
      msg <- sub("\n.*", "", as.character(f))
      ends[[paste0(arm, ".fit")]] <- fail_row(arm, "all_models", spec$models, msg)
      next
    }

    ## Stage 1 -- the whole candidate set, before any convergence screening.
    if (isTRUE(P10_REPORT_ALL_MODELS)) {
      r <- try(row_for(f, arm, "all_models", spec$models, f$elapsed,
                       NA_integer_), silent = TRUE)
      if (inherits(r, "try-error")) {
        ends[[paste0(arm, ".all")]] <- fail_row(
          arm, "all_models", spec$models, sub("\n.*", "", as.character(r)))
      } else {
        ends[[paste0(arm, ".all")]] <- r$ends
        wts[[paste0(arm, ".all")]] <- r$wts
      }
    }

    ## Stage 2 -- drop the models that did not converge and re-stack over the
    ## rest, which is the analyst's step and the study's primary result.
    dr <- try(p10_drop_nonconverged(f, P10_RHAT_CUTOFF), silent = TRUE)
    if (inherits(dr, "try-error")) {
      ends[[paste0(arm, ".conv")]] <- fail_row(
        arm, "converged", spec$models, sub("\n.*", "", as.character(dr)))
      next
    }
    rhs[[arm]] <- cbind(iteration = i, arm = arm, dr$rhat, row.names = NULL)
    if (is.null(dr$fit)) {
      ## Every model failed the cutoff. Recorded as an arm failure rather than
      ## silently falling back to the undropped average, which would reinstate
      ## exactly the components the rule exists to remove.
      ends[[paste0(arm, ".conv")]] <- fail_row(arm, "converged", spec$models,
                                               dr$error, n_kept = 0L)
      next
    }
    r <- try(row_for(dr$fit, arm, "converged", spec$models, f$elapsed,
                     length(dr$kept)), silent = TRUE)
    if (inherits(r, "try-error")) {
      ends[[paste0(arm, ".conv")]] <- fail_row(
        arm, "converged", spec$models, sub("\n.*", "", as.character(r)),
        n_kept = length(dr$kept))
      next
    }
    ends[[paste0(arm, ".conv")]] <- r$ends
    wts[[paste0(arm, ".conv")]] <- r$wts
    ## The curve is stored for the primary (converged) fit only. Storing both
    ## would double the largest artefact for a sensitivity the endpoints already
    ## carry.
    cv <- try(p10_curve(dr$fit, x_grid, scale = p10_zb_scale(sim, arm)),
              silent = TRUE)
    if (!inherits(cv, "try-error")) {
      crv[[arm]] <- cbind(iteration = i, arm = arm, cv, row.names = NULL)
    }
  }
  list(endpoints = cbind(iteration = i, f_neg = f_neg,
                         do.call(rbind, ends), row.names = NULL),
       weights = if (length(wts)) do.call(rbind, wts) else NULL,
       curves = if (length(crv)) do.call(rbind, crv) else NULL,
       rhat = if (length(rhs)) do.call(rbind, rhs) else NULL)
}

## ------------------------------------------------------------------ sweep ----
for (r in seq_len(nrow(cells))) {
  cel <- cells[r, ]
  k <- cel$cell_index                       # ORIGINAL index; never renumbered
  truth <- sim_truth(R = cel$R, delta = cel$delta, t = 7)
  design <- sim_design(truth, top_factor = cel$top_factor, n_conc = 12,
                       n_rep = 5, n_control = 10)
  x_grid <- seq(0, max(design$x), length.out = P10_GRID_N)

  blocks <- split(iters, ceiling(seq_along(iters) / CHUNK))
  todo <- Filter(function(b) !file.exists(file.path(
            OUT, sprintf("%s__%s_i%03d-%03d.rds", cel$cell, TAG,
                         min(b), max(b)))), blocks)
  if (length(todo) == 0) {
    cat("[skip] cell", k, cel$cell, "-- all blocks present\n"); next
  }

  ## The cell's prior sets, from the reference dataset seeded 6e5 + k, which no
  ## analysis iteration uses. Fixed within the cell so the Stan programs are
  ## compile-cache constants (Trap 5).
  ref <- simulate_dataset(truth, design, cv_control = cal$cv_control,
                          sigma_ratio = cal$sigma_ratio, seed = 6e5 + k,
                          sigma_mode = "absolute", sigma_0_abs = SIGMA_0)
  specs <- p10_cell_priors(ref)
  cat("=== cell", k, cel$cell, "| models:",
      paste(sprintf("%s=%d", names(specs),
                    vapply(specs, function(z) length(z$models), integer(1))),
            collapse = " "), "\n")

  ## Warm the compile cache SERIALLY. Every model recompiles when the cell
  ## changes, because the prior constants are literals in the Stan program, so
  ## this is ~46 compiles per cell and roughly three hours. Eighteen workers each
  ## meeting an uncompiled model at once is eighteen concurrent 1.5 GB C++
  ## compiles, which is how the first Phase 5 sweep exhausted a 31 GB box with no
  ## R-level error.
  cat("  warming compile cache ... ")
  tw <- Sys.time()
  warm <- try(run_iteration(0L, truth, design, specs, x_grid), silent = TRUE)
  cat(sprintf("%.1f min%s\n", as.numeric(difftime(Sys.time(), tw, units = "mins")),
              if (inherits(warm, "try-error")) " [FAILED]" else ""))
  if (inherits(warm, "try-error")) {
    cat("[abort] cell", k, "-- warm-up failed:\n", as.character(warm), "\n")
    next
  }

  for (b in todo) {
    path <- file.path(OUT, sprintf("%s__%s_i%03d-%03d.rds", cel$cell, TAG,
                                   min(b), max(b)))
    t0 <- Sys.time()
    ## Seed the master immediately before forking, from the cell index and the
    ## block's first iteration, so a block re-run in a fresh process after an
    ## interrupt is identical to the one it replaces. `mclapply` derives each
    ## worker's RNG stream from the master's state and bnec()'s initial values
    ## come from it; Phase 8 measured inits alone moving an endpoint by up to 7%.
    ## With 46 models in play that matters more here, not less.
    set.seed(9e5 + 1000 * k + min(b))
    res <- parallel::mclapply(b, function(i)
             try(run_iteration(i, truth, design, specs, x_grid), silent = TRUE),
           mc.cores = WORKERS, mc.preschedule = FALSE)
    bad <- vapply(res, function(z) inherits(z, "try-error"), logical(1))
    if (any(bad)) {
      msgs <- unique(vapply(res[bad], as.character, character(1)))
      cat(sprintf("  %d/%d worker failures; distinct messages:\n",
                  sum(bad), length(b)))
      for (m in utils::head(msgs, 5)) cat("    ", trimws(m), "\n")
    }
    good <- res[!bad]
    if (!length(good)) {
      cat("[abort]", basename(path), "-- every iteration failed; not written\n")
      next
    }
    stamp <- function(d) {
      if (is.null(d)) return(NULL)
      cbind(cell = cel$cell, cell_index = k, delta = cel$delta,
            top_factor = cel$top_factor, R = cel$R,
            sigma_ratio = cel$sigma_ratio, d, row.names = NULL)
    }
    ends <- stamp(do.call(rbind, lapply(good, `[[`, "endpoints")))
    ends <- merge(ends, true_endpoints(truth), by = "endpoint", all.x = TRUE)
    out <- list(endpoints = ends,
                weights = stamp(do.call(rbind, lapply(good, `[[`, "weights"))),
                curves  = stamp(do.call(rbind, lapply(good, `[[`, "curves"))),
                rhat    = stamp(do.call(rbind, lapply(good, `[[`, "rhat"))),
                x_grid  = x_grid)
    tmp <- paste0(path, ".tmp"); saveRDS(out, tmp); file.rename(tmp, path)
    cat(sprintf("[done] %s  %.1f min  %d/%d iterations  %d worker failures\n",
                basename(path), as.numeric(difftime(Sys.time(), t0, units = "mins")),
                length(unique(ends$iteration)), length(b), sum(bad)))
  }
}
cat("\nPhase 10 complete at", format(Sys.time()), "\n")
