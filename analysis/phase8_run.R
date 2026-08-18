## Phase 8, stage 2 -- top every arm up to 500 iterations.
##
## Iterations 241-500 for all EIGHT arms in all twelve scenarios:
## 12 x 260 x 8 = 24,960 fits, roughly 3.5 days at 18 workers. This buys Monte
## Carlo precision and nothing else -- the differences the study reports run from
## several percentage points to the gap between 0.95 and 0.00, so no ordering can
## move -- which is why it is sequenced last and why the vignette did not wait.
##
## WHY THIS SCRIPT EXISTS RATHER THAN A FLAG ON THE OTHERS. The six Gaussian arms
## were produced by phase5_run.R and the two family-floored arms by phase7_run.R,
## and both are the provenance of artefacts already reported. Adding an iteration
## range to them would edit the scripts that produced those artefacts. Instead
## the two iteration bodies are reproduced here unchanged, and
## `PHASE8_EQUIVALENCE=1` re-runs iterations 1-2 into a scratch directory so the
## reproduction can be checked against the sweep row by row before 3.5 days of
## compute are committed to it. It matches bit-for-bit; see the log.
##
## Writes are chunked and additive, exactly as in phase7_run.R: each block of
## CHUNK iterations goes to its own file under a temporary name and is renamed,
## so a file that exists is complete, a resume skips it by filename, and an
## interrupt costs at most one block. Nothing already on disk is ever rewritten.
##
## Run:    N_ITER_FROM=241 N_ITER_TO=500 CHUNK=40 WORKERS=18 Rscript analysis/phase8_run.R
## Resume: re-issue the identical command.
## Check:  PHASE8_EQUIVALENCE=1 Rscript analysis/phase8_run.R

source("R/setup.R"); source("R/data_prep.R"); source("R/diagnostics.R")
source("R/simulate.R"); load_bayesnec(); source("R/arms.R"); source("R/metrics.R")
use_compile_cache()
options(mc.cores = 1L)          # one whole iteration per worker, chains serial

EQUIV     <- identical(Sys.getenv("PHASE8_EQUIVALENCE", "0"), "1")
ITER_FROM <- as.integer(Sys.getenv("N_ITER_FROM", if (EQUIV) "1" else "241"))
ITER_TO   <- as.integer(Sys.getenv("N_ITER_TO",   if (EQUIV) "2" else "500"))
CHUNK     <- as.integer(Sys.getenv("CHUNK", if (EQUIV) "2" else "40"))
WORKERS   <- as.integer(Sys.getenv("WORKERS", as.character(STUDY_CORES)))
OUT       <- Sys.getenv("OUT", if (EQUIV) "analysis/phase8_equivalence" else
                               "analysis/phase5")
TAG       <- Sys.getenv("TAG", if (EQUIV) "equiv" else "topup")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

GAUSS_ARMS <- c("A", "B1", "B2", "B3", "C", "D")
ZB_ARMS    <- c("E", "F")
MCMC_SIM <- list(chains = 4L, iter = 2000L, warmup = 1000L, adapt_delta = 0.95,
                 max_treedepth = 12L, seed = 20260812L)

cal <- calibrate_sigma(read_sgr("c_proliferum"))
SIGMA_0 <- sigma_0_at(cal$cv_control, R_ref = 2.3, t = 7)

## Cell order is load-bearing, not cosmetic: the reference prior for cell k is
## drawn with seed 6e5 + k, so a different ordering would give the top-up
## iterations a different prior from iterations 1-240 in the same cell. This is
## the `rsep` grid of phase5_run.R and the grid of phase7_run.R, which are the
## same grid in the same order.
cells <- rbind(
  expand.grid(delta = c(2, 4, 8), top_factor = c(0.8, 1.0, 2.0), R = 2.3,
              sigma_ratio = cal$sigma_ratio, KEEP.OUT.ATTRS = FALSE),
  expand.grid(delta = 4, top_factor = 2.0, R = c(3.3, 17, 73),
              sigma_ratio = cal$sigma_ratio, KEEP.OUT.ATTRS = FALSE))
cells$cell <- sprintf("d%.1f_t%.1f_R%.1f_s%.1f", cells$delta, cells$top_factor,
                      cells$R, cells$sigma_ratio)

## Restrict to particular cells by index, e.g. CELL_INDEX=7,8. Used by the
## equivalence check to exercise a REACHING cell, where arm D is refitted at its
## own arm-A crossing, as well as a stops-short cell, where the crossing is not
## reached and arm D is reported as arm A. Those are different code paths and
## checking only one of them would check half the runner.
CELL_INDEX <- Sys.getenv("CELL_INDEX", "")
if (nzchar(CELL_INDEX)) {
  cells <- cells[as.integer(strsplit(CELL_INDEX, ",")[[1]]), , drop = FALSE]
}

iters <- ITER_FROM:ITER_TO
cat("Phase 8 stage 2 |", nrow(cells), "cells x", length(iters), "iterations x",
    length(c(GAUSS_ARMS, ZB_ARMS)), "arms =",
    nrow(cells) * length(iters) * 8, "fits\n")
cat("iterations", ITER_FROM, "-", ITER_TO, "| WORKERS =", WORKERS,
    "| CHUNK =", CHUNK, "| OUT =", OUT, "| TAG =", TAG, "\n\n")

## --------------------------------------------- the two arms sets, unchanged --
prep_zb <- function(sim, arm) {
  y0 <- pmax(sim$sgr, 0)
  if (identical(arm, "E")) data.frame(x = sim$x, y = y0 / max(y0))
  else                     data.frame(x = sim$x, y = y0)
}
fam_zb <- function(arm) {
  if (identical(arm, "E")) brms::Beta(link = "identity") else
    stats::Gamma(link = "identity")
}
fit_zb <- function(arm, sim, prior, mcmc) {
  bayesnec::bnec(bayesnec::bnf(y ~ crf(x, "nec3param")), data = prep_zb(sim, arm),
                 family = fam_zb(arm), prior = prior, chains = mcmc$chains,
                 iter = mcmc$iter, warmup = mcmc$warmup, seed = mcmc$seed,
                 backend = "cmdstanr",
                 control = list(adapt_delta = mcmc$adapt_delta,
                                max_treedepth = mcmc$max_treedepth))
}

run_iteration <- function(i, truth, design, prior_ref, priors_zb) {
  sim <- simulate_dataset(truth, design, cv_control = cal$cv_control,
                          sigma_ratio = cal$sigma_ratio, seed = 7e5 + i,
                          sigma_mode = "absolute", sigma_0_abs = SIGMA_0)
  ## `< 0`, matching phase5_run.R. phase7_run.R wrote `<= 0`; on a continuous
  ## response the two differ only on an event of probability zero, and this file
  ## has to carry ONE f_neg per iteration for all eight arms.
  f_neg <- mean(sim$sgr < 0)
  mcmc <- MCMC_SIM; mcmc$seed <- MCMC_SIM$seed + i
  pr <- prior_ref
  out <- list()

  ## -- the six Gaussian arms, as phase5_run.R fits them ---------------------
  fitA <- NULL
  for (arm in c("A", setdiff(GAUSS_ARMS, c("A", "D")))) {
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
  ## Arm D truncates at THIS iteration's arm-A crossing, never the true one.
  ## `zero_crossing()` returns Inf where the fitted curve never reaches zero
  ## inside the tested range, and arm D is then identical to arm A and is
  ## reported as such rather than refitted. Testing only for "try-error" here is
  ## what destroyed the first sweep: a bare Inf is not a try-error, and the
  ## fall-through killed every worker in the top_factor = 0.8 cells.
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

  ## -- the two family-floored arms, as phase7_run.R fits them ---------------
  for (arm in ZB_ARMS) {
    f <- try(suppressMessages(fit_zb(arm, sim, priors_zb[[arm]], mcmc)),
             silent = TRUE)
    if (inherits(f, "try-error")) {
      out[[arm]] <- data.frame(arm = arm, endpoint = NA, estimate = NA,
                               lower = NA, upper = NA, divergences = NA,
                               max_rhat = NA, ok = FALSE,
                               error = sub("\n.*", "", as.character(f)))
      next
    }
    dgn <- fit_diagnostics(f)
    et <- endpoint_table(f, arm, "sim")
    out[[arm]] <- data.frame(arm = arm, endpoint = et$endpoint,
                             estimate = et$estimate, lower = et$lower,
                             upper = et$upper, divergences = dgn$divergences,
                             max_rhat = dgn$max_rhat, ok = TRUE,
                             error = NA_character_)
  }
  cbind(iteration = i, f_neg = f_neg, do.call(rbind, out), row.names = NULL)
}

## ------------------------------------------------------------------ sweep ----
for (k in seq_len(nrow(cells))) {
  cel <- cells[k, ]
  truth <- sim_truth(R = cel$R, delta = cel$delta, t = 7)
  design <- sim_design(truth, top_factor = cel$top_factor, n_conc = 12,
                       n_rep = 5, n_control = 10)

  blocks <- split(iters, ceiling(seq_along(iters) / CHUNK))
  todo <- Filter(function(b) !file.exists(file.path(
            OUT, sprintf("%s__%s_i%03d-%03d.rds", cel$cell, TAG,
                         min(b), max(b)))), blocks)
  if (length(todo) == 0) { cat("[skip]", cel$cell, "-- all blocks present\n"); next }

  ## The cell's reference priors, drawn with seeds no analysis iteration uses,
  ## and identical to those used for iterations 1-240 in this cell.
  ref <- simulate_dataset(truth, design, cv_control = cal$cv_control,
                          sigma_ratio = cal$sigma_ratio, seed = 6e5 + k,
                          sigma_mode = "absolute", sigma_0_abs = SIGMA_0)
  ref_prep <- prepare_sgr(ref, "raw", meta = sim_meta())
  prior_ref <- arm_prior(ref_prep$x, ref_prep$y)
  priors_zb <- lapply(setNames(ZB_ARMS, ZB_ARMS), function(a) {
    d <- prep_zb(ref, a)
    bayesnec:::define_prior("nec3param", fam_zb(a), d$x, d$y, "uninformative")
  })

  ## Warm the compile cache SERIALLY. Eighteen workers each meeting an
  ## uncompiled model at once is eighteen concurrent 1.5 GB C++ compiles, which
  ## is how the first sweep exhausted a 31 GB box with no R-level error.
  cat("  warming compile cache for", cel$cell, "... ")
  tw <- Sys.time()
  warm <- try(run_iteration(0L, truth, design, prior_ref, priors_zb),
              silent = TRUE)
  cat(sprintf("%.1f min%s\n", as.numeric(difftime(Sys.time(), tw, units = "mins")),
              if (inherits(warm, "try-error")) " [FAILED]" else ""))
  if (inherits(warm, "try-error")) {
    cat("[abort]", cel$cell, "-- warm-up failed:\n", as.character(warm), "\n")
    next
  }

  for (b in todo) {
    path <- file.path(OUT, sprintf("%s__%s_i%03d-%03d.rds", cel$cell, TAG,
                                   min(b), max(b)))
    t0 <- Sys.time()
    ## Seed the master immediately before forking, from the cell index and the
    ## block's first iteration.
    ##
    ## `mclapply` derives each worker's RNG stream from the master's state, and
    ## the initial values `bnec()` draws come from that stream. Without this the
    ## inits a block gets depend on how much RNG the process happened to consume
    ## beforehand -- which differs between a block run mid-sweep and the same
    ## block re-run in a fresh process after an interrupt. Over a 3.5-day run
    ## that is not hypothetical. Fixing the seed per block makes a resumed block
    ## identical to the one it replaces, which is what "resume" has to mean.
    ##
    ## Measured while verifying this runner: initial values alone move an
    ## endpoint by up to 7% on arm B1 and change B2's divergence count, so this
    ## is a visible source of variation, not a formality. It is sampler noise
    ## rather than bias -- the inits are drawn from the same priors by the same
    ## code -- but it should be pinned rather than left to chance.
    set.seed(9e5 + 1000 * k + min(b))
    res <- parallel::mclapply(b, function(i)
             try(run_iteration(i, truth, design, prior_ref, priors_zb),
                 silent = TRUE),
           mc.cores = WORKERS, mc.preschedule = FALSE)
    bad <- vapply(res, function(z) inherits(z, "try-error"), logical(1))
    if (any(bad)) {
      msgs <- unique(vapply(res[bad], as.character, character(1)))
      cat(sprintf("  %d/%d worker failures in %s; distinct messages:\n",
                  sum(bad), length(b), cel$cell))
      for (m in utils::head(msgs, 5)) cat("    ", trimws(m), "\n")
    }
    out <- do.call(rbind, res[!bad])
    if (is.null(out) || nrow(out) == 0) {
      cat("[abort]", basename(path), "-- every iteration failed; not written\n")
      next
    }
    out <- cbind(endpoint = out$endpoint, delta = cel$delta,
                 top_factor = cel$top_factor, R = cel$R,
                 sigma_ratio = cel$sigma_ratio, cell = cel$cell,
                 out[, setdiff(names(out), "endpoint")])
    out <- merge(out, true_endpoints(truth), by = "endpoint", all.x = TRUE)
    tmp <- paste0(path, ".tmp"); saveRDS(out, tmp); file.rename(tmp, path)
    cat(sprintf("[done] %s  %.1f min  %d/%d iterations  %d worker failures\n",
                basename(path), as.numeric(difftime(Sys.time(), t0, units="mins")),
                length(unique(out$iteration)), length(b), sum(bad)))
  }
}
cat("\nPhase 8 stage 2 complete at", format(Sys.time()), "\n")
