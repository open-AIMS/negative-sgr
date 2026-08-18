## Phase 7, stage 1 -- arms E and F: the families that floor implicitly.
##
## E: floor negatives at zero, divide by the observed maximum, Beta, nec3param.
## F: floor negatives at zero, no scaling, Gamma, nec3param.
##
## Both take bnec()'s own default priors (that being the practice under
## examination) and let check_data() perform the boundary nudge. Negatives are
## floored explicitly because both families reject a negative outright and there
## would otherwise be no fit to score.
##
## WRITES ARE CHUNKED, not per-cell. A cell of 240 iterations takes over an hour
## and the Phase 5 runner only wrote at cell completion, so an interrupt cost the
## whole cell. Here each block of CHUNK iterations is written as its own file the
## moment it finishes, so an interrupt costs at most one block (~1-2 min).
## Nothing already on disk is ever rewritten: the twelve verified Phase 5 cell
## files are not touched, and a resumed run skips blocks whose file exists.
##
## Run:  N_ITER=240 CHUNK=40 WORKERS=18 Rscript analysis/phase7_run.R
## Resume: re-issue the identical command. Completed blocks are skipped.

source("R/setup.R"); source("R/data_prep.R"); source("R/diagnostics.R")
source("R/simulate.R"); load_bayesnec(); source("R/arms.R"); source("R/metrics.R")
use_compile_cache()
options(mc.cores = 1L)          # one whole iteration per worker, chains serial

N_ITER  <- as.integer(Sys.getenv("N_ITER", "240"))
CHUNK   <- as.integer(Sys.getenv("CHUNK", "40"))
WORKERS <- as.integer(Sys.getenv("WORKERS", as.character(STUDY_CORES)))
OUT     <- "analysis/phase5"
ARMS    <- c("E", "F")
MCMC_SIM <- list(chains = 4L, iter = 2000L, warmup = 1000L, adapt_delta = 0.95,
                 max_treedepth = 12L, seed = 20260812L)

cal <- calibrate_sigma(read_sgr("c_proliferum"))
SIGMA_0 <- sigma_0_at(cal$cv_control, R_ref = 2.3, t = 7)

cells <- rbind(
  expand.grid(delta = c(2, 4, 8), top_factor = c(0.8, 1.0, 2.0), R = 2.3,
              sigma_ratio = cal$sigma_ratio, KEEP.OUT.ATTRS = FALSE),
  expand.grid(delta = 4, top_factor = 2.0, R = c(3.3, 17, 73),
              sigma_ratio = cal$sigma_ratio, KEEP.OUT.ATTRS = FALSE))
cells$cell <- sprintf("d%.1f_t%.1f_R%.1f_s%.1f", cells$delta, cells$top_factor,
                      cells$R, cells$sigma_ratio)

cat("Phase 7 stage 1 |", nrow(cells), "cells x", N_ITER, "iterations x",
    length(ARMS), "arms =", nrow(cells) * N_ITER * length(ARMS), "fits\n")
cat("WORKERS =", WORKERS, "| CHUNK =", CHUNK, "iterations per file\n\n")

#' Data preparation for the two zero-bounded arms
#'
#' Flooring is explicit; the boundary nudge (exact zeros, exact ones) is left to
#' `check_data()` because that adjustment is part of the convention under test.
prep_zb <- function(sim, arm) {
  y0 <- pmax(sim$sgr, 0)
  if (identical(arm, "E")) data.frame(x = sim$x, y = y0 / max(y0))
  else                     data.frame(x = sim$x, y = y0)
}
fam_zb <- function(arm) {
  if (identical(arm, "E")) brms::Beta(link = "identity")
  else stats::Gamma(link = "identity")
}

fit_zb <- function(arm, sim, prior, mcmc = MCMC_SIM) {
  dat <- prep_zb(sim, arm)
  bayesnec::bnec(bayesnec::bnf(y ~ crf(x, "nec3param")), data = dat,
                 family = fam_zb(arm), prior = prior,
                 chains = mcmc$chains, iter = mcmc$iter, warmup = mcmc$warmup,
                 seed = mcmc$seed, backend = "cmdstanr",
                 control = list(adapt_delta = mcmc$adapt_delta,
                                max_treedepth = mcmc$max_treedepth))
}

run_iteration <- function(i, truth, design, priors) {
  sim <- simulate_dataset(truth, design, cv_control = cal$cv_control,
                          sigma_ratio = cal$sigma_ratio, seed = 7e5 + i,
                          sigma_mode = "absolute", sigma_0_abs = SIGMA_0)
  mu <- nec4param_curve(design$x, truth$top, truth$bot, truth$beta, truth$nec)
  f_neg <- mean(sim$sgr <= 0)
  truths <- c(ErC10 = true_ecx(truth, 10), ErC50 = true_ecx(truth, 50),
              NSEC = truth$nec)
  do.call(rbind, lapply(ARMS, function(a) {
    mc <- MCMC_SIM; mc$seed <- MCMC_SIM$seed + i
    f <- try(suppressMessages(fit_zb(a, sim, priors[[a]], mc)), silent = TRUE)
    if (inherits(f, "try-error")) {
      return(data.frame(endpoint = names(truths), iteration = i, f_neg = f_neg,
                        arm = a, estimate = NA_real_, lower = NA_real_,
                        upper = NA_real_, divergences = NA_real_,
                        max_rhat = NA_real_, ok = FALSE,
                        error = sub("\n.*", "", as.character(f)),
                        truth = unname(truths), stringsAsFactors = FALSE))
    }
    et <- endpoint_table(f, a, "sim"); dg <- fit_diagnostics(f)
    data.frame(endpoint = et$endpoint, iteration = i, f_neg = f_neg, arm = a,
               estimate = et$estimate, lower = et$lower, upper = et$upper,
               divergences = dg$divergences, max_rhat = dg$max_rhat, ok = TRUE,
               error = NA_character_, truth = unname(truths[et$endpoint]),
               stringsAsFactors = FALSE)
  }))
}

for (k in seq_len(nrow(cells))) {
  cel <- cells[k, ]
  truth <- sim_truth(R = cel$R, delta = cel$delta, t = 7)
  design <- sim_design(truth, top_factor = cel$top_factor)
  blocks <- split(seq_len(N_ITER), ceiling(seq_len(N_ITER) / CHUNK))
  todo <- Filter(function(b) !file.exists(file.path(
            OUT, sprintf("%s__ef_i%03d-%03d.rds", cel$cell, min(b), max(b)))),
          blocks)
  if (length(todo) == 0) { cat("[skip]", cel$cell, "-- all blocks present\n"); next }

  ## Reference prior, built once per cell from iteration 1's data and reused for
  ## every iteration in it. Rebuilding per dataset would write different
  ## constants into the Stan program and force a recompile on every fit (Trap 5).
  ref <- simulate_dataset(truth, design, cv_control = cal$cv_control,
                          sigma_ratio = cal$sigma_ratio, seed = 6e5 + k,
                          sigma_mode = "absolute", sigma_0_abs = SIGMA_0)
  priors <- lapply(setNames(ARMS, ARMS), function(a) {
    d <- prep_zb(ref, a)
    bayesnec:::define_prior("nec3param", fam_zb(a), d$x, d$y, "uninformative")
  })

  ## Warm the compile cache SERIALLY before fanning out. Eighteen workers each
  ## launching cc1plus on the same new model is what pinned this machine at
  ## 29 of 31 GB during Phase 5.
  cat("  warming compile cache for", cel$cell, "...")
  for (a in ARMS) invisible(try(suppressMessages(
        fit_zb(a, ref, priors[[a]],
               list(chains = 1L, iter = 200L, warmup = 100L,
                    adapt_delta = 0.8, max_treedepth = 10L,
                    seed = MCMC_SIM$seed))), silent = TRUE))
  cat(" done\n")

  for (b in todo) {
    path <- file.path(OUT, sprintf("%s__ef_i%03d-%03d.rds", cel$cell,
                                   min(b), max(b)))
    t0 <- Sys.time()
    res <- parallel::mclapply(b, function(i)
             try(run_iteration(i, truth, design, priors), silent = TRUE),
           mc.cores = WORKERS)
    bad <- vapply(res, function(z) inherits(z, "try-error"), logical(1))
    out <- do.call(rbind, res[!bad])
    if (is.null(out) || nrow(out) == 0) {
      cat("[abort]", basename(path), "-- every iteration failed; not written\n")
      next
    }
    out <- cbind(endpoint = out$endpoint, delta = cel$delta,
                 top_factor = cel$top_factor, R = cel$R,
                 sigma_ratio = cel$sigma_ratio, cell = cel$cell,
                 out[, setdiff(names(out), "endpoint")])
    ## Written to a temporary name and renamed, so a file that exists is
    ## always complete even if the process dies mid-write.
    tmp <- paste0(path, ".tmp"); saveRDS(out, tmp); file.rename(tmp, path)
    cat(sprintf("[done] %s  %.1f min  %d/%d iterations  %d worker failures\n",
                basename(path), as.numeric(difftime(Sys.time(), t0, units="mins")),
                length(unique(out$iteration)), length(b), sum(bad)))
  }
}
cat("\nPhase 7 stage 1 complete at", format(Sys.time()), "\n")
