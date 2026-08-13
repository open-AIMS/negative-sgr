## Session constants. Sourced by every entry point.
##
## bayesnec is loaded from a pinned worktree rather than from a library, so the
## commit under test is unambiguous and cannot drift between phases. The SHA is
## asserted at load time: if the worktree has moved, the run stops rather than
## silently producing results from a different branch state.

## `dev` rather than `cens-impl`: dev CONTAINS cens-impl (merged as PR #189) and
## adds three things this study needs -- the Censoring section in example1
## (commit ee76e815, which post-dates cens-impl and is the Phase 6 target), the
## normalisation detection in check_data.R, and the matching bnec() guidance.
## Pinning to cens-impl would have run the whole study against a superseded
## commit and written the vignette section into the wrong file.
##
## The pin is a DETACHED worktree at that commit, not the `dev` branch itself.
## `dev` now lives in the main working directory, where it will move as work
## continues; a detached checkout cannot drift, so the study keeps its own frozen
## copy. The SHA assertion below is what actually guarantees which code loads.
## (The directory is named for an earlier issue -- rename it with
## `git worktree move` once no run is in flight.)
BAYESNEC_WORKTREE <- "/mnt/c/Rworking/bayesnec-issue173"
BAYESNEC_BRANCH   <- "dev"
BAYESNEC_SHA      <- "374e511c665fee04ca6ca8e7b48a547a3928a28d"

## 22 cores. The competing study that had claimed about half of them has
## finished, so the cap is raised to 18, leaving 4 for the interactive Positron
## sessions. Note that the Phase 5 pilot TIMINGS were measured at 6 workers --
## per-fit costs are per-core and should carry over, but wall-clock projections
## made against the pilot must be rescaled.
STUDY_CORES <- 18L

## MCMC settings used everywhere except the pilot timing run.
MCMC <- list(chains = 4L, iter = 4000L, warmup = 2000L, adapt_delta = 0.99,
             max_treedepth = 12L, seed = 20260812L)

## Persistent Stan compile cache.
##
## cmdstanr writes each model's .stan file to a temporary directory and rebuilds
## the executable whenever it is missing, so a fresh R session recompiles every
## model from scratch -- 3-5 minutes each on this machine, which would dominate
## a simulation of thousands of fits. Pointing `cmdstanr_write_stan_file_dir` at
## a directory inside the project makes the file name a hash of the Stan code,
## so an unchanged model reuses its existing executable across sessions.
##
## The study uses only a handful of distinct Stan programs (arms A/B1/D share
## one, B2/B3 share one, C one, SQ one), so after the first run compilation
## disappears from the budget entirely.
CMDSTAN_CACHE <- file.path(getwd(), "cmdstan_cache")

use_compile_cache <- function() {
  dir.create(CMDSTAN_CACHE, showWarnings = FALSE, recursive = TRUE)
  options(cmdstanr_write_stan_file_dir = CMDSTAN_CACHE)
  invisible(CMDSTAN_CACHE)
}

load_bayesnec <- function(check_sha = TRUE) {
  if (check_sha) {
    sha <- system2("git", c("-C", BAYESNEC_WORKTREE, "rev-parse", "HEAD"),
                   stdout = TRUE)
    if (!identical(sha, BAYESNEC_SHA)) {
      stop("bayesnec worktree is at ", sha, " but this study is pinned to ",
           BAYESNEC_SHA, ". Check out the pinned commit or update SESSION.md ",
           "deliberately.")
    }
  }
  suppressMessages(pkgload::load_all(BAYESNEC_WORKTREE, quiet = TRUE,
                                     export_all = FALSE))
  invisible(TRUE)
}

#' Record a Phase 1 gate result
#'
#' Gate outcomes go to a machine-readable file as well as SESSION.md so that a
#' later phase can refuse to run against a failed gate.
gate_result <- function(id, description, passed, detail = "") {
  data.frame(id = id, description = description, passed = passed,
             detail = detail, stringsAsFactors = FALSE)
}
