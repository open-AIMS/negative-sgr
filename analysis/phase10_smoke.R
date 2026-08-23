## Phase 10 smoke test -- exercise every code path with the cheapest possible
## fits, before Gate 2 spends real time.
##
## The pilot fits 46 models per iteration and takes hours. This fits TWO models
## per arm with short chains, which is enough to prove the plumbing: a model set
## reaching `bnec()` as a `c(...)` formula, a `bayesmanecfit` coming back,
## `endpoint_table()` / `fit_diagnostics()` / `p10_weights()` / `p10_curve()`
## all reading it, and -- the two paths that had never run to completion -- a
## `cens()` aterm with a model set (arm C) and a `constant(0)` prior inside a
## named prior list (arm B3).
##
## The priors are cell 12's real ones, so the Stan programs compiled here are
## cache hits for the pilot rather than throwaway work.
##
## Run: Rscript analysis/phase10_smoke.R

source("R/setup.R"); source("R/setup_phase10.R")
source("R/data_prep.R"); source("R/diagnostics.R"); source("R/simulate.R")
load_bayesnec_p10(); source("R/arms.R"); source("R/metrics.R")
source("R/model_average.R")
p10_assert_arms_ready()
use_compile_cache()
options(mc.cores = 2L)

## Deliberately not P10_MCMC: this run is about whether the code executes, not
## about what it estimates. Nothing from it is reported.
SMOKE_MCMC <- list(chains = 2L, iter = 600L, warmup = 300L, adapt_delta = 0.9,
                   max_treedepth = 10L, seed = 1L)
PAIR <- list(gauss = c("nec4param", "ecxwb1"), zb = c("nec3param", "ecxexp"))

cal <- calibrate_sigma(read_sgr("c_proliferum"))
SIGMA_0 <- sigma_0_at(cal$cv_control, R_ref = 2.3, t = 7)
cells <- p10_cells(cal$sigma_ratio)
cel <- cells[cells$cell_index == 12L, ]
truth <- sim_truth(R = cel$R, delta = cel$delta, t = 7)
design <- sim_design(truth, top_factor = cel$top_factor, n_conc = 12,
                     n_rep = 5, n_control = 10)
x_grid <- seq(0, max(design$x), length.out = P10_GRID_N)

ref <- simulate_dataset(truth, design, cv_control = cal$cv_control,
                        sigma_ratio = cal$sigma_ratio, seed = 6e5 + 12,
                        sigma_mode = "absolute", sigma_0_abs = SIGMA_0)
sim <- simulate_dataset(truth, design, cv_control = cal$cv_control,
                        sigma_ratio = cal$sigma_ratio, seed = 7e5 + 1,
                        sigma_mode = "absolute", sigma_0_abs = SIGMA_0)
full <- p10_cell_priors(ref)

ok <- TRUE
for (arm in P10_ARMS) {
  want <- if (arm %in% c("E", "F")) PAIR$zb else PAIR$gauss
  mods <- intersect(full[[arm]]$models, want)
  spec <- list(models = mods, prior = full[[arm]]$prior[mods])
  cat("\n=== ", arm, " | ", paste(mods, collapse = " + "), " ===\n", sep = "")
  t0 <- Sys.time()
  f <- try(suppressMessages(p10_fit(arm, sim, spec, mcmc = SMOKE_MCMC)),
           silent = TRUE)
  if (inherits(f, "try-error")) {
    cat("  [FAIL] ", trimws(as.character(f)), "\n"); ok <- FALSE; next
  }
  cat("  class:", paste(class(f$fit), collapse = "/"),
      sprintf(" (%.1f min)\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  w <- try(p10_weights(f, spec$models), silent = TRUE)
  et <- try(endpoint_table(f$fit, arm, "sim"), silent = TRUE)
  dg <- try(fit_diagnostics(f$fit), silent = TRUE)
  cv <- try(p10_curve(f, x_grid, scale = p10_zb_scale(sim, arm)), silent = TRUE)
  ## The convergence screen. Short chains make this fire often, which is exactly
  ## what a smoke test wants: all three branches -- nothing dropped, some
  ## dropped, everything dropped -- get exercised cheaply here rather than
  ## discovered mid-sweep.
  dr <- try(p10_drop_nonconverged(f, P10_RHAT_CUTOFF), silent = TRUE)
  for (nm in c("w", "et", "dg", "cv", "dr")) {
    if (inherits(get(nm), "try-error")) {
      cat("  [FAIL]", nm, ":", trimws(as.character(get(nm))), "\n"); ok <- FALSE
    }
  }
  if (!inherits(w, "try-error")) {
    cat("  weights:", paste(sprintf("%s %.2f", w$model, w$wi), collapse = "  "),
        "| requested", w$n_requested[1], "fitted", w$n_fitted[1], "\n")
  }
  if (!inherits(et, "try-error")) {
    cat("  endpoints:",
        paste(sprintf("%s %.3f [%.3f, %.3f]", et$endpoint, et$estimate,
                      et$lower, et$upper), collapse = "  "), "\n")
    if (anyNA(et$estimate)) { cat("  [FAIL] NA endpoint\n"); ok <- FALSE }
  }
  if (!inherits(dr, "try-error")) {
    cat(sprintf("  rhat: %s | kept %d, dropped %d%s\n",
                paste(sprintf("%s %.3f", dr$rhat$model, dr$rhat$max_rhat),
                      collapse = "  "),
                length(dr$kept), length(dr$dropped),
                if (is.null(dr$fit)) "  -> NO SURVIVOR (recorded as arm failure)"
                else paste0("  -> ", paste(class(dr$fit)[1]))))
    ## Short chains usually put BOTH models over 1.01, which exercises the
    ## all-dropped branch but never reaches `amend()`. So the amend path is
    ## forced separately, at a cutoff placed between the two models' R-hats so
    ## that exactly one is dropped. Endpoints computed on the amended fit are
    ## the thing most likely to break and the thing the sweep depends on.
    forced <- dr
    if (is.null(dr$fit) && nrow(dr$rhat) > 1) {
      cut2 <- mean(sort(dr$rhat$max_rhat)[1:2])
      forced <- try(p10_drop_nonconverged(f, cut2), silent = TRUE)
      if (!inherits(forced, "try-error") && !is.null(forced$fit)) {
        cat(sprintf("  forced amend at cutoff %.4f: kept %s, dropped %s -> %s\n",
                    cut2, paste(forced$kept, collapse = " "),
                    paste(forced$dropped, collapse = " "),
                    class(forced$fit)[1]))
      } else {
        cat("  [FAIL] could not force an amend to exercise the drop path\n")
        ok <- FALSE; forced <- dr
      }
    }
    if (!is.null(forced$fit)) {
      et2 <- try(endpoint_table(p10_pull(forced$fit), arm, "sim"), silent = TRUE)
      w2 <- try(p10_weights(forced$fit, spec$models), silent = TRUE)
      if (inherits(et2, "try-error") || inherits(w2, "try-error")) {
        cat("  [FAIL] reading the amended fit:",
            trimws(as.character(if (inherits(et2, "try-error")) et2 else w2)), "\n")
        ok <- FALSE
      } else {
        cat("  after drop:",
            paste(sprintf("%s %.3f", et2$endpoint, et2$estimate),
                  collapse = "  "),
            "| weights", paste(sprintf("%s %.2f", w2$model, w2$wi),
                               collapse = " "), "\n")
        if (anyNA(et2$estimate)) {
          cat("  [FAIL] NA endpoint after drop\n"); ok <- FALSE
        }
      }
    }
  }
  if (!inherits(cv, "try-error")) {
    cat(sprintf("  curve: mu in [%.4f, %.4f] over %d grid points\n",
                min(cv$mu), max(cv$mu), nrow(cv)))
    ## B3 pins the lower asymptote at zero; a curve that dips below it means the
    ## constant(0) prior did not take, which is silent and would make B3 a
    ## different arm from the one described.
    if (arm == "B3" && min(cv$mu) < -1e-8) {
      cat("  [FAIL] B3 curve below zero\n"); ok <- FALSE
    }
  }
}
cat("\n", if (ok) "SMOKE PASSED -- run Gate 2.\n" else
             "SMOKE FAILED -- fix before Gate 2.\n", sep = "")
if (!ok) quit(status = 1)
