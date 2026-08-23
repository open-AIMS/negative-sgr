## Phase 10 helpers -- the simulation under bayesnec's model-averaged workflow.
##
## Additive, like R/setup_phase10.R: Phase 9 was in flight when this was written
## and owns R/arms.R and R/metrics.R, so nothing here edits them. Phase 9 defines
## its own `arm_prior_set()` and `weights_table()` inside
## `analysis/phase9_modelavg.R`; the versions here are named `p10_*` and are
## deliberately separate rather than shared, because the two phases differ in
## what a "prior set" means (Phase 9 builds one from the real data's arm A;
## Phase 10 builds one per cell from a reference simulated dataset, so that the
## prior is a compile-cache constant -- see Trap 5).

#' Assert the Phase 9 changes to R/arms.R and R/metrics.R are present
#'
#' Phase 10 depends on four things that were uncommitted working-tree changes
#' when it was written: `fit_arm()` accepting a model SET, `fix_bot_prior()`
#' recursing over a named list of priors, `bot_bearing_models()`, and
#' `fit_diagnostics()` handling a `bayesmanecfit`. (`zero_crossing()`'s
#' `bayesmanecfit` branch is a Phase 9 need, not a Phase 10 one: arm D is not in
#' this phase, so nothing here calls it.) If any is missing the sweep would not
#' error --
#' it would fit single models and quietly answer a different question, over
#' several days. Hence an assertion rather than a comment.
p10_assert_arms_ready <- function(env = parent.frame()) {
  missing <- character(0)
  ## Absence and wrong-version are separate failures and are reported
  ## separately: `formals()` on a name that does not exist raises "object not
  ## found", which points at this function rather than at the file that was
  ## never sourced.
  have <- function(nm) exists(nm, envir = env, mode = "function")
  for (nm in c("fit_arm", "fit_diagnostics", "bot_bearing_models")) {
    if (!have(nm)) missing <- c(missing, paste0(nm, "() [not found at all]"))
  }
  if (have("fit_arm")) {
    fa <- get("fit_arm", envir = env, mode = "function")
    if (!"model" %in% names(formals(fa))) {
      missing <- c(missing, "fit_arm(model=)")
    }
    ## `fit_arm()` must write a length-n model vector into the formula as
    ## c(...). Checked on the source rather than by fitting, which would cost
    ## minutes.
    if (!any(grepl("mod_arg", deparse(body(fa))))) {
      missing <- c(missing, "fit_arm() model-set formula construction (mod_arg)")
    }
  }
  if (have("fit_diagnostics")) {
    fd <- get("fit_diagnostics", envir = env, mode = "function")
    if (!any(grepl("bayesmanecfit", deparse(body(fd))))) {
      missing <- c(missing, "fit_diagnostics() bayesmanecfit branch")
    }
  }
  if (length(missing)) {
    stop("Phase 10 needs the Phase 9 changes to R/arms.R and R/metrics.R, but ",
         "these are absent: ", paste(missing, collapse = ", "),
         ". Commit or restore them before running.")
  }
  invisible(TRUE)
}

#' The twelve-cell simulation grid, carrying its original index
#'
#' Reproduces `phase8_run.R`'s grid in the same order, because the cell index is
#' load-bearing: the reference-prior dataset of cell k is drawn with seed
#' `6e5 + k`. Phase 8 subsets this grid by CELL_INDEX and *then* numbers `k` from
#' 1, so asking it for cell 12 alone gives k = 1 and a different reference prior.
#' Phase 10 keeps `cell_index` as a column and never renumbers, which is what
#' makes its simulated datasets identical to the ones the single-model arms saw.
p10_cells <- function(sigma_ratio) {
  cells <- rbind(
    expand.grid(delta = c(2, 4, 8), top_factor = c(0.8, 1.0, 2.0), R = 2.3,
                sigma_ratio = sigma_ratio, KEEP.OUT.ATTRS = FALSE),
    expand.grid(delta = 4, top_factor = 2.0, R = c(3.3, 17, 73),
                sigma_ratio = sigma_ratio, KEEP.OUT.ATTRS = FALSE))
  cells$cell_index <- seq_len(nrow(cells))
  cells$cell <- sprintf("d%.1f_t%.1f_R%.1f_s%.1f", cells$delta, cells$top_factor,
                        cells$R, cells$sigma_ratio)
  cells
}

#' Data preparation for the two family-floored arms
#'
#' Identical to `prep_zb()` in `phase8_run.R`, reproduced rather than sourced so
#' that a later edit to Phase 8 cannot silently change Phase 10. E scales to the
#' observed maximum, which is the usual preparation for a Beta; F does not.
p10_zb_prep <- function(sim, arm) {
  y0 <- pmax(sim$sgr, 0)
  if (identical(arm, "E")) {
    data.frame(x = sim$x, y = y0 / max(y0))
  } else {
    data.frame(x = sim$x, y = y0)
  }
}

#' The scale E's response was divided by, so its curve can be put back
#'
#' `ecx()` and `nsec()` are invariant to this (both sides of
#' `f(x) = max(f) * (1 - x/100)` carry the scaling), which is why no endpoint
#' needs correcting. The stored CURVE does, or E's would be plotted on a
#' different y axis from everything else.
p10_zb_scale <- function(sim, arm) {
  if (identical(arm, "E")) max(pmax(sim$sgr, 0)) else 1
}

p10_zb_family <- function(arm) {
  # "Beta", not "beta" -- base R's beta() is the beta function (Trap 7).
  if (identical(arm, "E")) brms::Beta(link = "identity") else
    stats::Gamma(link = "identity")
}

#' The family each arm fits under
p10_family <- function(arm) {
  if (arm %in% c("E", "F")) p10_zb_family(arm) else
    stats::gaussian(link = "identity")
}

#' The candidate set an arm actually gets
#'
#' Derived by asking `check_models()` rather than by hard-coding a list, because
#' the whole point of the phase is what `bayesnec` gives an analyst. Two
#' consequences, both reported rather than corrected:
#'
#' Under `gaussian` the eleven `zero_bounded` models are dropped, leaving 8 of
#' the 14 declining models. Under `Beta` and `Gamma` only the linear models go,
#' leaving 12 -- so the family-floored arms average over a RICHER set, six
#' members of which asymptote at zero by construction, which is exactly the
#' shape a floored dataset wants. That asymmetry is not an artefact to correct:
#' it is what an analyst choosing a Beta actually gets. But "E did better under
#' averaging" may be a statement about `check_models()` rather than about
#' averaging, and the report must say which.
#'
#' Arm B3 additionally loses `neclin` and `ecxlin`, which decline linearly
#' without bound and have no asymptote for `constant(0)` to act on. A pinned arm
#' therefore differs from arm A by its candidate set as well as by its
#' constraint.
p10_candidate_models <- function(arm, sim, model_set = P10_MODEL_SET) {
  requested <- bayesnec:::expand_model_set(model_set)
  dat <- if (arm %in% c("E", "F")) {
    p10_zb_prep(sim, arm)
  } else {
    data.frame(x = sim$x, y = sim$sgr)
  }
  ## `check_models()` reads the predictor off a bayesnec model frame via
  ## `retrieve_var()`, which needs the `bnec_pop` attribute -- a plain
  ## data.frame silently yields NULL there. Build the frame the way `bnec()`
  ## does rather than hand-rolling the attribute, so this exercises the real
  ## path. The `data` argument only ever drops models for a NEGATIVE predictor,
  ## which cannot arise in this design, but passing it keeps the gate honest.
  stopifnot(min(dat$x) >= 0)
  bdat <- stats::model.frame(bayesnec::bnf(y ~ crf(x, "nec4param")), data = dat)
  mods <- suppressMessages(
    bayesnec:::check_models(requested, p10_family(arm), bdat))
  if (identical(arm, "B3")) mods <- bot_bearing_models(mods)
  mods
}

#' One prior per candidate model, built once per cell
#'
#' Trap 5: `bnec()` writes prior constants into the Stan program as literals, so
#' a per-dataset prior forces a recompile per fit. The prior is therefore built
#' from the cell's REFERENCE dataset (seed 6e5 + k, which no analysis iteration
#' uses) and held fixed across that cell's iterations, exactly as Phases 5, 7 and
#' 8 do. `phase5_prior_check.csv` measured the substitution as immaterial --
#' 3-12% of an estimate's Monte Carlo spread.
#'
#' The Gaussian arms share ONE set, derived from arm A's response, so that an arm
#' contrast is not confounded with a prior change. E and F take `bnec()`'s own
#' defaults for their family and response, because taking the defaults is the
#' practice under examination (Phase 7's decision, carried over unchanged).
p10_prior_set <- function(x, y, models, family) {
  setNames(lapply(models, function(m) {
    bayesnec:::define_prior(model = m, family = family, predictor = x,
                            response = y, prior_type = "uninformative")
  }), models)
}

#' Build every prior set one cell needs, from that cell's reference dataset
p10_cell_priors <- function(ref_sim, arms = P10_ARMS, model_set = P10_MODEL_SET,
                            meta = sim_meta()) {
  ref_prep <- prepare_sgr(ref_sim, "raw", meta = meta)
  out <- list()
  for (arm in arms) {
    mods <- p10_candidate_models(arm, ref_sim, model_set)
    if (arm %in% c("E", "F")) {
      d <- p10_zb_prep(ref_sim, arm)
      pr <- p10_prior_set(d$x, d$y, mods, p10_zb_family(arm))
    } else {
      ## Arm A's response for every Gaussian arm, including the censored and
      ## floored ones. `fix_bot_prior()` pins `bot` inside B3's set at fit time.
      pr <- p10_prior_set(ref_prep$x, ref_prep$y, mods, p10_family(arm))
    }
    out[[arm]] <- list(models = mods, prior = pr)
  }
  out
}

#' Fit one arm as a model average
#'
#' A, C and B3 go through `fit_arm()`, which is the study's definition of those
#' arms and is model-set aware after the Phase 9 changes. E and F are fitted here
#' because `fit_arm()` has no zero-bounded-family path; the call mirrors
#' `fit_zb()` in `phase8_run.R` with a set in place of `"nec3param"`.
p10_fit <- function(arm, sim, spec, mcmc = P10_MCMC, meta = sim_meta()) {
  if (arm %in% c("E", "F")) {
    t0 <- Sys.time()
    d <- p10_zb_prep(sim, arm)
    mod_arg <- paste0("c(", paste0("\"", spec$models, "\"", collapse = ", "), ")")
    fit <- bayesnec::bnec(
      bayesnec::bnf(paste0("y ~ crf(x, ", mod_arg, ")")), data = d,
      family = p10_zb_family(arm), prior = spec$prior,
      chains = mcmc$chains, iter = mcmc$iter, warmup = mcmc$warmup,
      seed = mcmc$seed, backend = "cmdstanr", init = P10_INIT,
      control = list(adapt_delta = mcmc$adapt_delta,
                     max_treedepth = mcmc$max_treedepth))
    return(list(fit = fit, arm = arm, n = nrow(d), n_censored = 0L,
                escalated = FALSE,
                elapsed = as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
  ## `rhat_threshold = P10_RHAT_ESCALATE` (Inf) disables `fit_arm()`'s refit-on-
  ## non-convergence path. Under a model set it triggers on the maximum R-hat
  ## across every component and re-samples all of them, which is both expensive
  ## and not the analyst's response to one stuck model. See R/setup_phase10.R.
  ## `init = P10_INIT` for every arm, so that B3 -- which is forced onto Stan's
  ## random inits under a model set -- does not differ from the others by its
  ## initialisation as well as by its constraint. See R/setup_phase10.R.
  fit_arm(arm, sim, prior = spec$prior, mcmc = mcmc, model = spec$models,
          meta = meta, rhat_threshold = P10_RHAT_ESCALATE, init = P10_INIT)
}

## Match BOTH `<cell>.rds` and `<cell>__<tag>_iNNN-NNN.rds`. The Gaussian arms
## for iterations 1-240 live in the UNSUFFIXED phase5_run.R file; only E/F
## (phase7) and the 241-500 top-up (phase8) carry a `__tag`. A pattern requiring
## `__` therefore drops arms A, C and B3 for the first 240 iterations without
## erroring -- the merge simply returns fewer rows, and the comparison quietly
## becomes E-versus-F. Dots are escaped because the cell name is full of them.
sweep_files <- function(cl, dir = "analysis/phase5") {
  list.files(dir, full.names = TRUE,
             pattern = paste0("^", gsub(".", "\\.", cl, fixed = TRUE),
                              "(__.*)?\\.rds$"))
}

#' Accept either a `fit_arm()`-style wrapper or a bare bnecfit
#'
#' `p10_fit()` returns a wrapper carrying timings; `amend()` returns the fit
#' itself. Both are read by the same downstream functions, and a `bayesnecfit`
#' confusingly *has* a `$fit` of its own (the brmsfit), so the test is on class
#' rather than on the presence of that element.
p10_pull <- function(x) {
  if (inherits(x, "bnecfit")) x else x$fit
}

#' Maximum R-hat of every model in a candidate set
#'
#' Equivalent to `bayesnec::rhat()`'s per-model `failed` flag, computed directly
#' off each component's brmsfit rather than through `pull_out()`, which rebuilds
#' a `bayesnecfit` per model and is wasteful at 13,800 fits. The quantity is the
#' same: `rhat.bayesnecfit()` takes `any(brms::rhat(pull_brmsfit(x)) > cutoff)`
#' over all parameters, and `clean_rhat_names()` renames without reordering, so
#' the maximum is unaffected. Gate 1 asserts the equivalence against the
#' package's own `rhat()` on `manec_example` rather than leaving it as a claim;
#' it lives there rather than in the test suite because it needs bayesnec
#' loaded.
p10_model_rhat <- function(fit) {
  f <- p10_pull(fit)
  if (!inherits(f, "bayesmanecfit")) {
    return(data.frame(model = if (!is.null(f$model)) f$model else NA_character_,
                      max_rhat = max(brms::rhat(f$fit), na.rm = TRUE),
                      row.names = NULL))
  }
  do.call(rbind, lapply(names(f$mod_fits), function(m) {
    data.frame(model = m,
               max_rhat = max(brms::rhat(f$mod_fits[[m]]$fit), na.rm = TRUE),
               row.names = NULL)
  }))
}

#' Drop the models that did not converge, and re-stack over the rest
#'
#' The analyst's step, done the package's way: `rhat()` to find them,
#' `amend(drop = )` to remove them and recompute the stacking weights over the
#' survivors. `amend()` reuses the existing fits and only refits models being
#' ADDED, so a pure drop costs a re-weighting and not a re-sample.
#'
#' Three outcomes, all of which happen and none of which may be silent:
#'  * nothing dropped -- the fit is returned untouched, and no `amend()` is
#'    called at all, so the common case costs nothing;
#'  * some dropped -- `amend()` returns a `bayesmanecfit` over the survivors, or
#'    a `bayesnecfit` if exactly one survives, which is also what an analyst
#'    would get and is handled downstream by `p10_pull()`;
#'  * all dropped -- `expand_manec()` stops, so `fit` comes back NULL with a
#'    reason. That is recorded as an arm failure for that dataset rather than
#'    substituted with the undropped fit, which would quietly reintroduce the
#'    models the rule exists to remove.
p10_drop_nonconverged <- function(fit, cutoff = P10_RHAT_CUTOFF) {
  rh <- p10_model_rhat(fit)
  rh$dropped <- rh$max_rhat > cutoff
  bad <- rh$model[rh$dropped]
  kept <- rh$model[!rh$dropped]
  if (!length(bad)) {
    return(list(fit = p10_pull(fit), rhat = rh, dropped = character(0),
                kept = kept, error = NA_character_))
  }
  if (!length(kept)) {
    return(list(fit = NULL, rhat = rh, dropped = bad, kept = kept,
                error = paste0("all ", length(bad), " models exceeded rhat ",
                               cutoff)))
  }
  amended <- try(suppressMessages(
    bayesnec::amend(p10_pull(fit), drop = bad)), silent = TRUE)
  if (inherits(amended, "try-error")) {
    return(list(fit = NULL, rhat = rh, dropped = bad, kept = kept,
                error = sub("\n.*", "", as.character(amended))))
  }
  list(fit = amended, rhat = rh, dropped = bad, kept = kept,
       error = NA_character_)
}

#' Per-model stacking weights, and what the candidate set actually became
#'
#' Recorded because an averaged estimate is uninterpretable without it: two arms
#' can report the same ErC50 from entirely different mixtures, and an arm that
#' has moved its weight onto a different shape is a finding in itself. Models
#' `bnec()` asked for but could not fit are dropped from the average with a
#' message and no error, so a set of eight that averaged over six is a different
#' analysis and must be visible rather than absorbed.
p10_weights <- function(fit, requested) {
  f <- p10_pull(fit)
  if (inherits(f, "bayesmanecfit")) {
    ms <- f$mod_stats
    w <- data.frame(model = rownames(ms), wi = ms$wi, row.names = NULL)
  } else {
    ## A set that collapsed to one model. Reported with weight 1 so the two
    ## cases tabulate together.
    w <- data.frame(model = if (!is.null(f$model)) f$model else NA_character_,
                    wi = 1, row.names = NULL)
  }
  w$n_requested <- length(requested)
  w$n_fitted <- nrow(w)
  w$dropped <- paste(setdiff(requested, w$model), collapse = " ")
  w
}

#' The model-averaged curve on a fixed grid
#'
#' Stored so that the phase's actual question -- is a flexible shape adequate
#' over the region ErCx is read from -- can be answered as weighted lack-of-fit
#' between ErC10 and ErC50, rather than inferred from the endpoints. About 1.2 MB
#' for the whole phase; the refit needed to recover it later would cost days.
#'
#' E's curve is multiplied back by the scale its response was divided by, so all
#' five arms share a y axis. `nec` and `beta` are unaffected by that rescaling
#' and the endpoints are invariant to it.
p10_curve <- function(fit, x_grid, scale = 1) {
  pe <- brms::posterior_epred(p10_pull(fit), newdata = data.frame(x = x_grid),
                              re_formula = NA)
  data.frame(
    x = x_grid,
    mu = scale * apply(pe, 2, mean),
    lo = scale * apply(pe, 2, stats::quantile, 0.025),
    hi = scale * apply(pe, 2, stats::quantile, 0.975)
  )
}
