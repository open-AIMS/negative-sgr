## Phase 9 -- the case studies re-run under bayesnec's model-averaged workflow.
##
## Phase 3 fixed every arm at `nec4param`. In the simulation that is
## unambiguously right: it is the model the data are generated from, so the arm
## comparison carries no model misspecification. On the four real datasets the
## justification is much weaker. The functional form is unknown, and fitting a
## single model chosen in advance is not what an analyst using `bayesnec` would
## do -- the package's whole workflow is to fit a candidate set and average over
## it. This phase asks whether the Phase 3 conclusion survives that workflow.
##
## The question is specifically NOT "which model is right". It is whether the
## data-handling convention still changes the reported toxicity when the model
## is averaged rather than assumed, because averaging could in principle absorb
## the distortion: a floored dataset might simply shift weight onto a different
## model shape and land in the same place.
##
## Design, and the two places it is not a clean translation of Phase 3:
##
## 1. THE CANDIDATE SET IS `decline`. Under a Gaussian family `check_models()`
##    drops the six zero-bounded models, leaving eight. Hormesis models are not
##    included: they answer a different question (is there low-dose
##    stimulation), and admitting them would change what the arms are being
##    compared on. This is stated rather than defaulted into.
##
## 2. ARMS B2 AND B3 CANNOT USE THE WHOLE SET. They pin the lower asymptote at
##    zero, and `neclin` and `ecxlin` have no asymptote to pin -- they decline
##    linearly without bound. Those two are dropped from the pinned arms only.
##    The consequence must be reported with the result: a pinned arm differs
##    from arm A by its candidate set as well as by its constraint, so its
##    contrast is not purely a likelihood contrast the way Phase 3's was.
##
## The shared-prior principle carries over unchanged: one prior per model, all
## built from arm A's response, handed to every arm, so that an arm contrast is
## not confounded with a prior change. `bnec()` takes them as a named list.
##
## Run:   Rscript analysis/phase9_modelavg.R
## Pilot: PHASE9_PILOT=1 Rscript analysis/phase9_modelavg.R   (one dataset, arms A and B1)

source("R/setup.R"); source("R/data_prep.R"); source("R/diagnostics.R")
load_bayesnec(); source("R/arms.R"); source("R/metrics.R")
use_compile_cache()
options(width = 200, mc.cores = min(STUDY_CORES, MCMC$chains))

PILOT <- identical(Sys.getenv("PHASE9_PILOT", "0"), "1")
## Models that did not converge are dropped from the average before anything is
## reported; see the R-hat rule below.
RHAT_MAX <- as.numeric(Sys.getenv("PHASE9_RHAT_MAX", "1.01"))
FITDIR <- "analysis/phase9_fits"
dir.create(FITDIR, showWarnings = FALSE, recursive = TRUE)
MODEL_SET <- "decline"
OUT <- "analysis"

datasets <- if (PILOT) "c_proliferum" else dataset_names()
ARMS <- if (PILOT) c("A", "B1") else c("A", "B1", "B2", "B3", "C", "C2", "D", "SQ")

cat("Phase 9 | model set:", MODEL_SET, "| datasets:", length(datasets),
    "| arms:", length(ARMS), "| drop models with R-hat >", RHAT_MAX, "\n\n")

#' The priors every arm of one dataset shares, one per candidate model
#'
#' Phase 3's `arm_prior()` builds a single `brmsprior` for `nec4param`. The same
#' idea, extended: build one for each model in the candidate set, always from
#' arm A's response vector, and hand the whole named list to every arm.
#' `bnec()` selects by name per model (`validate_priors()`), and falls back to
#' its own default with a message for any model the list does not name.
arm_prior_set <- function(x, y, models) {
  setNames(lapply(models, function(m) {
    bayesnec:::define_prior(model = m, family = stats::gaussian(link = "identity"),
                            predictor = x, response = y,
                            prior_type = "uninformative")
  }), models)
}

#' Per-model weights and diagnostics from a model-averaged fit
#'
#' Recorded because the averaged estimate is uninterpretable without them: two
#' arms can report the same ErC50 from entirely different mixtures, and an arm
#' that has moved its weight onto a different model shape is a finding in
#' itself. `bayesmanecfit` carries these in `mod_stats`; a fit that collapsed to
#' a single model is returned with weight 1 so the two cases tabulate together.
weights_table <- function(f, dataset, arm) {
  if (inherits(f, "bayesmanecfit")) {
    ms <- f$mod_stats
    data.frame(dataset = dataset, arm = arm, model = rownames(ms),
               wi = ms$wi, row.names = NULL)
  } else {
    data.frame(dataset = dataset, arm = arm,
               model = if (!is.null(f$model)) f$model else NA_character_,
               wi = 1, row.names = NULL)
  }
}

ends <- list(); wts <- list(); dgn <- list(); errs <- list(); nodrop <- list()

sfx <- if (PILOT) "_pilot" else ""
write_if <- function(x, name) {
  if (length(x)) {
    utils::write.csv(do.call(rbind, x),
                     file.path(OUT, paste0("phase9_", name, sfx, ".csv")),
                     row.names = FALSE)
  }
}
flush_outputs <- function() {
  write_if(ends, "endpoints"); write_if(wts, "weights")
  write_if(nodrop, "endpoints_nodrop")
  write_if(dgn, "diagnostics"); write_if(errs, "failures")
}

for (ds in datasets) {
  dat <- read_sgr(ds)
  prep_a <- prepare_sgr(dat, if (any(dat$density == 0)) "bound" else "raw",
                        meta = dataset_meta())
  models_all <- bayesnec:::expand_model_set(MODEL_SET)
  priors <- arm_prior_set(prep_a$x, prep_a$y, models_all)
  cat("=== ", ds, " | ", nrow(prep_a), " rows ===\n", sep = "")

  crossing <- NULL
  for (arm in ARMS) {
    ## Arms B2/B3 lose the two models with no asymptote to pin (see header).
    ## The SQ benchmark takes bayesnec's own defaults on a Beta response, so it
    ## is given neither the shared priors nor a restricted set.
    mods <- if (arm %in% c("B2", "B3")) bot_bearing_models(models_all) else models_all
    pr <- if (identical(arm, "SQ")) NULL else priors[intersect(names(priors), mods)]

    ## INITIALISATION. Arms take bayesnec's own initial-value search, except
    ## the pinned arms, which cannot.
    ##
    ## An earlier version of this script used Stan's random inits everywhere,
    ## because bayesnec's search costs 612.8 s against 6.1 s on floored data --
    ## it rejects any draw whose predicted curve leaves the observed response
    ## range, and floored data has no negative range. That change was validated
    ## on ONE model and was wrong for the set. Measured on `c_proliferum` arm A
    ## over the whole candidate set, same priors and seed:
    ##
    ##   random inits    5 of 8 models with R-hat > 1.01, 85.2% of the weight
    ##   bayesnec inits  0 of 6 models with R-hat > 1.01, 0% of the weight
    ##
    ## Random inits converge badly on most of these curve shapes. The speed was
    ## real and the convergence cost was worse, so the search is back.
    ##
    ## B2 and B3 are the exception and have no choice: `constant(0)` has no
    ## sampling distribution for the search to draw from, and `bnec()` passes a
    ## single `init` to every model in a set with no per-model hook. They take
    ## Stan's inits and are therefore expected to lose models to the R-hat rule
    ## below -- which is itself a finding about pinning an asymptote inside a
    ## model-averaged workflow, not a nuisance to be hidden.
    t0 <- Sys.time()
    f <- try(suppressMessages(
      fit_arm(arm, dat, prior = pr, crossing = crossing,
              model = mods, meta = dataset_meta(),
              init = if (arm %in% c("B2", "B3")) "random" else NULL,
              rhat_threshold = Inf)),
      silent = TRUE)
    if (inherits(f, "try-error")) {
      cat(sprintf("  [fail] %-3s %s\n", arm, sub("\n.*", "", as.character(f))))
      errs[[length(errs) + 1]] <- data.frame(dataset = ds, arm = arm,
                                             error = as.character(f))
      next
    }
    ## THE R-HAT RULE. A model that did not converge is not evidence, and an
    ## average that leans on one is not a result. Models with R-hat > 1.01 are
    ## dropped and the stacking weights re-solved over what remains, which is
    ## what `amend()` does -- weights come from an optimisation over the set, so
    ## they cannot simply be renormalised by hand.
    ##
    ## Both versions are kept. `phase9_endpoints.csv` is the reported analysis,
    ## after dropping; `phase9_endpoints_nodrop.csv` is the same fit before it,
    ## so the size of the correction is visible rather than asserted. An arm
    ## whose numbers move a long way under this rule was resting on models that
    ## had not converged, and that is worth seeing.
    fit_use <- f$fit
    dropped_rhat <- character(0)
    dd <- manec_model_diagnostics(f$fit, ds, arm)
    bad <- dd$model[!is.na(dd$max_rhat) & dd$max_rhat > RHAT_MAX]
    if (length(bad) && inherits(f$fit, "bayesmanecfit")) {
      if (length(bad) >= length(dd$model)) {
        cat("  [warn]", arm, "-- every model exceeds R-hat", RHAT_MAX,
            "; reported without dropping\n")
      } else {
        amended <- try(suppressMessages(bayesnec::amend(f$fit, drop = bad)),
                       silent = TRUE)
        if (inherits(amended, "try-error")) {
          cat("  [warn]", arm, "-- amend() failed, reported without dropping\n")
        } else {
          fit_use <- amended
          dropped_rhat <- bad
        }
      }
    }
    saveRDS(f, file.path(FITDIR, sprintf("%s__%s.rds", ds, arm)))

    ## Arm D truncates at arm A's zero crossing, taken from the AVERAGED arm-A
    ## fit here rather than from a single model -- the same rule as Phase 3,
    ## applied to the workflow under test. Taken AFTER the R-hat rule, so the
    ## truncation point does not come from models the analysis discards.
    if (identical(arm, "A")) crossing <- zero_crossing(fit_use)

    ## `n_requested` is what was asked for; `n_averaged` is what the estimate
    ## actually rests on. They differ for two independent reasons and both
    ## matter: `check_models()` drops models the family cannot support (six of
    ## the fourteen declining models under a Gaussian), and any model whose
    ## chains fail is dropped from the average by `bnec()` with a message. A
    ## single "n_models" column would hide both.
    fitted_models <- if (inherits(fit_use, "bayesmanecfit"))
      names(fit_use$mod_fits) else fit_use$model
    et <- endpoint_table(fit_use, arm, ds)
    nodrop[[length(nodrop) + 1]] <- cbind(endpoint_table(f$fit, arm, ds),
      n_averaged = if (inherits(f$fit, "bayesmanecfit"))
        length(f$fit$mod_fits) else 1L)
    ends[[length(ends) + 1]] <- cbind(et, n = f$n, n_censored = f$n_censored,
                                      escalated = isTRUE(f$escalated),
                                      n_requested = length(mods),
                                      n_averaged = length(fitted_models),
                                      dropped_rhat = paste(dropped_rhat,
                                                           collapse = " "),
                                      elapsed = f$elapsed)
    w <- weights_table(fit_use, ds, arm)
    ## Models bnec asked for but could not fit are dropped from the average with
    ## a message and no error -- `ecxll4` loses a chain on some of these
    ## datasets. That is the package's own behaviour and what an analyst would
    ## get, so it is recorded rather than worked around: a candidate set of
    ## eight that averaged over six is a different analysis from one that
    ## averaged over eight.
    w$requested <- length(mods)
    w$dropped <- paste(setdiff(mods, w$model), collapse = " ")
    w$dropped_rhat <- paste(dropped_rhat, collapse = " ")
    wts[[length(wts) + 1]] <- w
    dgn[[length(dgn) + 1]] <- cbind(escalated = isTRUE(f$escalated), dd)
    cat(sprintf("  [ok]   %-3s %2d of %2d models averaged  %5.1f min\n", arm,
                length(fitted_models), length(mods),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    ## Written after every arm rather than once at the end. An earlier version
    ## of this script lost five hours of fits to a deliberate restart because
    ## nothing was on disk until the final line; the simulation runners learned
    ## the same lesson in Phase 7. Each file is rewritten whole from the
    ## accumulated list, which is cheap at this size and leaves a complete file
    ## at every instant.
    flush_outputs()
  }
}

flush_outputs()
if (length(ends)) print(do.call(rbind, ends), digits = 4, row.names = FALSE)
cat("\nPhase 9 complete at", format(Sys.time()), "\n")
