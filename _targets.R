## targets pipeline for the real-data phases (3 and 4).
##
## Phase 1 (branch verification) and Phase 2 (diagnostics) are standalone
## scripts in analysis/: they are cheap, they gate the pipeline rather than feed
## it, and their whole purpose is to be read by a human. Phase 5 (simulation) is
## its own pipeline in analysis/phase5_*.R because its unit of work is an
## iteration, not a dataset.
##
## Every fit is its own target, keyed by dataset and arm, so a re-run after a
## bayesnec commit change costs only the fits that actually changed.

library(targets)

source("R/setup.R")
source("R/data_prep.R")
source("R/diagnostics.R")
source("R/arms.R")
source("R/metrics.R")
source("R/figures.R")
load_bayesnec()

# Arms that did not apply to a dataset are carried as stubs with a NULL fit.
keep_fits <- function(fits) Filter(function(a) !is.null(a$fit), fits)

# Sequential: each brms fit already runs its own chains in parallel, and the
# machine is shared with another study, so there is nothing to gain from
# running targets in parallel on top of that.
tar_option_set(packages = c("brms", "ggplot2"), format = "rds")
options(mc.cores = min(STUDY_CORES, MCMC$chains))
use_compile_cache()

list(
  tar_target(ds_name, dataset_names()),
  tar_target(dat, read_sgr(ds_name), pattern = map(ds_name),
             iteration = "list"),

  ## The prior every arm of a dataset shares, built once from arm A's response.
  tar_target(
    prior_a,
    {
      a <- prepare_sgr(dat, if (any(dat$density == 0)) "bound" else "raw")
      arm_prior(a$x, a$y)
    },
    pattern = map(dat), iteration = "list"
  ),

  ## Arm A first: arm D's truncation point comes from it.
  tar_target(fit_A, fit_arm("A", dat, prior_a),
             pattern = map(dat, prior_a), iteration = "list"),
  tar_target(crossing_A, zero_crossing(fit_A$fit), pattern = map(fit_A)),

  tar_target(fit_B1, fit_arm("B1", dat, prior_a),
             pattern = map(dat, prior_a), iteration = "list"),
  tar_target(fit_B2, fit_arm("B2", dat, prior_a),
             pattern = map(dat, prior_a), iteration = "list"),
  tar_target(fit_B3, fit_arm("B3", dat, prior_a),
             pattern = map(dat, prior_a), iteration = "list"),
  tar_target(fit_C, fit_arm("C", dat, prior_a),
             pattern = map(dat, prior_a), iteration = "list"),
  tar_target(fit_C2, fit_arm("C2", dat, prior_a),
             pattern = map(dat, prior_a), iteration = "list"),
  tar_target(fit_D, fit_arm("D", dat, prior_a, crossing = crossing_A),
             pattern = map(dat, prior_a, crossing_A), iteration = "list"),
  tar_target(fit_SQ, fit_arm("SQ", dat, NULL),
             pattern = map(dat), iteration = "list"),

  # Sensitivity arm for the two datasets where "raw" means dropping rows.
  # Returns a stub rather than NULL on the datasets it does not apply to,
  # because a NULL target value is awkward for targets to store and branch over.
  tar_target(fit_A_dropped,
             if (any(dat$density == 0)) {
               fit_arm("A_raw_dropped", dat, prior_a)
             } else {
               list(fit = NULL, arm = "A_raw_dropped",
                    dataset = unique(dat$dataset), applicable = FALSE)
             },
             pattern = map(dat, prior_a), iteration = "list"),

  tar_target(
    all_fits,
    list(fit_A, fit_B1, fit_B2, fit_B3, fit_C, fit_C2, fit_D, fit_SQ,
         fit_A_dropped),
    pattern = map(fit_A, fit_B1, fit_B2, fit_B3, fit_C, fit_C2, fit_D, fit_SQ,
                  fit_A_dropped),
    iteration = "list"
  ),

  tar_target(
    endpoints,
    do.call(rbind, lapply(keep_fits(all_fits), function(a)
      cbind(endpoint_table(a$fit, a$arm, a$dataset),
            preparation = a$preparation, n = a$n, n_censored = a$n_censored,
            elapsed = a$elapsed))),
    pattern = map(all_fits)
  ),

  tar_target(
    parameters,
    do.call(rbind, lapply(keep_fits(all_fits), function(a)
      parameter_table(a$fit, a$arm, a$dataset))),
    pattern = map(all_fits)
  ),

  tar_target(
    diagnostics_tab,
    do.call(rbind, lapply(keep_fits(all_fits), function(a)
      cbind(dataset = a$dataset, arm = a$arm,
            escalated = isTRUE(a$escalated), fit_diagnostics(a$fit)))),
    pattern = map(all_fits)
  ),

  ## Phase 4: how much of the bot posterior is prior.
  tar_target(
    bot_contraction,
    do.call(rbind, lapply(
      Filter(function(a) !a$arm %in% c("B2", "B3", "SQ"), keep_fits(all_fits)),
      function(a) cbind(dataset = a$dataset, arm = a$arm,
                        contraction(a$fit, "bot")))),
    pattern = map(all_fits)
  ),

  tar_target(
    arm_curve_plot,
    plot_arm_curves(keep_fits(all_fits), dat),
    pattern = map(all_fits, dat), iteration = "list"
  )
)
