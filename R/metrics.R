## Endpoint extraction and fit diagnostics.
##
## Only ErC10, ErC50 and NSEC are extracted. Nothing above ErC50 is reported:
## the upper end of the effect axis is exactly the region the study argues is
## design-determined and uninterpretable, so reporting it would undercut the
## argument.

#' Sampler health for one fit
#'
#' Recorded for every fit in every arm and iteration. A study whose conclusion
#' rests on interval coverage cannot discard bad fits silently, so the exclusion
#' rule is applied downstream from these numbers rather than at fit time.
fit_diagnostics <- function(fit) {
  ## A model-averaged fit has no single `$fit`: a `bayesmanecfit` holds one
  ## brmsfit per candidate model under `$mod_fits`, and reading `$fit` off it
  ## returns NULL, which surfaces much later as an opaque error inside
  ## `array()`. Aggregate over the set instead -- worst R-hat, total divergent
  ## transitions, worst effective-sample ratio -- so that a single bad model in
  ## an average cannot be hidden by the good ones. Per-model numbers are the
  ## caller's business; see `manec_model_diagnostics()`.
  if (inherits(fit, "bayesmanecfit")) {
    per <- do.call(rbind, lapply(fit$mod_fits, function(m) fit_diagnostics(m)))
    return(data.frame(
      divergences = sum(per$divergences),
      max_treedepth_hit = sum(per$max_treedepth_hit),
      max_rhat = max(per$max_rhat, na.rm = TRUE),
      min_neff_ratio = min(per$min_neff_ratio, na.rm = TRUE),
      stringsAsFactors = FALSE
    ))
  }
  bf <- fit$fit
  np <- brms::nuts_params(bf)
  divergent <- sum(np$Value[np$Parameter == "divergent__"])
  treedepth <- np$Value[np$Parameter == "treedepth__"]
  rh <- brms::rhat(bf)
  ess <- brms::neff_ratio(bf)
  data.frame(
    divergences = divergent,
    max_treedepth_hit = sum(treedepth >= max(treedepth)) *
      as.integer(max(treedepth) >= 12),
    max_rhat = max(rh, na.rm = TRUE),
    min_neff_ratio = min(ess, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

#' Sampler health of every model inside a model-averaged fit
#'
#' The aggregate in `fit_diagnostics()` answers "is anything wrong here", which
#' is what an exclusion rule needs. This answers "where", which is what a reader
#' needs: a model that sampled badly but carries 1% of the weight is a different
#' situation from one that carries 60%.
manec_model_diagnostics <- function(fit, dataset = NA, arm = NA) {
  if (!inherits(fit, "bayesmanecfit")) {
    return(cbind(dataset = dataset, arm = arm,
                 model = if (!is.null(fit$model)) fit$model else NA_character_,
                 wi = 1, fit_diagnostics(fit)))
  }
  ms <- fit$mod_stats
  do.call(rbind, lapply(names(fit$mod_fits), function(m) {
    cbind(dataset = dataset, arm = arm, model = m,
          wi = ms$wi[match(m, rownames(ms))],
          fit_diagnostics(fit$mod_fits[[m]]))
  }))
}

#' The three reported endpoints, with 95% credible intervals
#'
#' ECx is `type = "absolute"`, which in bayesnec is the decline from `top` to
#' zero -- algebraically ErCx once the response is raw specific growth rate.
#' That is the OECD TG 201 definition and the reason no normalisation is applied
#' anywhere in this study.
endpoint_table <- function(fit, arm = NA, dataset = NA, xform = NULL) {
  get_one <- function(f, label) {
    v <- tryCatch(f(), error = function(e) c(NA_real_, NA_real_, NA_real_))
    data.frame(endpoint = label, estimate = unname(v[1]),
               lower = unname(v[2]), upper = unname(v[3]),
               stringsAsFactors = FALSE)
  }
  rows <- rbind(
    get_one(function() bayesnec::ecx(fit, ecx_val = 10, type = "absolute"),
            "ErC10"),
    get_one(function() bayesnec::ecx(fit, ecx_val = 50, type = "absolute"),
            "ErC50"),
    get_one(function() bayesnec::nsec(fit), "NSEC")
  )
  cbind(dataset = dataset, arm = arm, rows, row.names = NULL)
}

#' Posterior summaries of the curve parameters
#'
#' `bot` is included because Phase 4 needs its posterior, and because a fit with
#' `bot` held constant should show exactly zero posterior variation there --
#' a cheap check that the constant prior did what was intended.
parameter_table <- function(fit, arm = NA, dataset = NA) {
  dr <- brms::as_draws_df(fit$fit)
  pars <- grep("^b_(top|bot|beta|nec)_Intercept$", names(dr), value = TRUE)
  sig <- grep("^sigma$", names(dr), value = TRUE)
  out <- do.call(rbind, lapply(c(pars, sig), function(p) {
    v <- dr[[p]]
    data.frame(parameter = sub("^b_|_Intercept$", "", p),
               mean = mean(v), sd = stats::sd(v),
               q2.5 = unname(stats::quantile(v, 0.025)),
               q50 = unname(stats::quantile(v, 0.5)),
               q97.5 = unname(stats::quantile(v, 0.975)),
               stringsAsFactors = FALSE)
  }))
  cbind(dataset = dataset, arm = arm, out, row.names = NULL)
}

#' Prior-to-posterior contraction for one parameter (Phase 4)
#'
#' `1 - sd(posterior)/sd(prior)`: 0 means the data said nothing and the
#' posterior is the prior; near 1 means the parameter is data-determined.
#' `bnec()` draws prior samples by default (`sample_prior = "yes"`), so the
#' prior sd is taken from the fit rather than re-derived from the prior string.
contraction <- function(fit, parameter = "bot") {
  dr <- brms::as_draws_df(fit$fit)
  post <- dr[[paste0("b_", parameter, "_Intercept")]]
  prior <- dr[[paste0("prior_b_", parameter)]]
  if (is.null(prior)) {
    prior <- dr[[paste0("prior_b_", parameter, "_Intercept")]]
  }
  if (is.null(post) || is.null(prior)) {
    return(data.frame(parameter = parameter, prior_sd = NA_real_,
                      posterior_sd = NA_real_, contraction = NA_real_))
  }
  data.frame(parameter = parameter,
             prior_sd = stats::sd(prior),
             posterior_sd = stats::sd(post),
             contraction = 1 - stats::sd(post) / stats::sd(prior))
}

#' Monte Carlo standard errors for the simulation metrics
#'
#' Reported on every quantity in Phase 5. Coverage is a proportion, so its MCSE
#' is binomial; bias and RMSE use the usual delta-method forms.
mcse_summary <- function(est, truth, lower, upper) {
  ok <- !is.na(est) & !is.na(lower) & !is.na(upper)
  est <- est[ok]; lower <- lower[ok]; upper <- upper[ok]
  n <- length(est)
  if (n == 0) {
    return(data.frame(n = 0, bias = NA, bias_mcse = NA, rmse = NA,
                      rmse_mcse = NA, coverage = NA, coverage_mcse = NA))
  }
  err <- est - truth
  cov <- as.numeric(truth >= lower & truth <= upper)
  bias <- mean(err)
  rmse <- sqrt(mean(err^2))
  data.frame(
    n = n,
    bias = bias,
    bias_mcse = stats::sd(err) / sqrt(n),
    rmse = rmse,
    # Var(mean(err^2)) / (2 rmse)^2, the delta method on the square root.
    rmse_mcse = stats::sd(err^2) / (2 * rmse * sqrt(n)),
    coverage = mean(cov),
    coverage_mcse = sqrt(mean(cov) * (1 - mean(cov)) / n)
  )
}

#' Iterations needed for a target MCSE on coverage
#'
#' Used to turn the Phase 5 pilot into a budget rather than a guess.
iterations_for_coverage_mcse <- function(target_mcse, p = 0.95) {
  ceiling(p * (1 - p) / target_mcse^2)
}
