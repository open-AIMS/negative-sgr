## Per-dataset diagnostics. Nothing here fits a model.

#' Gap between zero and the smallest non-zero value, in recording units
#'
#' Reproduced from `bayesnec`'s example6 vignette. A ratio near 1 means the
#' smallest non-zero value is one rounding unit above zero, so zero is simply
#' the next value down (censoring or rounding). A large ratio means a real gap
#' between "no response" and "the smallest response measured" (a distinct
#' process, i.e. a hurdle case).
#'
#' Applied here to **cell density**, not to SGR. SGR is not zero-bounded, so its
#' zeros are not the kind of zero the diagnostic is about; the zeros in these
#' datasets arise on the density scale, where the diagnostic is meaningful.
gap_ratio <- function(y) {
  y <- y[!is.na(y)]
  nonzero <- y[y > 0]
  if (length(unique(nonzero)) < 2) return(NA_real_)
  resolution <- min(diff(sort(unique(nonzero))))
  min(nonzero) / resolution
}

#' Is the group-mean response monotone decreasing in concentration?
#'
#' Returns the number of increases, their largest size relative to the control
#' mean, and where they occur, rather than a bare TRUE/FALSE. Sampling noise
#' guarantees some increases in any real dataset; what matters is whether any of
#' them is large and consistent across replicates.
monotonicity <- function(dat, response = "sgr") {
  gm <- tapply(dat[[response]], dat$x, mean)
  concs <- as.numeric(names(gm))
  d <- diff(gm)
  ctrl <- gm[1]
  up <- which(d > 0)
  list(
    n_increases = length(up),
    max_increase = if (length(up)) max(d[up]) else 0,
    max_increase_rel_control = if (length(up)) max(d[up]) / abs(ctrl) else 0,
    at = if (length(up)) concs[up + 1] else numeric(0),
    group_means = gm
  )
}

#' Replicate-level reversals between adjacent concentrations
#'
#' `monotonicity()` works on group means, which cannot distinguish one aberrant
#' replicate from a whole treatment sitting above its neighbour. `dominance` is
#' the fraction of replicate pairs (i at concentration k, j at k+1) with
#' `y_j > y_i`: 1 means complete separation, 0.5 means the two treatments are
#' interleaved.
#'
#' A strict all-pairs rule was tried first and rejected: `c_proliferum`'s
#' 15-to-20 reversal has dominance 0.88 (22 of 25 pairs), not 1, because the
#' two treatments overlap slightly at their edges. Reporting the statistic
#' rather than a bare TRUE/FALSE keeps that visible; the default threshold is
#' set below 0.88 so the reversal is caught.
#'
#' `in_decline` marks reversals where both treatments have already fallen below
#' half the control mean, which separates a real reversal in the declining limb
#' from ordinary noise among the near-control treatments.
reversals <- function(dat, response = "sgr", min_dominance = 0.85) {
  concs <- sort(unique(dat$x))
  ctrl <- mean(dat[[response]][dat$x == concs[1]])
  out <- NULL
  for (i in seq_len(length(concs) - 1L)) {
    lower <- dat[[response]][dat$x == concs[i]]
    upper <- dat[[response]][dat$x == concs[i + 1L]]
    dom <- mean(outer(upper, lower, ">"))
    if (dom >= min_dominance) {
      dens_l <- mean(dat$density[dat$x == concs[i]])
      dens_u <- mean(dat$density[dat$x == concs[i + 1L]])
      out <- rbind(out, data.frame(
        from = concs[i], to = concs[i + 1L], dominance = dom,
        mean_lower = mean(lower), mean_upper = mean(upper),
        min_upper = min(upper), max_lower = max(lower),
        density_ratio = dens_u / dens_l,
        in_decline = mean(lower) < ctrl / 2 && mean(upper) < ctrl / 2
      ))
    }
  }
  if (is.null(out)) {
    out <- data.frame(from = numeric(0), to = numeric(0),
                      dominance = numeric(0), mean_lower = numeric(0),
                      mean_upper = numeric(0), min_upper = numeric(0),
                      max_lower = numeric(0), density_ratio = numeric(0),
                      in_decline = logical(0))
  }
  out
}

#' Everything about one dataset, in one row
#'
#' `d` is the deepest group-mean decline, so `delta = d / mu_0` is the
#' discarded effect fraction: how much of the effect axis a 100% cap throws
#' away. For a dataset whose deepest group mean is a censoring bound it is a
#' lower bound, flagged by `d_censored`.
dataset_summary <- function(dat, meta = dataset_meta()) {
  ds <- unique(dat$dataset)
  design <- recover_design(dat)
  ctrl <- dat$sgr[dat$x == min(dat$x)]
  mu_0 <- mean(ctrl)
  R <- exp(mu_0 * design$t)
  gm <- tapply(dat$sgr, dat$x, mean)
  # The floored zeros in the supplied SGR column make the deepest *supplied*
  # group mean an underestimate of the decline, so take the deepest mean over
  # the "bound" preparation, where undefined rows carry their LOD bound.
  bnd <- prepare_sgr(dat, "bound", meta = meta)
  gm_bound <- tapply(bnd$y, bnd$x, mean)
  d <- -min(gm_bound)
  lod_obs <- min_detected_density(dat)
  mono <- monotonicity(dat)
  data.frame(
    dataset = ds,
    n = nrow(dat),
    n_conc = length(unique(dat$x)),
    n_control = length(ctrl),
    n_negative = sum(dat$sgr < 0),
    n_zero_sgr = sum(dat$sgr == 0),
    n_zero_density = sum(dat$density == 0),
    min_positive_sgr = min(dat$sgr[dat$sgr > 0]),
    min_positive_density = min(dat$density[dat$density > 0]),
    mu_0 = mu_0,
    sd_control = stats::sd(ctrl),
    cv_control = stats::sd(ctrl) / mu_0,
    t = design$t,
    n_0 = design$n_0,
    design_r_squared = design$r_squared,
    R = R,
    effect_at_zero_density_scale = 1 - 1 / R,
    d = d,
    d_censored = any(dat$density == 0 &
                       dat$x == as.numeric(names(which.min(gm_bound)))),
    delta = d / mu_0,
    f_neg = mean(dat$sgr < 0 | dat$density == 0),
    lod_observed = lod_obs,
    lod_meta = meta$lod[meta$dataset == ds],
    lod_bound = attr(bnd, "lod_bound"),
    gap_ratio_density = gap_ratio(dat$density),
    n_increases = mono$n_increases,
    max_increase_rel_control = mono$max_increase_rel_control,
    n_reversals = nrow(reversals(dat)),
    n_reversals_in_decline = sum(reversals(dat)$in_decline),
    stringsAsFactors = FALSE
  )
}

#' The default bayesnec gaussian priors implied by a response vector
#'
#' `bnec(prior_type = "uninformative")`, the default, sets
#' `top ~ normal(quantile(y, 0.9), 2.5 * sd(y))` and
#' `bot ~ normal(quantile(y, 0.1), 2.5 * sd(y))`. Because both are functions of
#' the response vector, preparing the data differently also moves the prior --
#' which is why every arm of a dataset is given the arm-A prior explicitly
#' rather than the default.
default_gaussian_prior <- function(y) {
  list(
    top_mean = unname(stats::quantile(y, 0.9)),
    bot_mean = unname(stats::quantile(y, 0.1)),
    sd = 2.5 * stats::sd(y)
  )
}
