## Reading and preparing the four microalgal SGR datasets.
##
## The raw CSVs are inconsistent: the concentration column is "X" or "x", the
## density column is "CellDensity" or "CellDensity_Final", and the column order
## differs. Everything downstream sees the same three columns (x, sgr, density)
## plus a dataset label, so no arm function has to know which file it came from.

dataset_names <- function() {
  c("c_proliferum", "c_proliferum2", "r_salina", "r_salina2")
}

#' Static per-dataset facts
#'
#' `lod` is the counter's limit of detection in cells/mL. It is not recorded in
#' the CSVs; it is recovered arithmetically in `recover_design()` and confirmed
#' here. `species` and `duration_nominal` come from the test type, not the data.
dataset_meta <- function() {
  data.frame(
    dataset = dataset_names(),
    species = c("Cladocopium proliferum", "Cladocopium proliferum",
                "Rhodomonas salina", "Rhodomonas salina"),
    lod = c(NA, NA, 10, 100),
    # The density below which a population is extinct, used as the lower end of
    # the interval-censored preparation. 1 cell/mL is a physical floor, not a
    # measured one; it is deliberately conservative, since its only job is to
    # stop the lower asymptote running to minus infinity (see `prepare_sgr`,
    # preparation "interval").
    physical_floor = c(1, 1, 1, 1),
    duration_nominal = c(7, 7, 3, 3),
    stringsAsFactors = FALSE
  )
}

#' Read one raw CSV and normalise its columns
#'
#' @param dataset One of `dataset_names()`.
#' @param path Directory holding the raw CSVs.
#' @return A data.frame with columns dataset, x, sgr, density.
read_sgr <- function(dataset, path = "data-raw") {
  dataset <- match.arg(dataset, dataset_names())
  raw <- utils::read.csv(file.path(path, paste0(dataset, ".csv")),
                         check.names = FALSE)
  x_col <- grep("^[Xx]$", names(raw), value = TRUE)[1]
  d_col <- grep("^CellDensity", names(raw), value = TRUE)[1]
  s_col <- grep("^SGR$", names(raw), value = TRUE)[1]
  if (any(is.na(c(x_col, d_col, s_col)))) {
    stop("Could not identify x / density / SGR columns in ", dataset,
         ": found ", paste(names(raw), collapse = ", "))
  }
  # x is coerced to double because two of the four files record whole-number
  # concentrations, which read.csv types as integer; bayesnec rejects an integer
  # x_var outright ("does not currently support integer concentration data").
  data.frame(
    dataset = dataset,
    x = as.double(raw[[x_col]]),
    sgr = as.double(raw[[s_col]]),
    density = as.double(raw[[d_col]]),
    stringsAsFactors = FALSE
  )
}

#' Recover the test design (t and n_0) from the density/SGR pairs
#'
#' SGR = (log(N_t) - log(n_0)) / t, so log(N_t) = log(n_0) + t * SGR. With a
#' single inoculum density shared across the test, regressing log(density) on
#' SGR returns t as the slope and log(n_0) as the intercept. Zero-density rows
#' carry no information (log(0)) and are excluded.
#'
#' This is a regression rather than a two-point solve so that the residual SD
#' reports whether n_0 really was constant: an R^2 of exactly 1 means the SGR
#' column was computed deterministically from these densities and one n_0.
recover_design <- function(dat) {
  ok <- dat$density > 0
  fit <- stats::lm(log(density) ~ sgr, data = dat[ok, ])
  t_hat <- unname(stats::coef(fit)[2])
  list(
    t = t_hat,
    n_0 = unname(exp(stats::coef(fit)[1])),
    resid_sd = stats::sd(stats::residuals(fit)),
    # Computed directly rather than via summary.lm(), which warns
    # "essentially perfect fit" on exactly the datasets where the fit IS
    # perfect -- which is the informative case here, not a problem.
    r_squared = 1 - stats::var(stats::residuals(fit)) /
      stats::var(log(dat$density[ok])),
    n_used = sum(ok)
  )
}

#' The growth rate implied by a density sitting exactly at the counting limit
#'
#' A row recorded as zero density means "fewer than `lod` cells/mL", so its true
#' growth rate is at most this value. It is a left-censoring bound, not a
#' measurement.
lod_bound <- function(lod, n_0, t) {
  (log(lod) - log(n_0)) / t
}

#' Infer the LOD from the data alone
#'
#' The smallest non-zero density in a dataset that reached the counting limit is
#' the limit itself. Returned so the value asserted in `dataset_meta()` can be
#' checked rather than trusted.
infer_lod <- function(dat) {
  pos <- dat$density[dat$density > 0]
  if (length(pos) == 0) NA_real_ else min(pos)
}

#' Build one of the four data preparations
#'
#' The preparations differ in how they treat two kinds of row: genuinely
#' measured negative growth rates, and rows where the population fell below the
#' counting limit so that SGR is *undefined* rather than negative.
#'
#' | preparation | negatives   | undefined (density 0)                  |
#' |-------------|-------------|----------------------------------------|
#' | raw         | as measured | row dropped                            |
#' | bound       | as measured | set to the LOD bound, treated observed |
#' | supplied    | as measured | set to 0 (what the CSVs contain)       |
#' | floored     | set to 0    | set to 0                               |
#' | censored    | left-censored at 0 | left-censored at the LOD bound  |
#' | interval    | interval-censored [extinction, 0] | interval-censored [extinction, LOD bound] |
#'
#' `supplied` exists because the raw CSVs turned out **not** to be fully
#' floored: the labs substituted 0 for the undetected rows but left genuinely
#' measured negative growth rates intact (r_salina keeps -1.62 and -1.99;
#' r_salina2 keeps 17 negatives). The 100%-inhibition cap is therefore applied
#' downstream in their curve-fitting, not in the delivered data. `floored` is
#' our application of that convention, and the study must say so rather than
#' claiming the delivered data already embody it.
#'
#' `censored` additionally left-censors genuinely negative rows at 0 when
#' `cens_negatives = TRUE`, which is arm C's censoring-as-coarsening treatment.
#'
#' The returned `cens` column uses brms' character encoding ("none"/"left"). A
#' row sitting *exactly* at the LOD was detected, so it stays "none" even though
#' its response value equals the bound; only `density == 0` is censored.
prepare_sgr <- function(dat, preparation = c("raw", "bound", "supplied",
                                             "floored", "censored", "interval"),
                        cens_negatives = FALSE, meta = dataset_meta()) {
  preparation <- match.arg(preparation)
  ds <- unique(dat$dataset)
  stopifnot(length(ds) == 1)
  design <- recover_design(dat)
  lod <- meta$lod[meta$dataset == ds]
  floor_density <- meta$physical_floor[meta$dataset == ds]
  undefined <- dat$density == 0
  if (any(undefined) && is.na(lod)) {
    stop("Dataset ", ds, " has zero-density rows but no LOD in dataset_meta().")
  }
  bound <- if (is.na(lod)) NA_real_ else lod_bound(lod, design$n_0, design$t)

  out <- dat
  out$y <- out$sgr
  out$cens <- "none"
  out$undefined <- undefined

  if (preparation == "raw") {
    out <- out[!undefined, , drop = FALSE]
  } else if (preparation == "bound") {
    out$y[undefined] <- bound
  } else if (preparation == "supplied") {
    # Exactly what the CSV contains: undetected rows carry 0, measured
    # negatives are untouched. Asserted identical to the supplied column in
    # tests/testthat/test-data_prep.R.
    out$y[undefined] <- 0
  } else if (preparation == "floored") {
    # The convention under study, applied by us: substitute 0 for anything at
    # or below zero growth.
    out$y[undefined] <- 0
    out$y[out$y < 0] <- 0
  } else if (preparation == "censored") {
    out$y[undefined] <- bound
    out$cens[undefined] <- "left"
    if (cens_negatives) {
      neg <- !undefined & out$y < 0
      out$y[neg] <- 0
      out$cens[neg] <- "left"
    }
  } else if (preparation == "interval") {
    # Why this exists, measured rather than assumed (c_proliferum2, 4 chains,
    # 2000 iter): arm A puts bot at -1.27 [-1.64, -0.94]; arm C, which
    # left-censors the same rows at zero, puts it at -0.42 [-0.84, -0.16].
    #
    # The reason is saturation, not a push. A left-censored row contributes
    # Phi((0 - mu)/sigma), and sigma here is ~0.01-0.05, so that contribution is
    # already ~1 once mu is a few hundredths below zero. The likelihood is
    # therefore FLAT in `bot` across the whole region the data are actually in,
    # and the prior -- centred at quantile(y, 0.1) = -0.043 -- decides where in
    # that flat region the posterior sits. Censoring at zero does not coarsen
    # the information about `bot`; it removes the information that identified
    # it, and the reported asymptote becomes the prior's.
    #
    # Interval censoring does not restore identification either -- nothing can,
    # once the magnitudes are discarded -- but it bounds the flat region on both
    # sides with a statement that is true (the growth rate lies between
    # extinction and zero) instead of on one side only.
    # A row detected exactly AT the counting limit gets the wider interval
    # [floor, 0] while an undetected row gets the tighter [floor, LOD bound].
    # That inversion is the coarsening convention showing its cost, not a bug:
    # under a rule that reports declining replicates only as "not positive", a
    # measured value becomes less informative than a non-detect.
    floor_bound <- lod_bound(floor_density, design$n_0, design$t)
    out$y2 <- out$y
    coarse <- undefined | out$y < 0
    out$y[coarse] <- floor_bound
    out$y2[undefined] <- bound
    out$y2[!undefined & out$sgr < 0] <- 0
    out$cens[coarse] <- "interval"
    attr(out, "floor_bound") <- floor_bound
  }
  if (!"y2" %in% names(out)) {
    # brms ignores the upper bound on uncensored rows, but the column must be
    # present and finite for every row.
    out$y2 <- out$y
  }
  attr(out, "preparation") <- preparation
  attr(out, "lod_bound") <- bound
  attr(out, "design") <- design
  rownames(out) <- NULL
  out
}

#' Drop concentrations at or above a zero-crossing estimate (arm D)
#'
#' Stated as an algorithm rather than an outcome: the crossing must come from a
#' fit to the same dataset (or, in simulation, the same iteration), never from
#' the truth.
truncate_at_crossing <- function(dat, crossing) {
  keep <- dat$x < crossing
  if (sum(keep) < 10 || length(unique(dat$x[keep])) < 4) {
    warning("Truncation at ", signif(crossing, 3), " leaves ", sum(keep),
            " rows across ", length(unique(dat$x[keep])),
            " concentrations; arm D may not be estimable.")
  }
  out <- dat[keep, , drop = FALSE]
  attr(out, "crossing") <- crossing
  rownames(out) <- NULL
  out
}

#' Percent inhibition capped at 100 (the SQ benchmark's response)
#'
#' Returned on the 0-1 scale as a *proportion of the control retained*, so that
#' a bounded family is fitted to a declining curve in the usual orientation.
#' This is the literal lab practice and is deliberately not comparable to the
#' other arms: it changes the family, the likelihood and the data at once.
percent_inhibition <- function(dat) {
  control <- mean(dat$sgr[dat$x == min(dat$x)], na.rm = TRUE)
  y <- dat$sgr / control
  y[dat$density == 0] <- 0
  y[y < 0] <- 0
  # Between 15 and 28 rows per dataset exceed the control mean (ordinary
  # replicate variation, and in two datasets a treated maximum), which a
  # 0-1 bounded family cannot represent. Capping at 1 is the second half of the
  # same convention: inhibition is bounded to [0, 100]%. Both caps are counted
  # so the benchmark can report how much data the practice truncates.
  n_capped_high <- sum(y > 1)
  n_capped_low <- sum(y == 0)
  y[y > 1] <- 1
  out <- dat
  out$y <- y
  out$cens <- "none"
  attr(out, "control_mean") <- control
  attr(out, "n_capped_high") <- n_capped_high
  attr(out, "n_capped_low") <- n_capped_low
  rownames(out) <- NULL
  out
}
