## Phase 10 report -- did model averaging rescue the floored arms?
##
## Produces `phase10_metrics.csv` (the Phase 5 columns plus what averaging adds),
## `phase10_paired.csv` (the headline), `phase10_weights.csv` and a printed
## summary. Regenerate rather than editing numbers by hand, exactly as
## `vignette_tables.R` is used for the single-model results.
##
## THE PAIRED CONTRAST IS THE HEADLINE, not the two biases side by side. Because
## the sweep reuses Phase 8's seeds, iteration i of a cell is the SAME simulated
## dataset under both workflows, so `estimate_averaged - estimate_single` is a
## within-dataset difference and its standard error is far smaller than that of
## the difference of two independently estimated biases. Reporting the two biases
## alone would throw that away.
##
## Run: Rscript analysis/phase10_report.R > analysis/phase10_report.log

source("R/setup.R"); source("R/setup_phase10.R")
source("R/data_prep.R"); source("R/simulate.R"); source("R/metrics.R")
source("R/model_average.R")

IN  <- Sys.getenv("IN", "analysis/phase10")
OUT <- "analysis"

files <- list.files(IN, pattern = "\\.rds$", full.names = TRUE)
if (!length(files)) stop("No Phase 10 output in ", IN)
blocks <- lapply(files, readRDS)
ends_all <- do.call(rbind, lapply(blocks, `[[`, "endpoints"))
wts_all  <- do.call(rbind, lapply(blocks, `[[`, "weights"))
rhat     <- do.call(rbind, lapply(blocks, `[[`, "rhat"))

## The PRIMARY result is the converged stage: models that failed the R-hat
## cutoff are dropped and the weights re-stacked over the survivors, which is
## what the workflow does. The all-models stage is the sensitivity to that
## choice and is reported separately, never mixed in.
##
## The two stages also answer two different comparison questions against the
## single-model arms, and the distinction is worth stating rather than glossing.
## The single-model sweep applied NO R-hat screen -- with one model there is
## nothing to drop, and Trap 6 kept every fit. So:
##
##   STAGE = "converged"  -- screened average vs unscreened single model. The
##     realistic contrast: what each workflow actually delivers. It is not
##     like-for-like on screening, and the averaged arm has an advantage there.
##   STAGE = "all_models" -- unscreened average vs unscreened single model. The
##     like-for-like contrast, which isolates averaging from screening.
##
## Report the first as the headline and the second as the check. If they
## disagree, the difference IS the screening, and that is a finding about the
## workflow rather than a nuisance.
STAGE <- Sys.getenv("STAGE", "converged")
ends <- ends_all[ends_all$stage == STAGE, ]
wts  <- wts_all[wts_all$stage == STAGE, ]
cat("stage:", STAGE, "(set STAGE=all_models for the sensitivity)\n")

cat("Phase 10 report |", length(files), "block files |",
    length(unique(ends$cell_index)), "cells |",
    length(unique(ends$iteration)), "iterations\n\n")

## ---------------------------------------------------- convergence screening --
## How much work the R-hat rule actually did. This is a finding, not bookkeeping:
## arms E and F carry twelve candidates against the Gaussian arms' eight and six,
## so they have more opportunities to lose one, and an arm that routinely
## averages over half its set is a different analysis from one that keeps all of
## it.
if (!is.null(rhat)) {
  cat("\n-- models dropped at rhat >", P10_RHAT_CUTOFF, "--\n")
  drop_rate <- do.call(rbind, lapply(split(rhat, ~ cell + arm, drop = TRUE),
    function(d) data.frame(
      cell = d$cell[1], arm = d$arm[1],
      candidates = length(unique(d$model)),
      mean_dropped = round(mean(tapply(d$dropped, d$iteration, sum)), 2),
      pct_iterations_with_a_drop =
        round(100 * mean(tapply(d$dropped, d$iteration, any)), 1),
      row.names = NULL)))
  drop_rate <- drop_rate[order(drop_rate$cell,
                               match(drop_rate$arm, P10_ARMS)), ]
  print(drop_rate, row.names = FALSE)

  ## Which models fail, because "averaging dropped four models" and "averaging
  ## dropped the true model" are very different statements.
  cat("\n-- proportion of iterations each model was dropped --\n")
  pm <- aggregate(dropped ~ arm + model, data = rhat, FUN = mean)
  pm <- pm[pm$dropped > 0, ]
  pm <- pm[order(pm$arm, -pm$dropped), ]
  if (nrow(pm)) {
    print(transform(pm, dropped = round(dropped, 3)), row.names = FALSE)
  } else {
    cat("none\n")
  }
  n4 <- rhat[rhat$model == "nec4param", ]
  if (nrow(n4)) {
    cat("\nnec4param (the generating model) dropped in",
        round(100 * mean(n4$dropped), 1), "% of arm-iterations\n")
  }

  ## Arms where every model failed produce no estimate at all. That has to be
  ## visible: it is a silent loss of the arm, not a neutral event.
  lost <- ends_all[ends_all$stage == "converged" & !ends_all$ok &
                     !is.na(ends_all$n_kept) & ends_all$n_kept == 0, ]
  cat("arm-iterations with NO surviving model:", nrow(lost), "\n")
}

## ------------------------------------ exclusions, and common support --------
##
## THE ONE PLACE THE R-HAT RULE CAN REPRODUCE TRAP 6. Two cases, and only the
## second is a selection.
##
##  * PARTIAL DROP -- some models fail, at least one survives. The arm still
##    returns an estimate for that dataset, every dataset is scored, and what
##    changed is the estimator, uniformly. That is the treatment, not a filter.
##
##  * TOTAL DROP -- every model fails, so the arm returns NOTHING for that
##    dataset. If that happens at different rates across arms, each arm's
##    coverage is computed over a different subset of datasets, selected on
##    something correlated with how hard the dataset was. That is exactly the
##    mechanism Trap 6 forbids: the filter flatters whichever arm behaves worst.
##    Outright fit failures (bnec could not fit any model) do the same thing and
##    already existed in Phase 5/7/8, where arm C lost 5 of 240 in one cell.
##
## Handled in three layers rather than one:
##  1. the `all_models` stage never screens, so it is structurally free of
##     R-hat-induced missingness and is the check on all of this;
##  2. the headline metrics below are computed on COMMON SUPPORT -- only those
##     iterations where EVERY arm returned an estimate -- so no arm is scored on
##     a subset another arm was not;
##  3. the per-arm exclusion rate is reported next to coverage, always, so the
##     cost of layer 2 is visible rather than implied.
##
## What common support costs, stated plainly: if one arm fails on the hardest
## datasets, those datasets leave the comparison for every arm, so the ABSOLUTE
## bias and coverage of all arms improve slightly. The arm-to-arm COMPARISON,
## which is what this study is about, stays valid. `available` columns below give
## the unrestricted per-arm numbers so the difference can be seen.
classify <- function(d) {
  ifelse(d$ok, "ok",
         ifelse(!is.na(d$n_kept) & d$n_kept == 0, "all_models_dropped",
                "fit_or_endpoint_failure"))
}
ends_all$outcome <- classify(ends_all)
excl <- do.call(rbind, lapply(
  split(ends_all[ends_all$endpoint %in% "ErC50" | is.na(ends_all$endpoint), ],
        ~ cell + arm + stage, drop = TRUE),
  function(d) data.frame(
    cell = d$cell[1], arm = d$arm[1], stage = d$stage[1],
    n = length(unique(d$iteration)),
    n_all_dropped = sum(d$outcome == "all_models_dropped"),
    n_other_failure = sum(d$outcome == "fit_or_endpoint_failure"),
    pct_excluded = round(100 * mean(d$outcome != "ok"), 2),
    row.names = NULL)))
cat("\n-- per-arm exclusion (an arm returning NO estimate for a dataset) --\n")
print(excl[order(excl$stage, excl$cell, match(excl$arm, P10_ARMS)), ],
      row.names = FALSE)

## The number that decides whether any of this matters: how far apart the arms'
## exclusion rates are within a cell and stage. Equal rates are not a Trap 6
## problem however high they are; unequal rates are, however low.
spread <- do.call(rbind, lapply(split(excl, ~ cell + stage, drop = TRUE),
  function(d) data.frame(cell = d$cell[1], stage = d$stage[1],
                         min_pct = min(d$pct_excluded),
                         max_pct = max(d$pct_excluded),
                         spread_pts = round(max(d$pct_excluded) -
                                              min(d$pct_excluded), 2),
                         row.names = NULL)))
cat("\n-- spread in exclusion rate across arms (the Trap 6 quantity) --\n")
print(spread, row.names = FALSE)
if (any(spread$spread_pts > 2)) {
  cat("\n*** Exclusion rates differ across arms by more than 2 points in at\n",
      "*** least one cell. Report the COMMON-SUPPORT metrics, not the\n",
      "*** available-case ones, and say so in the text.\n")
}

## ------------------------------------------------------------- exclusions ----
## Recorded, never filtered. The Phase 5 rule holds and bites harder here: a set
## of twelve will show worse aggregate diagnostics than a set of six for
## arithmetic reasons alone, so a convergence filter would preferentially punish
## E and F -- the arms whose behaviour is under test.
fails <- ends[!ends$ok, ]
cat("-- failed fits (recorded, not excluded) --\n")
if (nrow(fails)) {
  print(table(fails$cell_index, fails$arm))
} else {
  cat("none\n")
}
## Two different losses, kept separate on purpose. `n_requested - n_fitted` is
## models bnec() could not fit at all (a lost chain, a failed initialisation);
## `n_fitted - n_kept` is models that fitted but failed the R-hat screen. The
## `dropped` column on the weights table conflates them, because it is just
## set-difference; these three counts are what disambiguates.
cat("\n-- candidate set: requested, fitted, kept after the rhat screen --\n")
setsz <- aggregate(cbind(requested = n_requested, fitted = n_fitted,
                         kept = n_kept) ~ arm + cell_index,
                   data = ends[ends$ok, ], FUN = function(z) round(mean(z), 2))
setsz$failed_to_fit <- round(setsz$requested - setsz$fitted, 2)
setsz$failed_rhat <- round(setsz$fitted - setsz$kept, 2)
print(setsz[order(setsz$cell_index, match(setsz$arm, P10_ARMS)), ],
      row.names = FALSE)

## --------------------------------------------------------------- metrics ----
avail <- ends[ends$ok & !is.na(ends$estimate), ]

## Common support: keep only (cell, endpoint, iteration) triples for which EVERY
## arm returned an estimate, so that no arm is scored on datasets another arm
## was not. See the note above.
key <- paste(avail$cell, avail$endpoint, avail$iteration, sep = "|")
n_arms_present <- tapply(avail$arm, key, function(z) length(unique(z)))
complete_keys <- names(n_arms_present)[n_arms_present == length(P10_ARMS)]
ok <- avail[key %in% complete_keys, ]
cat(sprintf("\ncommon support: %d of %d (cell, endpoint, iteration) triples ",
            length(complete_keys), length(unique(key))))
cat(sprintf("complete across all %d arms (%.1f%% dropped)\n", length(P10_ARMS),
            100 * (1 - length(complete_keys) / length(unique(key)))))

metrics <- do.call(rbind, lapply(split(ok, ~ cell_index + arm + endpoint, drop = TRUE),
  function(d) {
    m <- mcse_summary(d$estimate, d$truth[1], d$lower, d$upper)
    cbind(cell = d$cell[1], cell_index = d$cell_index[1], arm = d$arm[1],
          endpoint = d$endpoint[1], delta = d$delta[1], R = d$R[1],
          truth = d$truth[1], m, rel_bias = m$bias / d$truth[1],
          n_fitted = mean(d$n_fitted), row.names = NULL)
  }))
metrics$workflow <- "averaged"
metrics$support <- "common"

## The same quantities without the common-support restriction, so the cost of
## imposing it is visible rather than asserted.
metrics_avail <- do.call(rbind, lapply(
  split(avail, ~ cell_index + arm + endpoint, drop = TRUE),
  function(d) {
    m <- mcse_summary(d$estimate, d$truth[1], d$lower, d$upper)
    cbind(cell = d$cell[1], cell_index = d$cell_index[1], arm = d$arm[1],
          endpoint = d$endpoint[1], delta = d$delta[1], R = d$R[1],
          truth = d$truth[1], m, rel_bias = m$bias / d$truth[1],
          n_fitted = mean(d$n_fitted), workflow = "averaged",
          support = "available", row.names = NULL)
  }))
metrics <- rbind(metrics, metrics_avail)
utils::write.csv(metrics, file.path(OUT, "phase10_metrics.csv"), row.names = FALSE)
metrics <- metrics[metrics$support == "common", ]

## --------------------------------------------- the same cells, single model --
## Read from the Phase 5/7/8 store, restricted to the arms and cells Phase 10
## fitted, so the comparison is like with like.
single <- do.call(rbind, lapply(unique(ends$cell), function(cl) {
  fs <- sweep_files(cl)
  d <- do.call(rbind, lapply(fs, readRDS))
  d[d$arm %in% P10_ARMS & d$ok & !is.na(d$estimate), ]
}))

## Restrict BOTH workflows to the iterations Phase 10 actually ran, or the
## comparison would set 100 averaged iterations against 500 single-model ones and
## the difference in Monte Carlo error would masquerade as a difference in the
## workflow.
iters_run <- sort(unique(ok$iteration))
single <- single[single$iteration %in% iters_run, ]
cat("\ncomparing on", length(iters_run), "shared iterations per cell\n")

single_metrics <- do.call(rbind, lapply(split(single, ~ cell + arm + endpoint, drop = TRUE),
  function(d) {
    m <- mcse_summary(d$estimate, d$truth[1], d$lower, d$upper)
    cbind(cell = d$cell[1], arm = d$arm[1], endpoint = d$endpoint[1],
          truth = d$truth[1], m, rel_bias = m$bias / d$truth[1], row.names = NULL)
  }))
single_metrics$workflow <- "single"

## ------------------------------------------------------- the paired contrast -
## `estimate_averaged - estimate_single` on the SAME dataset. Summarised by the
## median of the per-iteration difference with a bootstrap interval, rather than
## the mean: an averaged fit that lands on a badly-fitting model occasionally
## returns an extreme endpoint, and the question is where the bulk moved.
set.seed(20260820)
boot_ci <- function(z, B = 2000) {
  if (length(z) < 5) return(c(NA_real_, NA_real_))
  bs <- replicate(B, stats::median(sample(z, length(z), replace = TRUE)))
  unname(stats::quantile(bs, c(0.025, 0.975)))
}
key <- c("cell", "arm", "endpoint", "iteration")
m <- merge(ok[, c(key, "estimate", "lower", "upper", "truth")],
           single[, c(key, "estimate", "lower", "upper")],
           by = key, suffixes = c("_avg", "_sgl"))
paired <- do.call(rbind, lapply(split(m, ~ cell + arm + endpoint, drop = TRUE),
  function(d) {
    dif <- d$estimate_avg - d$estimate_sgl
    ci <- boot_ci(dif)
    ## Coverage moves too, and it is the quantity the question is really about.
    cov_a <- mean(d$truth >= d$lower_avg & d$truth <= d$upper_avg)
    cov_s <- mean(d$truth >= d$lower_sgl & d$truth <= d$upper_sgl)
    ## Interval width, because coverage without it is unreadable: B2 in the
    ## single-model study reached coverage 1.000 with the worst RMSE in the set.
    data.frame(cell = d$cell[1], arm = d$arm[1], endpoint = d$endpoint[1],
               n = nrow(d), truth = d$truth[1],
               median_diff = stats::median(dif), lo = ci[1], hi = ci[2],
               rel_median_diff = stats::median(dif) / d$truth[1],
               abs_err_avg = stats::median(abs(d$estimate_avg - d$truth)),
               abs_err_sgl = stats::median(abs(d$estimate_sgl - d$truth)),
               coverage_avg = cov_a, coverage_sgl = cov_s,
               width_avg = stats::median(d$upper_avg - d$lower_avg),
               width_sgl = stats::median(d$upper_sgl - d$lower_sgl),
               row.names = NULL)
  }))
utils::write.csv(paired, file.path(OUT, "phase10_paired.csv"), row.names = FALSE)

for (ep in c("ErC10", "ErC50", "NSEC")) {
  cat("\n=== ", ep, " -- single vs averaged, paired on the same datasets ===\n", sep = "")
  z <- paired[paired$endpoint == ep, ]
  z <- z[order(z$cell, match(z$arm, P10_ARMS)), ]
  print(data.frame(
    cell = z$cell, arm = z$arm,
    `median shift` = round(z$rel_median_diff * 100, 1),
    `shift 95% CI` = paste0("[", round(z$lo, 3), ", ", round(z$hi, 3), "]"),
    `|err| sgl` = round(z$abs_err_sgl, 3),
    `|err| avg` = round(z$abs_err_avg, 3),
    `cov sgl` = z$coverage_sgl, `cov avg` = z$coverage_avg,
    `width sgl` = round(z$width_sgl, 2), `width avg` = round(z$width_avg, 2),
    check.names = FALSE), row.names = FALSE)
}

## ---------------------------------------------------------- bias/coverage ----
cat("\n=== relative bias and coverage, both workflows ===\n")
cmp <- merge(
  metrics[, c("cell", "arm", "endpoint", "rel_bias", "coverage", "rmse")],
  single_metrics[, c("cell", "arm", "endpoint", "rel_bias", "coverage", "rmse")],
  by = c("cell", "arm", "endpoint"), suffixes = c("_avg", "_sgl"))
## Exclusion rate travels WITH coverage, in the same table, never in a separate
## one that can be read past. A coverage figure computed after an arm has
## silently dropped some of its datasets is not comparable with one computed
## after it dropped none, and the reader has to be able to see that at a glance.
## Same discipline the study already applies to coverage and interval width:
## B2 reached coverage 1.000 in the single-model study with the worst RMSE in
## the set, which is why neither number is ever reported alone.
cmp <- merge(cmp, excl[excl$stage == STAGE,
                       c("cell", "arm", "pct_excluded", "n_all_dropped")],
             by = c("cell", "arm"), all.x = TRUE)
cmp <- cmp[order(cmp$endpoint, cmp$cell, match(cmp$arm, P10_ARMS)), ]
print(cmp, row.names = FALSE, digits = 3)
cat("\npct_excluded is the share of datasets for which this arm returned no\n",
    "estimate at all (every model dropped, or the fit failed). Coverage above\n",
    "is on COMMON SUPPORT, so it is computed on the same datasets for every\n",
    "arm; pct_excluded says how much each arm contributed to that restriction.\n")

## ---------------------------------------------------------------- weights ----
## The mechanism. An averaged estimate is uninterpretable without this: two arms
## can report the same ErC50 from entirely different mixtures.
## A model that was dropped in an iteration contributed weight ZERO to that
## iteration's average, so its mean must count those iterations rather than skip
## them. Averaging only over the rows where a model appears overstates exactly
## the models the R-hat rule removes -- which are the ones under discussion.
## Completed against the full (arm x model x iteration) grid before averaging.
grid <- do.call(rbind, lapply(split(wts, ~ cell + arm, drop = TRUE), function(d) {
  expand.grid(cell = d$cell[1], arm = d$arm[1],
              model = unique(d$model), iteration = unique(d$iteration),
              stringsAsFactors = FALSE)
}))
wts_full <- merge(grid, wts[, c("cell", "arm", "model", "iteration", "wi")],
                  by = c("cell", "arm", "model", "iteration"), all.x = TRUE)
wts_full$wi[is.na(wts_full$wi)] <- 0
mw <- aggregate(wi ~ cell + arm + model, data = wts_full, FUN = mean)
utils::write.csv(mw[order(mw$cell, mw$arm, -mw$wi), ],
                 file.path(OUT, "phase10_weights.csv"), row.names = FALSE)
cat("\n=== mean stacking weight by arm (top 4 models) ===\n")
for (cl in unique(mw$cell)) for (a in P10_ARMS) {
  z <- mw[mw$cell == cl & mw$arm == a, ]
  if (!nrow(z)) next
  z <- utils::head(z[order(-z$wi), ], 4)
  cat(sprintf("  %-22s %-2s  %s\n", cl, a,
              paste(sprintf("%s %.2f", z$model, z$wi), collapse = "  ")))
}

## The headline number the hypothesis names explicitly.
cat("\n=== weight on nec4param, the generating model ===\n")
n4 <- mw[mw$model == "nec4param", ]
print(n4[order(n4$cell, match(n4$arm, P10_ARMS)), c("cell", "arm", "wi")],
      row.names = FALSE, digits = 3)
cat("\nA small weight here is what the hypothesis predicted, and is NOT itself\n",
    "evidence either way: the question is whether the averaged ESTIMATE is\n",
    "closer to the truth, not whether the true equation was selected.\n")
