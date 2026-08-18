## Phase 5 report -- bias, RMSE and coverage with Monte Carlo standard errors.

source("R/setup.R")
source("R/data_prep.R")
source("R/metrics.R")

files <- list.files("analysis/phase5", pattern = "\\.rds$", full.names = TRUE)
if (length(files) == 0) stop("No sweep output in analysis/phase5/.")
res <- do.call(rbind, lapply(files, readRDS))
## Arm levels and grouping come from the shared palette file rather than being
## hardcoded here, so that adding the Phase 7 arms (E, F) does not require this
## script to be edited in two places. Levels are intersected with what is
## actually on disk, so the report works on a partial sweep: during Phase 7
## stage 1 only some cells carry E and F, and those arms simply appear with
## fewer cells rather than breaking the run.
source("analysis/arm_palette.R")
present <- intersect(arm_levels, unique(res$arm))
if (!setequal(present, unique(res$arm))) {
  stop("arms on disk that the palette does not name: ",
       paste(setdiff(unique(res$arm), arm_levels), collapse = ", "))
}
res$arm <- factor(res$arm, levels = present)
res$arm_group <- unname(arm_group[as.character(res$arm)])
cat("arms present:", paste(present, collapse = ", "), "\n")
cat("cells per arm:\n")
print(tapply(res$cell, res$arm, function(z) length(unique(z))))

options(width = 200)

## ----------------------------------------------------------- exclusions -----
## Reported per arm before anything is filtered. If one arm's exclusions are
## systematically different, they are not missing at random and that is itself a
## finding rather than something to clean away.
cat("\n===== EXCLUSION RATE PER ARM =====\n")
excl <- aggregate(cbind(failed = !ok) ~ arm + delta + top_factor + R,
                  res[!duplicated(res[, c("cell", "iteration", "arm")]), ],
                  mean)
print(aggregate(failed ~ arm, excl, mean), digits = 3, row.names = FALSE)
cat("\nby cell:\n")
print(excl[excl$failed > 0, ], digits = 3, row.names = FALSE)

converged <- res$ok & !is.na(res$estimate) & res$max_rhat <= 1.01 &
  res$divergences == 0
cat("\nfits excluded for divergences or Rhat > 1.01:",
    sum(res$ok & !converged, na.rm = TRUE), "of", sum(res$ok, na.rm = TRUE),
    "\n")
cat("per arm:\n")
print(round(tapply(res$ok & !converged, res$arm, mean, na.rm = TRUE), 4))

## Metrics are computed on ALL successful fits, not only the clean ones.
##
## Filtering on `divergences == 0` conditions on an OUTCOME that differs
## systematically by arm: arm A averages ~0 divergences while B2 averages 4-6,
## so the filter would retain nearly all of A and only the best-behaved B2 fits,
## flattering the arm the study finds worst. The exclusion rates above are the
## finding; applying them as a filter would hide it.
##
## Both versions are produced. `mets` is the headline (all successful fits);
## `mets_clean` is the sensitivity analysis. Where they disagree, the
## disagreement is reported rather than resolved by preferring one.
use <- res[which(res$ok & !is.na(res$estimate)), ]
use_clean <- res[which(converged), ]

## -------------------------------------------------------------- metrics -----
cat("\n===== BIAS, RMSE, COVERAGE (with MCSE) =====\n")
key <- list(use$cell, use$arm, use$endpoint)
mets <- do.call(rbind, lapply(split(use, key, drop = TRUE), function(s) {
  m <- mcse_summary(s$estimate, s$truth[1], s$lower, s$upper)
  cbind(cell = s$cell[1], arm = as.character(s$arm[1]),
        endpoint = s$endpoint[1], delta = s$delta[1],
        top_factor = s$top_factor[1], R = s$R[1],
        sigma_ratio = s$sigma_ratio[1], truth = s$truth[1], m)
}))
rownames(mets) <- NULL
mets$rel_bias <- mets$bias / mets$truth
## f_neg is the mediating quantity -- the fraction of responses the convention
## actually discards -- so every table carries it.
fn <- aggregate(f_neg ~ cell, use, mean)
mets <- merge(mets, fn, by = "cell", all.x = TRUE)
utils::write.csv(mets, "analysis/phase5_metrics.csv", row.names = FALSE)

## Sensitivity: the same metrics on clean fits only.
key_c <- list(use_clean$cell, use_clean$arm, use_clean$endpoint)
mets_clean <- do.call(rbind, lapply(split(use_clean, key_c, drop = TRUE),
  function(s) {
    m <- mcse_summary(s$estimate, s$truth[1], s$lower, s$upper)
    cbind(cell = s$cell[1], arm = as.character(s$arm[1]),
          endpoint = s$endpoint[1], m)
  }))
rownames(mets_clean) <- NULL
utils::write.csv(mets_clean, "analysis/phase5_metrics_cleanfits.csv",
                 row.names = FALSE)

cat("\n===== ALL-FITS vs CLEAN-FITS SENSITIVITY (ErC50) =====\n")
cmp <- merge(mets[mets$endpoint == "ErC50", c("cell", "arm", "bias", "n")],
             mets_clean[mets_clean$endpoint == "ErC50",
                        c("cell", "arm", "bias", "n")],
             by = c("cell", "arm"), suffixes = c(".all", ".clean"))
cmp$bias_shift <- cmp$bias.clean - cmp$bias.all
cmp$dropped <- cmp$n.all - cmp$n.clean
byarm <- aggregate(cbind(dropped, bias_shift) ~ arm, cmp, mean)
cat("Mean fits dropped per cell, and how much dropping them moves the bias:\n")
print(byarm, digits = 3, row.names = FALSE)
cat("A large `dropped` with a large `bias_shift` means the convergence filter\n",
    "is doing the arm a favour; report the all-fits number.\n")

for (en in c("ErC10", "ErC50", "NSEC")) {
  cat("\n--", en, "--\n")
  if (en == "NSEC") {
    cat("!! NSEC bias and coverage against the true `nec` are WITHDRAWN as\n",
        "   study outputs. NSEC is defined against the posterior spread of the\n",
        "   control response, so `nec` is its limit, not its expectation: at\n",
        "   negligible noise it returns 1.2985 against a true 1.3, but at\n",
        "   realistic noise coverage falls to 0 in the widest cells on EVERY\n",
        "   arm including A. Read the arm-to-arm CONTRAST only; the shared\n",
        "   design component cancels there. See SESSION.md.\n")
  }
  s <- mets[mets$endpoint == en, ]
  print(s[order(s$delta, s$top_factor, s$R, s$arm),
          c("delta", "top_factor", "R", "arm", "n", "rel_bias", "bias_mcse",
            "rmse", "coverage", "coverage_mcse")],
        digits = 3, row.names = FALSE)
}

## ----------------------------------------------------------- the regime -----
## The arm comparison is REGIME-DEPENDENT and must never be reported as a single
## pooled number.
##
## Where the design stops at or below the zero crossing, almost nothing is
## negative, `bot` is weakly identified in every arm alike, all six arms carry a
## shared bias of the same sign, and B1 can even come out ahead of A. That is
## not flooring helping: it is arm A carrying a baseline bias from an
## unidentified asymptote that flooring happens to offset. Only once the design
## identifies `bot` does A's bias vanish and the flooring bias stand alone.
## Pooling the two regimes averages a real effect against an artefact.
##
## Split on `f_neg` -- the fraction of responses the convention actually alters
## -- rather than on `top_factor`, which is only its design proxy. The threshold
## sits in an empty gap: cells run 0.003-0.042 below it and 0.148-0.152 above,
## with nothing in between, so the split is not sensitive to where in that gap
## the line is drawn.
F_NEG_REACHING <- 0.10
mets$regime <- ifelse(mets$f_neg > F_NEG_REACHING, "reaches", "stops short")

## `R` is held out of the headline. Under this parameterisation the generating
## model is exactly scale-equivariant in the growth rate, so the R cells vary
## signal-to-noise and nothing else (the sweep holds sigma absolute); averaging
## them into the arm comparison would silently weight it by noise level.
head_ln <- mets[mets$R == 2.3, ]

cat("\n===== HEADLINE: ARMS BY REGIME (R = 2.3 cells only) =====\n")
for (en in c("ErC10", "ErC50")) {
  cat("\n--", en, "--\n")
  s <- head_ln[head_ln$endpoint == en, ]
  tab <- aggregate(cbind(rel_bias, coverage, rmse, f_neg) ~ arm + regime, s,
                   mean)
  ## Merged rather than column-bound: `aggregate` on a cbind() formula drops
  ## rows where any response is NA, so a cell with one missing metric would
  ## shift the two results out of alignment and silently mislabel n_cells.
  nc <- aggregate(cell ~ arm + regime, s, function(z) length(unique(z)))
  names(nc)[names(nc) == "cell"] <- "n_cells"
  tab <- merge(tab, nc, by = c("arm", "regime"))
  tab <- tab[order(tab$regime, tab$rel_bias), ]
  print(tab, digits = 3, row.names = FALSE)
}
cat("\nRead the `reaches` block as the study's result and the `stops short`",
    "block\nas the null case. An arm ordering that only appears in `stops",
    "short` is a\nstatement about weak identification, not about the",
    "convention.\n")

## The mechanism claim: within the reaching regime the flooring penalty should
## grow with `delta`, because `delta` sets how wrong each altered point is.
## Reported per arm rather than pooled so a flat arm cannot be hidden by a
## steep one.
cat("\n===== DOES THE PENALTY GROW WITH delta, WHERE THE DESIGN REACHES? =====\n")
for (en in c("ErC10", "ErC50")) {
  cat("\n--", en, "--\n")
  s <- head_ln[head_ln$endpoint == en & head_ln$regime == "reaches", ]
  print(reshape(aggregate(rel_bias ~ delta + arm, s, mean),
                idvar = "arm", timevar = "delta", direction = "wide"),
        digits = 3, row.names = FALSE)
}

## ------------------------------------------------ the pre-registered claims --
cat("\n===== PRE-REGISTERED PREDICTIONS =====\n")

cat("\n1. Bias grows as the series extends past the zero-crossing.\n")
## "Floors" now means any arm that imposes the zero boundary, by whatever
## mechanism: by substituting data (B1), by constraining the mean (B2, B3), or
## by choosing a family whose support excludes negatives (E, F).
## unique(), not levels(): `mets$arm` is rebuilt as a character column by the
## cbind() above, so levels() is NULL here and the intersect would be empty.
FLOORING <- intersect(names(arm_group)[arm_group != "Measurement retained"],
                      unique(mets$arm))
p1 <- aggregate(rel_bias ~ top_factor + arm + endpoint,
                mets[mets$arm %in% FLOORING, ], mean)
print(reshape(p1, idvar = c("arm", "endpoint"), timevar = "top_factor",
              direction = "wide"), digits = 3, row.names = FALSE)

cat("\n2. Bias orders by Delta, not by R.\n")
byd <- aggregate(abs(rel_bias) ~ delta,
                 mets[mets$arm %in% FLOORING &
                        mets$endpoint != "NSEC", ], mean)
byr <- aggregate(abs(rel_bias) ~ R,
                 mets[mets$arm %in% FLOORING &
                        mets$endpoint != "NSEC", ], mean)
print(byd, digits = 3, row.names = FALSE)
print(byr, digits = 3, row.names = FALSE)
cat("If |bias| rises with Delta and is flat in R, the growth-rate-scale",
    "mechanism holds.\nIf it rises with R at fixed Delta, it does not.\n")
cat("\nNOTE: R is NOT an independent axis here. Under the study's\n",
    "parameterisation the generating model is exactly scale-equivariant in the\n",
    "growth rate, so R can only act through signal-to-noise, and only because\n",
    "the sweep holds the residual scale ABSOLUTE (sigma_mode = 'absolute',\n",
    "anchored at R = 2.3). Read the R rows as a noise sweep, not as a test of\n",
    "the control fold-change. See SESSION.md.\n")

cat("\n3. ErCx biased low under B2; NSEC low under B1 and high under B2.\n")
p3 <- aggregate(rel_bias ~ arm + endpoint, mets, mean)
print(reshape(p3, idvar = "endpoint", timevar = "arm", direction = "wide"),
      digits = 3, row.names = FALSE)

cat("\n4. Coverage degrades under B2/B3 more than the point estimates move.\n")
p4 <- aggregate(cbind(coverage, abs_rel_bias = abs(rel_bias)) ~ arm + endpoint,
                mets, mean)
print(p4, digits = 3, row.names = FALSE)

## The residual-scale mechanism the directions rest on.
cat("\n===== RESIDUAL SCALE BY ARM (the mechanism) =====\n")
cat("B1 compresses sigma by replacing extreme negatives with 0;\n",
    "B2 inflates it by making those points unreachable.\n",
    "Requires the parameter table from the sweep; skipped if absent.\n")
