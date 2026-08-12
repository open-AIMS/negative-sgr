## Phase 5 report -- bias, RMSE and coverage with Monte Carlo standard errors.

source("R/setup.R")
source("R/data_prep.R")
source("R/metrics.R")

files <- list.files("analysis/phase5", pattern = "\\.rds$", full.names = TRUE)
if (length(files) == 0) stop("No sweep output in analysis/phase5/.")
res <- do.call(rbind, lapply(files, readRDS))
res$arm <- factor(res$arm, levels = c("A", "B1", "B2", "B3", "C", "D"))

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

use <- res[which(converged), ]

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
utils::write.csv(mets, "analysis/phase5_metrics.csv", row.names = FALSE)

for (en in c("ErC10", "ErC50", "NSEC")) {
  cat("\n--", en, "--\n")
  s <- mets[mets$endpoint == en, ]
  print(s[order(s$delta, s$top_factor, s$R, s$arm),
          c("delta", "top_factor", "R", "arm", "n", "rel_bias", "bias_mcse",
            "rmse", "coverage", "coverage_mcse")],
        digits = 3, row.names = FALSE)
}

## ------------------------------------------------ the pre-registered claims --
cat("\n===== PRE-REGISTERED PREDICTIONS =====\n")

cat("\n1. Bias grows as the series extends past the zero-crossing.\n")
p1 <- aggregate(rel_bias ~ top_factor + arm + endpoint,
                mets[mets$arm %in% c("B1", "B2", "B3"), ], mean)
print(reshape(p1, idvar = c("arm", "endpoint"), timevar = "top_factor",
              direction = "wide"), digits = 3, row.names = FALSE)

cat("\n2. Bias orders by Delta, not by R.\n")
byd <- aggregate(abs(rel_bias) ~ delta,
                 mets[mets$arm %in% c("B1", "B2", "B3") &
                        mets$endpoint != "NSEC", ], mean)
byr <- aggregate(abs(rel_bias) ~ R,
                 mets[mets$arm %in% c("B1", "B2", "B3") &
                        mets$endpoint != "NSEC", ], mean)
print(byd, digits = 3, row.names = FALSE)
print(byr, digits = 3, row.names = FALSE)
cat("If |bias| rises with Delta and is flat in R, the growth-rate-scale",
    "mechanism holds.\nIf it rises with R at fixed Delta, it does not.\n")

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
