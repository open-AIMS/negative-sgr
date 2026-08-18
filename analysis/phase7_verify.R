## Phase 7 stage 1 -- verify the artefacts before anything is reported.
##
## Trap 9: check the artefact, not the exit code. A sweep has previously
## reported [done] on a cell where 239 of 240 iterations had failed, so every
## cell is inspected here and the report is regenerated only if all twelve pass.
##
## Two things are checked that the [done] lines cannot show:
##
## 1. The stage-1 blocks together hold 240 iterations x 2 arms x 3 endpoints =
##    1440 rows per cell, with both arms present in every block.
## 2. The twelve original cell files are UNTOUCHED -- 240 iterations, six arms,
##    and the row count implied by their recorded failures. Stage 1 is additive
##    by design; this is what proves it was.
##
## Run: Rscript analysis/phase7_verify.R  (add REPORT=0 to verify only)

options(width = 200)
REPORT <- !identical(Sys.getenv("REPORT", "1"), "0")
OUT <- "analysis/phase5"

cells <- c(sprintf("d%.1f_t%.1f_R2.3_s8.1",
                   rep(c(2, 4, 8), times = 3), rep(c(0.8, 1.0, 2.0), each = 3)),
           sprintf("d4.0_t2.0_R%.1f_s8.1", c(3.3, 17, 73)))

ok_all <- TRUE
summ <- list()
for (cl in cells) {
  blocks <- list.files(OUT, pattern = sprintf("^%s__ef_i.*\\.rds$", cl),
                       full.names = TRUE)
  base <- file.path(OUT, paste0(cl, ".rds"))

  ## -- stage 1 (arms E, F) --
  if (length(blocks) == 0) {
    cat("MISSING: no stage-1 blocks for", cl, "\n"); ok_all <- FALSE; next
  }
  r <- do.call(rbind, lapply(blocks, readRDS))
  n_it <- length(unique(r$iteration)); arms <- sort(unique(r$arm))
  bad_block <- vapply(blocks, function(p) {
    b <- readRDS(p); !setequal(unique(b$arm), c("E", "F"))
  }, logical(1))

  ## -- the untouched six-arm cell --
  if (!file.exists(base)) {
    cat("MISSING:", base, "\n"); ok_all <- FALSE; next
  }
  b6 <- readRDS(base)
  n_it6 <- length(unique(b6$iteration)); n_arm6 <- length(unique(b6$arm))
  ## A FAILED fit writes ONE row carrying the error, not three endpoint rows,
  ## so the expected row count is not a flat 4320: it is 240 x 6 fits, of which
  ## the failures contribute one row each. Only `d4.0_t0.8_R2.3` is affected --
  ## arm C, 5 of 240 -- and hardcoding 4320 flags that known-good cell as
  ## incomplete. Deriving the expectation from the recorded failures keeps the
  ## check strict: it still catches a cell that silently lost rows.
  n_fail6 <- sum(!b6$ok)
  exp_rows6 <- 3L * (240L * 6L - n_fail6) + n_fail6

  pass <- n_it == 240 && setequal(arms, c("E", "F")) && nrow(r) == 1440 &&
    !any(bad_block) && n_it6 == 240 && n_arm6 == 6 && nrow(b6) == exp_rows6
  ok_all <- ok_all && pass

  summ[[cl]] <- data.frame(
    cell = cl, blocks = length(blocks), ef_rows = nrow(r), ef_iter = n_it,
    ef_arms = paste(arms, collapse = ","),
    ef_failed = sum(!r$ok), ef_na = sum(is.na(r$estimate)),
    ef_div = round(mean(r$divergences, na.rm = TRUE), 2),
    base_rows = nrow(b6), base_exp = exp_rows6, base_fail = n_fail6,
    base_arms = n_arm6,
    f_neg = round(mean(r$f_neg), 4),
    pass = pass, stringsAsFactors = FALSE)
}
tab <- do.call(rbind, summ)
print(tab, row.names = FALSE)

cat("\nexclusions and divergences per arm, stage 1:\n")
ef <- do.call(rbind, lapply(list.files(OUT, pattern = "__ef_i.*\\.rds$",
                                       full.names = TRUE), readRDS))
print(data.frame(
  arm = c("E", "F"),
  fits = c(sum(ef$arm == "E" & ef$endpoint == "ErC50"),
           sum(ef$arm == "F" & ef$endpoint == "ErC50")),
  failed = tapply(!ef$ok, ef$arm, sum),
  na_estimate = tapply(is.na(ef$estimate), ef$arm, sum),
  mean_div = round(tapply(ef$divergences, ef$arm, mean, na.rm = TRUE), 3),
  max_rhat = round(tapply(ef$max_rhat, ef$arm, max, na.rm = TRUE), 4)),
  row.names = FALSE)

if (!ok_all) {
  cat("\n*** One or more cells are incomplete. phase5_metrics.csv NOT",
      "regenerated;\n    the existing six-arm file is left untouched. Re-issue",
      "the phase 7 run --\n    completed blocks are skipped.\n")
  quit(status = 1)
}
cat("\nAll 12 cells complete: 8 arms x 240 iterations everywhere.\n")

if (REPORT) {
  cat("\n=== regenerating the report over all eight arms ===\n")
  st <- system("Rscript analysis/phase5_report.R > analysis/phase5_report.log 2>&1")
  cat("phase5_report.R exit status:", st, "\n")
  if (st != 0) {
    cat("*** report FAILED; see analysis/phase5_report.log\n"); quit(status = 1)
  }

  ## The R-axis read-out, rebuilt from the new metrics for the same reason
  ## overnight_finish.R built it: it is the study's strongest result and should
  ## not have to be recomputed by hand in the next session.
  m <- read.csv("analysis/phase5_metrics.csv")
  s <- m[m$delta == 4 & m$top_factor == 2 &
           m$endpoint %in% c("ErC10", "ErC50"), ]
  w <- stats::reshape(s[, c("R", "arm", "endpoint", "rel_bias")],
                      idvar = c("arm", "endpoint"), timevar = "R",
                      direction = "wide")
  cov <- stats::reshape(s[, c("R", "arm", "endpoint", "coverage")],
                        idvar = c("arm", "endpoint"), timevar = "R",
                        direction = "wide")
  sink("analysis/phase5_r_axis.txt")
  cat("R axis at delta = 4, top_factor = 2, ALL EIGHT ARMS. R rises => control\n",
      "CV falls, so these columns are a NOISE sweep, not a test of control\n",
      "fold-change: the generating model is exactly scale-equivariant in the\n",
      "growth rate.\n\n")
  cat("--- relative bias ---\n")
  print(w[order(w$endpoint, w$arm), ], digits = 3, row.names = FALSE)
  cat("\n--- coverage ---\n")
  print(cov[order(cov$endpoint, cov$arm), ], digits = 3, row.names = FALSE)
  cat("\n--- f_neg and sigma_ratio (should be near-constant across R) ---\n")
  print(unique(s[, c("R", "f_neg", "sigma_ratio")]), digits = 4,
        row.names = FALSE)
  cat("\nWhat to look for: arms A/C/D collapsing toward zero bias as noise\n",
      "falls while the flooring arms approach a NON-ZERO asymptote. A quantity\n",
      "that does not vanish as noise vanishes is misspecification, not\n",
      "estimation error. E and F are new to this table: they floor by family\n",
      "rather than by substitution, so they belong with B1/B2/B3 in the\n",
      "reading, not with A/C/D.\n")
  sink()
  cat("wrote analysis/phase5_r_axis.txt\n")
  cat(readLines("analysis/phase5_r_axis.txt"), sep = "\n")
}
cat("\nDone", format(Sys.time()), "\n")
