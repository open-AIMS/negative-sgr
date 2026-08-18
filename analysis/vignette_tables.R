## The numbers the vignette embeds, generated rather than transcribed.
##
## `example7.Rmd.orig` carries its results as inline `read.csv(text = ...)`
## blocks, because the vignette must knit in seconds without the sweep. Those
## blocks were first written by hand, which is a transcription risk every time
## an arm is added -- and Phase 7 adds two. This script prints them in exactly
## the form the vignette expects, so updating the vignette is a copy of a
## generated block rather than eight rows of manual editing.
##
## Run: Rscript analysis/vignette_tables.R > analysis/vignette_tables.txt

source("R/setup.R"); source("R/simulate.R"); source("R/data_prep.R")
source("analysis/arm_palette.R")
options(width = 200)

m <- read.csv("analysis/phase5_metrics.csv")
m$regime <- ifelse(m$f_neg > 0.10, "reaches", "stops short")

## ------------------------------------------------------- sim_results block --
## Averaged over the delta levels within a regime, R = 2.3 only, exactly as
## phase5_report.R's headline table does it: R is a noise axis here, so pooling
## it into the arm comparison would weight that comparison by noise level.
h <- m[m$R == 2.3 & m$endpoint %in% c("ErC10", "ErC50"), ]
agg <- aggregate(cbind(rel_bias, coverage, rmse) ~ arm + endpoint + regime, h,
                 mean)
agg <- agg[order(agg$endpoint, agg$regime, match(agg$arm, arm_levels)), ]
cat("### sim_results\n")
cat("arm,estimate,regime,rel_bias,coverage,rmse\n")
for (i in seq_len(nrow(agg))) {
  cat(sprintf("%s,%s,%s,%.4f,%.4f,%.4f\n", agg$arm[i], agg$endpoint[i],
              agg$regime[i], agg$rel_bias[i], agg$coverage[i], agg$rmse[i]))
}

## --------------------------------------------------------- precision block --
## The control CV is derived from the generating model rather than looked up, so
## the axis label in the vignette cannot drift away from what was simulated.
cal <- calibrate_sigma(read_sgr("c_proliferum"))
SIGMA_0 <- sigma_0_at(cal$cv_control, R_ref = 2.3, t = 7)
cv_of <- function(R) {
  tr <- sim_truth(R = R, delta = 4, t = 7)
  100 * SIGMA_0 / tr$top
}
p <- m[m$delta == 4 & m$top_factor == 2 & m$endpoint %in% c("ErC10", "ErC50"), ]
p$cv <- round(vapply(p$R, cv_of, numeric(1)), 1)
p <- p[order(p$endpoint, match(p$arm, arm_levels), p$R), ]
cat("\n### precision\n")
cat("arm,estimate,R,cv,rel_bias,coverage\n")
for (i in seq_len(nrow(p))) {
  cat(sprintf("%s,%s,%s,%.1f,%.4f,%.3f\n", p$arm[i], p$endpoint[i],
              format(p$R[i], trim = TRUE), p$cv[i], p$rel_bias[i],
              p$coverage[i]))
}

## ------------------------------------------------- numbers used in the prose --
cat("\n### delta gradient, reaching regime, R = 2.3 (prose)\n")
g <- aggregate(rel_bias ~ arm + delta + endpoint,
               m[m$R == 2.3 & m$regime == "reaches" &
                   m$endpoint %in% c("ErC10", "ErC50"), ], mean)
print(stats::reshape(g, idvar = c("arm", "endpoint"), timevar = "delta",
                     direction = "wide"), digits = 3, row.names = FALSE)

cat("\n### RMSE in the reaching regime (prose; never quote coverage alone)\n")
r <- aggregate(cbind(rmse, coverage, rel_bias) ~ arm + endpoint,
               m[m$R == 2.3 & m$regime == "reaches" &
                   m$endpoint %in% c("ErC10", "ErC50"), ], mean)
print(r[order(r$endpoint, match(r$arm, arm_levels)), ], digits = 3,
      row.names = FALSE)

## Divergences per fit, reaching cells only: the vignette quotes these to make
## the point that the family-floored arms fail silently.
cat("\n### mean divergences per fit, reaching cells (prose)\n")
files <- list.files("analysis/phase5", pattern = "\\.rds$", full.names = TRUE)
res_all <- do.call(rbind, lapply(files, readRDS))
res <- res_all[res_all$R == 2.3 & res_all$top_factor == 2 &
                 res_all$endpoint == "ErC50", ]
d <- aggregate(cbind(divergences, max_rhat) ~ arm, res, mean, na.rm = TRUE)
print(d[order(match(d$arm, arm_levels)), ], digits = 3, row.names = FALSE)

## ---------------------------------------------------- nsec_precision block --
## NSEC is reported as a CONTRAST only: scored against the true `nec` it
## measures the metric's definition and the approach's error at once (see
## phase5_report.R), so this block carries no bias or coverage column.
##
## The two columns are deliberately computed differently, which is worth stating
## because it looks like an inconsistency. `nsec` is the MEAN across iterations,
## recovered as truth + bias so it cannot drift from the metrics file the rest of
## the vignette uses. `ratio_to_A` is the MEDIAN of the per-iteration ratio,
## which is a paired contrast -- each approach against the intact analysis of the
## SAME simulated dataset -- and therefore the quantity the vignette's prose
## about the direction of the contrast is actually about. A ratio of the two
## means would drop the pairing and moves the B1 figures by about a percentage
## point.
n <- m[m$delta == 4 & m$top_factor == 2 & m$endpoint == "NSEC", ]
n$cv <- round(vapply(n$R, cv_of, numeric(1)), 1)
n$nsec <- n$truth + n$bias

nres <- res_all[res_all$endpoint == "NSEC" & res_all$delta == 4 &
                  res_all$top_factor == 2 & res_all$ok, ]
ratio <- do.call(rbind, lapply(split(nres, nres$R), function(s) {
  w <- stats::reshape(s[, c("iteration", "arm", "estimate")],
                      idvar = "iteration", timevar = "arm", direction = "wide")
  a <- w[["estimate.A"]]
  do.call(rbind, lapply(setdiff(unique(s$arm), NA), function(ar) {
    v <- w[[paste0("estimate.", ar)]]
    data.frame(R = s$R[1], arm = ar,
               ratio_to_A = stats::median(v / a, na.rm = TRUE))
  }))
}))
n <- merge(n, ratio, by = c("R", "arm"))
n <- n[order(match(n$arm, arm_levels), n$R), ]
cat("\n### nsec_precision\n")
cat("arm,cv,nsec,ratio_to_A\n")
for (i in seq_len(nrow(n))) {
  cat(sprintf("%s,%.1f,%.3f,%.3f\n", n$arm[i], n$cv[i], n$nsec[i],
              n$ratio_to_A[i]))
}
