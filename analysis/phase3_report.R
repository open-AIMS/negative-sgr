## Phase 3 report -- the arm comparison, and the gate the plan defines.
##
## Reads the targets pipeline's outputs; runs no models. `targets::tar_make()`
## must have completed first.

source("R/setup.R")
source("R/data_prep.R")
source("R/diagnostics.R")
source("R/figures.R")
library(targets)

ep <- tar_read(endpoints)
pars <- tar_read(parameters)
dg <- tar_read(diagnostics_tab)
bc <- tar_read(bot_contraction)
summary_tab <- utils::read.csv("analysis/dataset_summary.csv")

## Two things have to be marked before any of these numbers are read.
##
## 1. `ecx()` returns max(x_vec) when the target decline is never reached inside
##    the fitted range -- it warns, but the value it returns still looks like an
##    estimate. Arm D truncates the design, so this is exactly where it bites.
## 2. A fit that did not converge is not a result. Flagged, not silently kept.
fits <- tar_read(all_fits)
xmax <- do.call(rbind, lapply(unlist(fits, recursive = FALSE), function(a) {
  if (is.null(a$fit)) return(NULL)
  data.frame(dataset = a$dataset, arm = a$arm,
             x_min = min(a$fit$fit$data$x),
             x_max = max(a$fit$fit$data$x), stringsAsFactors = FALSE)
}))
ep <- merge(ep, xmax, by = c("dataset", "arm"), all.x = TRUE)
## An estimate or interval endpoint sitting at either end of the evaluation
## grid is a boundary artefact, not a measurement: `ecx()` searches
## `seq(min(x_range), max(x_range), length = resolution)` and returns the
## closest grid point, so a curve that never reaches the target simply returns
## the nearest edge.
##
## The previous rule tested the UPPER end only, and additionally required a
## zero-width interval. Both conditions failed on r_salina2 arm B2, whose ErC10
## came back as the grid floor with an interval spanning most of the design
## (0.01 [0.01, 135]) and was reported as `usable`. Flag either end, and do not
## require the interval to be degenerate -- a boundary-pinned point estimate
## with a wide interval is precisely the pathological case.
##
## Tolerance is relative to the design span rather than absolute, so it does not
## depend on the grid's internal spacing.
ep$design_span <- ep$x_max - ep$x_min
btol <- 1e-4 * ep$design_span
ep$at_lower <- !is.na(ep$x_min) & ep$estimate <= ep$x_min + btol
ep$at_upper <- !is.na(ep$x_max) & ep$estimate >= ep$x_max - btol
ep$at_boundary <- ep$at_lower | ep$at_upper
## An interval running to a grid edge is truncated by the search range rather
## than by the data, so the reported uncertainty is a floor, not an estimate.
ep$interval_at_boundary <- !is.na(ep$x_min) &
  (ep$lower <= ep$x_min + btol | ep$upper >= ep$x_max - btol)
## And an interval covering most of the tested range means the data barely
## constrain the endpoint, whatever the point estimate looks like.
ep$interval_spans_design <- !is.na(ep$design_span) & ep$design_span > 0 &
  (ep$upper - ep$lower) / ep$design_span > 0.8
dg_key <- dg[, c("dataset", "arm", "max_rhat", "divergences")]
ep <- merge(ep, dg_key, by = c("dataset", "arm"), all.x = TRUE)
ep$not_converged <- !is.na(ep$max_rhat) & ep$max_rhat > 1.05
ep$usable <- !ep$at_boundary & !ep$not_converged &
  !ep$interval_at_boundary & !ep$interval_spans_design

arm_order <- c("A", "A_raw_dropped", "B1", "B2", "B3", "C", "C2", "D", "SQ")
ep$arm <- factor(ep$arm, levels = arm_order)
ep <- ep[order(ep$dataset, ep$endpoint, ep$arm), ]

utils::write.csv(ep, "analysis/phase3_endpoints.csv", row.names = FALSE)
utils::write.csv(pars, "analysis/phase3_parameters.csv", row.names = FALSE)
utils::write.csv(dg, "analysis/phase3_diagnostics.csv", row.names = FALSE)
utils::write.csv(bc, "analysis/phase4_bot_contraction.csv", row.names = FALSE)

options(width = 200)
cat("\n===== ENDPOINTS (95% credible intervals) =====\n")
for (ds in dataset_names()) {
  cat("\n--", ds, "--\n")
  sub <- ep[ep$dataset == ds, ]
  for (en in c("ErC10", "ErC50", "NSEC")) {
    s <- sub[sub$endpoint == en, ]
    cat(sprintf("  %-6s ", en))
    lab <- ifelse(s$at_boundary, " [NOT ESTIMABLE in range]",
           ifelse(s$not_converged, " [DID NOT CONVERGE]",
           ifelse(s$interval_at_boundary, " [INTERVAL TRUNCATED by range]",
           ifelse(s$interval_spans_design, " [INTERVAL SPANS DESIGN]", ""))))
    cat(paste0(sprintf("%s: %.4g [%.4g, %.4g]", s$arm, s$estimate, s$lower,
                       s$upper), lab, collapse = "  |  "), "\n")
  }
}

## ------------------------------------------------------- the arm contrast ---
## Divergence from arm A, expressed as a ratio, and whether arm A's point
## estimate falls inside the other arm's interval. A ratio alone can look large
## while both arms remain compatible; the containment flag is what "the arms
## agree" means operationally.
cat("\n===== UNUSABLE RESULTS (excluded from the gate) =====\n")
unusable <- ep[!ep$usable, c("dataset", "arm", "endpoint", "estimate",
                             "lower", "upper", "at_boundary",
                             "interval_at_boundary", "interval_spans_design",
                             "not_converged", "max_rhat")]
if (nrow(unusable)) print(unusable, digits = 4, row.names = FALSE) else
  cat("none\n")

cat("\n===== DIVERGENCE FROM ARM A =====\n")
ep_ok <- ep[ep$usable, ]
div <- do.call(rbind, lapply(split(ep_ok, list(ep_ok$dataset, ep_ok$endpoint),
                                   drop = TRUE), function(s) {
  a <- s[s$arm == "A", ]
  if (nrow(a) == 0) return(NULL)
  other <- s[s$arm != "A", ]
  data.frame(
    dataset = other$dataset, endpoint = other$endpoint, arm = other$arm,
    ratio_to_A = other$estimate / a$estimate,
    A_inside_arm_interval = a$estimate >= other$lower & a$estimate <= other$upper,
    arm_inside_A_interval = other$estimate >= a$lower & other$estimate <= a$upper,
    width_ratio = (other$upper - other$lower) / (a$upper - a$lower),
    stringsAsFactors = FALSE
  )
}))
rownames(div) <- NULL
utils::write.csv(div, "analysis/phase3_divergence.csv", row.names = FALSE)
print(div, digits = 3)

## ------------------------------------------------------------- the gate -----
## The plan's pre-registered prediction: divergence orders by Delta, not by R.
## Scored on the arms that embody the convention (B1, B2, B3), on ErC10 and
## ErC50, as the mean absolute log ratio to arm A.
cat("\n===== GATE: does divergence order by Delta or by R? =====\n")
score <- div[div$arm %in% c("B1", "B2", "B3") &
               div$endpoint %in% c("ErC10", "ErC50"), ]
by_ds <- tapply(abs(log(score$ratio_to_A)), score$dataset, mean)
gate <- data.frame(dataset = names(by_ds),
                   mean_abs_log_ratio = unname(by_ds))
gate <- merge(gate, summary_tab[, c("dataset", "delta", "R")], by = "dataset")
gate <- gate[order(-gate$mean_abs_log_ratio), ]
print(gate, digits = 4, row.names = FALSE)

rank_obs <- gate$dataset
rank_delta <- summary_tab$dataset[order(-summary_tab$delta)]
rank_R <- summary_tab$dataset[order(summary_tab$R)]
cat("\nobserved divergence order:", paste(rank_obs, collapse = " > "), "\n")
cat("predicted by Delta:        ", paste(rank_delta, collapse = " > "), "\n")
cat("predicted by R:            ", paste(rank_R, collapse = " > "), "\n")
cat("\nSpearman correlation with Delta:",
    round(stats::cor(gate$mean_abs_log_ratio, gate$delta, method = "spearman"), 3),
    "\nSpearman correlation with R:    ",
    round(stats::cor(gate$mean_abs_log_ratio, gate$R, method = "spearman"), 3),
    "\n(4 datasets, so |rho| = 1 is the only decisive value and even that is",
    "weak evidence -- the simulation is what separates the two.)\n")

n_disagree <- sum(!div$A_inside_arm_interval[
  div$arm %in% c("B1", "B2", "B3") & div$endpoint %in% c("ErC10", "ErC50")])
cat("\nArm-A estimate outside the convention arms' intervals in",
    n_disagree, "of",
    sum(div$arm %in% c("B1", "B2", "B3") & div$endpoint %in% c("ErC10", "ErC50")),
    "dataset x arm x endpoint combinations.\n")
cat(if (n_disagree == 0) {
  "GATE: all arms agree within their intervals -- the honest output is a\n  vignette section and a short note, not a paper.\n"
} else {
  "GATE: at least one arm disagrees with A -- Phase 5 quantifies the magnitude.\n"
})

## -------------------------------------------------------------- Phase 4 -----
cat("\n===== bot: PRIOR-TO-POSTERIOR CONTRACTION =====\n")
bc2 <- merge(bc, summary_tab[, c("dataset", "delta", "f_neg")], by = "dataset")
print(bc2[order(bc2$dataset, bc2$arm), ], digits = 3, row.names = FALSE)
cat("\nA contraction near 0 means bot is prior-driven: the curve never",
    "establishes a lower plateau,\nso the arm-A reference is itself a prior",
    "statement rather than a data statement.\n")

cat("\n===== SAMPLER DIAGNOSTICS =====\n")
print(dg[order(dg$dataset, dg$arm), ], digits = 4, row.names = FALSE)
bad <- dg[dg$divergences > 0 | dg$max_rhat > 1.01, ]
if (nrow(bad)) {
  cat("\nFits needing attention:\n"); print(bad, row.names = FALSE)
} else {
  cat("\nNo divergences and all Rhat <= 1.01.\n")
}

## -------------------------------------------------------------- figures -----
dir.create("analysis/figures", showWarnings = FALSE, recursive = TRUE)
plots <- tar_read(arm_curve_plot)
for (i in seq_along(plots)) {
  ds <- dataset_names()[i]
  ggplot2::ggsave(file.path("analysis/figures",
                            paste0("arm_curves_", ds, ".png")),
                  plots[[i]], width = 7, height = 5, dpi = 200)
}
cat("\nFigures written to analysis/figures/\n")
