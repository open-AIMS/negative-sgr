## Phase 9 report -- does the Phase 3 conclusion survive model averaging?
##
## Phase 3 fixed every arm at `nec4param`, which is right in the simulation (it
## is the generating model) and weak on real data, where the functional form is
## unknown and `bnec()`'s model-averaged workflow is what an analyst would use.
## Phase 9 refits the case studies that way. This report asks three questions,
## in order of what they would change:
##
## 1. Does the CONVENTION still change the answer? That is the study's claim,
##    and if averaging absorbed the distortion the claim would be weaker on real
##    data than the simulation implies.
## 2. Does averaging move the answer relative to the fixed model? A large move
##    would say the Phase 3 numbers were model-choice artefacts.
## 3. WHERE does flooring act under averaging -- on the estimate within a shape,
##    or by shifting weight onto different shapes? This question does not exist
##    in Phase 3 and is the one genuinely new thing here.
##
## Reads `phase9_*.csv`; runs no models. The usability rules are Phase 3's,
## reproduced so the two are scored the same way.

source("R/setup.R"); source("R/data_prep.R")
options(width = 200)

ep <- utils::read.csv("analysis/phase9_endpoints.csv")
wt <- utils::read.csv("analysis/phase9_weights.csv")
dg <- utils::read.csv("analysis/phase9_diagnostics.csv")
p3 <- utils::read.csv("analysis/phase3_endpoints.csv")

## ------------------------------------------------------------ usability -----
## Phase 3's rules, unchanged: `ecx()` returns a grid edge when the target
## decline is never reached, and a fit that did not converge is not a result.
## The design range is taken from the Phase 3 endpoints table for the same
## dataset, since the data a given arm sees are identical in both phases.
rng <- unique(p3[, c("dataset", "arm", "x_min", "x_max")])
ep <- merge(ep, rng, by = c("dataset", "arm"), all.x = TRUE)
ep$design_span <- ep$x_max - ep$x_min
btol <- 1e-4 * ep$design_span
ep$at_boundary <- (!is.na(ep$x_min) & ep$estimate <= ep$x_min + btol) |
                  (!is.na(ep$x_max) & ep$estimate >= ep$x_max - btol)
ep$interval_at_boundary <- !is.na(ep$x_min) &
  (ep$lower <= ep$x_min + btol | ep$upper >= ep$x_max - btol)
ep$interval_spans_design <- !is.na(ep$design_span) & ep$design_span > 0 &
  (ep$upper - ep$lower) / ep$design_span > 0.8

## Convergence, of the fit AS REPORTED.
##
## `phase9_diagnostics.csv` describes the fit BEFORE the R-hat rule -- it is the
## full picture, including the models the rule removed. The reported analysis is
## the fit after removal, so the two have to be joined rather than read
## interchangeably: take the models that survive (they are the ones carrying a
## weight in `phase9_weights.csv`) and attach each one's own R-hat. A model's
## R-hat does not change when a different model is dropped, so this is exact and
## needs no refitting.
##
## An earlier version of this report printed the pre-drop table under a heading
## claiming it was the reported fit, and so showed R-hat up to 2.72 in an
## analysis that had already excluded those models.
dgn_reported <- merge(wt[, c("dataset", "arm", "model", "wi")],
                      dg[, c("dataset", "arm", "model", "max_rhat",
                             "divergences")],
                      by = c("dataset", "arm", "model"), all.x = TRUE)
conv <- do.call(rbind, lapply(split(dgn_reported,
                                    list(dgn_reported$dataset, dgn_reported$arm),
                                    drop = TRUE),
  function(s) data.frame(dataset = s$dataset[1], arm = s$arm[1],
    n_models = nrow(s),
    max_rhat = max(s$max_rhat, na.rm = TRUE),
    wtd_rhat = sum(s$max_rhat * s$wi, na.rm = TRUE) / sum(s$wi, na.rm = TRUE),
    divergences = sum(s$divergences, na.rm = TRUE),
    over_rule = sum(s$max_rhat > 1.01, na.rm = TRUE),
    stringsAsFactors = FALSE)))
rownames(conv) <- NULL

## What the rule removed, and how much weight it had been carrying before it
## was removed. That is the number that says whether the rule mattered.
dropped <- do.call(rbind, lapply(split(dg, list(dg$dataset, dg$arm),
                                       drop = TRUE), function(s) {
  bad <- s[!is.na(s$max_rhat) & s$max_rhat > 1.01, ]
  if (!nrow(bad)) return(NULL)
  data.frame(dataset = s$dataset[1], arm = s$arm[1],
             dropped = paste(bad$model, collapse = " "),
             n_dropped = nrow(bad),
             weight_held = sum(bad$wi, na.rm = TRUE),
             worst_rhat = max(bad$max_rhat), stringsAsFactors = FALSE)
}))
rownames(dropped) <- NULL

ep <- merge(ep, conv[, c("dataset", "arm", "wtd_rhat", "max_rhat")],
            by = c("dataset", "arm"), all.x = TRUE)
ep$not_converged <- !is.na(ep$max_rhat) & ep$max_rhat > 1.01
ep$usable <- !ep$at_boundary & !ep$not_converged &
  !ep$interval_at_boundary & !ep$interval_spans_design

cat("===== CONVERGENCE OF THE FITS AS REPORTED (after the R-hat rule) =====\n")
cat("Every model still in an average must be at or below R-hat 1.01; anything\n",
    "in `over_rule` other than zero is a bug in the rule, not a result.\n")
print(conv[order(-conv$max_rhat), ], digits = 3, row.names = FALSE)

cat("\n===== WHAT THE R-HAT RULE REMOVED =====\n")
cat("`weight_held` is the stacking weight those models carried BEFORE removal.\n",
    "Where it is large the reported estimate would otherwise have rested on\n",
    "chains that had not mixed.\n")
if (!is.null(dropped)) {
  print(dropped[order(-dropped$weight_held), ], digits = 3, row.names = FALSE)
} else cat("nothing dropped\n")

cat("\n===== UNUSABLE RESULTS (excluded from the gate) =====\n")
un <- ep[!ep$usable, c("dataset", "arm", "endpoint", "estimate", "lower",
                       "upper", "at_boundary", "interval_at_boundary",
                       "interval_spans_design", "not_converged")]
if (nrow(un)) print(un, digits = 4, row.names = FALSE) else cat("none\n")

## ------------------------------------------------- 1. the convention gate ---
cat("\n===== 1. DOES THE CONVENTION STILL CHANGE THE ANSWER? =====\n")
ok <- ep[ep$usable, ]
div <- do.call(rbind, lapply(split(ok, list(ok$dataset, ok$endpoint),
                                   drop = TRUE), function(s) {
  a <- s[s$arm == "A", ]
  if (nrow(a) == 0) return(NULL)
  o <- s[s$arm != "A", ]
  if (nrow(o) == 0) return(NULL)
  data.frame(dataset = o$dataset, endpoint = o$endpoint, arm = o$arm,
             ratio_to_A = o$estimate / a$estimate,
             A_inside_arm_interval = a$estimate >= o$lower & a$estimate <= o$upper,
             width_ratio = (o$upper - o$lower) / (a$upper - a$lower),
             stringsAsFactors = FALSE)
}))
rownames(div) <- NULL
utils::write.csv(div, "analysis/phase9_divergence.csv", row.names = FALSE)
print(div[order(div$endpoint, div$arm, div$dataset), ], digits = 3,
      row.names = FALSE)

sel <- div$arm %in% c("B1", "B2", "B3") & div$endpoint %in% c("ErC10", "ErC50")
cat("\nGATE, scored exactly as Phase 3 scored it -- arm A's estimate outside\n",
    "the convention arms' intervals in", sum(!div$A_inside_arm_interval[sel]),
    "of", sum(sel), "usable combinations.\n")
cat("Phase 3, fixed nec4param: 13 of 18.\n")

## ------------------------------------------- 2. averaging against fixing ----
cat("\n===== 2. WHAT DOES AVERAGING ITSELF MOVE? =====\n")
cat("Same data, same arm, same prior; the only change is a candidate set in\n",
    "place of an assumed model. A large move here would mean the Phase 3\n",
    "numbers were partly an artefact of choosing nec4param.\n")
cmp <- merge(ep[, c("dataset", "arm", "endpoint", "estimate", "lower", "upper")],
             p3[, c("dataset", "arm", "endpoint", "estimate", "lower", "upper")],
             by = c("dataset", "arm", "endpoint"),
             suffixes = c(".manec", ".nec4"))
cmp$ratio <- cmp$estimate.manec / cmp$estimate.nec4
cmp$nec4_inside_manec <- cmp$estimate.nec4 >= cmp$lower.manec &
  cmp$estimate.nec4 <= cmp$upper.manec
utils::write.csv(cmp, "analysis/phase9_vs_phase3.csv", row.names = FALSE)
for (en in c("ErC10", "ErC50", "NSEC")) {
  s <- cmp[cmp$endpoint == en, ]
  cat("\n--", en, "--\n")
  print(s[order(s$dataset, s$arm),
          c("dataset", "arm", "estimate.nec4", "estimate.manec", "ratio",
            "nec4_inside_manec")], digits = 4, row.names = FALSE)
}
cat("\nmedian |log ratio| by arm (how far averaging moved each arm):\n")
agg <- aggregate(abs(log(ratio)) ~ arm, cmp[cmp$endpoint != "NSEC", ], median)
names(agg)[2] <- "median_abs_log_ratio"
print(agg[order(agg$median_abs_log_ratio), ], digits = 3, row.names = FALSE)

## ---------------------------------------------- 3. where flooring acts -------
cat("\n===== 3. DOES FLOORING SHIFT THE WEIGHT ONTO DIFFERENT SHAPES? =====\n")
cat("The question Phase 3 could not ask. If a floored dataset simply moves\n",
    "weight from one curve shape to another, the convention is acting on model\n",
    "SELECTION and not only on the estimate within a model.\n")
w <- reshape(wt[, c("dataset", "arm", "model", "wi")],
             idvar = c("dataset", "model"), timevar = "arm", direction = "wide")
names(w) <- sub("^wi[.]", "", names(w))
w[is.na(w)] <- 0
for (ds in unique(w$dataset)) {
  s <- w[w$dataset == ds, ]
  s <- s[order(-s[[if ("A" %in% names(s)) "A" else 3]]), ]
  cat("\n--", ds, "-- stacking weight per model\n")
  print(s[, setdiff(names(s), "dataset")], digits = 2, row.names = FALSE)
}

cat("\ntotal weight on the top-weighted model of arm A, by arm:\n")
top <- do.call(rbind, lapply(unique(wt$dataset), function(ds) {
  a <- wt[wt$dataset == ds & wt$arm == "A", ]
  if (!nrow(a)) return(NULL)
  best <- a$model[which.max(a$wi)]
  s <- wt[wt$dataset == ds, ]
  do.call(rbind, lapply(split(s, s$arm), function(z)
    data.frame(dataset = ds, arm = z$arm[1], A_best_model = best,
               weight_on_it = sum(z$wi[z$model == best]),
               stringsAsFactors = FALSE)))
}))
print(reshape(top, idvar = c("dataset", "A_best_model"), timevar = "arm",
              direction = "wide"), digits = 2, row.names = FALSE)

cat("\n===== MODELS DROPPED BY bnec (fitting failures, not exclusions) =====\n")
drop <- unique(wt[wt$dropped != "" & !is.na(wt$dropped),
                  c("dataset", "arm", "requested", "dropped")])
if (nrow(drop)) print(drop, row.names = FALSE) else cat("none\n")

cat("\nPhase 9 report complete", format(Sys.time()), "\n")
