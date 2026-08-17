## Phase 4 -- is the A-vs-B2 gap an artefact of the default `bot` prior?
##
## Phase 3 reports that arm A and arm B2 disagree on every dataset, and Phase 4's
## contraction table shows `bot` is substantially prior-informed in arm A. Those
## two facts together admit an alternative reading: that the gap is not the cost
## of pinning `bot` at zero but the cost of arm A's `bot` prior -- a prior nobody
## chose deliberately, since `bnec()` derives it from the response vector. This
## script rules that reading in or out.
##
## Design. One dataset, arm A fitted three times under `bot` priors that differ
## only in SCALE -- the package default, four times wider, four times tighter --
## against a single arm B2 fit. Location is deliberately held at the default
## `quantile(y, 0.1)`: moving it toward zero would smuggle B2's assumption into
## arm A's prior and the two arms would converge by construction, which answers
## a different question. Scale is the axis on which "a prior nobody chose" is
## actually a worry.
##
## The wider prior is the decisive one. At 10 * sd(y) it carries essentially no
## information about `bot` over the range the data occupy, so if arm A's
## endpoints hold there, whatever separates A from B2 came from the likelihood.
##
## Dataset: c_proliferum. It is the only one of the four where arm A and arm B2
## are both usable on ErC50 and NSEC (Phase 3's `usable` flag), B2 samples
## cleanly enough to be worth comparing (2 divergences, vs 8-12 elsewhere), and
## arm A's `bot` is the most prior-driven of the datasets whose asymptote is
## actually observed rather than substituted (contraction 0.387). That last
## point is what makes it the hard case: if the prior were going to matter
## anywhere, it would matter here.
##
## Run: Rscript analysis/phase4_prior_sweep.R

source("R/setup.R")
source("R/data_prep.R")
source("R/diagnostics.R")
load_bayesnec()
source("R/arms.R")
source("R/metrics.R")

use_compile_cache()
## Deliberately modest: the Phase 5 sweep may be holding 18 of 22 cores. Four
## fits at 4 sequential chains each is minutes, and contending with the sweep
## would slow both.
options(mc.cores = 4L)

DATASET <- Sys.getenv("PRIOR_SWEEP_DATASET", "c_proliferum")
SCALES <- c(default = 1, wider = 4, tighter = 0.25)

dat <- read_sgr(DATASET)
adat <- prepare_sgr(dat, if (any(dat$density == 0)) "bound" else "raw",
                    meta = dataset_meta())
prior_default <- arm_prior(adat$x, adat$y)

#' Rescale the `bot` prior's sd, holding its location
#'
#' Edited as a data.frame and re-classed for the same reason `fix_bot_prior()`
#' is: `brms::validate_prior()` would need a formula and data this function does
#' not have, and brms reads the prior as a data.frame downstream.
scale_bot_prior <- function(prior, factor) {
  cls <- class(prior)
  pr <- as.data.frame(prior)
  i <- which(pr$nlpar == "bot")
  parts <- as.numeric(strsplit(gsub("^normal\\(|\\)$", "", pr$prior[i]),
                               ",\\s*")[[1]])
  if (length(parts) != 2 || anyNA(parts)) {
    stop("bot prior is not a two-argument normal(): ", pr$prior[i])
  }
  pr$prior[i] <- sprintf("normal(%.10g, %.10g)", parts[1], parts[2] * factor)
  class(pr) <- cls
  pr
}

bot_prior_string <- function(prior) {
  pr <- as.data.frame(prior)
  pr$prior[pr$nlpar == "bot"]
}

#' Pull one summary column for `bot` out of a `parameter_table()`
#'
#' `parameter_table()` names its rows with `sub("^b_|_Intercept$", "", p)`, and
#' `sub()` replaces only the FIRST match, so the rows come back as
#' "bot_Intercept" rather than "bot". Matching on the stem rather than on
#' equality keeps this working whichever way that is later fixed.
#'
#' Returns NA rather than a zero-length vector when `bot` is absent: under arms
#' B2/B3 the prior is `constant(0)`, and Stan does not declare a constant as a
#' parameter, so there is no draws column to summarise.
bot_stat <- function(pt, column) {
  i <- which(sub("_Intercept$", "", pt$parameter) == "bot")
  if (length(i) != 1) NA_real_ else pt[[column]][i]
}

cat("dataset:", DATASET, "| preparation:",
    if (any(dat$density == 0)) "bound" else "raw", "| n =", nrow(adat), "\n")
cat("default bot prior:", bot_prior_string(prior_default),
    " (sd = 2.5 * sd(y), sd(y) =", signif(stats::sd(adat$y), 4), ")\n\n")

## ------------------------------------------------------------- the fits -----
## Arm B2 is fitted here rather than read from the Phase 3 store so the whole
## comparison comes out of one session with one bayesnec load, and so a stale
## target cannot silently supply the reference (Trap 9).
plan <- rbind(
  data.frame(arm = "A", label = names(SCALES), factor = unname(SCALES),
             stringsAsFactors = FALSE),
  data.frame(arm = "B2", label = "fixed_bot", factor = NA_real_,
             stringsAsFactors = FALSE)
)

rows_ep <- list()
rows_par <- list()
for (k in seq_len(nrow(plan))) {
  arm <- plan$arm[k]
  label <- plan$label[k]
  pr <- if (arm == "B2") prior_default else
    scale_bot_prior(prior_default, plan$factor[k])

  cat("fitting arm", arm, "/", label, ":", bot_prior_string(pr), "\n")
  f <- fit_arm(arm, dat, prior = pr, meta = dataset_meta())

  dg <- fit_diagnostics(f$fit)
  et <- endpoint_table(f$fit, arm, DATASET)
  ct <- contraction(f$fit, "bot")
  pt <- parameter_table(f$fit, arm, DATASET)

  rows_ep[[k]] <- cbind(
    prior_label = label, bot_prior = bot_prior_string(pr), et,
    max_rhat = dg$max_rhat, divergences = dg$divergences,
    escalated = f$escalated,
    bot_mean = bot_stat(pt, "mean"),
    bot_sd = bot_stat(pt, "sd"),
    prior_sd = ct$prior_sd, contraction = ct$contraction,
    row.names = NULL)
  rows_par[[k]] <- cbind(prior_label = label, pt, row.names = NULL)

  cat("  bot =", signif(bot_stat(pt, "mean"), 4),
      "| contraction =", signif(ct$contraction, 3),
      "| max Rhat =", signif(dg$max_rhat, 5),
      "| divergences =", dg$divergences, "\n")
}

ep <- do.call(rbind, rows_ep)
par <- do.call(rbind, rows_par)
utils::write.csv(ep, "analysis/phase4_prior_sweep.csv", row.names = FALSE)
utils::write.csv(par, "analysis/phase4_prior_sweep_parameters.csv",
                 row.names = FALSE)

options(width = 200)
cat("\n===== ENDPOINTS UNDER EACH bot PRIOR =====\n")
print(ep[order(ep$endpoint, ep$prior_label),
         c("endpoint", "arm", "prior_label", "bot_prior", "estimate", "lower",
           "upper", "bot_mean", "bot_sd", "contraction", "divergences")],
      digits = 4, row.names = FALSE)

## ------------------------------------------------------------- the gap ------
## The quantity under test. Two comparators are reported for each prior:
##   * how far arm A moves from its own default-prior fit (prior sensitivity),
##   * how far arm A sits from arm B2 (the gap the study reports).
## The gap is an artefact only if the first is comparable to the second.
cat("\n===== IS THE GAP A PRIOR ARTEFACT? =====\n")
b2 <- ep[ep$arm == "B2", c("endpoint", "estimate", "lower", "upper")]
names(b2) <- c("endpoint", "b2_est", "b2_lower", "b2_upper")
a <- ep[ep$arm == "A", ]
a_def <- a[a$prior_label == "default", c("endpoint", "estimate")]
names(a_def) <- c("endpoint", "a_default")

gap <- merge(merge(a, b2, by = "endpoint"), a_def, by = "endpoint")
gap$shift_from_default <- gap$estimate / gap$a_default - 1
gap$gap_to_B2 <- gap$estimate / gap$b2_est - 1
## Phase 3's operational test of "the arms disagree", reapplied under each prior.
## Reported in BOTH directions because containment is asymmetric and the
## informative direction depends on which interval is wider: on c_proliferum
## B2's ErC50 interval is 2.85 times arm A's, so it swallows A's estimate while
## its own falls outside A's. Reporting only `A_outside_B2` would score that as
## agreement when what it shows is that B2 cannot resolve the endpoint.
gap$A_outside_B2_interval <- gap$estimate < gap$b2_lower |
  gap$estimate > gap$b2_upper
gap$B2_outside_A_interval <- gap$b2_est < gap$lower | gap$b2_est > gap$upper
gap$b2_width_ratio <- (gap$b2_upper - gap$b2_lower) / (gap$upper - gap$lower)
gap$ratio_shift_to_gap <- abs(gap$shift_from_default) / abs(gap$gap_to_B2)

print(gap[order(gap$endpoint, gap$prior_label),
          c("endpoint", "prior_label", "estimate", "a_default", "b2_est",
            "shift_from_default", "gap_to_B2", "ratio_shift_to_gap",
            "A_outside_B2_interval", "B2_outside_A_interval",
            "b2_width_ratio")],
      digits = 3, row.names = FALSE)

cat("\nRead `ratio_shift_to_gap`: the prior effect as a fraction of the A-vs-B2\n",
    "gap. Well below 1 on every endpoint, with `A_outside_B2_interval` holding\n",
    "under all three priors, means the gap is a likelihood effect and the\n",
    "default prior is not carrying it. Near or above 1 would mean the study\n",
    "cannot distinguish the two and Phase 3 must be reported with that caveat.\n")

cat("\n===== bot POSTERIOR vs PRIOR SCALE =====\n")
cat("If `bot_sd` tracks `prior_sd` the parameter is prior-driven; if it is",
    "flat\nacross a 16-fold change in prior scale it is data-driven.\n")
print(a[order(a$prior_label), c("endpoint", "prior_label", "prior_sd",
                                "bot_sd", "bot_mean", "contraction")],
      digits = 4, row.names = FALSE)
