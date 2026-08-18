## Phase 4 extension -- is the plateau prior-sensitive in EVERY method, or only
## in the intact analysis?
##
## The original sweep varied the `bot` prior for arm A alone, which answers
## "is the A-vs-B2 gap a prior artefact" but not "which methods depend on the
## prior". Those are different questions and the second is the one a reader
## asks. Arms B2 and B3 are excluded because their `bot` is `constant(0)`:
## there is no prior to vary, which is itself the point.
##
## Run: Rscript analysis/phase4_prior_sweep_allarms.R

source("R/setup.R"); source("R/data_prep.R"); source("R/diagnostics.R")
load_bayesnec(); source("R/arms.R"); source("R/metrics.R")
use_compile_cache()
options(mc.cores = 4L)

DATASET <- "c_proliferum"
SCALES <- c(tighter = 0.25, default = 1, wider = 4)
FREE_BOT <- c("A", "B1", "C", "D")

dat <- read_sgr(DATASET)
adat <- prepare_sgr(dat, "raw", meta = dataset_meta())
prior_default <- arm_prior(adat$x, adat$y)

scale_bot_prior <- function(prior, factor) {
  cls <- class(prior); pr <- as.data.frame(prior)
  i <- which(pr$nlpar == "bot")
  parts <- as.numeric(strsplit(gsub("^normal\\(|\\)$", "", pr$prior[i]), ",\\s*")[[1]])
  pr$prior[i] <- sprintf("normal(%.10g, %.10g)", parts[1], parts[2] * factor)
  class(pr) <- cls; pr
}
bot_stat <- function(pt, col) {
  i <- which(sub("_Intercept$", "", pt$parameter) == "bot")
  if (length(i) != 1) NA_real_ else pt[[col]][i]
}

## Arm D needs a zero-crossing from this dataset's own arm-A fit, so arm A at
## the default prior is fitted first and its crossing reused for every D fit.
## Recomputing the crossing per prior would confound the prior effect on D with
## a change in which data D even sees.
a_ref <- fit_arm("A", dat, prior = prior_default, meta = dataset_meta())
crossing <- zero_crossing(a_ref$fit)
cat("zero crossing from arm A (default prior):", signif(crossing, 4), "\n\n")

rows <- list()
for (arm in FREE_BOT) {
  for (lab in names(SCALES)) {
    pr <- scale_bot_prior(prior_default, SCALES[[lab]])
    cat("fitting", arm, "/", lab, "\n")
    f <- try(fit_arm(arm, dat, prior = pr, crossing = crossing,
                     meta = dataset_meta()), silent = TRUE)
    if (inherits(f, "try-error")) { cat("  FAILED\n"); next }
    et <- endpoint_table(f$fit, arm, DATASET)
    pt <- parameter_table(f$fit, arm, DATASET)
    dg <- fit_diagnostics(f$fit)
    rows[[length(rows) + 1]] <- cbind(
      arm = arm, prior_label = lab, prior_scale = unname(SCALES[[lab]]), et,
      bot_mean = bot_stat(pt, "mean"), bot_sd = bot_stat(pt, "sd"),
      contraction = contraction(f$fit, "bot")$contraction,
      divergences = dg$divergences, row.names = NULL)
    cat("  bot =", signif(bot_stat(pt, "mean"), 4), "\n")
  }
}
ep <- do.call(rbind, rows)
utils::write.csv(ep, "analysis/phase4_prior_sweep_allarms.csv", row.names = FALSE)

options(width = 200)
cat("\n===== HOW FAR DOES A 16-FOLD PRIOR CHANGE MOVE EACH METHOD? =====\n")
for (en in c("ErC10", "ErC50", "NSEC")) {
  cat("\n--", en, "--\n")
  s <- ep[ep$endpoint == en, ]
  o <- do.call(rbind, lapply(split(s, s$arm), function(d) {
    d <- d[order(d$prior_scale), ]
    base <- d$estimate[d$prior_label == "default"]
    data.frame(arm = d$arm[1],
               tighter = d$estimate[d$prior_label == "tighter"],
               default = base,
               wider = d$estimate[d$prior_label == "wider"],
               span_pct = 100 * (max(d$estimate) - min(d$estimate)) / base,
               bot_tighter = d$bot_mean[d$prior_label == "tighter"],
               bot_default = d$bot_mean[d$prior_label == "default"],
               bot_wider = d$bot_mean[d$prior_label == "wider"],
               bot_span_x = max(abs(d$bot_mean)) / min(abs(d$bot_mean)))
  }))
  print(o, digits = 3, row.names = FALSE)
}
cat("\n`span_pct` is how far the endpoint moves across the whole 16-fold prior\n",
    "range, as a percentage of its default value. `bot_span_x` is the same for\n",
    "the plateau, as a multiple. A method where the plateau moves a lot but the\n",
    "endpoint barely moves is reporting a robust endpoint from a parameter the\n",
    "data do not determine.\n")
