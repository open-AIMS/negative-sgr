## Phase 1 -- verify the pinned bayesnec branch before any study fits.
##
## Run on c_proliferum2: 70 rows, 10 negatives, no zero densities, minimum
## density 990, so no detection-limit complications to confound the checks.
##
## Writes analysis/phase1_gates.csv and prints a summary. Later phases refuse to
## run if the gating checks (1, 3, 5, 6) did not pass.

source("R/setup.R")
source("R/data_prep.R")
source("R/diagnostics.R")
load_bayesnec()
source("R/arms.R")

options(mc.cores = min(STUDY_CORES, 4L))
use_compile_cache()
PILOT <- list(chains = 2L, iter = 2000L, warmup = 1000L, adapt_delta = 0.95,
              max_treedepth = 12L, seed = 20260812L)

## bnec() has no `brm_args` argument -- brms arguments travel through `...`.
## A list passed as brm_args reaches brm() as one unused named argument and
## chains/backend/init silently revert to defaults, which is how the first run
## of this script ended up on the rstan backend with 4 chains. Helper keeps the
## call sites honest.
bnec_args <- function(...) {
  c(list(family = "gaussian", chains = PILOT$chains, iter = PILOT$iter,
         warmup = PILOT$warmup, seed = PILOT$seed, backend = "cmdstanr",
         control = list(adapt_delta = PILOT$adapt_delta)), list(...))
}

gates <- list()
dat <- read_sgr("c_proliferum2")
raw <- prepare_sgr(dat, "raw")
pr <- arm_prior(raw$x, raw$y)

## ---------------------------------------------------------------- gate 1 ----
## cens() survives the whole bayesnec path. c_proliferum2 has no censored rows
## of its own, so a synthetic left-censoring declaration is used: the 10 rows
## with negative SGR are declared "truth <= 0", which is arm C's treatment and
## is exactly what has to work.
cat("\n== gate 1: cens() through the whole path ==\n")
cdat <- prepare_sgr(dat, "raw")
cdat$y[cdat$y < 0] <- 0
cdat$cens[cdat$sgr < 0] <- "left"
g1 <- try({
  f_c <- do.call(bayesnec::bnec, bnec_args(
    formula = bayesnec::bnf(y | cens(cens) ~ crf(x, "nec4param")),
    data = cdat, prior = pr))
  list(ecx = bayesnec::ecx(f_c, ecx_val = 50),
       nsec = bayesnec::nsec(f_c),
       stancode = brms::stancode(f_c$fit))
}, silent = TRUE)
ok1 <- !inherits(g1, "try-error")
detail1 <- if (ok1) {
  paste0("ecx50 = ", signif(g1$ecx[1], 4),
         "; nsec = ", signif(g1$nsec[1], 4),
         "; stan code contains cens block: ",
         grepl("_lcdf|cens", g1$stancode))
} else as.character(g1)
gates$g1 <- gate_result(1, "cens() through formula, brms, ecx and nsec",
                        ok1, detail1)
cat(detail1, "\n")

## ---------------------------------------------------------------- gate 2 ----
## Diagnostic only (the main results use nec4param alone): does censoring change
## which of the decline models survive?
cat("\n== gate 2: candidate set, censored vs uncensored (diagnostic) ==\n")
decline <- bayesnec::models()$decline
set_unc <- bayesnec:::check_models(decline, stats::gaussian(link = "identity"))
## check_models() selects on family, link and the sign of x only -- it never
## inspects the censoring declaration -- so the *valid* set cannot differ
## between a censored and an uncensored fit of the same data. What censoring
## could still change is which of those models converge, and that is not
## testable without fitting all eight twice. The main results use nec4param
## alone, which removes the question, so this is recorded as a source fact
## rather than run as an experiment.
cens_args <- names(formals(bayesnec:::check_models))
ok2 <- identical(cens_args, c("model", "family", "data"))
detail2 <- paste0(length(set_unc), " decline models valid under gaussian: ",
                  paste(set_unc, collapse = ", "),
                  "; check_models() args are (",
                  paste(cens_args, collapse = ", "),
                  ") -- no censoring input, so the valid set is invariant")
gates$g2 <- gate_result(2, "candidate set unchanged by censoring", ok2, detail2)
cat(detail2, "\n")

## ---------------------------------------------------------------- gate 3 ----
## Per-row bounds: the bound is carried in the response column and a separate
## column takes "none"/"left". Verified by checking that brms received the
## censoring vector row-for-row.
cat("\n== gate 3: per-row censoring bounds ==\n")
g3 <- try({
  sd_c <- brms::standata(f_c$fit)
  list(cens_vec = sd_c$cens, y = sd_c$Y)
}, silent = TRUE)
ok3 <- !inherits(g3, "try-error") &&
  identical(as.integer(g3$cens_vec), as.integer(cdat$cens == "left") * -1L)
detail3 <- if (inherits(g3, "try-error")) as.character(g3) else
  paste0(sum(g3$cens_vec == -1), " of ", length(g3$cens_vec),
         " rows passed to brms as left-censored; expected ",
         sum(cdat$cens == "left"))
gates$g3 <- gate_result(3, "per-row censoring bounds reach brms", ok3, detail3)
cat(detail3, "\n")

## ---------------------------------------------------------------- gate 4 ----
## Extract, do not assume, the default bot prior.
cat("\n== gate 4: default bot prior ==\n")
g4 <- try({
  f_default <- do.call(bayesnec::bnec, bnec_args(
    formula = bayesnec::bnf(y ~ crf(x, "nec4param")), data = raw))
  bayesnec::pull_prior(f_default)
}, silent = TRUE)
ok4 <- !inherits(g4, "try-error")
detail4 <- if (ok4) {
  bp <- as.data.frame(g4)
  paste0("bot: ", bp$prior[bp$nlpar == "bot"],
         " | top: ", bp$prior[bp$nlpar == "top"])
} else as.character(g4)
gates$g4 <- gate_result(4, "default bot prior extracted from a fit", ok4,
                        detail4)
cat(detail4, "\n")

## ---------------------------------------------------------------- gate 5 ----
## The gaussian candidate set drops every zero-bounded model, so arm A is
## available off the shelf -- and, as a consequence, nec3param is NOT available
## as the fixed-bot arm.
cat("\n== gate 5: gaussian drops zero-bounded models ==\n")
dropped <- intersect(bayesnec::models()$zero_bounded, decline)
ok5 <- length(intersect(set_unc, bayesnec::models()$zero_bounded)) == 0 &&
  "nec4param" %in% set_unc
nec3_blocked <- inherits(
  try(bayesnec:::check_models("nec3param", stats::gaussian(link = "identity")),
      silent = TRUE), "try-error")
detail5 <- paste0("zero-bounded decline models dropped: ",
                  paste(dropped, collapse = ", "),
                  "; nec4param retained: ", "nec4param" %in% set_unc,
                  "; nec3param alone under gaussian errors: ", nec3_blocked)
gates$g5 <- gate_result(5, "gaussian drops zero-bounded models", ok5, detail5)
cat(detail5, "\n")

## ---------------------------------------------------------------- gate 6 ----
## NEW. Arms B2/B3 need bot held at exactly 0. nec3param cannot do it (gate 5),
## so the only route is a constant(0) prior on a nec4param. This gate decides
## whether half the 2x2 is buildable.
cat("\n== gate 6: constant(0) prior on bot ==\n")
g6 <- try({
  prc <- arm_prior(raw$x, raw$y, fix_bot = TRUE)
  ini <- fixed_bot_inits(raw$x, raw$y, PILOT$chains, seed = PILOT$seed)
  f_b2 <- do.call(bayesnec::bnec, bnec_args(
    formula = bayesnec::bnf(y ~ crf(x, "nec4param")), data = raw,
    prior = prc, init = ini))
  post_bot <- brms::as_draws_df(f_b2$fit)[["b_bot_Intercept"]]
  list(fit = f_b2,
       bot_range = if (is.null(post_bot)) NA else range(post_bot),
       ecx = bayesnec::ecx(f_b2, ecx_val = 50),
       nsec = bayesnec::nsec(f_b2))
}, silent = TRUE)
ok6 <- !inherits(g6, "try-error") &&
  (all(is.na(g6$bot_range)) || all(g6$bot_range == 0))
detail6 <- if (inherits(g6, "try-error")) as.character(g6) else
  paste0("bot posterior range: ",
         paste(signif(g6$bot_range, 4), collapse = " to "),
         "; ecx50 = ", signif(g6$ecx[1], 4),
         "; nsec = ", signif(g6$nsec[1], 4))
gates$g6 <- gate_result(6, "constant(0) prior on bot (arm B2/B3 feasibility)",
                        ok6, detail6)
cat(detail6, "\n")

## ------------------------------------------------------------------ out -----
res <- do.call(rbind, gates)
dir.create("analysis", showWarnings = FALSE)
utils::write.csv(res, "analysis/phase1_gates.csv", row.names = FALSE)
cat("\n===== PHASE 1 SUMMARY =====\n")
print(res[, c("id", "passed", "description")])
blocking <- res$id %in% c(1, 3, 5, 6)
cat("\nBlocking gates passed:", all(res$passed[blocking]), "\n")
if (ok1 && ok6) {
  saveRDS(list(cens_fit = g1, b2_fit = g6$fit),
          "analysis/phase1_fits.rds")
}
