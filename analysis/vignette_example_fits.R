## Eight approaches fitted to one simulated dataset, for the vignette figure.
##
## Scenario 8 (delta = 4, top_factor = 2.0, R = 2.3) is the central case in the
## reaching regime: the design runs well past the zero crossing, so every
## approach has something to act on, and the true ErC10/ErC50 are known.
##
## The point of the figure is the fitted CURVE, not the endpoint: it shows
## directly how flooring and pinning bend the fit in the region the endpoints
## are read from.

source("R/setup.R"); source("R/data_prep.R"); source("R/diagnostics.R")
source("R/simulate.R"); load_bayesnec(); source("R/arms.R"); source("R/metrics.R")
use_compile_cache()
options(mc.cores = 4L)

cal <- calibrate_sigma(read_sgr("c_proliferum"))
SIGMA_0 <- sigma_0_at(cal$cv_control, R_ref = 2.3, t = 7)
truth <- sim_truth(R = 2.3, delta = 4, t = 7)
design <- sim_design(truth, top_factor = 2.0)
sim <- simulate_dataset(truth, design, cv_control = cal$cv_control,
                        sigma_ratio = cal$sigma_ratio, seed = 20260818,
                        sigma_mode = "absolute", sigma_0_abs = SIGMA_0)
saveRDS(list(truth = truth, sim = sim), "analysis/vignette_example_data.rds")

ap <- prepare_sgr(sim, "raw", meta = sim_meta())
prior <- arm_prior(ap$x, ap$y)
cat("true top", signif(truth$top,4), "bot", signif(truth$bot,4),
    "ErC10", signif(true_ecx(truth,10),4), "ErC50", signif(true_ecx(truth,50),4), "\n")

## Arm D needs the crossing from arm A, so A is fitted first.
fits <- list()
fits$A <- fit_arm("A", sim, prior = prior, meta = sim_meta())
crossing <- zero_crossing(fits$A$fit)
cat("estimated zero crossing:", signif(crossing, 4), "\n")
for (a in c("B1", "B2", "B3", "C", "D")) {
  cat("fitting", a, "\n")
  fits[[a]] <- fit_arm(a, sim, prior = prior, crossing = crossing,
                       meta = sim_meta())
}

## ------------------------------------------------ the two family-floored arms
## E and F are not `fit_arm()` arms: they change the family, not the data
## preparation, and they take `bnec()`'s own default priors for that family
## rather than the shared arm-A prior, because taking the defaults is the
## practice under examination (plan revision 6, Phase 7).
prep_zb <- function(sim, arm) {
  y0 <- pmax(sim$sgr, 0)
  if (identical(arm, "E")) data.frame(x = sim$x, y = y0 / max(y0))
  else                     data.frame(x = sim$x, y = y0)
}
fam_zb <- function(arm) {
  if (identical(arm, "E")) brms::Beta(link = "identity") else
    stats::Gamma(link = "identity")
}
## E's response is divided by its maximum, so its fitted curve and its `top`
## come back on the scaled response and would sit near 1 on a figure whose other
## seven curves are growth rates. Multiplying by that same maximum puts it back
## on the growth-rate axis; it is a pure rescaling and leaves `nec`, `beta` and
## every reported endpoint unchanged (see the vignette's note on ecx(type =
## "absolute") carrying any positive scaling of the response).
Y_MAX <- max(pmax(sim$sgr, 0))
zb_scale <- function(arm) if (identical(arm, "E")) Y_MAX else 1

for (a in c("E", "F")) {
  cat("fitting", a, "\n")
  d <- prep_zb(sim, a)
  pr <- bayesnec:::define_prior("nec3param", fam_zb(a), d$x, d$y,
                                "uninformative")
  fits[[a]] <- list(fit = bayesnec::bnec(
    bayesnec::bnf(y ~ crf(x, "nec3param")), data = d, family = fam_zb(a),
    prior = pr, chains = MCMC$chains, iter = MCMC$iter, warmup = MCMC$warmup,
    seed = MCMC$seed, backend = "cmdstanr",
    control = list(adapt_delta = MCMC$adapt_delta,
                   max_treedepth = MCMC$max_treedepth)))
}

## Fitted curves on a common grid, plus endpoints, for the vignette.
xx <- seq(0, max(sim$x), length.out = 200)
curves <- do.call(rbind, lapply(names(fits), function(a) {
  p <- fitted(fits[[a]]$fit, newdata = data.frame(x = xx), re_formula = NA)
  k <- zb_scale(a)
  data.frame(arm = a, x = xx, est = p[, "Estimate"] * k,
             lo = p[, "Q2.5"] * k, hi = p[, "Q97.5"] * k)
}))
ends <- do.call(rbind, lapply(names(fits), function(a)
  endpoint_table(fits[[a]]$fit, a, "sim")))
## `parameter_table()` greps for `sigma`, which a Gamma (`shape`) and a Beta
## (`phi`) do not have, so the dispersion is pulled separately and labelled with
## the family's own name for it. A Gaussian sigma and a Gamma shape are not the
## same quantity and must never be read down a single column.
zb_pars <- function(a) {
  dr <- brms::as_draws_df(fits[[a]]$fit$fit)
  k <- zb_scale(a)
  loc <- grep("^b_(top|bot|beta|nec)_Intercept$", names(dr), value = TRUE)
  disp <- intersect(c("sigma", "shape", "phi"), names(dr))
  do.call(rbind, lapply(c(loc, disp), function(p) {
    nm <- sub("_Intercept$", "", sub("^b_", "", p))     # Trap 10
    ## Only `top` and `bot` live on the response scale; `nec` is a
    ## concentration, `beta` a rate, and a dispersion parameter is not a
    ## rescaling of either, so none of the three is multiplied by k.
    v <- dr[[p]] * if (nm %in% c("top", "bot")) k else 1
    data.frame(dataset = "sim", arm = a, parameter = nm, mean = mean(v),
               sd = stats::sd(v), q2.5 = unname(stats::quantile(v, 0.025)),
               q50 = unname(stats::quantile(v, 0.5)),
               q97.5 = unname(stats::quantile(v, 0.975)))
  }))
}
pars <- do.call(rbind, lapply(names(fits), function(a)
  if (a %in% c("E", "F")) zb_pars(a)
  else parameter_table(fits[[a]]$fit, a, "sim")))
## Trap 10: `parameter_table()` strips with sub("^b_|_Intercept$", ...), and
## sub() replaces only the FIRST match, so its rows come back named
## "bot_Intercept" while the E/F path above names them "bot". Left alone the
## artefact carries two naming conventions in one column and every downstream
## filter silently drops half the arms. Normalised here rather than in
## `parameter_table()`, which the Phase 3 and Phase 4 outputs already depend on.
pars$parameter <- sub("_Intercept$", "", pars$parameter)

write.csv(curves, "analysis/vignette_example_curves.csv", row.names = FALSE)
write.csv(ends, "analysis/vignette_example_endpoints.csv", row.names = FALSE)
write.csv(pars, "analysis/vignette_example_parameters.csv", row.names = FALSE)

options(width = 200)
cat("\n=== endpoints, this dataset ===\n"); print(ends, digits = 4, row.names = FALSE)
cat("\n=== bot ===\n")
print(pars[grepl("^bot", pars$parameter), c("arm","mean","sd","q2.5","q97.5")],
      digits = 4, row.names = FALSE)
cat("\n=== dispersion (sigma for the Gaussian arms, shape/phi for F/E) ===\n")
print(pars[pars$parameter %in% c("sigma", "shape", "phi"),
           c("arm", "parameter", "mean", "sd")], digits = 4, row.names = FALSE)
cat("\n=== curve parameters, all arms (E on the growth-rate scale) ===\n")
print(stats::reshape(pars[pars$parameter %in% c("top","bot","beta","nec"),
                          c("arm","parameter","mean")],
                     idvar = "arm", timevar = "parameter", direction = "wide"),
      digits = 5, row.names = FALSE)
cat("\nE's y was divided by", signif(Y_MAX, 5),
    "before fitting; its top/bot above are back on the growth-rate scale.\n")
