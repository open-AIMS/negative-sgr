## Phase 7 -- what actually drives arm F's bias, and why it changes sign.
##
## RESUME.md records a variance-structure explanation for arm F that was
## PROPOSED AND WITHDRAWN: a Gamma's mean-squared variance overweighting the low
## tail predicts a bias of one sign at every `delta`, and F's ErC50 bias is
## +13.0% at delta 2 and -5.6% at delta 4. Until something explains the sign
## change the defensible statement is that F is badly biased and poorly covered,
## without saying why. This script is the minimum evidence needed to say more.
##
## The comparison that isolates the family is F against B3, NOT F against A.
## B3 floors negatives at zero and pins the lower asymptote at zero, so it sees
## EXACTLY the data F sees and imposes exactly the same zero boundary; the only
## thing left different is the likelihood. F-vs-A confounds three changes at
## once (data, asymptote, family), which is why it could not settle the
## question the first time. E is carried alongside because it is the other
## family arm and its scaling makes it the same contrast for Beta.
##
## Seeds are `7e5 + i`, identical to the sweep, so these are the sweep's own
## simulated datasets and the endpoints here must reproduce the sweep's.
##
## Run: N_ITER=20 Rscript analysis/phase7_f_mechanism.R

source("R/setup.R"); source("R/data_prep.R"); source("R/diagnostics.R")
source("R/simulate.R"); load_bayesnec(); source("R/arms.R"); source("R/metrics.R")
use_compile_cache()
options(mc.cores = 1L, width = 200)

N_ITER  <- as.integer(Sys.getenv("N_ITER", "20"))
WORKERS <- as.integer(Sys.getenv("WORKERS", as.character(STUDY_CORES)))
MCMC_SIM <- list(chains = 4L, iter = 2000L, warmup = 1000L, adapt_delta = 0.95,
                 max_treedepth = 12L, seed = 20260812L)

cal <- calibrate_sigma(read_sgr("c_proliferum"))
SIGMA_0 <- sigma_0_at(cal$cv_control, R_ref = 2.3, t = 7)

## The three reaching cells, at the sweep's own cell indices. `k` matters: the
## reference prior is drawn with seed 6e5 + k, so reusing the index reproduces
## the exact prior each sweep used rather than a new one that happens to look
## similar.
cells <- data.frame(k = 7:9, delta = c(2, 4, 8), top_factor = 2.0, R = 2.3,
                    sigma_ratio = cal$sigma_ratio)
cells$cell <- sprintf("d%.1f_t%.1f_R%.1f_s%.1f", cells$delta, cells$top_factor,
                      cells$R, cells$sigma_ratio)

prep_zb <- function(sim, arm) {
  y0 <- pmax(sim$sgr, 0)
  if (identical(arm, "E")) data.frame(x = sim$x, y = y0 / max(y0))
  else                     data.frame(x = sim$x, y = y0)
}
fam_zb <- function(arm) {
  if (identical(arm, "E")) brms::Beta(link = "identity") else
    stats::Gamma(link = "identity")
}
fit_zb <- function(arm, sim, prior, mcmc) {
  bayesnec::bnec(bayesnec::bnf(y ~ crf(x, "nec3param")), data = prep_zb(sim, arm),
                 family = fam_zb(arm), prior = prior, chains = mcmc$chains,
                 iter = mcmc$iter, warmup = mcmc$warmup, seed = mcmc$seed,
                 backend = "cmdstanr",
                 control = list(adapt_delta = mcmc$adapt_delta,
                                max_treedepth = mcmc$max_treedepth))
}

#' Curve parameters plus the family's dispersion parameter
#'
#' `parameter_table()` in R/metrics.R greps for `sigma` only, which is right for
#' the six Gaussian arms and returns nothing for a Gamma (`shape`) or a Beta
#' (`phi`). Dispersion is the whole point of the variance-structure hypothesis,
#' so it is extracted here rather than inferred from its absence.
#'
#' Trap 10: `sub()` replaces only the first match, so `parameter_table()` yields
#' "bot_Intercept". Names are stripped with an anchored suffix pattern instead.
par_row <- function(fit) {
  dr <- brms::as_draws_df(fit$fit)
  pars <- grep("^b_(top|bot|beta|nec)_Intercept$", names(dr), value = TRUE)
  disp <- intersect(c("sigma", "shape", "phi"), names(dr))
  v <- lapply(c(pars, disp), function(p) mean(dr[[p]]))
  nm <- c(sub("_Intercept$", "", sub("^b_", "", pars)), disp)
  ## Dispersion is renamed to a common column so the arms line up in one table;
  ## the family-specific name is kept alongside because a Gamma `shape` and a
  ## Gaussian `sigma` are not the same quantity and must never be averaged.
  out <- as.data.frame(setNames(v, nm))
  out$dispersion_par <- if (length(disp)) disp else NA_character_
  out$dispersion <- if (length(disp)) out[[disp]] else NA_real_
  out
}

run_iteration <- function(i, truth, design, prior_ref, priors_zb) {
  sim <- simulate_dataset(truth, design, cv_control = cal$cv_control,
                          sigma_ratio = cal$sigma_ratio, seed = 7e5 + i,
                          sigma_mode = "absolute", sigma_0_abs = SIGMA_0)
  f_neg <- mean(sim$sgr < 0)
  truths <- c(ErC10 = true_ecx(truth, 10), ErC50 = true_ecx(truth, 50))
  mcmc <- MCMC_SIM; mcmc$seed <- MCMC_SIM$seed + i

  one <- function(a) {
    f <- try(suppressMessages(
      if (a %in% c("E", "F")) fit_zb(a, sim, priors_zb[[a]], mcmc)
      else fit_arm(a, sim, prior = prior_ref, mcmc = mcmc, meta = sim_meta())$fit),
      silent = TRUE)
    if (inherits(f, "try-error")) return(NULL)
    et <- endpoint_table(f, a, "sim")
    et <- et[et$endpoint %in% names(truths), ]
    dg <- fit_diagnostics(f)
    cbind(arm = a, iteration = i, f_neg = f_neg,
          par_row(f), divergences = dg$divergences, max_rhat = dg$max_rhat,
          ErC10 = et$estimate[et$endpoint == "ErC10"],
          ErC50 = et$estimate[et$endpoint == "ErC50"],
          ErC10_lwr = et$lower[et$endpoint == "ErC10"],
          ErC10_upr = et$upper[et$endpoint == "ErC10"],
          ErC50_lwr = et$lower[et$endpoint == "ErC50"],
          ErC50_upr = et$upper[et$endpoint == "ErC50"],
          row.names = NULL)
  }
  ## bind_rows-style fill: A carries a `bot` column and B3/E/F do not (Stan does
  ## not declare a parameter whose prior is constant(0), and nec3param has none),
  ## so rbind() on the raw frames would fail on differing column sets.
  rows <- Filter(Negate(is.null), lapply(c("A", "B3", "E", "F"), one))
  if (!length(rows)) return(NULL)
  nms <- unique(unlist(lapply(rows, names)))
  do.call(rbind, lapply(rows, function(r) {
    for (n in setdiff(nms, names(r))) r[[n]] <- NA
    r[, nms, drop = FALSE]
  }))
}

all_res <- list()
for (j in seq_len(nrow(cells))) {
  cel <- cells[j, ]
  truth <- sim_truth(R = cel$R, delta = cel$delta, t = 7)
  design <- sim_design(truth, top_factor = cel$top_factor)
  ref <- simulate_dataset(truth, design, cv_control = cal$cv_control,
                          sigma_ratio = cal$sigma_ratio, seed = 6e5 + cel$k,
                          sigma_mode = "absolute", sigma_0_abs = SIGMA_0)
  prior_ref <- arm_prior(prepare_sgr(ref, "raw", meta = sim_meta())$x,
                         prepare_sgr(ref, "raw", meta = sim_meta())$y)
  priors_zb <- lapply(setNames(c("E", "F"), c("E", "F")), function(a) {
    d <- prep_zb(ref, a)
    bayesnec:::define_prior("nec3param", fam_zb(a), d$x, d$y, "uninformative")
  })

  cat("\n=== cell", cel$cell, "| warming compile cache ...")
  invisible(try(suppressMessages(run_iteration(0L, truth, design, prior_ref,
                                               priors_zb)), silent = TRUE))
  cat(" done\n")
  t0 <- Sys.time()
  res <- parallel::mclapply(seq_len(N_ITER), function(i)
           try(run_iteration(i, truth, design, prior_ref, priors_zb),
               silent = TRUE), mc.cores = WORKERS)
  res <- Filter(function(z) !inherits(z, "try-error") && !is.null(z), res)
  r <- do.call(rbind, res)
  r <- cbind(cell = cel$cell, delta = cel$delta, top_factor = cel$top_factor,
             R = cel$R, r, row.names = NULL)
  r$true_ErC10 <- true_ecx(truth, 10)
  r$true_ErC50 <- true_ecx(truth, 50)
  r$true_nec <- truth$nec
  r$true_beta <- truth$beta
  r$true_top <- truth$top
  r$true_bot <- truth$bot
  all_res[[j]] <- r
  cat(sprintf("  %d/%d iterations in %.1f min\n", length(res), N_ITER,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

res <- do.call(rbind, all_res)
saveRDS(res, "analysis/phase7_f_mechanism.rds")
utils::write.csv(res, "analysis/phase7_f_mechanism.csv", row.names = FALSE)

## --------------------------------------------------------------- read-out ---
cat("\n===== ENDPOINT BIAS, THESE ITERATIONS ONLY =====\n")
cat("Reproduces the sweep's ordering on", N_ITER, "iterations; it is not a",
    "replacement\nfor the sweep's", 240, "and its MCSE is correspondingly wide.\n")
bias <- do.call(rbind, lapply(split(res, list(res$delta, res$arm), drop = TRUE),
  function(s) data.frame(delta = s$delta[1], arm = s$arm[1], n = nrow(s),
    ErC10_rel = mean(s$ErC10 - s$true_ErC10, na.rm = TRUE) / s$true_ErC10[1],
    ErC50_rel = mean(s$ErC50 - s$true_ErC50, na.rm = TRUE) / s$true_ErC50[1])))
print(reshape(bias[, c("delta", "arm", "ErC50_rel")], idvar = "arm",
              timevar = "delta", direction = "wide"), digits = 3,
      row.names = FALSE)
cat("\nErC10:\n")
print(reshape(bias[, c("delta", "arm", "ErC10_rel")], idvar = "arm",
              timevar = "delta", direction = "wide"), digits = 3,
      row.names = FALSE)

cat("\n===== THE CURVE PARAMETERS THAT PRODUCE THEM =====\n")
cat("nec and beta are what ErCx is computed from; `top` sets the reference the",
    "\ndecline is measured against. Truth is printed per cell for comparison.\n")
for (d in unique(res$delta)) {
  s <- res[res$delta == d, ]
  cat("\n-- delta =", d, "| true nec", round(s$true_nec[1], 4),
      "| true beta", round(s$true_beta[1], 4),
      "| true top", round(s$true_top[1], 4),
      "| true bot", round(s$true_bot[1], 4), "\n")
  agg <- do.call(rbind, lapply(split(s, s$arm, drop = TRUE), function(z)
    data.frame(arm = z$arm[1], n = nrow(z),
               nec = mean(z$nec, na.rm = TRUE),
               beta = mean(z$beta, na.rm = TRUE),
               top = mean(z$top, na.rm = TRUE),
               bot = mean(z$bot, na.rm = TRUE),
               disp_par = z$dispersion_par[1],
               dispersion = mean(z$dispersion, na.rm = TRUE),
               divergences = mean(z$divergences, na.rm = TRUE))))
  print(agg, digits = 4, row.names = FALSE)
}

cat("\n===== THE ISOLATING CONTRAST: F MINUS B3, SAME DATA =====\n")
cat("Both floor negatives at zero and both hold the lower asymptote at zero,\n",
    "so any difference here is the LIKELIHOOD and nothing else. If the sign\n",
    "change lives in `nec`, the variance-structure story is about where the\n",
    "breakpoint lands; if it lives in `beta`, it is about curve shape.\n")
wide <- reshape(res[res$arm %in% c("B3", "F"),
                    c("delta", "iteration", "arm", "nec", "beta", "top",
                      "ErC10", "ErC50")],
                idvar = c("delta", "iteration"), timevar = "arm",
                direction = "wide")
cmp <- do.call(rbind, lapply(split(wide, wide$delta), function(s)
  data.frame(delta = s$delta[1], n = nrow(s),
             d_nec = mean(s$nec.F - s$nec.B3, na.rm = TRUE),
             d_beta = mean(s$beta.F - s$beta.B3, na.rm = TRUE),
             d_top = mean(s$top.F - s$top.B3, na.rm = TRUE),
             d_ErC10 = mean(s$ErC10.F - s$ErC10.B3, na.rm = TRUE),
             d_ErC50 = mean(s$ErC50.F - s$ErC50.B3, na.rm = TRUE))))
print(cmp, digits = 4, row.names = FALSE)

cat("\nDone", format(Sys.time()), "\n")
