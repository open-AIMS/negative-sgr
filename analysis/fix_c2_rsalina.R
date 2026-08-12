## Arm C2 on r_salina: one stuck chain, not a global failure.
##
## Chains 1/2/4 agree (beta ~0.43, nec ~4.86, bot ~-2.43); chain 3 sits at
## beta 2.59, nec 4.99. beta is the log decay rate, so exp(0.43) = 1.5 against
## exp(2.59) = 13.3 -- the chain wandered up a ridge where the curve is simply
## "very steep".
##
## The cause is the design, not the code. r_salina's entire transition happens
## inside one dose step: mean SGR is 0.80 at 5 ug/L and below the counting limit
## at 7.5. Under arm A the LOD-implied values are treated as observed, which
## pins the curve at -1.986 and with it beta. Under C2 those rows are only known
## to lie in [-2.75, -1.99], so the curve has slack and beta is bounded below
## but not above. That is a real limitation of interval censoring on this
## design and belongs in the results.
##
## The repair is sampler settings only -- longer warmup and a tighter step size.
## The prior and the model are untouched, so the arm stays comparable.
source("R/setup.R"); source("R/data_prep.R"); source("R/diagnostics.R")
load_bayesnec(); source("R/arms.R"); source("R/metrics.R")
use_compile_cache(); options(mc.cores = 4L)

d <- read_sgr("r_salina")
a <- prepare_sgr(d, "bound")
pr <- arm_prior(a$x, a$y)
m <- list(chains = 4L, iter = 12000L, warmup = 8000L, adapt_delta = 0.999,
          max_treedepth = 15L, seed = 20260812L)
f <- fit_arm("C2", d, pr, mcmc = m)
dgn <- fit_diagnostics(f$fit)
cat("\n== C2 r_salina, repaired settings ==\n"); print(dgn)
dr <- brms::as_draws_df(f$fit$fit)
for (p in c("b_top_Intercept","b_bot_Intercept","b_beta_Intercept","b_nec_Intercept","sigma")) {
  cat(sprintf("%-18s mean %10.4g | chains: %s\n", p, mean(dr[[p]]),
      paste(sprintf("%.4g", tapply(dr[[p]], dr$.chain, mean)), collapse=", ")))
}
print(endpoint_table(f$fit, "C2", "r_salina"))
saveRDS(f, "analysis/c2_rsalina_repaired.rds")
