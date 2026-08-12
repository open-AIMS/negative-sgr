## Simulation machinery. No models are fitted here; these guard the algebra the
## metrics depend on.

test_that("the curve matches bayesnec's nec4param equation", {
  # bayesnec: bot + (top - bot) * exp(-exp(beta) * (x - nec) * step(x - nec)).
  # Below nec the curve is flat at top; at nec it equals top exactly.
  tr <- sim_truth(R = 2.3, delta = 3.7, t = 7, nec = 5, beta = log(0.4))
  expect_equal(nec4param_curve(0, tr$top, tr$bot, tr$beta, tr$nec), tr$top)
  expect_equal(nec4param_curve(tr$nec, tr$top, tr$bot, tr$beta, tr$nec), tr$top)
  expect_equal(nec4param_curve(tr$nec - 1, tr$top, tr$bot, tr$beta, tr$nec),
               tr$top)
  # Asymptotes to bot from above, never below it. Far enough out the
  # exponential underflows and the curve equals bot exactly, so the strict
  # inequality only holds where the exponential is representable.
  expect_gt(nec4param_curve(50, tr$top, tr$bot, tr$beta, tr$nec), tr$bot)
  expect_lt(nec4param_curve(50, tr$top, tr$bot, tr$beta, tr$nec),
            tr$bot + 1e-6)
  expect_equal(nec4param_curve(1e6, tr$top, tr$bot, tr$beta, tr$nec), tr$bot)
})

test_that("R and Delta parameterise top and bot as intended", {
  tr <- sim_truth(R = 16, delta = 2, t = 3)
  expect_equal(tr$mu_0, log(16) / 3)
  expect_equal(tr$top, tr$mu_0)
  expect_equal(tr$bot, -2 * tr$mu_0)
  expect_equal(tr$bot / tr$top, -2)
  # R and Delta move independently -- the whole point of the simulation, since
  # the four real datasets confound them.
  a <- sim_truth(R = 2.3, delta = 2)
  b <- sim_truth(R = 73, delta = 2)
  expect_equal(a$bot / a$top, b$bot / b$top)
  expect_false(isTRUE(all.equal(a$top, b$top)))
})

test_that("the zero-crossing is where the curve reaches zero growth", {
  for (delta in c(0.8, 1.8, 3.7)) {
    tr <- sim_truth(R = 2.3, delta = delta, nec = 5, beta = log(0.4))
    expect_equal(nec4param_curve(tr$zero_crossing, tr$top, tr$bot, tr$beta,
                                 tr$nec), 0)
    # And it depends on delta and the decay rate only, not on the growth rate.
    tr2 <- sim_truth(R = 73, delta = delta, nec = 5, beta = log(0.4))
    expect_equal(tr$zero_crossing, tr2$zero_crossing)
  }
})

test_that("analytic ErCx agrees with the package's grid estimator", {
  # If these disagree, the simulation would be measuring discretisation rather
  # than bias.
  for (delta in c(0.8, 3.7)) {
    tr <- sim_truth(R = 2.3, delta = delta, nec = 5, beta = log(0.4))
    chk <- check_ecx_estimator(tr, resolution = 1000)
    expect_true(all(abs(chk$analytic - chk$grid) <= chk$grid_spacing))
  }
})

test_that("ErCx sits above nec and below the zero-crossing", {
  # Neither ErC10 nor ErC50 ever lies in the region the 100% cap alters, which
  # is why the cap can only bias them indirectly.
  tr <- sim_truth(R = 2.3, delta = 3.7, nec = 5, beta = log(0.4))
  e10 <- true_ecx(tr, 10); e50 <- true_ecx(tr, 50)
  expect_gt(e10, tr$nec)
  expect_gt(e50, e10)
  expect_lt(e50, tr$zero_crossing)
})

test_that("the design places the top concentration relative to the crossing", {
  tr <- sim_truth(R = 2.3, delta = 1.8, nec = 5, beta = log(0.4))
  for (f in c(0.8, 1, 2)) {
    d <- sim_design(tr, top_factor = f, n_conc = 12, n_rep = 5, n_control = 10)
    expect_equal(max(d$x), f * tr$zero_crossing)
    expect_equal(sum(d$x == 0), 10)
    expect_equal(length(unique(d$x[d$x > 0])), 12)
    expect_equal(nrow(d), 10 + 12 * 5)
  }
})

test_that("sigma is flat when the ratio is 1 and rises when it is not", {
  tr <- sim_truth(R = 2.3, delta = 1.8, nec = 5, beta = log(0.4))
  x <- c(0, 5, 10, 50)
  flat <- sim_sigma(x, tr, cv_control = 0.096, sigma_ratio = 1)
  expect_equal(flat, rep(0.096 * tr$top, length(x)))
  rising <- sim_sigma(x, tr, cv_control = 0.096, sigma_ratio = 8)
  expect_equal(rising[1], 0.096 * tr$top)
  expect_true(all(diff(rising) >= 0))
  expect_lte(max(rising), 8 * 0.096 * tr$top)
})

test_that("simulated data round-trip through the real preparation code", {
  tr <- sim_truth(R = 2.3, delta = 1.8, nec = 5, beta = log(0.4))
  d <- sim_design(tr, top_factor = 2)
  sim <- simulate_dataset(tr, d, seed = 1)
  expect_equal(nrow(sim), nrow(d))
  # The density column exists so prepare_sgr() works unchanged. With no
  # detection limit, no row is ever undefined and raw/bound/supplied coincide.
  expect_true(all(sim$density > 0))
  raw <- prepare_sgr(sim, "raw", meta = sim_meta())
  sup <- prepare_sgr(sim, "supplied", meta = sim_meta())
  expect_equal(nrow(raw), nrow(sim))
  expect_equal(raw$y, sup$y)
  fl <- prepare_sgr(sim, "floored", meta = sim_meta())
  expect_gte(min(fl$y), 0)
})

test_that("the same seed gives the same data", {
  tr <- sim_truth(); d <- sim_design(tr)
  expect_equal(simulate_dataset(tr, d, seed = 42)$sgr,
               simulate_dataset(tr, d, seed = 42)$sgr)
  expect_false(isTRUE(all.equal(simulate_dataset(tr, d, seed = 42)$sgr,
                                simulate_dataset(tr, d, seed = 43)$sgr)))
})

test_that("MCSE formulas behave", {
  set.seed(1)
  est <- rnorm(500, 10, 2)
  m <- mcse_summary(est, 10, est - 4, est + 4)
  expect_equal(m$n, 500)
  expect_lt(abs(m$bias), 3 * m$bias_mcse + 0.01)
  expect_gt(m$coverage, 0.9)
  # Binomial MCSE, and the inverse used to size the sweep.
  expect_equal(m$coverage_mcse,
               sqrt(m$coverage * (1 - m$coverage) / 500))
  expect_equal(iterations_for_coverage_mcse(0.02), ceiling(0.95 * 0.05 / 0.0004))
})
