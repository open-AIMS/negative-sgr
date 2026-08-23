## Tests for the Phase 10 helpers that do not need bayesnec loaded. The ones that
## do -- p10_candidate_models(), p10_fit(), p10_curve() -- are exercised by
## analysis/phase10_gate1_models.R and Gate 2 instead, which is where a check
## that needs a fitted object belongs.

test_that("p10_cells reproduces the phase 8 grid, in order", {
  # The cell index is load-bearing: the reference-prior dataset of cell k is
  # drawn with seed 6e5 + k, so a grid in a different order silently gives every
  # cell a different prior from the one the single-model sweep used. This asserts
  # the literal grid from phase8_run.R rather than a property of it.
  cells <- p10_cells(8.09)
  expect_equal(nrow(cells), 12L)
  expect_equal(cells$delta, c(2, 4, 8, 2, 4, 8, 2, 4, 8, 4, 4, 4))
  expect_equal(cells$top_factor, c(rep(0.8, 3), rep(1.0, 3), rep(2.0, 3),
                                   rep(2.0, 3)))
  expect_equal(cells$R, c(rep(2.3, 9), 3.3, 17, 73))
  expect_equal(cells$cell_index, 1:12)
  # The three Phase 10 cells are the ones the plan names.
  expect_equal(cells$cell[cells$cell_index == 12], "d4.0_t2.0_R73.0_s8.1")
  expect_equal(cells$cell[cells$cell_index == 8],  "d4.0_t2.0_R2.3_s8.1")
  expect_equal(cells$cell[cells$cell_index == 9],  "d8.0_t2.0_R2.3_s8.1")
})

test_that("subsetting p10_cells does not renumber the cell index", {
  # This is the specific defect in phase8_run.R that Phase 10 must not inherit:
  # it subsets `cells` and then iterates k over seq_len(nrow(cells)), so asking
  # for cell 12 alone yields k = 1 and therefore seed 6e5 + 1.
  cells <- p10_cells(8.09)
  sub <- cells[match(c(12L, 8L, 9L), cells$cell_index), ]
  expect_equal(sub$cell_index, c(12L, 8L, 9L))
  expect_equal(sub$R, c(73, 2.3, 2.3))
  expect_false(identical(sub$cell_index, seq_len(nrow(sub))))
})

test_that("p10_zb_prep floors and scales the way each family needs", {
  sim <- data.frame(x = c(0, 1, 2, 3), sgr = c(0.12, 0.05, -0.02, -0.30))
  e <- p10_zb_prep(sim, "E")
  f <- p10_zb_prep(sim, "F")
  # Both floor: neither family has support below zero.
  expect_true(all(e$y >= 0))
  expect_true(all(f$y >= 0))
  # E scales to the observed maximum, which is the usual preparation for a Beta;
  # F does not, because a Gamma accommodates the positive response directly.
  expect_equal(max(e$y), 1)
  expect_equal(f$y, c(0.12, 0.05, 0, 0))
  expect_equal(e$y, c(1, 0.05 / 0.12, 0, 0))
  # The x column is carried through untouched -- an arm must never move the
  # concentration axis.
  expect_equal(e$x, sim$x)
  expect_equal(f$x, sim$x)
})

test_that("p10_zb_scale returns what puts E's curve back on the response scale", {
  sim <- data.frame(x = c(0, 1, 2), sgr = c(0.12, 0.05, -0.30))
  # E was divided by the maximum of the FLOORED response, not of the raw one.
  expect_equal(p10_zb_scale(sim, "E"), 0.12)
  expect_equal(p10_zb_scale(sim, "F"), 1)
  # Gaussian arms are never rescaled.
  expect_equal(p10_zb_scale(sim, "A"), 1)
  expect_equal(p10_zb_scale(sim, "B3"), 1)
})

test_that("p10_weights records what the candidate set actually became", {
  requested <- c("nec4param", "ecx4param", "ecxwb1", "ecxll4")
  ms <- data.frame(wi = c(0.6, 0.3, 0.1))
  rownames(ms) <- c("nec4param", "ecx4param", "ecxwb1")
  fake <- list(fit = structure(list(mod_stats = ms), class = "bayesmanecfit"))
  w <- p10_weights(fake, requested)
  expect_equal(nrow(w), 3L)
  expect_equal(w$wi, c(0.6, 0.3, 0.1))
  # A set of four that averaged over three is a different analysis and has to be
  # visible in the output rather than silently absorbed.
  expect_equal(unique(w$n_requested), 4L)
  expect_equal(unique(w$n_fitted), 3L)
  expect_equal(unique(w$dropped), "ecxll4")
})

test_that("p10_weights handles a set that collapsed to a single model", {
  # bnec() returns a bayesnecfit, not a bayesmanecfit, when only one model
  # survives. Reported with weight 1 so the two cases tabulate together.
  fake <- list(fit = structure(list(model = "nec4param"), class = "bayesnecfit"))
  w <- p10_weights(fake, c("nec4param", "ecxwb1"))
  expect_equal(nrow(w), 1L)
  expect_equal(w$wi, 1)
  expect_equal(w$model, "nec4param")
  expect_equal(w$dropped, "ecxwb1")
})

test_that("p10_assert_arms_ready refuses to run without the Phase 9 changes", {
  # The failure this guards against is silent: without the model-set changes the
  # sweep would fit single models and answer a different question over days.
  e <- new.env()
  # A reverted arms.R: fit_arm() exists but takes no model set.
  assign("fit_arm", function(arm, dat, prior) NULL, envir = e)
  assign("fit_diagnostics", function(fit) NULL, envir = e)
  assign("bot_bearing_models", function(models) models, envir = e)
  expect_error(p10_assert_arms_ready(e), "fit_arm\\(model=\\)")
  expect_error(p10_assert_arms_ready(e), "fit_diagnostics\\(\\) bayesmanecfit")

  # Nothing sourced at all: the message must name the file, not raise an
  # "object not found" that points back into this function.
  expect_error(p10_assert_arms_ready(new.env()), "not found at all")
})

test_that("p10_pull accepts a fit wrapper or a bare bnecfit", {
  # A bayesnecfit confusingly HAS a $fit of its own (the brmsfit), so the test
  # has to be on class rather than on the presence of that element.
  manec <- structure(list(mod_fits = list()), class = c("bayesmanecfit", "bnecfit"))
  nec <- structure(list(fit = "brmsfit-here", model = "nec4param"),
                   class = c("bayesnecfit", "bnecfit"))
  expect_identical(p10_pull(list(fit = manec)), manec)
  expect_identical(p10_pull(manec), manec)
  expect_identical(p10_pull(nec), nec)          # not nec$fit
  expect_identical(p10_pull(list(fit = nec)), nec)
})

test_that("the rhat drop rule keeps, drops and empties correctly", {
  # p10_drop_nonconverged() is driven entirely by p10_model_rhat(), so stubbing
  # that is enough to exercise all three outcomes without a fitted object.
  fake <- structure(list(mod_fits = list()), class = c("bayesmanecfit", "bnecfit"))
  with_rhat <- function(vals) {
    rlang_env <- environment()
    p10_model_rhat <<- function(fit) data.frame(model = names(vals),
                                                max_rhat = unname(vals),
                                                row.names = NULL)
  }
  keep <- p10_model_rhat                          # restore afterwards
  on.exit(p10_model_rhat <<- keep, add = TRUE)

  # nothing over the cutoff: returned untouched, and amend() is never called --
  # the common case must cost nothing.
  with_rhat(c(nec4param = 1.001, ecxwb1 = 1.004))
  r <- p10_drop_nonconverged(fake, cutoff = 1.01)
  expect_identical(r$fit, fake)
  expect_length(r$dropped, 0)
  expect_setequal(r$kept, c("nec4param", "ecxwb1"))
  expect_true(is.na(r$error))

  # all over the cutoff: no fit, a reason, and NOT a silent fallback to the
  # undropped average -- that would reinstate the models the rule removes.
  with_rhat(c(nec4param = 1.2, ecxwb1 = 1.3))
  r <- p10_drop_nonconverged(fake, cutoff = 1.01)
  expect_null(r$fit)
  expect_setequal(r$dropped, c("nec4param", "ecxwb1"))
  expect_length(r$kept, 0)
  expect_match(r$error, "all 2 models exceeded rhat")

  # the cutoff is strict: exactly 1.01 is kept, just above it is dropped.
  with_rhat(c(a = 1.01, b = 1.0100001))
  r <- p10_drop_nonconverged(fake, cutoff = 1.01)
  expect_equal(r$kept, "a")
  expect_equal(r$dropped, "b")
})

test_that("the study's two convergence rules are distinct and both recorded", {
  # Guards a real confusion: Trap 6 forbids filtering fitted DATASETS out of the
  # summaries on convergence grounds, while this phase deliberately drops
  # non-converged MODELS from inside an average. The constants encode the
  # difference -- a cutoff that screens models, and escalation switched off so
  # one bad component cannot re-sample the whole set.
  expect_equal(P10_RHAT_CUTOFF, 1.01)
  expect_true(is.infinite(P10_RHAT_ESCALATE))
  expect_true(P10_REPORT_ALL_MODELS)
})

test_that("every arm is initialised the same way", {
  # Under a model set, fit_arm() forces `init = "random"` for B3 because bnec()
  # passes one init to every model and they do not share a parameter vector.
  # If the other arms used bayesnec's own search, B3 would differ from them by
  # its initialisation as well as by its pinned asymptote -- a nuisance sitting
  # directly on the contrast the study is built on.
  expect_equal(P10_INIT, "random")
})

test_that("partial and total model drops are different kinds of event", {
  # The distinction the whole exclusion question turns on. A PARTIAL drop leaves
  # the arm with an estimate, so every dataset is still scored and what changed
  # is the estimator, applied uniformly -- that is the treatment. A TOTAL drop
  # leaves the arm with nothing for that dataset, and if the rate differs by arm
  # the arms end up scored on differently selected subsets, which is Trap 6.
  # p10_drop_nonconverged() must make the two distinguishable downstream: a
  # partial drop returns a fit, a total drop returns NULL with a reason.
  fake <- structure(list(mod_fits = list()), class = c("bayesmanecfit", "bnecfit"))
  keep <- p10_model_rhat
  on.exit(p10_model_rhat <<- keep, add = TRUE)

  # The partial branch calls bayesnec::amend(), which cannot run against a stub
  # object, so the assertion is on the DECISION rather than on the returned fit:
  # something survived, and the failure reason is not the all-dropped rule.
  # Whether amend() itself works on a real fit is checked by the smoke test,
  # which is where a claim needing a fitted object belongs.
  p10_model_rhat <<- function(fit) data.frame(model = c("a", "b", "c"),
                                              max_rhat = c(1.001, 1.2, 1.3))
  partial <- p10_drop_nonconverged(fake, cutoff = 1.01)
  expect_length(partial$kept, 1)
  expect_equal(partial$kept, "a")
  expect_setequal(partial$dropped, c("b", "c"))
  expect_false(grepl("exceeded rhat", partial$error))   # not the total-drop path

  p10_model_rhat <<- function(fit) data.frame(model = c("a", "b", "c"),
                                              max_rhat = c(1.2, 1.3, 1.4))
  total <- p10_drop_nonconverged(fake, cutoff = 1.01)
  expect_null(total$fit)                   # NOT scored -- the Trap 6 case
  expect_length(total$kept, 0)
  expect_false(is.na(total$error))

  # The two must not be confusable by inspecting `kept` alone being short.
  expect_true(length(partial$kept) > 0 && length(total$kept) == 0)
})
