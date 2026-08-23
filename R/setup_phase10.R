## Phase 10 constants. Sourced AFTER R/setup.R by every Phase 10 entry point.
##
## This file is deliberately additive and overrides nothing in R/setup.R. Phase 9
## (the case studies under model averaging) was in flight when Phase 10 was
## written and shares R/setup.R, R/arms.R and R/metrics.R; changing any of those
## would have altered a running analysis. Everything Phase 10 needs that does not
## already exist lives here and in R/model_average.R.

## ---------------------------------------------------------------- the pin ----
##
## Phase 10 runs on a DIFFERENT bayesnec build from Phases 1-9, and the reason is
## not convenience.
##
## `ecx()` and `nsec()` on a `bayesmanecfit` resampled the component draw index
## with an unseeded `sample()`, independently at every call site (bayesnec #216).
## The instability lands almost entirely on the interval -- on `manec_example`
## the model-averaged NSEC lower bound spanned 0.702-0.954 over six calls.
## Coverage is an interval property and is the quantity this phase exists to
## measure, so running on the old pin would have added bayesnec's own resampling
## noise to exactly the arms under test.
##
## The pin is therefore `374e511c` -- the Phase 1-9 pin, unchanged -- with the
## two #216 commits cherry-picked onto it and nothing else. Conflicts arose only
## where #216 was written on top of #180 (the dropped posterior cache), which is
## not in this pin; the resolutions keep the old storage and change only where
## the randomness is drawn. See NEWS.md on that branch.
##
## Why not re-pin to current `dev`: it has moved 30+ commits since 374e511c,
## touching define_prior.R, check_models.R, inits_functions.R, bnec.R and
## sysdata.rda. Any of those could shift a single-model number, which would put
## the averaged results and the 500-iteration single-model results on different
## packages -- the exact confound this phase exists to avoid.
##
## What makes that safe is that #216 touches the `bayesmanecfit` path only. Gate
## 0 asserts it: model-averaged output becomes reproducible AND a stored
## single-model fit returns bit-identical endpoints. Both halves must pass.
BAYESNEC_WORKTREE_P10 <- "/mnt/c/Rworking/bayesnec-negsgr-p10"
BAYESNEC_BRANCH_P10   <- "study-pin-216"
BAYESNEC_SHA_P10      <- "ef3954cf"      # asserted by prefix; see below

load_bayesnec_p10 <- function(check_sha = TRUE) {
  if (check_sha) {
    sha <- system2("git", c("-C", BAYESNEC_WORKTREE_P10, "rev-parse", "HEAD"),
                   stdout = TRUE)
    ## Prefix match, because this branch is study-owned and may be amended
    ## (a NEWS wording fix, say) without the analysis changing. A full-SHA
    ## assertion would then fail for a reason that is not about the analysis.
    ## The branch is never rebased onto a different base; if it is, change the
    ## constant deliberately.
    if (!startsWith(sha, BAYESNEC_SHA_P10)) {
      stop("Phase 10 bayesnec worktree is at ", sha, " but is pinned to ",
           BAYESNEC_SHA_P10, ". Check out the pinned commit or update ",
           "R/setup_phase10.R deliberately.")
    }
  }
  suppressMessages(pkgload::load_all(BAYESNEC_WORKTREE_P10, quiet = TRUE,
                                     export_all = FALSE))
  invisible(TRUE)
}

## ------------------------------------------------------------- the design ----

## The five arms. A is the reference and is NOT optional: averaging over a set of
## eight of which seven are wrong costs arm A some bias and some coverage too,
## and without measuring that, a movement in B3 cannot be attributed to the
## flooring rather than to averaging itself.
P10_ARMS <- c("A", "C", "B3", "E", "F")

## `decline`, not `all`. The hormesis models answer a different question (is
## there low-dose stimulation), the truth here is monotone, and admitting them
## would change what the arms are compared on rather than sharpening it.
## Consistent with Phase 9. The cost of the choice is that this phase measures
## model averaging over DECLINING shapes, not what bnec() does out of the box,
## and the report says so.
P10_MODEL_SET <- "decline"

## Cell indices into the Phase 5/7/8 grid -- NOT a renumbering. Phase 8's runner
## subsets `cells` by CELL_INDEX and then iterates `k` over seq_len(nrow(cells)),
## so running cell 12 alone gives k = 1 and therefore a different reference-prior
## dataset (seed 6e5 + k) from the one the full sweep used. Phase 10 carries the
## original index instead, which is what keeps its simulated datasets identical
## to those the single-model arms saw.
##
## Order is the run order and it is load-bearing: cell 12 is where the
## single-model failure is most extreme AND where averaging has its best chance,
## because precise data resolve the stacking weights sharply. If averaging
## cannot rescue cell 12 the hypothesis is dead and 8 and 9 are confirmation.
P10_CELL_ORDER <- c(12L, 8L, 9L)

## MCMC settings. Identical to phase8_run.R's MCMC_SIM -- not the MCMC list in
## R/setup.R, which is the real-data setting. Copied rather than sourced so that
## a later edit to Phase 8 cannot silently change Phase 10.
P10_MCMC <- list(chains = 4L, iter = 2000L, warmup = 1000L, adapt_delta = 0.95,
                 max_treedepth = 12L, seed = 20260812L)

## Convergence handling. THIS IS PART OF THE TREATMENT, NOT A FILTER ON RESULTS.
##
## The distinction matters because the study already carries an opposite-looking
## rule (Trap 6): simulation fits are never excluded from the summaries on
## convergence grounds, because arm A averages ~0 divergences and B2 averages
## 4-6, so a filter would preferentially discard the worst fits of the worst arm
## and flatter it. That rule governs whether a fitted DATASET enters the bias and
## coverage calculation, and it stands unchanged.
##
## This is a different thing. Dropping a non-converged MODEL from inside a
## candidate set, before averaging over the rest, is a step the analyst performs
## -- `bayesnec` documents it, via `rhat()` then `amend(drop = ...)` -- so it
## belongs inside the definition of the arm. A study of "what does model
## averaging give an analyst" that averaged over visibly unconverged components
## would be studying a workflow nobody uses.
##
## 1.01 rather than bayesnec's default of 1.05: it is the Vehtari et al. (2021)
## recommendation and it is what the study's own workflow uses. Because the
## per-model R-hat is recorded for every fit, the drop RATE at any other cutoff
## is reportable after the fact without refitting.
P10_RHAT_CUTOFF <- 1.01

## Endpoints are computed twice per arm: once over the whole candidate set and
## once over the converged subset. The pair brackets the sensitivity to the
## cutoff at the cost of extra `ecx()`/`nsec()` evaluations and no extra
## sampling. Set FALSE if the pilot shows the second evaluation is expensive.
P10_REPORT_ALL_MODELS <- TRUE

## Initial values: Stan's own random inits, for every arm.
##
## Two reasons, and the first is about comparability rather than speed.
##
## UNDER A MODEL SET, ARM B3 HAS NO CHOICE. `bnec()` passes one `init` argument
## to every model it fits, and the models do not share a parameter vector --
## inits drawn for `nec4param` are meaningless for `ecxll5` -- so `fit_arm()`
## hands initialisation to Stan whenever `bot` is pinned and a set is being
## fitted. If the other four arms used bayesnec's own initial-value search, B3
## would differ from them by its inits as well as by its constraint, and the
## contrast this study is built on would be contaminated by a nuisance.
##
## AND THE SEARCH IS EXPENSIVE ON EXACTLY THE DATA THIS STUDY IS ABOUT.
## `make_good_inits()` rejects any draw whose predicted curve leaves the observed
## response range. Floored data has no negative range, so nearly every draw is
## rejected and it grinds to its 10,000-trial fallback: measured by the Phase 9
## work at **612.8 s against 6.1 s** for random inits on the same data, priors
## and seed, with the sampling itself taking about a second. Arms B3, E and F all
## fit floored or non-negative responses, and arm C's response column carries
## zeros on the censored rows, so four of the five arms are exposed to it -- at
## 46 model fits per iteration.
##
## Legitimate because inits move the sampler's starting point, not its target.
## Verified rather than assumed, by the Phase 9 work: same data and seed, ErC10
## 2.7836 both ways, ErC50 6.347 against 6.367, NSEC 2.323 against 2.312, zero
## divergences and R-hat 1.005/1.007 either way.
##
## THE COST: Phase 5/7/8 used bayesnec's search for the single-model arms, so the
## paired contrast now differs in inits as well as in workflow. Phase 8 measured
## inits alone moving an endpoint by up to 7% on arm B1 -- sampler noise, not
## bias, and it averages out across iterations, but it is not nothing and it is
## why this is written down rather than defaulted into. Set to NULL to revert to
## bayesnec's search for the four arms that can use it; B3 will still take random
## inits regardless, which is the inconsistency this exists to remove.
P10_INIT <- "random"

## Escalation is OFF for Phase 10, and that is a consequence of the rule above.
##
## `fit_arm()` refits once with 3x iterations and 4x warmup when any R-hat
## exceeds a threshold. Under a model SET, `fit_diagnostics()` reports the
## maximum across every model, so one bad component escalates all twelve. Two
## reasons not to: the max of twelve trips a fixed threshold far more often than
## the max of one, so the budget would inflate unpredictably; and re-sampling
## every model is not what an analyst does about one stuck component -- dropping
## it is, which is what P10_RHAT_CUTOFF now handles.
##
## This does not meaningfully confound the paired contrast with Phase 8, where
## escalation was live for the Gaussian arms: mean maximum R-hat across that
## whole sweep was 1.004-1.012, so escalation fired rarely there. It also makes
## Phase 10 internally consistent, since arms E and F never had an escalation
## path at all.
P10_RHAT_ESCALATE <- Inf

## Grid on which each averaged fit's curve is stored, so that "is a flexible
## shape adequate over the region ErCx is read from" can be answered directly
## rather than inferred from the endpoints. ~1.2 MB for the whole phase; the
## refit needed to recover it later would cost days.
P10_GRID_N <- 100L
