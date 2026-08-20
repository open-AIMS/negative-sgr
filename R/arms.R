## The experimental arms. One function per arm, common signature:
##
##   fit_arm(arm, dat, prior, crossing = NULL, ...) -> list(fit, arm, data, ...)
##
## Arm C2 was added after arm C turned out to leave `bot` unidentified. A
## left-censored row contributes Phi((0 - mu)/sigma), which saturates at 1 once
## mu is a few sigma below zero, so the likelihood is flat in `bot` over the
## whole region the data occupy and the prior decides where the posterior lands.
## Measured on c_proliferum2: arm A gives bot = -1.27, arm C gives -0.42, with
## the same data and the same prior. C2 bounds that flat region from below with
## a statement that is true -- growth lies between extinction and zero -- rather
## than leaving it open. C is kept because its behaviour is a result: censoring
## negatives at zero is a decision not to estimate the asymptote, and an analyst
## should know that the number reported for it is then the prior's.
##
## Every arm uses the gaussian family, the nec4param model and the SAME prior,
## so that the only thing varying between arms is the data or the treatment of
## `bot`. The prior is built once per dataset from the arm-A response vector --
## see `arm_prior()` for why that matters.

arm_names <- function() c("A", "B1", "B2", "B3", "C", "C2", "D")

#' Build the prior shared by every arm of one dataset
#'
#' `bnec()`'s default gaussian priors are functions of the response vector:
#' `top ~ normal(quantile(y, 0.9), 2.5 * sd(y))` and
#' `bot ~ normal(quantile(y, 0.1), 2.5 * sd(y))`. The arms differ in their
#' response vector, so taking the default in each arm would confound the
#' likelihood change under study with a prior change. On `r_salina` the default
#' `bot` prior is centred at exactly 0 *because of the substituted zeros* -- the
#' convention being tested leaks into the prior.
#'
#' So: call the package's own prior constructor once, on arm A's data, and hand
#' the result to every arm. This is the default a practitioner fitting arm A
#' would get, not an invention.
#'
#' @param x,y Predictor and response from the arm-A preparation.
#' @param fix_bot Replace the `bot` prior with `constant(0)` (arms B2, B3).
arm_prior <- function(x, y, fix_bot = FALSE, model = "nec4param") {
  pr <- bayesnec:::define_prior(
    model = model,
    family = stats::gaussian(link = "identity"),
    predictor = x,
    response = y,
    prior_type = "uninformative"
  )
  if (fix_bot) pr <- fix_bot_prior(pr)
  pr
}

#' Replace an existing prior's `bot` with a point mass at zero
#'
#' Split out of `arm_prior()` so that arms B2/B3 can take *the prior they were
#' given* and pin `bot`, rather than deriving a fresh one from their own data.
#'
#' That distinction is not cosmetic. brms writes prior constants into the Stan
#' program as literals, so a prior rebuilt per dataset makes the program unique
#' and forces a 3-5 minute recompile for every fit. `fit_arm()` used to call
#' `arm_prior(adat$x, adat$y, fix_bot = TRUE)`, which meant B2 and B3 kept
#' recompiling through the Phase 5 sweep even with the prior otherwise held
#' fixed -- 18 workers each running a 1.2 GB `cc1plus` is what pinned the
#' machine at 29 of 31 GB.
#'
#' Behaviour is unchanged for Phase 3: arms A and B2 take the same preparation
#' there, so the prior derived from B2's data was already identical to the one
#' passed in.
#'
#' Edited as a data.frame and re-classed rather than passed through
#' `brms::validate_prior()`, which needs a formula and data this function does
#' not have. brms reads the prior as a data.frame downstream, so the class
#' attribute is all that has to survive.
fix_bot_prior <- function(prior) {
  ## A model SET arrives as a named list of brmsprior objects, one per model
  ## (see `?bnec`). Pin `bot` in each of them. Models with no `bot` parameter --
  ## `neclin` and `ecxlin` are the two in the Gaussian declining set -- pass
  ## through unchanged, because the subsetting expression matches no rows.
  ## Whether such a model belongs in a pinned-asymptote arm at all is a question
  ## for the caller, not for this function; see `bot_bearing_models()`.
  if (is.list(prior) && !inherits(prior, "brmsprior")) {
    return(lapply(prior, fix_bot_prior))
  }
  cls <- class(prior)
  pr <- as.data.frame(prior)
  pr$prior[pr$nlpar == "bot"] <- "constant(0)"
  pr$lb[pr$nlpar == "bot"] <- NA
  pr$ub[pr$nlpar == "bot"] <- NA
  class(pr) <- cls
  pr
}

#' The models in a set that have a lower asymptote to pin
#'
#' Arms B2 and B3 constrain the lower asymptote to zero. Two of the eight
#' declining models bayesnec retains under a Gaussian family -- `neclin` and
#' `ecxlin` -- have no asymptote at all: they decline linearly without bound, so
#' there is no parameter the constraint can act on and no sense in which they
#' obey it. They are therefore dropped from the candidate set of the pinned
#' arms, which is a consequence of the convention rather than a choice about it.
#'
#' The cost has to be reported: dropping them changes the candidate set, so a
#' pinned arm differs from arm A by its model set as well as by the constraint.
bot_bearing_models <- function(models) {
  ## `show_params()` returns a LIST of brmsformula objects, one per model, so
  ## the formula has to be unwrapped before its parameter names are visible.
  ## Reading `$pforms` off the outer list silently yields NULL and would drop
  ## every model, which looks like "no model has a bot" rather than an error.
  has_bot <- function(m) {
    "bot" %in% names(bayesnec:::show_params(m)[[1]]$pforms)
  }
  Filter(has_bot, models)
}

#' Initial values, needed only when `bot` is held constant
#'
#' `bayesnec:::make_inits()` draws inits by parsing each prior string and
#' calling the matching `r<dist>()`. Its lookup table holds gamma, normal, beta
#' and uniform only, so `constant(0)` has no draw function and the default init
#' path errors before `brm()` is reached. `bnec()` skips that path entirely when
#' `init` is supplied in `brm_args`, so arms B2/B3 supply their own: drawn from
#' the *free* prior for the sampled parameters, with `b_bot` omitted because
#' Stan does not declare a constant parameter.
fixed_bot_inits <- function(x, y, chains, model = "nec4param", seed = NULL,
                            prior_free = NULL) {
  # Inits are drawn against a prior in which `bot` is a point mass at 0 rather
  # than the free prior, so that the curve whose predictions get accepted or
  # rejected is the one actually being fitted. `make_good_inits()` rejects draws
  # whose predicted curve leaves the range of the response, which is what stops
  # a negative `top` reaching the sampler.
  pseudo <- as.data.frame(
    if (is.null(prior_free)) arm_prior(x, y, fix_bot = FALSE, model = model)
    else prior_free)
  pseudo$prior[pseudo$nlpar == "bot"] <- "normal(0, 1e-8)"
  pseudo$lb[pseudo$nlpar == "bot"] <- NA
  pseudo$ub[pseudo$nlpar == "bot"] <- NA
  class(pseudo) <- c("brmsprior", "data.frame")
  inits <- bayesnec:::make_good_inits(model, x, y, seed = seed,
                                      priors = pseudo, chains = chains)
  if (identical(names(inits), "random")) return(inits)
  # Stan does not declare a parameter whose prior is constant, so an init for
  # b_bot would be an unknown variable at initialisation.
  lapply(inits, function(z) z[setdiff(names(z), "b_bot")])
}

#' Fit one arm
#'
#' @param arm One of `arm_names()`, or "SQ" for the lab-practice benchmark.
#' @param dat A raw dataset from `read_sgr()`.
#' @param prior The shared prior from `arm_prior()`; rebuilt with
#'   `fix_bot = TRUE` internally for B2/B3.
#' @param crossing For arm D only: the zero-crossing estimate to truncate at,
#'   which must come from this dataset's own arm-A fit.
#' @param mcmc A list of MCMC settings (see `R/setup.R`).
#' @param rhat_threshold Refit once with longer warmup and a smaller step size
#'   if any Rhat exceeds this. Recorded in the returned `escalated` flag.
fit_arm <- function(arm, dat, prior, crossing = NULL, mcmc = MCMC,
                    model = "nec4param", meta = dataset_meta(),
                    rhat_threshold = 1.05, init = NULL) {
  ds <- unique(dat$dataset)
  has_undefined <- any(dat$density == 0)
  fix_bot <- arm %in% c("B2", "B3")

  ## Which data preparation each arm sees. On the Rhodomonas sets there is no
  ## "raw" SGR for rows where the population fell below the counting limit, so
  ## the arms that would use `raw` use `bound` instead -- the other common lab
  ## convention, substitute the detection limit. That substitution is declared
  ## here rather than hidden inside "arm A".
  prep <- switch(
    arm,
    A  = if (has_undefined) "bound" else "raw",
    B1 = "floored",
    B2 = if (has_undefined) "bound" else "raw",
    B3 = "floored",
    C  = "censored",
    C2 = "interval",
    D  = if (has_undefined) "bound" else "raw",
    A_raw_dropped = "raw",
    SQ = NA_character_,
    stop("Unknown arm: ", arm)
  )

  if (identical(arm, "SQ")) {
    adat <- percent_inhibition(dat)
    # "Beta", not "beta": validate_family() does get(family)(link = "identity"),
    # and base R's beta() is the beta *function* of two arguments, which fails
    # with "unused argument (link = 'identity')". brms::Beta() is the family.
    fam <- "Beta"
  } else {
    adat <- prepare_sgr(dat, prep, cens_negatives = identical(arm, "C"),
                        meta = meta)
    fam <- "gaussian"
  }

  if (identical(arm, "D")) {
    if (is.null(crossing)) stop("Arm D needs a zero-crossing estimate.")
    adat <- truncate_at_crossing(adat, crossing)
  }

  censored <- any(adat$cens != "none")
  interval <- any(adat$cens == "interval")
  # brms' interval form is cens(indicator, upper): the response column carries
  # the lower bound and the second argument the upper. The one-argument form
  # cannot express an interval, so the two arms use different formulas.
  ## `model` may name one model ("nec4param"), a bayesnec model GROUP
  ## ("decline"), or an explicit set of several. The first two are one string
  ## and go in quoted; a set has to be written into the formula as a `c(...)`
  ## call, because pasting a length-n vector into one quoted slot would
  ## silently produce n formulas and fail downstream.
  mod_arg <- if (length(model) > 1) {
    paste0("c(", paste0("\"", model, "\"", collapse = ", "), ")")
  } else {
    paste0("\"", model, "\"")
  }
  form <- if (interval) {
    bayesnec::bnf(paste0("y | cens(cens, y2) ~ crf(x, ", mod_arg, ")"))
  } else if (censored) {
    bayesnec::bnf(paste0("y | cens(cens) ~ crf(x, ", mod_arg, ")"))
  } else {
    bayesnec::bnf(paste0("y ~ crf(x, ", mod_arg, ")"))
  }

  # bnec() has no `brm_args` argument: everything destined for brm() travels
  # through `...`. Passing a list(brm_args = ...) silently reaches brm() as one
  # unused named argument, so chains, backend and init all revert to defaults
  # without any warning. Hence do.call over a flat list.
  prior_free <- prior
  args <- list(formula = form, data = adat, family = fam,
               chains = mcmc$chains, iter = mcmc$iter, warmup = mcmc$warmup,
               seed = mcmc$seed, backend = "cmdstanr",
               control = list(adapt_delta = mcmc$adapt_delta,
                              max_treedepth = mcmc$max_treedepth))
  manec <- is.list(prior) && !inherits(prior, "brmsprior")
  ## An explicit `init` overrides bayesnec's initial-value search for every arm.
  ##
  ## Measured on floored data (`c_proliferum`, B1, `nec4param`, identical priors
  ## and seed): bayesnec's search takes **612.8 s** against **6.1 s** for Stan's
  ## own random inits -- a hundredfold difference, and the whole run time, since
  ## the sampling itself takes about a second. The cause is documented in
  ## `make_good_inits()`: it rejects any draw whose predicted curve leaves the
  ## observed response range, and floored data has no negative range, so almost
  ## every draw is rejected and it grinds to the 10,000-trial fallback.
  ##
  ## Checked before adopting rather than assumed: the two agree to within Monte
  ## Carlo error on the same data and seed -- ErC10 2.7836 both ways, ErC50
  ## 6.347 against 6.367, NSEC 2.323 against 2.312, zero divergences and R-hat
  ## 1.005/1.007 either way. Initial values move the sampler's starting point,
  ## not its target.
  if (!is.null(init)) args$init <- init
  if (fix_bot) {
    # Pin `bot` in the prior we were GIVEN. Deriving a fresh one here would make
    # the Stan program unique to this dataset and force a recompile -- see
    # `fix_bot_prior()`. Inits may stay data-derived: they are passed as values,
    # not compiled into the program, so they cost nothing in cache terms.
    prior <- fix_bot_prior(prior)
    if (!is.null(init)) {
      # caller has been explicit; leave their choice alone
    } else if (manec) {
      # Under a model SET, `bnec()` passes one `init` argument to every model it
      # fits, and the models do not share a parameter vector -- inits drawn for
      # `nec4param` are meaningless for `ecxll5`. There is no per-model init
      # hook. So hand initialisation to Stan instead, which is legitimate here
      # for the reason `fixed_bot_inits()` exists at all: that helper is only
      # needed because bayesnec's own init path cannot draw from `constant(0)`,
      # and Stan's random inits never consult the prior. Inits change the
      # sampler's path, not its target.
      args$init <- "random"
    } else {
      args$init <- fixed_bot_inits(adat$x, adat$y, mcmc$chains, model,
                                   seed = mcmc$seed, prior_free = prior_free)
    }
  }
  if (!identical(arm, "SQ")) args$prior <- prior

  t0 <- Sys.time()
  fit <- do.call(bayesnec::bnec, args)

  ## Escalate once on non-convergence, rather than reporting a stuck chain as a
  ## result or hand-patching the one fit that failed.
  ##
  ## Arm C2 on r_salina is the case this exists for. Three chains agreed
  ## (beta ~0.41, nec 4.86) and one sat at beta 2.59, giving Rhat 1.47. The
  ## cause is the design: that test's whole transition happens inside one dose
  ## step (mean SGR 0.80 at 5 ug/L, below the counting limit at 7.5). Under the
  ## arms that treat the LOD-implied value as observed, the curve is pinned
  ## there and `beta` with it; under interval censoring those rows are only
  ## known to lie in [-2.75, -1.99], so the curve has slack and `beta` is
  ## bounded below but not above. A chain can wander up that ridge.
  ##
  ## Longer warmup and a smaller step size fix it (Rhat 1.001, 11 divergences in
  ## 16000). Only the sampler settings change -- model, data and prior are
  ## untouched -- so the arm stays comparable with the rest. The escalation is
  ## recorded so it can be reported rather than hidden.
  ## `fit_diagnostics()` knows how to read both a single fit and a model
  ## average; `brms::rhat(fit$fit)` does not, and on a `bayesmanecfit` fails
  ## inside `array()` with nothing pointing back to here.
  escalated <- FALSE
  if (fit_diagnostics(fit)$max_rhat > rhat_threshold) {
    escalated <- TRUE
    args$iter <- mcmc$iter * 3L
    args$warmup <- mcmc$warmup * 4L
    args$control <- list(adapt_delta = 0.999, max_treedepth = 15L)
    if (fix_bot && !manec && is.null(init)) {
      args$init <- fixed_bot_inits(adat$x, adat$y, mcmc$chains, model,
                                   seed = mcmc$seed)
    }
    fit <- do.call(bayesnec::bnec, args)
  }
  list(
    fit = fit,
    arm = arm,
    dataset = ds,
    preparation = prep,
    family = fam,
    censored = censored,
    n = nrow(adat),
    n_censored = sum(adat$cens != "none"),
    crossing = crossing,
    prior = prior,
    escalated = escalated,
    elapsed = as.numeric(difftime(Sys.time(), t0, units = "secs"))
  )
}

#' Posterior median concentration at which the fitted curve crosses zero growth
#'
#' Arm D's truncation point. Stated as an algorithm over the arm-A posterior so
#' that it never touches the truth in simulation. Returns `Inf` when the curve
#' does not cross zero inside the tested range, in which case arm D is the same
#' data as arm A and should be reported as such.
zero_crossing <- function(fit, resolution = 1000) {
  ## The predictor range has to be read differently from a single fit and a
  ## model average: a `bayesmanecfit` has no `$fit`, it has one brmsfit per
  ## candidate model under `$mod_fits`. Reading `$fit$data$x` off one returns
  ## NULL, and `range(NULL)` is `c(Inf, -Inf)`, so the failure surfaces as
  ## "'from' must be a finite number" from `seq()` with nothing naming the
  ## cause. Every model in the set was fitted to the same data, so the first
  ## one answers for all of them.
  bf <- if (inherits(fit, "bayesmanecfit")) fit$mod_fits[[1]]$fit else fit$fit
  x_rng <- range(bf$data$x)
  x_seq <- seq(x_rng[1], x_rng[2], length.out = resolution)
  pred <- stats::fitted(fit, newdata = data.frame(x = x_seq),
                        re_formula = NA, summary = TRUE)
  med <- pred[, "Estimate"]
  below <- which(med <= 0)
  if (length(below) == 0) Inf else x_seq[min(below)]
}
