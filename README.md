# Negative growth rates in microalgal concentration–response analysis

A simulation study of what happens to reported toxicity estimates when a
response that can legitimately be negative is prevented from being so.

## Background

Algal growth inhibition tests are among the most standardised assays in
ecotoxicology, and OECD Test Guideline 201 defines their primary endpoint as the
average specific growth rate,

*µ* = (ln *X<sub>j</sub>* − ln *X<sub>i</sub>*) / (*t<sub>j</sub>* − *t<sub>i</sub>*),

which is a rate of change and therefore lives on the real line. A population that
declines over the course of a test has a negative growth rate, and percentage
inhibition defined against a control accordingly exceeds 100% when it does.

Negative values are inconvenient. They cannot be expressed as a proportion of the
control without that proportion exceeding one, several widely used
concentration–response equations cannot generate them, and they look anomalous on
a plot. It is common practice to remove them: by substituting zero for the
measured value before fitting, by constraining the fitted curve so that its lower
asymptote cannot fall below zero, by doing both, or — most often of all — by
adopting a response distribution such as a Beta or a Gamma that has no support
below zero, so that the flooring occurs as a side effect of a choice made for
some other reason. TG 201 Annex 5 advises against the substitution explicitly, on
the grounds that it distorts the error distribution, but the practice persists
because the alternatives are less convenient.

This repository asks what those conventions cost. The question is not whether
they are tidy but whether the toxicity estimates subsequently reported — the
concentrations at which 10% and 50% effects occur, and the no-significant-effect
concentration — are biased, and whether their credible intervals retain their
nominal coverage.

## Design

The comparison is made on simulated data, because the truth must be known before
a convention can be scored rather than merely compared against another
convention. Data are generated from the same four-parameter threshold equation
that is subsequently fitted, so that the only misspecification under study is the
one each convention deliberately introduces.

Eight approaches are compared, falling into three groups.

| | data as measured | negatives floored at zero |
|---|---|---|
| lower asymptote free | **A** | **B1** |
| lower asymptote fixed at zero | **B2** | **B3** |

Three approaches retain the measurement: **A**, the intact analysis and the
reference; **C**, which declares the negative values left-censored at zero, using
the magnitude of neither; and **D**, which truncates the concentration series
below the fitted zero crossing and discards the remainder. Three impose the
boundary deliberately, by substituting data (**B1**), constraining the curve
(**B2**), or both (**B3**). Two impose it implicitly through the choice of
response distribution: **E** scales the floored response to a proportion and
fits a Beta, and **F** fits a Gamma to the floored response unscaled.

The first six share a single prior, derived once from the intact response, so
that a contrast between them isolates the likelihood rather than confounding it
with a change of prior. Approaches E and F take the package's own defaults for
their family, because taking the defaults is the practice under examination.

Twelve scenarios cross the depth of the lower asymptote with the position of the
highest tested concentration relative to the concentration at which the true
curve crosses zero growth, and then vary measurement precision at a fixed design.
The precision axis is the informative one: estimation error shrinks as an
experiment becomes more precise, and model misspecification does not, so the two
can be separated only by varying the noise and observing which discrepancies
persist.

Each scenario was simulated 500 times and every simulated dataset fitted by all
eight approaches — 48,000 model fits. A further phase re-ran five of the
approaches on three scenarios with the equation no longer fixed but averaged over
the candidate set the software would ordinarily fit, 200 iterations per scenario
and 27,600 additional fits, paired to the single-model results by shared random
seeds.

Four real algal growth tests, two species against two anonymised contaminants,
are analysed alongside. Only one of the two species is driven below the counting
resolution by the exposure, which is what makes the undefined-value problem
described under *Data* a live one for half the datasets and not the other. They establish whether the choice makes a practical
difference; they cannot establish which approach is right, because the true
toxicity of a real substance is unknown.

## Principal findings

**Where the concentration series stops matters more than which convention is
used.** Where the design does not reach the concentration at which the true curve
crosses zero, the lower asymptote is barely identified by any approach and all
eight carry a similar bias of the same sign. Only where the series runs past the
crossing do the approaches separate, and only those scenarios support a
comparison.

**The conventions that impose the zero boundary bias the effect concentration
downward, making a substance appear more toxic than it is.** Averaged over the
scenarios that reach past the crossing, the relative bias in ErC50 is
approximately −5% for flooring the data alone, −6% for flooring with the
asymptote pinned, −9% for pinning alone, and −9% for the Beta, against
approximately +1% for the intact analysis and −0.3% for left-censoring. Interval
coverage falls correspondingly, to 0.62 for the floored and pinned approach
against a nominal 0.95.

**The no-significant-effect concentration is affected more severely, and in the
opposite direction**, being biased upward by tens of per cent — that is, toward
declaring a substance safer than it is, which for a protective metric is the more
serious failure.

**The failures do not diminish as measurement improves; several worsen.** At the
finest precision tested, three of the approaches contained the true ErC50 in none
of 500 simulated datasets. A laboratory that improves its technique and reduces
its control variability will, if it floors its negative values, produce an
analysis that is more confident and no less wrong. The two approaches whose
distributional choice does the flooring produce the largest biases at high
precision while also producing the cleanest sampler diagnostics, so the failure
is silent.

**Fitting a candidate set rather than a single equation recovers much, but not
all, of the loss.** Averaging over the candidate set reduces the ErC50 bias of the
floored approaches by roughly two-thirds and restores much of their coverage,
with root mean squared error falling alongside, so the gain is accuracy rather
than intervals widened until they contain the truth. The no-significant-effect
concentration behaves differently in a way that depends on measurement precision:
averaging repairs it at realistic noise and reverses it at high precision, where
one approach moves from a bias of +53% to −63% with coverage of 0.01. This is the
same distinction the precision axis draws throughout — averaging removes the
component of the error that is estimation and exposes the component that is
misspecification. Notably, the equation that generated the data is almost never
selected in the constrained approaches, carrying a thousandth of the stacking
weight, and the effect concentration is recovered regardless.

## Repository organisation

`R/` holds functions and no side effects: the pinned commit and compute settings
(`setup.R`, `setup_phase10.R`), reading and preparing the source data
(`data_prep.R`), per-dataset diagnostics (`diagnostics.R`), one function per
approach with a common signature (`arms.R`), the simulation truth, designs and
data generation (`simulate.R`), endpoint extraction and Monte Carlo standard
errors (`metrics.R`), model-averaging helpers (`model_average.R`) and plotting
(`figures.R`).

`analysis/` holds one script per phase, each named for the phase it implements,
together with the results it produced as `.csv`. Phase 1 verifies the pinned
software against six gates before anything is fitted; phases 2 to 4 characterise
the real datasets and establish that the findings do not rest on the prior;
phases 5, 7 and 8 are the simulation sweep and its extension to the
distributional approaches; phase 9 re-analyses the real datasets under model
averaging; and phase 10 is the model-averaged simulation, with its own gates for
the software pin, the candidate sets and a pilot.

`data-raw/` holds the four source datasets, read-only. `tests/` encodes the known
data irregularities and the helper behaviour as `testthat` tests.

The raw simulation output of phase 10 is retained under `analysis/phase10/` as 36
block files. The corresponding output of the earlier sweeps is not retained, being
substantially larger; the derived results are, as `.csv`.

## Reproduction

Analyses are run against a frozen build of the modelling package rather than
whichever version happens to be installed, and the commit is asserted at load
time so that a run stops rather than silently producing results from a different
state. Phase 10 uses a second pin, identical to the first except for a fix to the
reproducibility of model-averaged output, without which interval coverage cannot
be measured because the intervals themselves vary between calls. `renv.lock`
records the full dependency state.

```sh
Rscript analysis/phase1_verify.R                       # software gates
Rscript analysis/phase2_diagnostics.R                  # dataset diagnostics
Rscript -e 'targets::tar_make()'                       # fit the real datasets
Rscript analysis/phase3_report.R
DESIGN=rsep N_ITER=240 WORKERS=18 Rscript analysis/phase5_run.R    # the sweep
N_ITER=240 CHUNK=40 WORKERS=18 Rscript analysis/phase7_run.R
N_ITER_FROM=241 N_ITER_TO=500 CHUNK=40 WORKERS=18 Rscript analysis/phase8_run.R
Rscript analysis/phase5_report.R
CELLS=12,8,9 N_ITER_FROM=1 N_ITER_TO=200 CHUNK=16 WORKERS=16 \
  Rscript analysis/phase10_run.R                       # model-averaged sweep
Rscript analysis/phase10_report.R
```

Every sweep skips work already on disk, so an interrupted run is resumed by
reissuing the same command. The full simulation is several days of computation on
eighteen cores. Tests run with `cd tests && Rscript testthat.R`.

Design decisions, the reasoning behind them, and the observed results of each
phase are recorded in the study plan; `SESSION.md` records provenance, the pinned
commit, and the package behaviours established before fitting that shaped the
design.

## Data

The four source datasets are laboratory growth tests of two microalgae — the
symbiotic dinoflagellate *Cladocopium proliferum* and the cryptophyte
*Rhodomonas salina* — against two contaminants identified only as A and B. The
two species grow on very different scales, which is part of why both are
informative here. The contaminant concentrations are on undisclosed
scales, preserve relative spacing only, and are not comparable between
contaminants. The same data are distributed as the `alga` dataset of the
`bayesnec` package.

Values recorded as zero density are not measurements of zero but records that the
count fell below the counting limit of 10 cells/mL, for which the growth rate is
undefined rather than negative. Analyses that require a value there substitute the
rate implied by the counting limit, and this substitution is declared in the
approach rather than concealed within it.
