# Negative growth rates in microalgal concentration–response analysis

Does forcing the lower asymptote of a concentration–response curve to zero growth
— the Australian laboratory convention of substituting 0 for negative specific
growth rate, equivalently capping percent inhibition at 100% — bias ErC10, ErC50
and NSEC, and does it degrade their interval coverage?

The plan this repository implements is
`../bayesnec/ignore/negative-sgr-study-plan.md` (revision 2). Provenance,
pinned commit and the package-behaviour facts that shaped the design are in
[`SESSION.md`](SESSION.md).

## Layout

```
R/                    functions only, no side effects
  setup.R             pinned commit, core cap, Stan compile cache
  data_prep.R         reading the CSVs and the five data preparations
  diagnostics.R       per-dataset summaries, gap ratio, reversals
  arms.R              one function per arm, common signature
  simulate.R          simulation truth, designs, data generation
  metrics.R           endpoint extraction, diagnostics, MCSE
  figures.R           plots
_targets.R            Phase 3/4 pipeline: one target per fit
analysis/
  phase1_verify.R     branch verification (six gates)
  phase2_diagnostics.R  dataset summaries and figures
  phase3_report.R     arm comparison and the paper gate
  phase5_pilot.R      timing and budget, run before the sweep
  phase5_run.R        the simulation sweep
  phase5_report.R     bias, RMSE, coverage with MCSE
data-raw/             the four CSVs, read-only
tests/                testthat, encoding the known data issues
```

## Running it

```r
# 1. verify the pinned bayesnec branch (six gates; ~10 min first time,
#    ~2 min once the Stan compile cache is warm)
Rscript analysis/phase1_verify.R

# 2. dataset diagnostics and figures (seconds, no fitting)
Rscript analysis/phase2_diagnostics.R

# 3. and 4. fit every arm on every dataset, then report
Rscript -e 'targets::tar_make()'
Rscript analysis/phase3_report.R

# 5. measure before committing, then sweep
Rscript analysis/phase5_pilot.R
N_ITER=50 DESIGN=core Rscript analysis/phase5_run.R
Rscript analysis/phase5_report.R
```

Tests: `cd tests && Rscript testthat.R`.

## The arms

|                | `bot` free | `bot` fixed at 0 |
|----------------|------------|------------------|
| raw SGR        | **A**      | **B2**           |
| floored SGR    | **B1**     | **B3**           |

plus **C** (left-censored at zero — censoring as coarsening), **D** (truncate the
design below the zero-crossing) and **SQ**, the literal lab practice reported as a
benchmark and never as an experimental arm.

Every arm of a dataset uses the Gaussian family, the `nec4param` model, and the
**same prior**, built once from that dataset's arm-A response. This matters:
`bnec()`'s default `bot` prior is `normal(quantile(y, 0.1), 2.5 * sd(y))`, a
function of the response vector, so letting each arm take its own default would
confound the likelihood contrast under study with a prior contrast.

`bot` is held at zero with `prior(constant(0), nlpar = "bot")`, not by switching
to `nec3param`: `check_models()` drops every zero-bounded model under `gaussian`,
and `ecx()` refuses `type = "absolute"` for a gaussian fit with no `bot`
parameter. Because `bayesnec:::make_inits()` cannot draw from a constant prior,
those arms supply their own initial values.

## What "raw" means per dataset

`r_salina` and `r_salina2` contain rows where the population fell below the
counting limit, so SGR is **undefined** rather than negative. There is no raw SGR
for those rows. Arms A, B2 and D therefore use the `bound` preparation on those
two datasets — substitute the detection limit, the other common lab convention —
and the substitution is declared rather than hidden inside "arm A". An
`A_raw_dropped` sensitivity arm bounds the effect of that choice.
