# Session and provenance record

Everything in this repository is produced against a single frozen `bayesnec`
commit. A development branch shifting mid-study is the most likely source of an
unreproducible result, so the SHA is asserted at load time by
`load_bayesnec()` in `R/setup.R` and the run stops if the worktree has moved.

## Pinned dependency

| item | value |
|---|---|
| package | `bayesnec` |
| repository | https://github.com/open-AIMS/bayesnec.git |
| branch | `dev` |
| commit | `374e511c665fee04ca6ca8e7b48a547a3928a28d` |
| commit date | 2026-08-10 |
| local worktree | `/mnt/c/Rworking/bayesnec-issue173`, **detached** at that commit |

**Why `dev` and not `cens-impl`.** The study was first pinned to `cens-impl` on the
reasoning that it was the branch adding `cens()`. That was wrong: `dev` **contains**
`cens-impl` (merged as PR #189) along with seven other feature branches, and adds three
things this study needs —

- commit `ee76e815` "add a Censoring section to example1", which post-dates `cens-impl`
  and is the Phase 6 target. On `cens-impl` that section does not exist, which led to a
  withdrawn conclusion that the vignette material lived in `example6`;
- the normalisation detection in `check_data.R` (`check_normalisation()`,
  `on_rational_grid()`), which is what makes the Phase 2 rule "bayesnec will emit messages
  if it detects percent-of-control" true — the SQ benchmark divides by the control mean
  and should trip it;
- the matching `bnec()` guidance on supplying the raw response.

Only two branches in the stack are **not** in `dev`: `beta-ub-impl` (parked, PR #187
closed) and `issue-191-dispersion` (the distributional-sigma work, issue #191, which this
study needs only in the negative — Phase 5's heteroscedastic cells are generated
heteroscedastic and fitted homoscedastic precisely because `bnec()` has no route to a
distributional sigma). At the time of pinning `dev` equalled `origin/dev` with no
modified tracked files.

The study loads from a *detached* worktree rather than from the `dev` branch. `dev` was
moved into the main working directory on 2026-08-12 and will advance as work continues;
a detached checkout cannot, so the study keeps a frozen copy of the pinned tree and the
SHA assertion in `load_bayesnec()` cannot be satisfied by a drifted branch.

Re-running Phase 1 against `dev` reproduced all six gate results to the digit
(ErC50 352.4, NSEC 187.5, `bot` posterior range 0 to 0), so the repin changes no fitted
result — only which vignette the Phase 6 edit belongs in, and which package messages fire.

## Environment

| item | value |
|---|---|
| OS | Debian on WSL2 (Linux 6.18.33.2-microsoft-standard-WSL2) |
| R | 4.6.1 (2026-06-24) |
| backend | `cmdstanr` |
| cores available | 22 |
| cores used by this study | 18 from 2026-08-14; 6 before that (`STUDY_CORES` in `R/setup.R`) |

**Shared machine, until 2026-08-14.** Another study
(`R/reanalysis_fit_worker.R`, 4 worker processes each running cmdstan chains)
occupied roughly half this box through Phases 1-3 and the Phase 5 pilot.
Parallelism was capped at 6 so that timing figures stayed reproducible and the
other run was not starved, so **every wall-clock number recorded for those
phases was measured under contention and is an upper bound**. That study has
since finished and `STUDY_CORES` is 18, leaving 4 cores for interactive
sessions. Per-fit costs are per-core and carry over; wall-clock projections made
against the 6-worker pilot must be rescaled.

**Stan compile cache.** `cmdstanr` rebuilds every model in a fresh session,
which costs 3-5 minutes per model here. `use_compile_cache()` points
`cmdstanr_write_stan_file_dir` at `cmdstan_cache/`, where the file name is a
hash of the Stan code, so an unchanged model reuses its executable across
sessions. The study has only a handful of distinct Stan programs (arms A/B1/D
share one, B2/B3 share one, C one, SQ one).

## Package-behaviour facts established before fitting

Recorded here because each one changed the design of the study.

1. **`ecx(type = "absolute")` measures the decline from `top` to zero.**
   `ecx_x_absolute()` sets `range_y <- c(0, max(y))`. On the raw SGR scale this
   is exactly ErCx as OECD TG 201 defines it, and it means **zero growth is 100%
   effect by construction, for every dataset, whatever the control growth rate**.
   The `1 - 1/R` result is on the final-density scale and does not describe the
   reported endpoints. See the revised organising hypothesis in the plan.

2. **`check_models()` drops every `zero_bounded` model under `gaussian`** and
   errors if nothing survives, so `nec3param` is unavailable. **`ecx()`
   independently refuses `type = "absolute"` for a gaussian fit with no `bot`
   parameter.** Arms B2/B3 therefore cannot be built as `nec3param`; they need a
   `constant(0)` prior on `bot` in a `nec4param`.

3. **The default `bot` prior is a function of the response vector**:
   `normal(quantile(y, 0.1), 2.5 * sd(y))` under `prior_type = "uninformative"`
   (the `bnec()` default). It therefore differs between arms, because the arms
   differ in their data. On `r_salina` the default is `normal(0.0000, 2.1209)` --
   centred at exactly zero *because of the substituted zeros*. Every arm of a
   dataset is given the arm-A prior explicitly so the contrast is a likelihood
   contrast and not a prior contrast.

4. **`bnec()` has no `brm_args` argument.** Arguments destined for `brm()`
   travel through `...`. A list passed as `brm_args = list(...)` reaches `brm()`
   as one unused named argument, and `chains`, `backend` and `init` silently
   revert to their defaults with no warning. The first Phase 1 run was invalid
   for this reason: it ran on the rstan backend with 4 chains instead of
   cmdstanr with 2, and arm B2's hand-built inits were ignored, which produced
   the `fcts[[fct_i]](1, v1, v2) : attempt to apply non-function` failure when
   bayesnec built its own inits from a `constant(0)` prior it cannot sample.

5. **`bayesnec:::make_inits()` can only draw from gamma, normal, beta and
   uniform priors.** A `constant()` prior has no draw function, so arms B2/B3
   must supply `init` themselves; `bnec()` skips its init search when `init` is
   supplied.

## Data facts established before fitting

Full table in `analysis/dataset_summary.csv`; tests in
`tests/testthat/test-data_prep.R`.

- `t` and `n_0` recover from the density/SGR pairs with `R^2 = 1` on
  `c_proliferum` and `r_salina` and `> 0.99999` on the other two, so the supplied
  SGR column was computed deterministically from those densities and a single
  inoculum density. `t` = 7, 7, 3, 3 d; `n_0` = 7968, 8684, 3871, 3123 cells/mL.
- **The two LODs were taken to differ**: 10 cells/mL for `r_salina`, 100 for
  `r_salina2`. Each bound (`(ln LOD - ln n_0)/t` = -1.9862 and -1.1470) equals
  the SGR of the lowest detected row to all printed digits. **The `r_salina2`
  value is now in doubt** -- that arithmetic only shows the lowest row sits at
  the assumed limit, which is true of any assumed limit chosen this way. See
  "Open: `r_salina2`'s detection limit is probably 10, not 100" below.
- **Every zero in the two *Rhodomonas* SGR columns is a zero cell density**,
  where SGR is undefined rather than zero. The *Cladocopium* sets have none.
- **The supplied CSVs are not fully floored.** The labs substituted 0 for the
  undetected rows but left genuinely measured negatives intact (`r_salina` keeps
  -1.62 and -1.99; `r_salina2` keeps 17 negatives). The 100% cap is applied
  downstream in their curve fitting, not in the delivered data. The `floored`
  preparation is therefore *our* application of the convention, and the study
  says so rather than claiming the delivered data embody it.
- `c_proliferum` reverses in the declining limb: 22 of 25 replicate pairs at
  20 ug/L sit above those at 15 ug/L (dominance 0.88) with mean final density
  2.28x higher. Flagged, not modelled around.

## Phase 1 gate results

All six pass. `analysis/phase1_gates.csv` is the machine-readable record and
`analysis/phase1_run.log` the full session output.

| gate | result |
|---|---|
| 1 `cens()` end to end | PASS — ErC50 352.4, NSEC 187.5, censoring block present in the Stan code |
| 2 candidate set invariant | PASS — `check_models()` takes (model, family, data) only, so censoring cannot change the valid set; 8 decline models valid under gaussian |
| 3 per-row bounds | PASS — 10 of 70 rows reached brms as left-censored, matching the declaration |
| 4 default `bot` prior | PASS — `normal(-0.0426, 0.3137)` on `c_proliferum2` |
| 5 gaussian drops zero-bounded | PASS — `nec3param, ecxexp, ecxsigm, ecxwb1p3, ecxwb2p3, ecxll3` dropped; `nec3param` alone errors |
| 6 `constant(0)` on `bot` | PASS — `bot` posterior range 0 to 0; `ecx()` and `nsec()` both work on the result |

Gate 6 is the one that decided the design: it makes arms B2 and B3 buildable.

## Findings that changed the design after Phase 1

**Arm C does not identify `bot`.** Measured on `c_proliferum2` (4 chains, 2000
iterations, same data and same prior in both arms):

| arm | `bot` posterior | contraction |
|---|---|---|
| A | −1.27 [−1.64, −0.94] | 0.43 |
| C | −0.42 [−0.84, −0.16] | 0.44 |

The cause is saturation. A left-censored row contributes `Phi((0 - mu)/sigma)`,
and `sigma` here is 0.01–0.05, so that contribution is already ≈1 once `mu` is a
few hundredths below zero. The likelihood is flat in `bot` across the whole
region the data occupy, and the prior — centred at `quantile(y, 0.1) = -0.043` —
decides where the posterior sits. Censoring negatives at zero is a decision
**not to estimate** the asymptote; whatever is reported for it afterwards is the
prior's.

An earlier reading, that the censored likelihood rises without bound as `mu`
falls and therefore pushes `bot` downwards, is wrong and was withdrawn once the
fit was run. The contribution is monotone in `mu` but saturates, so there is no
push — only flatness. The observable symptom is speed: arm C takes about 1.7x
arm A at matched settings, and far worse at `adapt_delta = 0.99`.

**Saturation itself is not a new finding** — `example1`'s Censoring section on
`dev` already states it, and contrasts it with substitution, which "actively
pulls the curve back up toward the value that was substituted". What the
measurement adds is the consequence for a four-parameter model: saturation
leaves `bot` unidentified, so the fitted lower asymptote is reported from the
prior. The Phase 6 subsection continues that paragraph rather than restating it.

**Arm C2 added.** The same coarsening stated as an interval — the response
carries the growth rate implied by extinction (1 cell/mL) and `y2` carries the
upper bound (0 for a measured negative, the LOD bound for an undetected row).
This does not restore identification, which nothing can once the magnitudes are
discarded, but it bounds the flat region on both sides with a statement that is
true rather than on one side only. Arm C is retained: its behaviour is a result.

## The control fold-change R is not an axis on the growth-rate scale

Found on 2026-08-14, before the Phase 5 sweep was launched, by checking the cells
rather than trusting them.

Parameterised as the study parameterises it -- `mu_0 = log(R)/t`, `top = mu_0`,
`bot = -Delta * mu_0`, and `sigma_0 = cv_control * top` -- the generating model is
**exactly equivariant under rescaling the growth rate**. `nec`, `beta` and the
concentration design do not involve `R` at all, so changing `R` multiplies every
simulated response by a constant and leaves every x-axis quantity untouched.

Verified numerically at `delta = 4`, `top_factor = 2`, common seed:

| quantity | R = 2.3 | R = 73 |
|---|---|---|
| `y` | reference | `5.1512 x` reference, max deviation 3.3e-16 |
| design `x` | identical | identical |
| rows floored (`y < 0`) | identical | identical |
| `f_neg` | 0.1714 | 0.1714 |
| `top` prior | `normal(0.1266, 0.1921)` | `normal(0.6520, 0.9896)` = `5.1512 x` |
| `nec`, `beta` priors | identical | identical |

Every estimand is an x-axis quantity (ErC10, ErC50, NSEC, the zero crossing), and
every arm is scale-equivariant too -- flooring is `max(y, 0)`, censoring is
`y <= 0`, a fixed `bot` is `0`, and arm D truncates at an x-axis crossing. So
**every arm returns the same answer at every R**, and the `rsep` R sweep as
originally written would have produced three copies of one cell.

Two consequences.

**It settles the Phase 3 ordering failure analytically.** The pre-registered
prediction that divergence should order with `R` did not fail for want of power
at n = 4. Under this model it is structurally impossible: `R` has no path to any
endpoint. The `Delta` ordering remains an open empirical question; the `R`
ordering does not.

**The equivariance is an artefact of the sigma model, not the physics.** Tying
`sigma` to `mu_0` holds the *coefficient of variation of the growth rate* fixed.
But counting error lives on the log-density scale: with SD `s` on `ln N` at each
time point, `SGR = (ln N_t - ln N_0)/t` has SD about `s*sqrt(2)/t`, independent of
the growth rate. Two tests with the same counting precision and different control
growth rates have the *same* absolute SGR noise and *different* signal. Holding
`sigma_0` fixed is therefore what lets `R` enter, and it enters as
signal-to-noise -- which is plausibly what OECD TG 201's `R >= 16` validity
criterion is implicitly buying, and is worth saying in the paper.

`sim_sigma()` and `simulate_dataset()` gained `sigma_mode` (`"cv"`, the original,
kept so the equivariance stays demonstrable; `"absolute"`) and `sigma_0_at()`
anchors the absolute mode. Phase 5 anchors at `R = 2.3`, so the nine core cells
are numerically identical to the pilot -- checked, difference exactly 0 -- and
only the three R cells change. Across the sweep the control CV now falls 9.6% ->
6.7% -> 2.8% -> 1.9% as `R` goes 2.3 -> 3.3 -> 17 -> 73, and `f_neg` falls 0.171
-> 0.143. `f_neg` is now recorded per iteration as the mediating quantity.

**The R axis is a noise axis, and should be reported as one.** Under the
absolute-sigma parameterisation, varying `R` at fixed `sigma_0` is exactly the
same manipulation as varying `sigma_0` at fixed `R`, up to an overall rescaling
of the response -- the model only ever sees the ratio. So the sweep does not
test "does the control fold-change matter?" in any richer sense; it tests
signal-to-noise, and `R` is one of the two things that sets it. That is the
honest reading, and it is the interesting one: it means a claim like "high-R
tests resist flooring bias" is the same claim as "precise tests resist flooring
bias", and TG 201's `R >= 16` validity criterion is buying noise control rather
than anything specific to growth. Anything stronger than that would need `R` to
enter the mean function independently, which on the growth-rate scale it does
not.

Note that `f_neg` moves much less than the CV does. Most negative observations
come from the region where the true mean is genuinely below zero, and that region
is fixed by the design, which is scale-invariant. What `R` changes is how
*precisely* those negative values are measured. The mechanism to look for in the
results is therefore not "high R produces fewer negatives" but "high R produces
well-determined negatives, so discarding them destroys more information".

## Change log

| date | change |
|---|---|
| 2026-08-12 | Repository created. Plan revised to revision 2 after review. |
| 2026-08-12 | Phase 1 (all gates pass) and Phase 2 complete. |
| 2026-08-12 | Arm C measured, arm C2 added; plan revision 3. Phase 3 pipeline rerun from scratch. |
| 2026-08-12 | Repinned from `cens-impl` to `dev` after finding the Censoring section had moved into `example1` in a commit `cens-impl` predates. Phase 1 re-verified (identical), Phase 3 restarted. |
| 2026-08-13 | Arm C2 non-convergence on `r_salina` traced to one stuck chain and fixed by a general refit-once-on-Rhat rule in `fit_arm()`; Phase 3 re-run (1h25m) so every fit comes from one code path. Simulation truth recalibrated from the arm-A posterior after the hand-picked values proved to be a near-step curve. Phase 5 pilot and budget complete. Session paused; see `RESUME.md`. |
| 2026-08-14 | `r_salina2` LOD queried (10 vs 100), pending protocol confirmation. Scale-equivariance in `R` found and fixed before the sweep (`sigma_mode`, `sigma_0_at()`); `STUDY_CORES` raised to 18; Phase 5 `rsep` sweep launched at 240 iterations per cell. |

## Open: `r_salina2`'s detection limit is probably 10, not 100

Found 2026-08-14 while cross-checking against the package's shipped `alga`
dataset, which is these same four tests. **This changes a Phase 3 result and
needs protocol confirmation before it is either adopted or dismissed.**

`dataset_meta()` asserts `lod = c(NA, NA, 10, 100)`, and `infer_lod()` "checks"
it as the smallest positive recorded density. That is not an independent check --
it is the same inference twice, and it conflates the smallest density *observed*
with the smallest density *observable*.

The densities themselves argue against 100:

| test | lowest recorded densities | reading |
|---|---|---|
| `r_salina` x A | 0, 10, 30, 36210, ... | 3 rows sit at exactly 10, on a grid of 10. The limit is 10. |
| `r_salina` x B (`r_salina2`) | 0, 100, 230, 290, 330, 450, ... | 230 and 290 put the reporting grid at 10, not 100. One row at 100 is just the sample minimum. |

Same species, same duration, same laboratory, and test A demonstrably reports
densities of 10 and 30. There is no positive evidence that test B could not have
reported a density below 100; nothing simply fell in that window among its four
below-limit rows.

If the limit is 10, `r_salina2`'s bound moves from
`(log 100 - log 3123)/3 = -1.147` to `(log 10 - log 3123)/3 = -1.915`. That is
not cosmetic: arm A's `bot` on `r_salina2` is -1.188 [-1.335, -1.057] with an SD
of 0.071, pinned against the assumed bound, so the assumed value is most of the
answer. Arms A, B2 and D take the `bound` preparation and C2 takes it as the
interval's upper end, so four of the seven arms move on that dataset. The nine
`c_proliferum` fits and the whole simulation are unaffected -- the simulation
uses `sim_meta()`, where `lod` is `NA`.

Not changed unilaterally, because it is a fact about the counting protocol rather
than about the data, and the protocols have not yet been consulted (see the
standing protocol item in `RESUME.md`). The cost of adopting it is a re-run of
the `r_salina2` column of Phase 3, not of the whole pipeline.

Note the package's own `notes/alga_dataset.md` states a single counting
resolution of 10 for the dataset and gives one bound formula using it, which is
consistent with the reading above and inconsistent with this study's 100.

## Corrections made during the work, for the record

Three claims of mine were withdrawn after being tested rather than assumed. They
are listed because each one would have shaped a result if it had gone unchecked.

1. **"Arm C pushes `bot` down without bound."** Wrong. The censored contribution
   is monotone in `mu` but *saturates*, so the likelihood is flat rather than
   sloped and the prior decides. Withdrawn when the fit was run.
2. **"The Censoring section lives in `example6`, not `example1`."** Wrong — read
   off `cens-impl`, which predates commit `ee76e815`. The plan was right.
3. **"Arm C is pathologically slow."** Wrong; that was `adapt_delta = 0.99` on
   `c_proliferum` under the superseded pin. In the final run arm C took 8m33s
   against B2's 16m and D's 30m.

And one design error caught before it consumed compute: the simulation truth was
hand-set to `nec = 5, beta = log(0.4)`, which puts the whole transition from
control to zero growth inside 0.6 concentration units. No design resolves that,
every arm would have failed identically, and the simulation would have returned a
confident null. `resolves_transition()` now refuses any cell with fewer than three
tested concentrations between ErC10 and the zero crossing.
