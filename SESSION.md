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
| cores used by this study | 6 (`STUDY_CORES` in `R/setup.R`) |

**Shared machine.** Another study (`R/reanalysis_fit_worker.R`, 4 worker
processes each running cmdstan chains) occupies roughly half this box for the
duration. Parallelism is capped so that timing figures are reproducible and the
other run is not starved. Any wall-clock number recorded here was measured under
that contention and is an upper bound.

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
- **The two LODs differ**: 10 cells/mL for `r_salina`, **100** for `r_salina2`.
  Each is confirmed arithmetically -- the bound it implies
  (`(ln LOD - ln n_0)/t` = -1.9862 and -1.1470) equals the SGR of the rows
  detected exactly at that limit, to all printed digits.
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

## Change log

| date | change |
|---|---|
| 2026-08-12 | Repository created. Plan revised to revision 2 after review. |
| 2026-08-12 | Phase 1 (all gates pass) and Phase 2 complete. |
| 2026-08-12 | Arm C measured, arm C2 added; plan revision 3. Phase 3 pipeline rerun from scratch. |
| 2026-08-12 | Repinned from `cens-impl` to `dev` after finding the Censoring section had moved into `example1` in a commit `cens-impl` predates. Phase 1 re-verified (identical), Phase 3 restarted. |
| 2026-08-13 | Arm C2 non-convergence on `r_salina` traced to one stuck chain and fixed by a general refit-once-on-Rhat rule in `fit_arm()`; Phase 3 re-run (1h25m) so every fit comes from one code path. Simulation truth recalibrated from the arm-A posterior after the hand-picked values proved to be a near-step curve. Phase 5 pilot and budget complete. Session paused; see `RESUME.md`. |

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
