# Resume here

Updated 2026-08-18. **Nothing is running.** Safe to power off.

**The simulation is complete: 12 of 12 cells, 17,280 fits, 0 worker failures.**
The last two cells landed 2026-08-17 (287.1 min and 321.5 min); both verified
against their expected shape, not just their `[done]` line. `phase5_report.R`
has been re-run over all 12 cells, `renv.lock` is written, and the noise-axis
read-out is in `analysis/phase5_r_axis.txt`.

**Phases 1-5 are finished.** Plan revision 6 (2026-08-18) adds two more:
Phase 7 brings in the zero-bounded families that floor implicitly (arms E and
F), and Phase 8 sequences the remaining compute so the vignette is not blocked
behind it. See "What is left".

## Starting a fresh Claude session

Point it at this file and say something like *"read RESUME.md in
/mnt/c/Rworking/negative-sgr and let's finish this"*. Everything needed is
below. Two things it will not know unless told:

- The prompt log required by `/mnt/c/Rworking/CLAUDE.md` §10 lives at
  `../bayesnec/prompts/negative-sgr-study.md` — append, don't start a new file.
- **Read "Traps" below before running anything.** Each one cost this project
  real time.

## Re-running things (nothing needs re-running as of 2026-08-17)

The sweep is complete. Should a cell ever need rebuilding, delete its `.rds` from
`analysis/phase5/` and re-issue the command below — completed cells are skipped
automatically, so it is safe to run at any time.

```bash
cd /mnt/c/Rworking/negative-sgr
DESIGN=rsep N_ITER=240 WORKERS=18 nohup Rscript analysis/phase5_run.R \
  > analysis/phase5_run_rcells.log 2>&1 &
```

**Verify it is actually running** — do not trust a launch message:

```bash
pgrep -fc "phase5_ru[n].R"                 # expect ~20
stat -c %y analysis/phase5_run_rcells.log  # should be seconds old
```

`analysis/overnight_finish.sh` waits for the sweep to exit and then re-runs the
report, writes `phase5_r_axis.txt` and `renv.lock`, and verifies every cell's
shape before it will overwrite `phase5_metrics.csv`. Launch it detached with
`setsid nohup ./analysis/overnight_finish.sh > /dev/null 2>&1 < /dev/null &`.

## State

| phase | status |
|---|---|
| Plan | **revision 6**, `../bayesnec/ignore/negative-sgr-study-plan.md` (rev 5 backed up in `../bayesnec/superceded/`). §Phase 4 and §Phase 5 carry "Observed result" subsections; §Phase 7 and §Phase 8 are the outstanding work |
| 1 branch verification | **done** — six gates pass, `analysis/phase1_gates.csv` |
| 2 diagnostics | **done** — `analysis/dataset_summary.csv` regenerated post-LOD |
| 3 arms on real data | **done** — rebuilt post-LOD, no errors, boundary flags corrected |
| 4 `bot` prior sensitivity | **done** — contraction table plus the prior sweep, `analysis/phase4_prior_sweep.{R,csv,log}` |
| 5 simulation | **done — 12 of 12 cells, 17,280 fits, 0 worker failures.** Only non-zero exclusion anywhere: arm C, 5 of 240 iterations in `d4.0_t0.8_R2.3` |
| 5 report | **done over all 12 cells** — `phase5_metrics.csv`, `phase5_metrics_cleanfits.csv`, `phase5_r_axis.txt`, `phase5_report.log` |
| 6 vignette | `example7.Rmd` written and precompiled on `negsgr-cens-vignette` (`0122ba5f`), **not** for `dev`. Awaiting arms E and F before it is final. The `example1` censoring edits are commit `4df470ea` |
| 6 paper artefacts | **not started** |
| 7 zero-bounded families | **not started** — arms E and F, §Phase 7 |
| 8 iteration top-up | **not started** — 240 to 500, §Phase 8 |
| environment | `renv.lock` written (106 packages, R 4.6.1) plus `analysis/session_info.txt` |

125 tests pass (`tests/testthat/`) — re-run 2026-08-17 via `cd tests && Rscript
testthat.R`. Running `testthat::test_dir()` from the project root fails: the
runner sources `../R/*.R` and the tests do not load them themselves.

## What is left

**Plan revision 6 (2026-08-18) adds Phases 7 and 8.** Read
`../bayesnec/ignore/negative-sgr-study-plan.md` §Phase 7 and §Phase 8 before
starting; the summary below is only a pointer.

1. **Phase 7, stage 1 — arms E and F at 240 iterations** (~20-24 h, 5,760 fits).
   E = floor negatives, divide by max, `Beta(link = "identity")`, `nec3param`.
   F = floor negatives, no scaling, `Gamma(link = "identity")`, `nec3param`.
   Both take `bnec()`'s own default priors, held fixed within a cell, and let
   `check_data()` do the boundary nudge. Write to
   `analysis/phase5/<cell>__ef.rds` — **additive, never rewrite an existing
   cell file.**
2. **Finalise the vignette** once stage 1 lands. At that point all eight arms
   have 240 iterations everywhere and the comparison is complete and balanced.
3. **Phase 8, stage 2 — top every arm up to 500 iterations** (~3.5 days,
   24,960 fits). Iterations 241-500, all eight arms, to
   `analysis/phase5/<cell>__topup.rds`. Seeds `7e5 + i` keep these distinct from
   1-240 with no further thought. This buys Monte Carlo precision only and
   cannot change an ordering, which is why it is last.
4. **Paper artefacts** (§Phase 6). Lead with the noise axis, not the `delta`
   gradient.
5. **Re-run the case studies under model averaging.** Phase 3 fixes `nec4param`
   for cross-arm comparability, which is right for the simulation (it is the
   generating model) and not defensible on real data. Deferred, not forgotten.
6. **The vignette branch has not been reviewed by anyone but Claude.** It is
   parked deliberately and nothing goes to `dev` without your say.

**Reporting changes stage 1 forces:** the six-colour arm palette needs two more
colours and re-validation, not extension by eye; the two-panel split becomes
three ("measurement retained", "zero boundary imposed", "zero bounded by the
family"); and the divergence and exclusion tables must cover E and F, whose
sampling behaviour is not predictable in advance.

## Open decisions

- **`renv.lock` is a record, not an renv project.** Written with
  `renv::lockfile_create()` / `lockfile_write()`, so it captures the dependency
  state without creating `renv/`, a project library, or an `.Rprofile`. Running
  `renv::init()` properly would change how every future R session in this
  project behaves and was deliberately left to you.
- Note `/mnt/c/Rworking/CLAUDE.md` records WSL R as 4.5.2; it is actually 4.6.1,
  which is what the lockfile pins.

## Findings

### Phase 5 — the core result

`phase5_report.R` now performs the regime split itself (§HEADLINE: ARMS BY
REGIME), so these numbers come straight out of its log rather than being
recomputed by hand.

**ErC50** where the design reaches the negative region — mean over the three
`R = 2.3, top_factor = 2.0` cells (`f_neg` 0.151):

| arm | what it does | bias | coverage | note |
|---|---|---|---|---|
| C | left-censor negatives | **-0.1%** | 0.829 | matches A — censoring recovers what flooring loses |
| A | raw, `bot` free | **+0.7%** | 0.915 | reference; the only arm near nominal coverage |
| D | truncate at crossing | -1.6% | 0.853 | |
| B1 | floor to 0, `bot` free | -4.8% | 0.703 | |
| B3 | floor + `bot` = 0 | -6.4% | 0.608 | worst coverage |
| B2 | raw + `bot` = 0 | -9.3% | 0.967 | worst bias; 4-6 divergences/fit. Its coverage is bought with width — RMSE 0.35 against A's 0.13 |

**[2026-08-17] Arm A's entry was previously recorded as -0.6%. Recomputed from
`phase5_metrics.csv` it is +0.7%; every other entry reproduces exactly. The
conclusion is unchanged — A is the only arm inside +/-1% — but do not requote
the old figure.**

The flooring penalty **grows with `delta`** (B1: -2.5 -> -5.5 -> -6.3, B3: -4.1
-> -7.4 -> -7.6, B2: -8.2 -> -10.0 -> -9.8 as delta goes 2 -> 4 -> 8) while A
and C stay inside +/-3.6% and show no trend. That gradient is the mechanism the
four real datasets were too few to show, and it is the strongest single piece of
evidence in the study.

**The effect is regime-dependent.** Where the design stops at or below the zero
crossing (`f_neg` <= 4%), every arm carries a shared -3 to -7% bias from weak
identification of `bot`, the arms are indistinguishable, and B1 can look
*better* than A. That is not flooring helping — it is arm A carrying a baseline
bias that flooring happens to offset. Once the design identifies `bot`, A's bias
vanishes and the flooring bias stands out. Do not quote the narrow cells alone.

### The noise axis is the strongest result — read `analysis/phase5_r_axis.txt`

Because the sweep holds the residual scale absolute, rising `R` is a pure
**precision** sweep at constant `f_neg` (0.143-0.152). Reading
`R` = 2.3 -> 3.3 -> 17 -> 73, i.e. falling noise:

| arm | ErC10 bias | ErC50 bias | ErC50 coverage |
|---|---|---|---|
| A | +19.0 -> +12.0 -> +3.4 -> **+1.4%** | -0.3 -> -0.8 -> -0.6 -> **-0.4%** | 0.94 -> 0.93 -> 0.88 -> **0.87** |
| C | +11.1 -> +7.1 -> +2.4 -> **+1.0%** | -0.6 -> -1.4 -> -1.3 -> **-1.0%** | 0.83 -> 0.83 -> 0.79 -> **0.78** |
| D | +13.6 -> +8.1 -> +2.5 -> **+1.0%** | -2.0 -> -2.0 -> -1.3 -> **-0.9%** | 0.84 -> 0.83 -> 0.80 -> **0.80** |
| B1 | +17.9 -> +15.9 -> +13.3 -> **+13.2%** | -5.5 -> -7.0 -> -8.1 -> **-8.2%** | 0.70 -> 0.49 -> 0.03 -> **0.000** |
| B3 | +22.6 -> +20.8 -> +16.3 -> **+15.0%** | -7.4 -> -8.8 -> -10.0 -> **-10.2%** | 0.59 -> 0.33 -> 0.004 -> **0.000** |
| B2 | +54.7 -> +52.0 -> +50.0 -> **+49.6%** | -10.0 -> -10.6 -> -11.2 -> **-11.3%** | 0.97 -> 0.99 -> 1.00 -> **1.000** |

**A/C/D converge to zero bias; B1/B2/B3 converge to non-zero asymptotes.** A
quantity that does not vanish as noise vanishes is misspecification, not
estimation error. This is a stronger argument than the `delta` gradient because
it needs no comparison across cells — each arm is its own control.

**Coverage under flooring collapses to 0 of 240.** The interval shrinks with
precision while its centre stays displaced, so a better experiment makes a
floored analysis *more confident and no less wrong*. That is the paper's
sentence. B2 fails in mirror image — coverage 1.000 with an 11%-low point
estimate — so **never report coverage without RMSE**, or the two worst arms
score as the two best.

### The ErC10 "reversal" was a high-noise artefact — resolved

At `R` = 2.3 alone every arm is biased high on ErC10 and arm A places fifth of
six (C +15%, D +17%, B1 +20%, B3 +24%, A +24%, B2 +57%). That is a shared
estimation difficulty swamping the arm effect, not a reordering: ErC10 needs a
10% decline from `top` resolved against the residual, and once precision is
adequate the ErC10 ordering *is* the ErC50 ordering, C ~ D ~ A << B1 < B3 << B2.
**Do not quote the `R` = 2.3 ErC10 row alone.**

### What to recommend

**Arms A and C are close and both sound; the claim is that B1/B2/B3 are not** —
not that one of A or C wins. At adequate precision they are indistinguishable on
bias (~1% on both endpoints). On coverage they trade: A is better on ErC50
(0.87-0.94 vs C's 0.78-0.83) and C is better on ErC10 (0.90-0.95 vs A's
0.83-0.98). **Neither reaches nominal coverage on ErC50 — say so plainly.** Arm
C's advantage is practical, not statistical: it still works when the negative
values were never recorded, which arm A does not. Its cost is the identification
failure in Phase 3 and the vignette — where the asymptote is entirely censored,
`bot` is the prior's and must be reported as such.

### Phase 4 — the A-vs-B2 gap is not a prior artefact

Arm A on `c_proliferum` under `bot` priors differing only in scale, against a
fixed arm B2 (`analysis/phase4_prior_sweep.csv`):

| `bot` prior | prior SD | `bot` posterior | ErC10 | ErC50 |
|---|---|---|---|---|
| tighter | 0.101 | -0.504 (SD 0.070) | 1.643 | 3.164 |
| **default** | 0.404 | **-1.022** (SD 0.239) | **1.622** | **3.324** |
| wider | 1.617 | -2.393 (SD 0.915) | 1.562 | 3.424 |
| arm B2 | — | 0 (fixed) | 0.922 | 2.503 |

A sixteen-fold change in the prior scale moves `bot` by a factor of 4.8 — its
posterior SD tracks the prior's almost exactly, contraction 0.32-0.44 throughout
— yet moves ErC50 by at most 4.8% against a 33% gap to B2. **The endpoints are
robust where the parameter that generates them is not.**

**But the binary containment test is not robust.** B2's ErC50 estimate falls
outside arm A's interval under the default and wider priors and *inside* it
under the tighter one, sitting only 0.12 below A's default lower bound. Report
ratios and their uncertainty, not containment flags — the flag compares two
intervals whose widths differ by 2.4-3.1x and inherits that asymmetry. The
Phase 3 gate (13 of 18) is unaffected, but no single combination is decisive.

### Two outputs withdrawn on evidence

- **The `R` hypothesis is untestable** on the growth-rate scale, not merely
  unsupported. `top`, `bot` and residual scale are all proportional to
  `mu_0 = log(R)/t` while `nec`, `beta` and the design do not involve `R`, so
  every reported endpoint is invariant. Verified to 3e-16. This explains Phase
  3's failed `R` ordering without appealing to n = 4. **The `R` cells were still
  worth running** — because the sweep holds the residual scale absolute they
  became a precision sweep, which is where the study's strongest result came
  from. Describe them as a noise axis, never as a fold-change axis.
- **NSEC coverage against the true `nec` is withdrawn.** Right in the limit
  (1.2985 vs 1.3 at negligible noise) but at realistic noise coverage falls to
  **0 of 240** in the widest cells on *every* arm including A, bias +120% to
  +266%. NSEC is defined against the posterior spread of the control response,
  so `nec` is its limit, not its expectation. Report **arm contrasts only**.

ErC10/ErC50 are unaffected: `ecx()` matches `true_ecx()` within 0.006 at
negligible noise (`analysis/phase5_estimand_check.R`).

### Phase 3 — real data

- **Gate passes**: arm A's estimate falls outside the convention arms' intervals
  in **13 of 18** usable combinations.
- **Arm B2 is unusable on all four datasets** for at least one endpoint.
- **`bot` contraction is the cleanest supporting result** — flooring pins the
  asymptote (0.966-0.997) where arm A has 0.39-0.42.
- **`Delta` does not order the four datasets** (Spearman 0.40; `R` gives -0.20),
  even after the LOD correction moved `r_salina2` from last to second. n = 4 is
  close to powerless; the simulation settles it.
- **Arm C leaves `bot` unidentified where the asymptote is entirely censored.**
  On `r_salina`: substitution gives -1.985 [-2.013, -1.958] SD 0.014;
  left-censoring -5.64 [-11.0, -2.51] SD 2.23 (the prior's tail; 5 divergences
  in 8000, Rhat 1.00 — nothing a routine check would stop on); interval
  censoring -2.45 [-2.68, -2.10], from a fit that needed a longer warmup and
  smaller step size. **This is the vignette.**

## Traps

Each of these cost real time. Read before running anything.

1. **`pkill`/`pgrep` patterns match your own shell.** `pkill -f phase5_run.R`
   kills the invoking bash too (exit 144). Use a bracket class:
   `pkill -f "phase5_ru[n].R"`. This has bitten the project three times,
   including killing a monitor.
2. **`tail -f` does not work on `/mnt/c`** — it is a 9p mount with no inotify
   and errors with "No data available". Poll instead.
3. **The `targets` store is on ext4 at `/home/rfisher/negsgr_targets`**, not in
   the project. `_targets.yaml` points at it. It was moved because 9p failed
   read-back verification on large objects under load ("Error storing output:
   file read error"), twice. Do not move it back.
4. **`bnec()` has no `brm_args` argument.** Everything for `brm()` goes through
   `...`; a `brm_args = list(...)` silently becomes one unused argument and
   `chains`/`backend`/`init` revert to defaults with no warning.
5. **Priors are compiled into the Stan program as literals.** `bnec()`'s
   gaussian defaults are functions of the response, so a per-dataset prior means
   a recompile per fit (3-5 min each). The sweep holds the prior fixed within a
   cell and warms the cache serially. Measured immaterial: switching the prior
   moves an endpoint by 3-12% of its Monte Carlo spread.
6. **Do not filter simulation fits on `divergences == 0`.** Arm A averages ~0
   and B2 averages 4-6, so the filter flatters the worst arm. `phase5_report.R`
   reports all successful fits with a clean-fits sensitivity table alongside.
7. **`family = "Beta"`, not `"beta"`** — base R's `beta()` is the beta function
   and fails with "unused argument (link = 'identity')".
8. **A detection limit is a protocol fact, not a data fact.** The smallest
   observed positive density bounds it from above; it does not identify it. That
   error put `r_salina2`'s LOD at 100 instead of 10. Hence
   `min_detected_density()`, deliberately not named `infer_lod()`.
9. **Check artefact timestamps, not exit codes.** A `targets` run reported "42
   completed" while leaving a stale pre-correction object behind, and a sweep
   reported `[done]` on a cell where 239 of 240 iterations had failed.
10. **`parameter_table()` returns `"bot_Intercept"`, not `"bot"`.** It names
    rows with `sub("^b_|_Intercept$", "", p)`, and `sub()` replaces only the
    *first* match, so the suffix survives. `pt$mean[pt$parameter == "bot"]` is
    therefore always length zero and fails downstream with "arguments imply
    differing number of rows". Match on `sub("_Intercept$", "", parameter)`.
    Note also that arms B2/B3 have no `bot` column at all — Stan does not
    declare a parameter whose prior is `constant(0)` — so handle length 0
    rather than assuming length 1.
11. **Run the tests as `cd tests && Rscript testthat.R`.** `testthat.R` sources
    `../R/*.R`; the test files do not, so `testthat::test_dir("tests/testthat")`
    from the project root fails with "could not find function `sim_truth`" and
    looks like a broken suite.

## Environment

- Pinned to `bayesnec` `dev` @ `374e511c`, loaded from
  `/mnt/c/Rworking/bayesnec-issue173`, **detached** so it cannot drift;
  `load_bayesnec()` asserts the SHA and stops if it moved.
- `STUDY_CORES = 18` of 22.
- Worktrees: main `/mnt/c/Rworking/bayesnec` on `dev`;
  `/mnt/c/Rworking/bayesnec-negsgr` on `negsgr-cens-vignette` (parked, HEAD
  `4df470ea`, committed locally and **not pushed**).
  `beta-ub-impl` and `issue-191-dispersion` belong to another session —
  **leave alone**, including the untracked `example7` drafts in
  `/mnt/c/Rworking/bayesnec-dispersion`.
- **Nothing goes to `dev`** without the user saying so. That is a standing
  decision from 2026-08-15, not an oversight.
