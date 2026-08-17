# Resume here

Updated 2026-08-17 08:30. **Safe to power off at any time.** Cells are written
atomically at completion, so a finished cell cannot be damaged and an
interrupted one is simply recomputed.

**The sweep is running as of 2026-08-17 08:09** (`analysis/phase5_run_rcells.log`,
18 workers) on the last two cells, `d4.0_t2.0_R17.0_s8.1` and
`d4.0_t2.0_R73.0_s8.1`. Expect ~7 h total. Nothing else is waiting on it except
a re-run of `phase5_report.R` and `renv.lock`.

## Starting a fresh Claude session

Point it at this file and say something like *"read RESUME.md in
/mnt/c/Rworking/negative-sgr and let's finish this"*. Everything needed is
below. Two things it will not know unless told:

- The prompt log required by `/mnt/c/Rworking/CLAUDE.md` §10 lives at
  `../bayesnec/prompts/negative-sgr-study.md` — append, don't start a new file.
- **Read "Traps" below before running anything.** Each one cost this project
  real time.

## Restart the sweep, if it was interrupted

```bash
cd /mnt/c/Rworking/negative-sgr
DESIGN=rsep N_ITER=240 WORKERS=18 nohup Rscript analysis/phase5_run.R \
  > analysis/phase5_run_rcells.log 2>&1 &
```

Completed cells are skipped automatically, so this is safe to re-issue at any
time. Remaining: `d4.0_t2.0_R17.0_s8.1`, `d4.0_t2.0_R73.0_s8.1`, ~3.5 h each on
18 workers.

**No primary finding depends on these two cells.** They extend the `R` axis
only, and `R` is a pure signal-to-noise sweep rather than a test of control
fold-change, because the generating model is exactly scale-equivariant in the
growth rate. The headline reads the `R = 2.3` cells alone and is already final.

**Verify it is actually running** — do not trust a launch message:

```bash
pgrep -fc "phase5_ru[n].R"                 # expect ~20
stat -c %y analysis/phase5_run_rcells.log  # should be seconds old
```

## State

| phase | status |
|---|---|
| Plan | **revision 5**, `../bayesnec/ignore/negative-sgr-study-plan.md` (rev 4 backed up in `../bayesnec/superceded/`). §Phase 4 and §Phase 5 now carry "Observed result" subsections scoring the pre-registered predictions |
| 1 branch verification | **done** — six gates pass, `analysis/phase1_gates.csv` |
| 2 diagnostics | **done** — `analysis/dataset_summary.csv` regenerated post-LOD |
| 3 arms on real data | **done** — rebuilt post-LOD, no errors, boundary flags corrected |
| 4 `bot` prior sensitivity | **done** — contraction table plus the prior sweep, `analysis/phase4_prior_sweep.{R,csv,log}` |
| 5 simulation | **10 of 12 cells written; the last two running since 2026-08-17 08:09.** 0 worker failures across ~14,400 fits |
| 5 report | **done for the 10 cells** — `phase5_metrics.csv`, `phase5_metrics_cleanfits.csv`, `phase5_report.log`. **Re-run when the sweep lands** |
| 6 vignette | **done, parked** on `negsgr-cens-vignette` — **not** for `dev`. Numbers re-checked 2026-08-17; three corrected, commit `4df470ea` |
| 6 paper artefacts | **not started** (out of scope for the 2026-08-17 session by agreement) |

125 tests pass (`tests/testthat/`) — re-run 2026-08-17 via `cd tests && Rscript
testthat.R`. Running `testthat::test_dir()` from the project root fails: the
runner sources `../R/*.R` and the tests do not load them themselves.

## What is left, in order

1. **Re-run `analysis/phase5_report.R`** once the two R cells land. The headline
   is unaffected — it reads the `R = 2.3` cells only — so only the R/noise rows
   change.
2. **Read the R rows** and decide whether the noise axis earns a place in the
   paper. Two of four R points were already in; this run adds R = 17 and 73.
3. **Paper artefacts** (§Phase 6). Not started.
4. **`renv.lock` not written.** Safe once nothing is running.

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

### ErC10 tells a different story — do not claim arm A is uniformly best

In the same reaching regime, **every** arm is biased *high* on ErC10 and arm A
places fifth of six:

| arm | C | D | B1 | B3 | A | B2 |
|---|---|---|---|---|---|---|
| ErC10 rel. bias | +15% | +17% | +20% | +24% | +24% | +57% |

ErC10 asks for a 10% decline from `top` to be resolved against the noise, and
that shared estimation difficulty dominates the arm differences rather than the
reverse. **The claim that survives on both endpoints is that arm C is best or
joint-best and arm B2 is worst.** So the recommendation is arm C — left-censor
the negatives — not arm A. Arm C also remains available when the negative
values were never recorded, which arm A does not. Its one cost is the
identification failure in Phase 3 and the vignette: where the asymptote is
entirely censored, `bot` is the prior's and must be reported as such.

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
  3's failed `R` ordering without appealing to n = 4.
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
