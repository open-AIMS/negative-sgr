# Resume here

Updated 2026-08-16 11:25. **Safe to power off at any time.** Cells are written
atomically at completion, so a finished cell cannot be damaged and an
interrupted one is simply recomputed.

## Starting a fresh Claude session

Point it at this file and say something like *"read RESUME.md in
/mnt/c/Rworking/negative-sgr and let's finish this"*. Everything needed is
below. Two things it will not know unless told:

- The prompt log required by `/mnt/c/Rworking/CLAUDE.md` §10 lives at
  `../bayesnec/prompts/negative-sgr-study.md` — append, don't start a new file.
- **Read "Traps" below before running anything.** Each one cost this project
  real time.

## Restart the sweep (optional — see below)

```bash
cd /mnt/c/Rworking/negative-sgr
DESIGN=rsep N_ITER=240 WORKERS=18 nohup Rscript analysis/phase5_run.R \
  > analysis/phase5_run.log 2>&1 &
```

Completed cells are skipped automatically. **10 of 12 cells are done.**
Remaining: `d4.0_t2.0_R17.0_s8.1`, `d4.0_t2.0_R73.0_s8.1`, ~3.5 h each on 18
workers.

**These two cells are optional and no primary finding depends on them.** They
extend the `R` axis only, and `R` is a pure signal-to-noise sweep rather than a
test of control fold-change, because the generating model is exactly
scale-equivariant in the growth rate. Two of the four R points are already in.
Skip them unless the noise axis is wanted for the paper.

**Verify it is actually running** — do not trust a launch message:

```bash
pgrep -fc "phase5_run[.]R"          # expect ~20
stat -c %y analysis/phase5_run.log  # should be seconds old
```

## State

| phase | status |
|---|---|
| Plan | **revision 4**, `../bayesnec/ignore/negative-sgr-study-plan.md` (rev 3 backed up in `../bayesnec/superceded/`) |
| 1 branch verification | **done** — six gates pass, `analysis/phase1_gates.csv` |
| 2 diagnostics | **done** — `analysis/dataset_summary.csv` regenerated post-LOD |
| 3 arms on real data | **done** — rebuilt post-LOD, no errors, boundary flags corrected |
| 4 `bot` prior sensitivity | **partly done** — contraction table done; the prior *sweep* (default/wider/tighter, 3 fits) **not run** |
| 5 simulation | **10 of 12 cells**, 0 worker failures across ~14,400 fits |
| 6 vignette | **done, parked** on `negsgr-cens-vignette` — **not** for `dev` |
| 6 paper artefacts | **not started** |

125 tests pass (`tests/testthat/`).

## What is left, in order

1. **`analysis/phase5_report.R`** — run it once the sweep is final. It writes
   `phase5_metrics.csv` plus a clean-fits sensitivity table.
2. **Phase 4 prior sweep** — 3 fits, never run. Compare `bot` default/wider/
   tighter on one dataset to confirm the A-vs-B2 gap is not a prior artefact.
3. **Fold the Phase 5 results into the plan's conclusions** (§Phase 5) and write
   the paper artefacts.
4. **Re-check the vignette's quoted numbers** against the rebuilt
   `analysis/phase3_parameters.csv` and `phase4_bot_contraction.csv` before that
   branch goes anywhere. `r_salina`'s LOD did not change and seeds are fixed, so
   they should reproduce — but they were read off the pre-rebuild table.
5. **`renv.lock` not written.** Safe once nothing is running.

## Findings

### Phase 5 — the core result

Where the design reaches the negative region (`f_neg` ~15%, `top_factor = 2.0`),
ErC50 relative bias:

| arm | what it does | bias | note |
|---|---|---|---|
| A | raw, `bot` free | **-0.6%** | reference; coverage 0.938 |
| C | left-censor negatives | **-0.1%** | matches A — censoring recovers what flooring loses |
| D | truncate at crossing | -1.6% | |
| B1 | floor to 0, `bot` free | -4.8% | |
| B3 | floor + `bot` = 0 | -6.4% | |
| B2 | raw + `bot` = 0 | -9.3% | worst; 4-6 divergences/fit |

The flooring penalty **grows with `delta`** (B1: -2.5 -> -5.5 -> -6.3 as delta
goes 2 -> 4 -> 8), which is the mechanism the four real datasets were too few to
show.

**The effect is regime-dependent.** Where the design stops at or below the zero
crossing (`f_neg` <= 4%), every arm carries a shared -3 to -5% bias from weak
identification of `bot`, the arms are indistinguishable, and B1 can look
*better* than A. That is not flooring helping — it is arm A carrying a baseline
bias that flooring happens to offset. Once the design identifies `bot`, A's bias
vanishes and the flooring bias stands out. Do not quote the narrow cells alone.

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
  left-censoring -5.64 [-11.0, -2.51] SD 2.23 (the prior's tail, chains clean);
  interval censoring -2.45 [-2.68, -2.10]. **This is the vignette.**

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

## Environment

- Pinned to `bayesnec` `dev` @ `374e511c`, loaded from
  `/mnt/c/Rworking/bayesnec-issue173`, **detached** so it cannot drift;
  `load_bayesnec()` asserts the SHA and stops if it moved.
- `STUDY_CORES = 18` of 22.
- Worktrees: main `/mnt/c/Rworking/bayesnec` on `dev`;
  `/mnt/c/Rworking/bayesnec-negsgr` on `negsgr-cens-vignette` (parked).
  `beta-ub-impl` and `issue-191-dispersion` belong to another session —
  **leave alone**, including the untracked `example7` drafts in
  `/mnt/c/Rworking/bayesnec-dispersion`.
- **Nothing goes to `dev`** without the user saying so. That is a standing
  decision from 2026-08-15, not an oversight.
