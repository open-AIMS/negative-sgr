# Resume here

Updated 2026-08-16 11:20. **Safe to power off at any time.** The sweep writes one
`.rds` per cell atomically when that cell finishes, so a completed cell cannot be
damaged and an interrupted cell is simply recomputed. Nothing else is running.

## To restart after power-on

```bash
cd /mnt/c/Rworking/negative-sgr
DESIGN=rsep N_ITER=240 WORKERS=18 nohup Rscript analysis/phase5_run.R \
  > analysis/phase5_run.log 2>&1 &
```

Completed cells are skipped automatically. As of shutdown **10 of 12 cells are
done**; cells `d4.0_t2.0_R17.0_s8.1` and `d4.0_t2.0_R73.0_s8.1` remain, about
3.5 h each on 18 workers.

**These two cells are optional.** They extend the R axis only, and `R` is now
known to be a pure signal-to-noise sweep rather than a test of control
fold-change (the generating model is exactly scale-equivariant in the growth
rate). Two of the four R points are already in. Every primary finding is
complete without them. Skip them if the machine is needed for something else.

Verify the sweep is genuinely running, rather than trusting a launch message:

```bash
pgrep -fc "phase5_run[.]R"          # expect ~20
stat -c %y analysis/phase5_run.log  # should be seconds old
```

## State

| phase | status |
|---|---|
| Plan review | **done** — plan at revision 3, `../bayesnec/ignore/negative-sgr-study-plan.md` |
| 1 branch verification | **done** — six gates pass on `dev`, `analysis/phase1_gates.csv` |
| 2 diagnostics | **done** — `analysis/dataset_summary.csv`, three figures |
| 3 arms on real data | **done** — rebuilt post-LOD, no errors, boundary flags corrected |
| 4 `bot` prior sensitivity | **partly done** — contraction table done; the prior *sweep* is **not** run |
| 5 simulation | **10 of 12 cells done**, 0 worker failures across ~14,400 fits; 2 optional R cells left |
| 6 vignette | **done, parked** — branch `negsgr-cens-vignette`, not for `dev` yet |
| 6 paper artefacts | **not started** |

125 tests pass (`tests/testthat/`).

## What is running

```
DESIGN=rsep N_ITER=240 WORKERS=18 Rscript analysis/phase5_run.R
```

12 cells x 240 iterations x 6 arms = **17,280 fits**, writing one `.rds` per cell
to `analysis/phase5/` as each finishes. Existing cells are skipped on restart, so
an interruption resumes rather than restarts. Progress:

```bash
grep -E "^\[done\]|^\[skip\]" analysis/phase5_run.log
ls analysis/phase5/ | wc -l          # cells finished, out of 12
```

Projection is 12-48 h wall depending on whether the pilot's 270 s/iteration was
measured per-core or per-4-chains; the pilot ran under contention at 6 workers,
so treat it as a range until the first cell lands. Summarise with
`analysis/phase5_report.R`.

To stop it: `pkill -f "analysis/phase5_run.R"` — note the pattern must not match
your own shell, which has bitten this project twice.

## The headline finding from this session

**The control fold-change `R` is not an axis on the growth-rate scale.** With
`sigma` proportional to `mu_0`, the generating model is exactly equivariant under
rescaling the growth rate, so every arm returns the same answer at every `R`
(verified: `y(R=73) = 5.1512 * y(R=2.3)` to 3e-16, identical rows floored,
identical `f_neg`). The `rsep` R sweep would have produced three copies of one
cell.

This settles Phase 3's failed R-ordering prediction analytically — not
underpowered at n = 4, structurally impossible — and it is fixed for the sweep by
`sigma_mode = "absolute"` (counting error lives on the log-density scale, so SGR
noise does not scale with the growth rate). Anchored at `R = 2.3`, so the nine
core cells are numerically unchanged. Full argument in `SESSION.md`.

Interpret the R axis as a **signal-to-noise axis**: varying `R` at fixed `sigma_0`
is the same manipulation as varying `sigma_0` at fixed `R`.

## Decisions taken 2026-08-14

1. **Both Rhodomonas detection limits are 10** (protocol-confirmed). `r_salina2`
   was recorded at 100 on the smallest-observed-density reasoning, which bounds
   the limit rather than identifying it. Its bound moves -1.147 -> **-1.915**.
   `infer_lod()` is renamed `min_detected_density()` to stop the mistake
   recurring. **Phase 3 is being rebuilt in full** (`dataset_meta()` feeds every
   target, so `targets` invalidated all of them; that also keeps the table on one
   code path). Log: `analysis/phase3_run_lod10.log`.
2. **Nothing goes to `dev` for now.** The vignette work stays parked on
   `negsgr-cens-vignette` (worktree `/mnt/c/Rworking/bayesnec-negsgr`), unmerged
   and un-PRed, to be revisited alongside the `example7` rewrite (#193).

**Check after the rebuild:** the vignette quotes `r_salina` `bot` values
(-1.985 / -5.64 / -2.45 and contractions 0.996 / 0.379 / 0.953). `r_salina`'s LOD
is unchanged and the seeds are fixed, so these should reproduce -- but they came
from the pre-rebuild table and must be re-checked against the new
`analysis/phase3_parameters.csv` and `phase4_bot_contraction.csv` before that
branch goes anywhere.

## Results so far (Phase 3, unchanged)

- **The gate passes.** Arm A's estimate falls outside the convention arms'
  intervals in 14 of 24 dataset x arm x endpoint combinations.
- **The pre-registered ordering prediction fails** — and the `R` half of it is
  now known to be untestable by construction. The `Delta` half remains open and
  is what the sweep is for.
- **Arm B2 is broken, not merely biased** — divergences on all four datasets,
  treedepth saturation, `c_proliferum2` ErC50 342 [0.01, 429].
- **B1/B3 sample cleanly and move the answer materially**, non-overlapping
  intervals: `c_proliferum` ErC50 3.36 [2.64, 4.41] -> 6.35 [5.77, 7.01].
- **`bot` contraction is the cleanest supporting result** — flooring pins the
  asymptote (0.966-0.997) where arm A has 0.39-0.42.
- **Two results flagged unusable** in `analysis/phase3_endpoints.csv` — arm D's
  ErC50 on `c_proliferum` and `r_salina`, where `ecx()` returns `max(x)`.
- **Arm C leaves `bot` unidentified where the asymptote is entirely censored.**
  On `r_salina`: substitution gives -1.985 [-2.013, -1.958] with SD 0.014,
  left-censoring -5.64 [-11.0, -2.51] with SD 2.23 (the prior's tail, chains
  clean), interval censoring -2.45 [-2.68, -2.10]. This is the vignette.

## Open items, roughly in order

1. **Re-run `analysis/phase3_report.R`** once the rebuild finishes, and diff the
   tables against the pre-LOD versions — `r_salina2` should move, nothing else
   should.
2. **Phase 5 report** once the sweep lands — `analysis/phase5_report.R`.
3. **Phase 4 prior sweep** — 3 fits, deferred only because the machine is busy.
4. **`renv.lock` not written** — deferred again, since `renv::init()` writes a
   `.Rprofile` and a sweep is running. Do it when the machine is idle.
5. **Protocol confirmation** — both LODs are now confirmed (decision 1).
   `mu_0`, `t` and `n_0` are still recovered arithmetically rather than read from
   the protocols; they reproduce the supplied SGR column exactly, so this is a
   provenance formality rather than an open risk.
6. **Paper artefacts** — not started.

## Plan file needs revision 4

`../bayesnec/ignore/negative-sgr-study-plan.md` is at revision 3 and now carries
three statements known to be wrong. Revision 3 is backed up at
`../bayesnec/superceded/negative-sgr-study-plan-rev3.md`. To change:

1. **`r_salina2`'s LOD is 10, not 100** (the revision-3 header states 100). Its
   `d` lower bound moves from >= 1.147 to >= 1.915, so its `Delta` rises and the
   pre-registered `Delta` ordering may no longer read
   `c_proliferum > c_proliferum2 > r_salina > r_salina2`. **Recompute from the
   rebuilt `analysis/dataset_summary.csv` before rewriting the ordering** --
   do not hand-edit the number.
2. **"If divergence instead orders by `R`, the mechanism above is wrong."** That
   test cannot be run: on the growth-rate scale `R` has no path to any endpoint.
   The plan's requirement that "the simulation must carry `R` and `Delta` as
   separate factors" is only satisfiable by making `R` a signal-to-noise axis,
   which is what `sigma_mode = "absolute"` now does. State this as a result.
3. **The revision-3 header still says the vignette target is `example6`, not
   `example1`.** That claim was withdrawn on 2026-08-12 (it was read off
   `cens-impl`, which predates the commit adding the section to `example1`) but
   the header was never corrected.

Also worth folding in: arm D is degenerate wherever the design does not reach
the zero crossing (all `top_factor = 0.8` cells), and the simulation now holds
the prior fixed within a cell.

## Environment notes

- Pinned to `bayesnec` `dev` @ `374e511c`, loaded from
  `/mnt/c/Rworking/bayesnec-issue173`, **detached** so it cannot drift.
  `load_bayesnec()` asserts the SHA.
- `STUDY_CORES = 18` (was 6; the competing study has finished). 22 cores total.
- Worktrees: main `/mnt/c/Rworking/bayesnec` on `dev`;
  `/mnt/c/Rworking/bayesnec-negsgr` on `negsgr-cens-vignette`. `beta-ub-impl` and
  `issue-191-dispersion` are the other session's — leave alone.
- `tail -f` does not work on `/mnt/c` (DrvFs, no inotify). Poll instead.
