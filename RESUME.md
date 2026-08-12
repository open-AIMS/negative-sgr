# Resume here

Written 2026-08-13. Nothing was running when the machine was shut down; every
result below is on disk. `SESSION.md` holds provenance and the package-behaviour
findings; `README.md` holds the layout and how to run each phase.

## State

| phase | status |
|---|---|
| Plan review | **done** — plan at revision 3, `../bayesnec/ignore/negative-sgr-study-plan.md`, original in `../bayesnec/superceded/` |
| 1 branch verification | **done** — all six gates pass on `dev`, `analysis/phase1_gates.csv` |
| 2 diagnostics | **done** — `analysis/dataset_summary.csv`, three figures, 109 tests pass |
| 3 arms on real data | **done** — 36 fits, `analysis/phase3_*.csv` |
| 4 `bot` prior sensitivity | **partly done** — contraction table done (`analysis/phase4_bot_contraction.csv`); the prior *sweep* (default/wider/tighter on one scenario) is **not** run |
| 5 simulation | **pilot done**, sweep **not run** — needs your budget decision, below |
| 6 vignette + paper | **not started** |

To pick up:

```r
setwd("/mnt/c/Rworking/negative-sgr")
targets::tar_make()                      # no-op, everything is cached
Rscript analysis/phase3_report.R         # regenerates every table and figure
```

The `_targets/` cache and `cmdstan_cache/` are both on disk, so nothing recompiles
and nothing refits.

## The decision waiting for you: Phase 5 budget

Measured per-fit cost, 4 pilot iterations, seconds (lean = 4 chains x 2000 iter,
`adapt_delta` 0.95; full = 4000 iter, 0.99):

| arm | lean | full |
|---|---|---|
| A | 15.6 | 42.7 |
| B1 | 39.5 | 73.1 |
| B2 | 7.9 | 69.8 |
| B3 | 10.4 | 28.7 |
| C | 42.2 | 76.4 |
| **D** | **154.4** | **171.0** |
| all six | 270 | 462 |

Iterations needed: 30 for MCSE +/-4%, 119 for +/-2%, 476 for +/-1% on 95% coverage.

| option | cells | iter | core-hours | wall on 6 cores | coverage MCSE |
|---|---|---|---|---|---|
| A | 9 (core) | 50 | 135 | 22 h | +/-3.1% |
| B | 12 (`rsep`, separates R from Delta) | 50 | 180 | 30 h | +/-3.1% |
| C | 12 (`rsep`) | 240 | 864 | 144 h | +/-1.4% |
| D | 72 (full factorial) | 240 | 8864 | 1477 h | +/-1.4% |

**Recommendation: option B, `DESIGN=rsep N_ITER=50`.** Phase 3 showed the four
real datasets cannot separate `R` from `Delta`, so the `rsep` cells are the only
ones that can settle the mechanism, and +/-3% MCSE is enough to see an ordering.
Then extend the cells that matter rather than running the full factorial.

```bash
cd /mnt/c/Rworking/negative-sgr
DESIGN=rsep N_ITER=50 WORKERS=6 nohup Rscript analysis/phase5_run.R \
  > analysis/phase5_run.log 2>&1 &
```

Results are written per cell to `analysis/phase5/<cell>.rds` **as each cell
finishes**, and existing cells are skipped, so the sweep resumes rather than
restarts after any interruption. Summarise with `analysis/phase5_report.R`.

**One lever worth considering first.** Arm D costs 154 s against arm A's 16 s on
the same data — ten times, and more than half the total. Dropping arm D would cut
the per-iteration cost from 270 s to 116 s, a 57% saving. The cost is almost
certainly the same thing that made arm D's ErC50 unestimable on real data:
truncating at the zero crossing removes the points that identify the steep part,
leaving a weakly identified fit that the sampler labours over. That is a finding
worth keeping, but it may not need 50 iterations per cell to establish.

## Results so far

**The gate passes.** Arm A's estimate falls outside the convention arms'
intervals in 14 of 24 dataset x arm x endpoint combinations.

**The pre-registered ordering prediction FAILS.** Observed divergence order is
`c_proliferum > r_salina > r_salina2 > c_proliferum2`; Spearman 0.4 with `Delta`
and -0.2 with `R`. Neither the revised (`Delta`) nor the original (`R`) hypothesis
orders the four datasets. `c_proliferum2` has the second-highest `Delta` and the
least divergence. With n = 4 this is close to powerless as a test, but it is not
support, and the simulation is now the only thing that can settle it.

**Arm B2 is broken, not merely biased.** Fixing `bot` at zero gives divergences on
all four datasets, treedepth saturation, and intervals spanning nearly the whole
design (`c_proliferum2` ErC50 342 [0.01, 429]). The plan's premise that fixing
`bot` gives "a clean one-degree-of-freedom question" does not survive contact.

**B1/B3 sample cleanly and move the answer materially**, with non-overlapping
intervals: `c_proliferum` ErC50 3.36 [2.64, 4.41] -> 6.35 [5.77, 7.01];
`r_salina2` ErC10 56.1 [44.1, 68.1] -> 87.1 [81.1, 92.1].

**The NSEC direction prediction also fails** — B1 moves NSEC up on two datasets
and down on two.

**`bot` contraction is the cleanest supporting result.** Flooring pins the
asymptote (contraction 0.966-0.997) where arm A has 0.39-0.42: the convention
manufactures precision about the parameter it fabricated.

**Two results are flagged unusable** in `analysis/phase3_endpoints.csv`
(`at_boundary`, `not_converged`, `usable` columns): arm D's ErC50 on
`c_proliferum` and `r_salina`, where the 50% decline is never reached inside the
truncated range and `ecx()` returns `max(x)` — it warns, but the value still looks
like an estimate.

## Open items, roughly in order

1. **Phase 5 sweep** — your budget call above.
2. **Phase 4 prior sweep** — not run. Compare `bot` default/wider/tighter on one
   dataset to confirm the A-vs-B2 gap is not a prior artefact. Cheap: 3 fits.
3. **Phase 6 vignette** — new subsection in `example1`'s Censoring section on
   `dev`. The section already states the saturation mechanism, so the subsection
   adds the *consequence*: saturation leaves `bot` unidentified, so the fitted
   asymptote is reported from the prior (arm A -1.27 vs arm C -0.42 on
   `c_proliferum2`, same data, same prior). Also replace `example6`'s stale
   disclaimer that `cens()` "is not available in bayesnec today - see issue #181",
   which `dev` has implemented.
4. **Issue #193** overlaps this study — "Rewrite `example7` as a growth-data case
   study using the `alga` dataset". The vignette deliverable may belong there
   rather than only as an `example1` subsection.
5. **`renv.lock` not yet written.** Deferred to avoid `renv::init()` writing a
   `.Rprofile` under running jobs. Safe to do now, nothing is running.
6. **Protocol confirmation still outstanding** — `mu_0`, `t`, `n_0` and both LODs
   are recovered arithmetically from the CSVs and reproduce exactly, but have not
   been checked against the actual test protocols. Needed before publication.
7. **`R/nsec.R` on `dev` has a roxygen typo** — `#'#'` on the references line
   (introduced between `cens-impl` and `dev`). Cosmetic, unrelated to this study,
   worth fixing next time that file is touched.

## Environment notes

- Study pinned to `bayesnec` `dev` @ `374e511c`, loaded from
  `/mnt/c/Rworking/bayesnec-issue173`, which is **detached** at that commit so it
  cannot drift. `load_bayesnec()` asserts the SHA and stops if it has moved.
- The main worktree `/mnt/c/Rworking/bayesnec` is now on `dev` (was
  `hurdle-vignette`).
- `STUDY_CORES = 6` in `R/setup.R`. Raise it if the machine is free — the other
  study that was competing for cores had finished by the end of this session.
- Branch cleanup was offered and **not** done: 8 branches are fully merged into
  `dev` (`Issue-157`, `Issue-162`, `cens-impl`, `hurdle-gamma-impl`,
  `hurdle-vignette`, `issue-173-guidance`, `issue-178-vignette-figs`,
  `zi-beta-impl`) and 4 worktrees hold nothing unmerged. `beta-ub-impl` and
  `issue-191-dispersion` were deliberately left alone.
