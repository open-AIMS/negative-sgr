# Resume here

Updated 2026-08-21 06:30. **Nothing is running. The simulation and the
model-averaged case studies are both finished.**

> **TWO SESSIONS ARE WORKING IN THIS CHECKOUT.** A second agent added
> `R/model_average.R`, `R/setup_phase10.R`, six `analysis/phase10_*.R` scripts
> and `tests/testthat/test-model_average.R` on 2026-08-20, edited
> `tests/testthat.R` to source them, and edited the "What is left" section
> below. None of that has been touched, run or committed from this session, and
> Phase 9's outputs are all named `phase9_*`, so nothing collides on disk. The
> test count going 125 -> 159 is theirs, not a regression. **Phase 9 and their
> Phase 10 overlap and need reconciling by you.** This is the hazard the
> worktree memory note warns about: two agents, one working tree.

**Phase 8 stage 2 completed 2026-08-20 14:00: 12 of 12 cells, 24,960 fits, 0
failures, 0 NA estimates.** Every arm now has 500 iterations in every one of the
twelve scenarios. Verified by shape before anything was reported: each cell holds
260 x 8 x 3 = 6,240 top-up rows and the 240-iteration files are untouched. The
run survived a WSL crash on the 19th and resumed with the loss of one block.

**No conclusion moved, which is what stage 2 was for.** Mean change in relative
bias 0.005 and in coverage 0.011 against the 240-iteration numbers; bias MCSE
fell from about 0.0046 to 0.0030. The ErC10 ordering is identical. The ErC50
ordering differs in one place -- D and F trade positions -- on a gap of 0.0008
against a combined MCSE of 0.0073, i.e. 0.12 MCSE, so they were never
distinguishable and this is not a reordering in any meaningful sense. **Do not
report D and F as ranked against each other on ErC50.**

`phase5_metrics.csv`, `phase5_metrics_cleanfits.csv`, `phase5_r_axis.txt`, the
vignette's three embedded result blocks and every figure quoted in its prose have
all been regenerated at 500 and committed.

## Re-running the simulation (nothing needs re-running as of 2026-08-20)

All three sweeps are complete and every one skips finished work, so any of these
is safe to issue at any time.

```bash
cd /mnt/c/Rworking/negative-sgr
# iterations 1-240, arms A-D and B1-B3
DESIGN=rsep N_ITER=240 WORKERS=18 nohup Rscript analysis/phase5_run.R \
  > analysis/phase5_run_rcells.log 2>&1 &
# iterations 1-240, arms E and F
N_ITER=240 CHUNK=40 WORKERS=18 nohup Rscript analysis/phase7_run.R \
  >> analysis/phase7_run.log 2>&1 &
# iterations 241-500, all eight arms
N_ITER_FROM=241 N_ITER_TO=500 CHUNK=40 WORKERS=18 nohup Rscript \
  analysis/phase8_run.R >> analysis/phase8_run.log 2>&1 &
```

Then `Rscript analysis/phase7_verify.R`, which checks every cell's shape and
regenerates the report only if all twelve pass, and
`Rscript analysis/vignette_tables.R > analysis/vignette_tables.txt` to rebuild
the vignette's blocks.

**Stopping a run:** `pkill -f "[-]-file=analysis/phase8"`. Match on the
`--file=` argument, never on the script name -- see Trap 1.

**Reading progress:** blocks are files. `ls analysis/phase5/*__topup*.rds | wc -l`
(84 is full), `ls analysis/phase5/*__ef*.rds | wc -l` (72 is full). Timing is
very uneven: the `stops short` cells run 8-11 min per block and the reaching and
precision cells 40-120. The long ones are `fit_arm()`'s escalation path -- an
iteration trips R-hat > 1.05 and is refitted with `num_warmup=4000`,
`adapt_delta=0.999`, `max_depth=15`, four chains sequentially, which took 119
minutes once. **A worker at 8% CPU for over an hour is that, not a hang**; check
the cmdstan arguments of its child process before concluding anything is stuck.

## Running Phase 10 (the simulation under model averaging)

**Built 2026-08-20, not yet launched.** Gates 0 and 1 pass. The launch is
sequenced deliberately; do not skip a step, and do not start before Phase 9 has
finished and `R/arms.R` / `R/metrics.R` are committed -- Phase 10 depends on
those working-tree changes and `p10_assert_arms_ready()` will refuse to run
without them.

```bash
cd /mnt/c/Rworking/negative-sgr
Rscript analysis/phase10_gate0_pin.R        # done, PASSED -- the #216 pin
Rscript analysis/phase10_gate1_models.R     # done, PASSED -- candidate sets
Rscript analysis/phase10_smoke.R            # 2 models per arm, cheap plumbing check
PHASE10_PILOT=1 Rscript analysis/phase10_run.R > analysis/phase10_pilot.log 2>&1
Rscript analysis/phase10_gate2_pilot.R      # verifies the pilot; do not skip
# then, and only then:
CELLS=12,8,9 N_ITER_FROM=1 N_ITER_TO=100 CHUNK=20 WORKERS=18 nohup Rscript \
  analysis/phase10_run.R > analysis/phase10_run.log 2>&1 &
Rscript analysis/phase10_report.R > analysis/phase10_report.log
```

**Stopping:** `pkill -f "[-]-file=analysis/phase10"`. Match on `--file=`, never
on the script name -- see Trap 1, which has cost this project five occurrences.

**Progress:** blocks are files. `ls analysis/phase10/*.rds | wc -l`; 15 is full
at 3 cells x 100 iterations x CHUNK 20.

**The pin is different from every other phase.** Phase 10 loads bayesnec from
`/mnt/c/Rworking/bayesnec-negsgr-p10` (branch `study-pin-216`), which is
`374e511c` plus the two #216 commits and nothing else. `R/setup_phase10.R` holds
the constants and `load_bayesnec_p10()` asserts the SHA. Phases 1-9 are
untouched and still load `/mnt/c/Rworking/bayesnec-issue173`. Gate 0 asserts both
halves of what makes this safe: model-averaged output is now reproducible, and a
stored single-model fit returns bit-identical endpoints, so the 48,000 existing
fits stand.

**Phase 9 has the same #216 problem and is NOT protected.** It reports
model-averaged `ecx()` and `nsec()` on the real datasets from the old pin, so
those numbers are not reproducible between calls. Worth re-running on the
Phase 10 pin, or at minimum flagging in whatever it reports.

**All five arms use Stan's random inits** (`P10_INIT = "random"`). Not a
performance tweak: under a model set `fit_arm()` forces random inits for B3,
because `bnec()` passes one `init` to every model and they do not share a
parameter vector — so if the other arms used bayesnec's own search, B3 would
differ by its initialisation as well as by its constraint. The Phase 9 work also
measured that search at **612.8 s against 6.1 s** on floored data, and four of
the five arms fit floored or non-negative responses. Verified equivalent on the
same data and seed (ErC10 identical, ErC50 6.347 vs 6.367, zero divergences
either way). Cost: Phase 5/7/8 used the search, so the paired contrast now
differs in inits too — sampler noise rather than bias, but recorded.

**Non-convergence is handled inside the arm, not by filtering results.** Models
with R-hat > 1.01 are dropped from the candidate set and the weights re-stacked
over the survivors (`rhat()` then `amend(drop=)`), because that is the workflow.
This is NOT Trap 6, which forbids excluding fitted datasets from the summaries
and still holds. Endpoints are recorded at two stages -- `all_models` and
`converged` -- and `converged` is primary. `fit_arm()`'s escalation is disabled
(`P10_RHAT_ESCALATE = Inf`): under a model set it re-samples all twelve models
when any one exceeds the threshold, which would inflate the budget
unpredictably and is not what an analyst does about one stuck component.

**Budget.** 46 model fits per iteration (A 8, C 8, B3 6, E 12, F 12 -- confirmed
by Gate 1, not assumed). 100 iterations x 3 cells is 13,800 fits, about 2.5 days
at 18 workers, **plus roughly 3 hours per cell of compilation** -- every model
recompiles when the cell changes because the prior constants are Stan literals
(Trap 5), and the warm-up is serial on purpose. Re-derive from the pilot's
measured per-arm timings before committing.

**The pairing is verified, and it is the phase's sharpest instrument.** The
simulated datasets are bit-identical to those the single-model arms saw: 500/500
iterations match on `f_neg` in all three cells. So `phase10_report.R` reports the
per-iteration paired difference `averaged - single`, not two independent biases.
Note `phase8_run.R` renumbers `k` when it subsets cells and would break this;
`p10_cells()` carries the original index and a test asserts it does.

## Starting a fresh Claude session

Point it at this file and say something like *"read RESUME.md in
/mnt/c/Rworking/negative-sgr and let's finish this"*. Everything needed is
below. Two things it will not know unless told:

- The prompt log required by `/mnt/c/Rworking/CLAUDE.md` §10 lives at
  `../bayesnec/prompts/negative-sgr-study.md` — append, don't start a new file.
- **Read "Traps" below before running anything.** Each one cost this project
  real time.

## State

| phase | status |
|---|---|
| Plan | **revision 6**, `../bayesnec/ignore/negative-sgr-study-plan.md` (revs 3-6 backed up in `../bayesnec/superceded/`). §Phase 4, §Phase 5 and §Phase 7 carry "Observed result" subsections; §Phase 8 stage 2 is the only outstanding compute |
| 1 branch verification | **done** — six gates pass, `analysis/phase1_gates.csv` |
| 2 diagnostics | **done** — `analysis/dataset_summary.csv` regenerated post-LOD |
| 3 arms on real data | **done** — rebuilt post-LOD, no errors, boundary flags corrected |
| 4 `bot` prior sensitivity | **done** — contraction table plus the prior sweep, `analysis/phase4_prior_sweep.{R,csv,log}` |
| 5 simulation | **done — 12 of 12 cells, 17,280 fits, 0 worker failures.** Only non-zero exclusion anywhere: arm C, 5 of 240 iterations in `d4.0_t0.8_R2.3` |
| 5 report | **done over all 12 cells** — `phase5_metrics.csv`, `phase5_metrics_cleanfits.csv`, `phase5_r_axis.txt`, `phase5_report.log` |
| 6 vignette | **done at 500 iterations** — `negsgr-cens-vignette` @ `8dc898ba`, committed locally and **not pushed**, **not** for `dev`. The `example1` censoring edits are commit `4df470ea` |
| 6 paper artefacts | **not started** |
| 7 zero-bounded families | **done 2026-08-18 — 12 of 12 cells, 5,760 fits, 0 failures, 0 NA estimates.** Verified by `analysis/phase7_verify.R`; report regenerated over all eight arms. The F mechanism was tested and is recorded below |
| 9 case studies, model-averaged | **done 2026-08-21** — 32 arm-fits, 0 failures, `analysis/phase9_{modelavg,report}.R`. Answers the question §Phase 3 deferred; see Findings |
| 8 iteration top-up | **done 2026-08-20** — `analysis/phase8_run.R`, iterations 241-500 for all eight arms, 24,960 fits, 0 failures. Every arm now has 500 iterations everywhere |
| 9 case studies, averaged | **run 2026-08-20 20:01** — `phase9_endpoints.csv`, `phase9_weights.csv`, `phase9_diagnostics.csv` (plus a `_randominit` variant). `R/arms.R` and `R/metrics.R` still **uncommitted**; commit before anything else touches them |
| 10 simulation, averaged | **built 2026-08-20, not launched** — runner, gates and report written; Gates 0 and 1 and the smoke test PASS. Runs on its own pin (`bayesnec-negsgr-p10`, `374e511c` + #216). Waiting on Phase 9 finishing and `R/arms.R` / `R/metrics.R` being committed |
| environment | `renv.lock` written (106 packages, R 4.6.1) plus `analysis/session_info.txt` |

125 tests pass (`tests/testthat/`) — re-run 2026-08-18 via `cd tests && Rscript
testthat.R`. Running `testthat::test_dir()` from the project root fails: the
runner sources `../R/*.R` and the tests do not load them themselves.

## What is left

1. **Paper artefacts** (§Phase 6) are now the only substantive work left. Lead
   with the noise axis, not the `delta` gradient. Every number they need is in
   `phase5_metrics.csv` and `phase5_r_axis.txt` at 500 iterations.
2. ~~Re-run the case studies under model averaging (Phase 9).~~ **DONE
   2026-08-21**, commit `ca1b205`. The `R/arms.R` and `R/metrics.R` changes it
   needed are committed. Results in the Findings section below.
2b. **Run the simulation under model averaging** (Phase 10, planned
   2026-08-20). Does averaging rescue the floored arms' bias and coverage, or is
   the distortion a property of the convention rather than of fixing one
   equation? Five arms (A, C, B3, E, F) over a `decline` set, three scenarios
   (cells 12, 8, 9 — in that order), paired to the existing single-model fits by
   seed reuse. ~13,800 fits, ~2.5 days plus ~9 h of compilation.
   **Blocked on bayesnec #216** — model-averaged `ecx()`/`nsec()` resample with
   an unseeded `sample()` and the noise lands on the interval, which is what
   coverage measures. Recommended route: cherry-pick `c7718866` and `d4b95783`
   from `issue-216-deterministic-model-averaging` onto the pin `374e511c`
   (tested: one conflict hunk each in `R/helpers.R` and `R/expand_classes.R`),
   which leaves every existing single-model number valid. Full plan in §Phase 10
   of `../bayesnec/ignore/negative-sgr-study-plan.md`.
3. **The vignette branch has not been reviewed by anyone but Claude.** It is
   parked deliberately and nothing goes to `dev` without your say.
4. **Why the `nec` displacement under arm F reverses with `delta` is not
   explained.** The vignette says so explicitly rather than guessing. If it is
   worth chasing, `analysis/phase7_f_mechanism.R` is the harness to extend.

**The reporting changes stage 1 forced are all discharged:** the palette is eight
colours and re-validated (`analysis/arm_palette.R`), the two-panel split is three
("measurement retained", "zero boundary imposed", "zero bounded by the family"),
and the exclusion and divergence tables cover E and F.

## Open decisions

- **`renv.lock` is a record, not an renv project.** Written with
  `renv::lockfile_create()` / `lockfile_write()`, so it captures the dependency
  state without creating `renv/`, a project library, or an `.Rprofile`. Running
  `renv::init()` properly would change how every future R session in this
  project behaves and was deliberately left to you.
- Note `/mnt/c/Rworking/CLAUDE.md` records WSL R as 4.5.2; it is actually 4.6.1,
  which is what the lockfile pins.

## Findings

### Phase 7 — the family-floored arms, complete

**Both arms are worse than every Gaussian convention at high precision, and both
sample cleanly while being so.** Mean divergent transitions per fit in the
reaching cells: E 0.00, F 0.01, against B2's 5.81 and D's 1.04. Exclusion rate
0.000 for both. Over the whole 500-iteration sweep E produced **not one divergent
transition in any of its 6,000 fits**. This strengthens rather than weakens the
study's existing point that flooring failures are silent -- the two approaches
with the worst high-precision bias are also the two best behaved by diagnostics.

ErC50 down the noise axis (`R` = 2.3 -> 3.3 -> 17 -> 73, i.e. falling noise):

| arm | ErC50 bias | ErC50 coverage |
|---|---|---|
| E | -10.6 -> -13.2 -> -15.0 -> **-14.8%** | 0.58 -> 0.24 -> 0.00 -> **0.000** |
| F | -5.6 -> -10.6 -> -14.7 -> **-15.6%** | 0.89 -> 0.80 -> 0.37 -> **0.116** |

At the finest precision these are the two largest ErC50 biases in the study,
ahead of B2's -11.3%. **E is the more important of the two**: it is the family
analogue of B3, behaves like it only worse, and choosing a Beta on a scaled
response is the default thing to do with these data.

**F's pooled bias is a trap -- never quote it.** Its reaching-regime ErC50 bias
averages -1.6%, which reads as respectable, but it is a mean over a quantity that
changes sign: **+14.0% at `delta` 2, -5.6% at 4, -13.4% at 8**. Its RMSE is 1.11
against arm A's 0.38, the largest of any arm, and that is what gives it away.

**The variance-structure explanation is now positively contradicted, not merely
unsupported.** `analysis/phase7_f_mechanism.R` refits F against **B3** -- same
floored data, same zero asymptote, so the contrast isolates the likelihood and
nothing else -- on the sweep's own first 20 iterations per reaching cell. It
reproduces the sweep bit-for-bit (paired difference exactly 0.00), so these are
the sweep's fits with the parameter table retained. The displacement is carried
by **`nec`**, not `beta`:

| `delta` | `nec` F-B3 | sign | `beta` F-B3 | ErC50 F-B3 |
|---|---|---|---|---|
| 2 | +1.50 (SE 0.33) | 19/20 positive | +0.16 (SE 0.08) | +0.80 (SE 0.13) |
| 4 | +0.78 (SE 0.43) | 15/20 positive | +0.27 (SE 0.10) | +0.23 (SE 0.30) |
| 8 | -0.16 (SE 0.10) | 3/20 positive | +0.01 (SE 0.06) | -0.18 (SE 0.05) |

`nec` and ErC50 displacements correlate at 0.56-0.96. The old story -- a Gamma's
mean-squared variance overweighting the low tail and flattening the descent --
would live in `beta`, and `beta` is the parameter that does not move. **Do not
revive it.** The Gamma relocates the breakpoint; it does not bend the curve.

The one mechanism statement the data support is about dispersion. A Gamma with an
identity link has a constant coefficient of variation, and fitting one to floored
data forces that CV to **47-64%**, against 10-13% for the Gaussian on the same
floored data and a true control CV of 9.6%. **Why the displacement reverses with
`delta` is still not established**, and the vignette says so rather than guessing.

E and F were **not** run on the four real datasets, deliberately: there is no
true value to score them against there, and what they contribute -- that the bias
survives to the noise-free limit -- exists only where the truth is known.

### Phase 5 — the core result

`phase5_report.R` now performs the regime split itself (§HEADLINE: ARMS BY
REGIME), so these numbers come straight out of its log rather than being
recomputed by hand.

**ErC50** where the design reaches the negative region — mean over the three
`R = 2.3, top_factor = 2.0` cells (`f_neg` 0.152), all eight arms at 500
iterations:

| arm | what it does | bias | coverage | RMSE | note |
|---|---|---|---|---|---|
| C | left-censor negatives | **-0.3%** | 0.845 | 0.46 | matches A — censoring recovers what flooring loses |
| A | raw, `bot` free | **+0.6%** | 0.929 | 0.38 | reference; the only arm near nominal coverage |
| F | floor + Gamma, `nec3param` | -1.6% | 0.770 | **1.11** | **do not read as accurate** — a mean over +14.0/-5.6/-13.4 across `delta`, and the largest RMSE of any arm |
| D | truncate at crossing | -1.7% | 0.850 | 0.48 | |
| B1 | floor to 0, `bot` free | -4.9% | 0.705 | 0.52 | |
| B3 | floor + `bot` = 0 | -6.5% | 0.617 | 0.53 | worst coverage of the Gaussian arms |
| E | scale + Beta, `nec3param` | -8.7% | 0.641 | 0.59 | the default thing to do with these data |
| B2 | raw + `bot` = 0 | -9.5% | 0.970 | 0.63 | worst bias; ~5.8 divergences/fit. Its coverage is bought with width |

**D and F are not ranked against each other here.** Their gap is 0.0008 against a
combined MCSE of 0.0073 — 0.12 MCSE. They swapped places between the 240- and
500-iteration runs, which is what two indistinguishable numbers do.

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
| A | +19.1 -> +12.1 -> +3.5 -> **+1.5%** | -0.4 -> -0.9 -> -0.6 -> **-0.4%** | 0.95 -> 0.94 -> 0.90 -> **0.89** |
| C | +10.4 -> +6.7 -> +2.3 -> **+1.0%** | -0.7 -> -1.4 -> -1.3 -> **-1.0%** | 0.85 -> 0.82 -> 0.79 -> **0.77** |
| D | +12.6 -> +7.6 -> +2.4 -> **+1.0%** | -2.1 -> -2.0 -> -1.3 -> **-1.0%** | 0.84 -> 0.82 -> 0.79 -> **0.79** |
| B1 | +16.8 -> +15.0 -> +13.1 -> **+13.1%** | -5.6 -> -7.0 -> -8.1 -> **-8.2%** | 0.71 -> 0.51 -> 0.02 -> **0.002** |
| B3 | +21.4 -> +19.9 -> +15.7 -> **+14.9%** | -7.6 -> -8.9 -> -10.1 -> **-10.2%** | 0.60 -> 0.34 -> 0.004 -> **0.000** |
| B2 | +54.0 -> +51.5 -> +49.8 -> **+49.5%** | -10.2 -> -10.8 -> -11.3 -> **-11.4%** | 0.97 -> 0.99 -> 1.00 -> **1.000** |
| E | +10.7 -> +12.5 -> +12.6 -> **+12.4%** | -10.6 -> -13.2 -> -15.0 -> **-14.8%** | 0.58 -> 0.24 -> 0.00 -> **0.000** |
| F | +49.3 -> +39.5 -> +23.8 -> **+15.7%** | -5.6 -> -10.6 -> -14.7 -> **-15.6%** | 0.89 -> 0.80 -> 0.37 -> **0.116** |

**A/C/D converge to zero bias; B1/B2/B3 and E/F converge to non-zero
asymptotes**, and E and F end furthest from the truth on ErC50 of any arm. A
quantity that does not vanish as noise vanishes is misspecification, not
estimation error. This is a stronger argument than the `delta` gradient because
it needs no comparison across cells — each arm is its own control.

**Coverage under flooring collapses to 0 of 500.** The interval shrinks with
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
**Do not quote the `R` = 2.3 ErC10 row alone.** F is the one arm whose ErC10
bias falls steeply with precision (+49.7% -> +15.8%) without converging: it
flattens an order of magnitude away from A's +1.4%, so a reader who sees only
that it improves will draw the wrong conclusion.

### What to recommend

**Arms A and C are close and both sound; the claim is that B1/B2/B3, E and F
are not** —
not that one of A or C wins. At adequate precision they are indistinguishable on
bias (~1% on both endpoints). On coverage they trade: A is better on ErC50
(0.87-0.94 vs C's 0.78-0.83) and C is better on ErC10 (0.90-0.95 vs A's
0.83-0.98). **Neither reaches nominal coverage on ErC50 — say so plainly.** Arm
C's advantage is practical, not statistical: it still works when the negative
values were never recorded, which arm A does not. Its cost is the identification
failure in Phase 3 and the vignette — where the asymptote is entirely censored,
`bot` is the prior's and must be reported as such.

**The strongest single recommendation to come out of Phase 7** is the one about
families, because it is the choice most analysts actually make: moving to a Beta
or a Gamma to accommodate a response that will not go negative is the worst of
the eight options tested, and it is the only one whose flooring is invisible in
the code.

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

1. **`pkill`/`pgrep` patterns match your own shell -- and the bracket class is
    not enough.** The usual fix is a bracket class, `pkill -f "phase5_ru[n].R"`,
    so the pattern cannot match the command that contains it. That fails the
    moment the SAME command line mentions the target file anywhere else: a
    `pkill -f "phase8_ru[n].R"` issued alongside a script that contained the
    literal string `analysis/phase8_run.R` matched its own shell and killed it
    (exit 144), twice in one session. Either issue the `pkill` on a line that
    mentions nothing else, or match something only the process has --
    `pkill -f "[-]-file=analysis/phase8"` matches the R process and not the
    shell that launched it. The same trap breaks `until ! pgrep -f ...` waiting
    loops, which then never exit because they always find themselves. The bare
    form, `pkill -f phase5_run.R`, has killed the invoking bash three times in
    this project including once killing a monitor; the bracket-class form has
    now done it twice more by the route above. Five occurrences, one trap.
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
11. **A cell file's row count is not `iterations x arms x 3`.** A *failed* fit
    writes one row carrying the error message, not three endpoint rows, so
    `d4.0_t0.8_R2.3_s8.1.rds` holds 4310 rows rather than 4320 -- the five arm-C
    failures contribute one row each. A verifier that hardcodes 4320 flags that
    known-good cell as incomplete and refuses to regenerate the report.
    `analysis/phase7_verify.R` derives the expectation from the recorded
    failures instead.
12. **Run the tests as `cd tests && Rscript testthat.R`.** `testthat.R` sources
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

## Where the vignette numbers come from

`analysis/vignette_tables.R` generates every result block the vignette embeds --
`sim_results`, `precision`, `nsec_precision`, plus the figures quoted in prose --
in exactly the form `example7.Rmd.orig` expects. It reproduces every pre-Phase-7
number exactly, which is the check that it is generating the same quantities the
hand-written blocks held. Two definitions in it are not obvious and are recorded
in the script: the NSEC `nsec` column is a MEAN across iterations while its
`ratio_to_A` column is the MEDIAN of the per-iteration PAIRED ratio, and the
headline tables average over `delta` within a regime at `R` = 2.3 only, because
`R` is a noise axis and pooling it would weight the arm comparison by noise
level.

Regenerate with `Rscript analysis/vignette_tables.R > analysis/vignette_tables.txt`
and paste the blocks; do not edit the numbers in the vignette by hand.

