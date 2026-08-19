# Resume here

Updated 2026-08-19 10:35. **Phase 8 stage 2 is RUNNING** -- launched 2026-08-18
22:58, interrupted by a WSL crash at about 10:00 on the 19th and resumed at
10:29 with 46 of 84 blocks already on disk. The crash cost the one block in
flight and nothing else: no stray `.tmp` files, every completed block readable,
41,280 rows and 0 failures across the blocks that survived. Expect completion
late on 20 August.

**Timing is very uneven and the block count misleads.** The six `stops short`
cells ran 8-11 min per block; the reaching and precision cells run 40-120. The
long ones are `fit_arm()`'s escalation path doing its job -- an iteration trips
R-hat > 1.05 and is refitted with `num_warmup=4000`, `adapt_delta=0.999`,
`max_depth=15`, four chains sequentially. A worker sitting at 8% CPU for over an
hour is that, not a hang; check with
`ps -eo etime,pcpu,cmd | grep "[-]-file=analysis/phase8"` and look at the
cmdstan arguments of its child.

```bash
pkill -f "[-]-file=analysis/phase8"   # stop it -- see Trap 1, and note the
                                      # pattern must not appear anywhere else
                                      # on the command line, filenames included
# resume, from /mnt/c/Rworking/negative-sgr:
N_ITER_FROM=241 N_ITER_TO=500 CHUNK=40 WORKERS=18 nohup Rscript \
  analysis/phase8_run.R > analysis/phase8_run.log 2>&1 &
```

Safe to interrupt at any instant. Writes are chunked at 40 iterations to
`analysis/phase5/<cell>__topup_i<from>-<to>.rds`, each written under a temporary
name and renamed, so a file that exists is complete and a resume skips it by
filename. An interrupt costs at most one block. Progress:
`ls analysis/phase5/*__topup*.rds | wc -l` -- **84 blocks** is the full stage
(12 cells x 6 blocks of 40 plus a final 20).

Nothing in stage 2 can change an ordering or a conclusion; it narrows intervals
only. The vignette and every reported number stand without it.

**Phase 7 stage 1 is complete: 12 of 12 cells, 5,760 fits, 0 failures, 0 NA
estimates.** All eight arms now have 240 iterations in every one of the twelve
scenarios. `phase5_metrics.csv`, `phase5_metrics_cleanfits.csv` and
`phase5_r_axis.txt` have been regenerated over all eight, verified first by
`analysis/phase7_verify.R` (shape, not `[done]` lines -- Trap 9).

**Phases 1-5 and 7 are finished.** What remains is Phase 6 (the vignette and the
paper artefacts) and Phase 8 stage 2, the iteration top-up. See "What is left".

## Re-running the simulation (nothing needs re-running as of 2026-08-18)

Both sweeps are complete and both skip finished work, so either command is safe
to issue at any time. Phase 7 writes per 40-iteration block, Phase 5 per cell.

```bash
cd /mnt/c/Rworking/negative-sgr
# arms E and F (phase 7 stage 1)
N_ITER=240 CHUNK=40 WORKERS=18 nohup Rscript analysis/phase7_run.R \
  >> analysis/phase7_run.log 2>&1 &
# arms A-D, B1-B3 (phase 5)
DESIGN=rsep N_ITER=240 WORKERS=18 nohup Rscript analysis/phase5_run.R \
  > analysis/phase5_run_rcells.log 2>&1 &
```

**Verify it is actually running** -- do not trust a launch message:

```bash
pgrep -fc "phase7_ru[n].R"                 # expect ~20
ls analysis/phase5/*__ef*.rds | wc -l      # 72 blocks is the full stage 1
```

After any re-run, `Rscript analysis/phase7_verify.R` checks every cell's shape
and regenerates the report only if all twelve pass.

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
| 6 vignette | **done** — `example7` extended to all eight approaches and re-knitted, `negsgr-cens-vignette` @ `6320d936`, committed locally and **not pushed**, **not** for `dev`. The `example1` censoring edits are commit `4df470ea` |
| 6 paper artefacts | **not started** |
| 7 zero-bounded families | **done 2026-08-18 — 12 of 12 cells, 5,760 fits, 0 failures, 0 NA estimates.** Verified by `analysis/phase7_verify.R`; report regenerated over all eight arms. The F mechanism was tested and is recorded below |
| 8 iteration top-up | **running** since 2026-08-18 22:58 — `analysis/phase8_run.R`, iterations 241-500 for all eight arms, 24,960 fits, ~3.5 days |
| environment | `renv.lock` written (106 packages, R 4.6.1) plus `analysis/session_info.txt` |

125 tests pass (`tests/testthat/`) — re-run 2026-08-18 via `cd tests && Rscript
testthat.R`. Running `testthat::test_dir()` from the project root fails: the
runner sources `../R/*.R` and the tests do not load them themselves.

## What is left

1. **Phase 8, stage 2 is running** (launched 2026-08-18 22:58; ~3.5 days,
   24,960 fits). `analysis/phase8_run.R`, iterations 241-500, all eight arms.
   When it finishes: re-run `analysis/phase5_report.R` (or `phase7_verify.R`,
   which calls it) and regenerate the vignette blocks with
   `analysis/vignette_tables.R`. Expect the numbers to move in the third
   decimal; if any ordering changes, something is wrong, not interesting.
2. **Paper artefacts** (§Phase 6). Lead with the noise axis, not the `delta`
   gradient.
3. **Re-run the case studies under model averaging.** Phase 3 fixes `nec4param`
   for cross-arm comparability, which is right for the simulation (it is the
   generating model) and not defensible on real data. Deferred, not forgotten.
4. **The vignette branch has not been reviewed by anyone but Claude.** It is
   parked deliberately and nothing goes to `dev` without your say.
5. **Why the `nec` displacement under arm F reverses with `delta` is not
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
reaching cells: E 0.00, F 0.03, against B2's 5.67 and D's 1.21. Exclusion rate
0.000 for both, over 2,880 fits each. This strengthens rather than weakens the
study's existing point that flooring failures are silent -- the two approaches
with the worst high-precision bias are also the two best behaved by diagnostics.

ErC50 down the noise axis (`R` = 2.3 -> 3.3 -> 17 -> 73, i.e. falling noise):

| arm | ErC50 bias | ErC50 coverage |
|---|---|---|
| E | -10.3 -> -13.2 -> -15.0 -> **-14.8%** | 0.60 -> 0.23 -> 0.00 -> **0.000** |
| F | -5.6 -> -10.5 -> -14.9 -> **-15.6%** | 0.89 -> 0.80 -> 0.36 -> **0.129** |

At the finest precision these are the two largest ErC50 biases in the study,
ahead of B2's -11.3%. **E is the more important of the two**: it is the family
analogue of B3, behaves like it only worse, and choosing a Beta on a scaled
response is the default thing to do with these data.

**F's pooled bias is a trap -- never quote it.** Its reaching-regime ErC50 bias
averages -2.0%, which reads as respectable, but it is a mean over a quantity that
changes sign: **+13.0% at `delta` 2, -5.6% at 4, -13.4% at 8**. Its RMSE is 1.07
against arm A's 0.40, the largest of any arm, and that is what gives it away.

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
| E | +11.6 -> +13.0 -> +12.7 -> **+12.5%** | -10.3 -> -13.2 -> -15.0 -> **-14.8%** | 0.60 -> 0.23 -> 0.00 -> **0.000** |
| F | +49.7 -> +40.6 -> +22.8 -> **+15.8%** | -5.6 -> -10.5 -> -14.9 -> **-15.6%** | 0.89 -> 0.80 -> 0.36 -> **0.129** |

**A/C/D converge to zero bias; B1/B2/B3 and E/F converge to non-zero
asymptotes**, and E and F end furthest from the truth on ErC50 of any arm. A
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

