## Phase 10, Gate 2 -- verify the pilot before committing days of compute.
##
## Run `PHASE10_PILOT=1 Rscript analysis/phase10_run.R` first; this reads what it
## wrote. Kept separate from the runner so that the checks can be re-run, and
## extended, without re-fitting.
##
## Four things are checked, and none of them is optional.
##
## 1. THE PAIRING. Every simulated dataset must be the one the single-model arms
##    saw, or the averaged-versus-single contrast is not paired and the phase
##    loses its sharpest comparison. Checked against `f_neg` stored by the
##    Phase 5/7/8 sweep for the same cell and iteration. `f_neg` is discrete
##    (k/70), so one match proves little and forty prove a lot; the check uses
##    every iteration the pilot produced and additionally regenerates the
##    remaining iterations of the cell without fitting, which is free.
##
## 2. THE MACHINERY WORKS ON A `bayesmanecfit` FOR EVERY ARM. Two paths had never
##    been run to completion when this was written: arm C, a `cens()` aterm with
##    a model SET, and arm B3, a `constant(0)` prior inside a NAMED PRIOR LIST
##    with `init = "random"` -- the Phase 9 change whose own pilot crashed before
##    reaching it. A missing endpoint here is a bug, not a hard dataset.
##
## 3. THE BUDGET. Replace the 4.5-worker-minute-per-model-fit estimate with
##    measured per-arm cost and re-derive the wall-clock before committing.
##
## 4. A FIRST LOOK AT THE WEIGHTS, which may answer the question before the sweep
##    runs.
##
## Run: Rscript analysis/phase10_gate2_pilot.R

source("R/setup.R"); source("R/setup_phase10.R")
source("R/data_prep.R"); source("R/simulate.R")
source("R/model_average.R")

IN  <- Sys.getenv("OUT", "analysis/phase10_pilot")
OUT <- "analysis/phase10_gate2.csv"

files <- list.files(IN, pattern = "\\.rds$", full.names = TRUE)
if (!length(files)) {
  stop("No pilot output in ", IN, ". Run PHASE10_PILOT=1 Rscript ",
       "analysis/phase10_run.R first.")
}
blocks <- lapply(files, readRDS)
ends_all <- do.call(rbind, lapply(blocks, `[[`, "endpoints"))
wts_all  <- do.call(rbind, lapply(blocks, `[[`, "weights"))
crv  <- do.call(rbind, lapply(blocks, `[[`, "curves"))
rhat <- do.call(rbind, lapply(blocks, `[[`, "rhat"))
## The converged stage is the primary result; gates are asserted on it.
ends <- ends_all[ends_all$stage == "converged", ]
wts  <- wts_all[wts_all$stage == "converged", ]

gates <- list()

## ------------------------------------------------------- 1. the pairing ------
cal <- calibrate_sigma(read_sgr("c_proliferum"))
SIGMA_0 <- sigma_0_at(cal$cv_control, R_ref = 2.3, t = 7)
all_cells <- p10_cells(cal$sigma_ratio)

## Over ALL three cells, not only the one the pilot fitted. Cell 12's `f_neg` is
## degenerate -- at a 1.9% control CV exactly ten of seventy points are negative
## in every dataset -- so it cannot corroborate anything on its own. Cells 8 and
## 9 carry 9 and 6 distinct values across 500 iterations and are what actually
## test the seeding; the pairing in cell 12 then follows from the code path being
## identical and only `R` differing. The gate below requires the informative
## cells to be present rather than counting matches blindly.
pair_rows <- list()
for (k in P10_CELL_ORDER) {
  cel <- all_cells[all_cells$cell_index == k, ]
  truth <- sim_truth(R = cel$R, delta = cel$delta, t = 7)
  design <- sim_design(truth, top_factor = cel$top_factor, n_conc = 12,
                       n_rep = 5, n_control = 10)
  old_files <- sweep_files(cel$cell)
  if (!length(old_files)) next
  old <- do.call(rbind, lapply(old_files, function(f) {
    d <- readRDS(f); unique(d[, c("iteration", "f_neg")])
  }))
  old <- unique(old)
  ## Regenerate every iteration the stored sweep holds -- no fitting, so this is
  ## free and far stronger than checking the two the pilot fitted.
  mine <- data.frame(iteration = old$iteration, f_neg = vapply(old$iteration,
    function(i) {
      s <- simulate_dataset(truth, design, cv_control = cal$cv_control,
                            sigma_ratio = cal$sigma_ratio, seed = 7e5 + i,
                            sigma_mode = "absolute", sigma_0_abs = SIGMA_0)
      mean(s$sgr < 0)
    }, numeric(1)))
  m <- merge(old, mine, by = "iteration", suffixes = c("_stored", "_regen"))
  pair_rows[[as.character(k)]] <- data.frame(
    cell_index = k, n = nrow(m),
    n_match = sum(m$f_neg_stored == m$f_neg_regen),
    distinct_f_neg = length(unique(m$f_neg_stored)))
}
pair <- do.call(rbind, pair_rows)
print(pair, row.names = FALSE)
informative <- !is.null(pair) & pair$distinct_f_neg >= 5
gates[[1]] <- gate_result(
  "P10-G2-1", "simulated datasets identical to those the single-model arms saw",
  !is.null(pair) && all(pair$n_match == pair$n) && sum(pair$n) > 20 &&
    sum(informative) >= 2,
  paste0(sum(pair$n_match), "/", sum(pair$n), " match on f_neg; ",
         sum(informative), " cells with a discriminating f_neg"))

## ------------------------------- 2. every arm produced every endpoint --------
want <- expand.grid(arm = P10_ARMS, endpoint = c("ErC10", "ErC50", "NSEC"),
                    stringsAsFactors = FALSE)
got <- unique(ends[ends$ok & !is.na(ends$estimate), c("arm", "endpoint")])
have <- merge(want, cbind(got, present = TRUE), all.x = TRUE)
have$present[is.na(have$present)] <- FALSE
gates[[2]] <- gate_result(
  "P10-G2-2", "all five arms returned all three endpoints on a bayesmanecfit",
  all(have$present),
  if (all(have$present)) "" else
    paste("missing:", paste(apply(have[!have$present, 1:2], 1, paste,
                                  collapse = "/"), collapse = " ")))

## Arms C and B3 specifically -- the two untested code paths.
gates[[3]] <- gate_result(
  "P10-G2-3", "arm C: cens() aterm with a model set averaged over >1 model",
  any(wts$arm == "C") && max(wts$n_fitted[wts$arm == "C"]) > 1,
  paste0("models fitted: ", paste(unique(wts$n_fitted[wts$arm == "C"]),
                                  collapse = ",")))
gates[[4]] <- gate_result(
  "P10-G2-4", "arm B3: constant(0) prior in a named list averaged over >1 model",
  any(wts$arm == "B3") && max(wts$n_fitted[wts$arm == "B3"]) > 1,
  paste0("models fitted: ", paste(unique(wts$n_fitted[wts$arm == "B3"]),
                                  collapse = ",")))

## B3 pins `bot` at zero, so its averaged curve must not go below it by more than
## rounding. A curve that dips is a constant(0) prior that did not take -- which
## is silent, and would make B3 a different arm from the one described.
b3 <- crv[crv$arm == "B3", ]
gates[[5]] <- gate_result(
  "P10-G2-5", "arm B3's averaged curve respects the zero lower asymptote",
  nrow(b3) > 0 && min(b3$mu) > -1e-8,
  paste0("min fitted mu = ", signif(min(b3$mu), 4)))

## The R-hat rule is part of the treatment, so the pilot must show it working
## rather than merely not crashing: both stages present, and a drop rate that is
## neither impossible (never fires) nor pathological (empties an arm).
gates[[6]] <- gate_result(
  "P10-G2-6", "both stages recorded (all_models and converged)",
  setequal(unique(ends_all$stage), c("all_models", "converged")),
  paste("stages:", paste(unique(ends_all$stage), collapse = ", ")))
## An arm that loses EVERY model returns no estimate for that dataset. Equal
## rates across arms are harmless however high; UNEQUAL rates reproduce Trap 6,
## because each arm's coverage would then be computed over a differently
## selected subset. The report handles it with common support, but the pilot
## should show whether it is going to bite at all -- and two iterations cannot,
## so this gate is a tripwire, not a measurement.
lost <- !ends_all$ok & !is.na(ends_all$n_kept) & ends_all$n_kept == 0
gates[[7]] <- gate_result(
  "P10-G2-7", "no arm-iteration lost every model to the rhat cutoff",
  !any(lost),
  if (!any(lost)) "none" else
    paste("arms affected:",
          paste(sort(unique(ends_all$arm[lost])), collapse = ", "),
          "-- CHECK the per-arm exclusion spread in the full run"))
if (!is.null(rhat)) {
  cat("\n-- convergence screening at rhat >", P10_RHAT_CUTOFF, "--\n")
  print(aggregate(cbind(max_rhat = max_rhat, dropped = dropped) ~ arm,
                  data = rhat, FUN = function(z) round(mean(z), 3)),
        row.names = FALSE)
  cat("\nmodels dropped per arm-iteration:\n")
  print(tapply(rhat$dropped, list(rhat$arm, rhat$iteration), sum))
}

## ------------------------------------------------- 3. the measured budget ----
cost <- merge(
  aggregate(cbind(models = n_requested) ~ arm, data = ends, FUN = max),
  aggregate(cbind(mins = elapsed) ~ arm,
            data = ends[ends$ok & !duplicated(ends[, c("iteration", "arm")]), ],
            FUN = function(z) round(mean(z) / 60, 2)), by = "arm")
cost <- cost[match(P10_ARMS, cost$arm), ]
cost$mins_per_model <- round(cost$mins / cost$models, 2)
cat("\n-- measured cost per arm --\n"); print(cost, row.names = FALSE)

## `elapsed` MEASURES THE FITTING ONLY, and the fitting is the minority of the
## cost. It is recorded around the `p10_fit()` call and therefore excludes
## everything that follows: endpoint_table() -- three ecx()/nsec() evaluations on
## a model average at resolution 1000, across 8-12 models -- plus p10_weights(),
## fit_diagnostics(), the R-hat drop with its amend() re-stack, and p10_curve().
## All of that runs TWICE, once per stage. On the pilot the summed `elapsed` came
## to 8.9 minutes per iteration while the block itself took 31.5, so the fitting
## was 28% of the work and a budget built on `elapsed` alone understates the run
## by a factor of about 3.5.
##
## So the budget is taken from the BLOCK WALL-CLOCK, which measures everything.
## Chains run sequentially (mc.cores = 1L in the runner) and one whole iteration
## occupies one worker, so wall-clock per iteration IS worker-minutes per
## iteration. Compilation is excluded from both -- it lands in iteration 0, which
## is not written out -- and has to be added back separately.
fit_only <- sum(cost$mins, na.rm = TRUE)
block_lines <- tryCatch(
  grep("^\\[done\\]", readLines("analysis/phase10_pilot.log", warn = FALSE),
       value = TRUE), error = function(e) character(0))
per_iter <- if (length(block_lines)) {
  m <- regmatches(block_lines[1],
                  regexpr("[0-9.]+(?= min)", block_lines[1], perl = TRUE))
  n <- regmatches(block_lines[1],
                  regexpr("[0-9]+(?=/[0-9]+ iterations)", block_lines[1],
                          perl = TRUE))
  w <- as.integer(Sys.getenv("PILOT_WORKERS", "2"))
  as.numeric(m) * min(as.integer(n), w) / as.integer(n)
} else {
  NA_real_
}
cat(sprintf("\nfitting only: %.1f worker-min/iteration\n", fit_only))
cat(sprintf("BLOCK WALL-CLOCK: %.1f worker-min/iteration -- use this one; the\n",
            per_iter))
cat("  difference is endpoint/curve computation, done twice per arm.\n")
cat("\nmodel fits per iteration:", sum(cost$models),
    "| worker-minutes per iteration:", round(per_iter, 1), "\n")
## Compilation is per CELL, because the prior constants are Stan literals and the
## reference dataset changes with the cell (Trap 5). Measured on the pilot rather
## than guessed: cell 12's warm-up ran ~15 min with the programs already cached
## and no C++ compile at all. Budget ~40 min for a cold cell, not the 3 h the
## plan originally assumed.
for (w in c(16, 18)) {
  for (n in c(100, 200, 500)) {
    cat(sprintf("  %3d iterations x 3 cells at %d workers: %4.1f h + ~2 h warm-up\n",
                n, w, n * 3 * per_iter / w / 60))
  }
}
cat("The plan budgeted 207 worker-minutes per iteration from the Phase 8 record.\n",
    "If the measured figure differs materially, re-scale before launching.\n")

## ------------------------------------------------ 4. a first look ------------
if (!is.null(wts)) {
  mw <- aggregate(wi ~ arm + model, data = wts, FUN = mean)
  mw <- mw[order(mw$arm, -mw$wi), ]
  cat("\n-- mean pseudo-BMA weight by arm (pilot, n =",
      length(unique(wts$iteration)), "iterations) --\n")
  for (a in P10_ARMS) {
    z <- mw[mw$arm == a, ]
    if (!nrow(z)) next
    cat(sprintf("  %-2s  %s\n", a,
                paste(sprintf("%s %.2f", z$model, z$wi), collapse = "  ")))
  }
  nec4 <- mw[mw$model == "nec4param", c("arm", "wi")]
  cat("\nweight on nec4param (the generating model):\n")
  print(nec4, row.names = FALSE)
}

cat("\n-- endpoints --\n")
print(ends[ends$ok, c("iteration", "arm", "endpoint", "estimate", "lower",
                      "upper", "truth", "n_fitted")], row.names = FALSE,
      digits = 4)

## ------------------------------------------------------------- report --------
g <- do.call(rbind, gates)
utils::write.csv(g, OUT, row.names = FALSE)
cat("\n"); print(g, row.names = FALSE)
if (all(g$passed)) {
  cat("\nGATE 2 PASSED -- re-derive the budget from the log, then launch.\n")
} else {
  cat("\nGATE 2 FAILED on:", paste(g$id[!g$passed], collapse = ", "),
      "\nDo not launch the sweep.\n")
  quit(status = 1)
}
