## Phase 10, Gate 1 -- the candidate sets are what the plan says they are.
##
## The plan derives them by reading `check_models()`. This runs it. The
## distinction is not pedantry: the LOD error recorded in the plan (§Phase 2,
## [R4]) came from reading a quantity off the data instead of establishing it,
## and cost a rebuild of every real-data result.
##
## The counts asserted here are load-bearing for the budget as well as the
## interpretation: 8 + 8 + 6 + 12 + 12 = 46 model fits per iteration is what
## turns 100 iterations x 3 cells into 2.5 days.
##
## Run: Rscript analysis/phase10_gate1_models.R

source("R/setup.R"); source("R/setup_phase10.R")
source("R/data_prep.R"); source("R/simulate.R")
load_bayesnec_p10(); source("R/arms.R"); source("R/metrics.R")
source("R/model_average.R")
p10_assert_arms_ready()

OUT <- "analysis/phase10_gate1.csv"

cal <- calibrate_sigma(read_sgr("c_proliferum"))
SIGMA_0 <- sigma_0_at(cal$cv_control, R_ref = 2.3, t = 7)
cells <- p10_cells(cal$sigma_ratio)

## Cell 8 -- delta 4, top_factor 2.0, R 2.3 -- the cell the vignette's
## eight-curve figure is drawn from. Iteration 1's dataset, i.e. the one the
## single-model sweep actually saw.
cel <- cells[cells$cell_index == 8L, ]
truth <- sim_truth(R = cel$R, delta = cel$delta, t = 7)
design <- sim_design(truth, top_factor = cel$top_factor, n_conc = 12,
                     n_rep = 5, n_control = 10)
sim <- simulate_dataset(truth, design, cv_control = cal$cv_control,
                        sigma_ratio = cal$sigma_ratio, seed = 7e5 + 1,
                        sigma_mode = "absolute", sigma_0_abs = SIGMA_0)

expected <- list(
  A  = c("nec4param", "neclin", "ecxlin", "ecx4param", "ecxwb1", "ecxwb2",
         "ecxll5", "ecxll4"),
  C  = c("nec4param", "neclin", "ecxlin", "ecx4param", "ecxwb1", "ecxwb2",
         "ecxll5", "ecxll4"),
  B3 = c("nec4param", "ecx4param", "ecxwb1", "ecxwb2", "ecxll5", "ecxll4"),
  E  = c("nec3param", "nec4param", "ecxexp", "ecxsigm", "ecx4param", "ecxwb1",
         "ecxwb2", "ecxwb1p3", "ecxwb2p3", "ecxll5", "ecxll4", "ecxll3"),
  F  = c("nec3param", "nec4param", "ecxexp", "ecxsigm", "ecx4param", "ecxwb1",
         "ecxwb2", "ecxwb1p3", "ecxwb2p3", "ecxll5", "ecxll4", "ecxll3")
)

rows <- list()
for (arm in P10_ARMS) {
  got <- p10_candidate_models(arm, sim)
  exp_a <- expected[[arm]]
  ok <- setequal(got, exp_a)
  rows[[arm]] <- data.frame(
    arm = arm,
    family = p10_family(arm)$family,
    n = length(got),
    n_expected = length(exp_a),
    matches_plan = ok,
    unexpected = paste(setdiff(got, exp_a), collapse = " "),
    missing = paste(setdiff(exp_a, got), collapse = " "),
    models = paste(got, collapse = " "),
    stringsAsFactors = FALSE
  )
}
g <- do.call(rbind, rows)
utils::write.csv(g, OUT, row.names = FALSE)
print(g[, c("arm", "family", "n", "n_expected", "matches_plan", "unexpected",
            "missing")], row.names = FALSE)

cat("\nmodel fits per iteration:", sum(g$n), "\n")
cat("nec4param (the generating model) present in every arm's set:",
    all(vapply(P10_ARMS, function(a)
      "nec4param" %in% p10_candidate_models(a, sim), logical(1))), "\n\n")

## The asymmetry, stated as a number so the report cannot omit it.
zb <- bayesnec:::mod_groups$zero_bounded
cat("zero-bounded shapes available to each arm (the shapes a floored dataset ",
    "wants):\n", sep = "")
for (arm in P10_ARMS) {
  m <- p10_candidate_models(arm, sim)
  cat(sprintf("  %-2s  %d of %d\n", arm, sum(m %in% zb), length(m)))
}

## ------------------------------------- the R-hat screen matches the package --
## `p10_model_rhat()` reads each component's brmsfit directly instead of going
## through `pull_out()`, which rebuilds a bayesnecfit per model and is wasteful
## at this scale. That is only legitimate if it gives the same answer, so it is
## checked against `bayesnec::rhat()` rather than asserted in a comment.
data(manec_example, package = "bayesnec")
## `rhat()` is an S3 method registered against the brms generic, not a bayesnec
## export -- `bayesnec::rhat` does not resolve. Dispatch through brms.
pkg <- suppressMessages(brms::rhat(manec_example,
                                   rhat_cutoff = P10_RHAT_CUTOFF))
pkg_failed <- setNames(vapply(pkg, `[[`, logical(1), "failed"),
                       manec_example$success_models)
mine <- p10_model_rhat(manec_example)
mine_failed <- setNames(mine$max_rhat > P10_RHAT_CUTOFF, mine$model)
same <- identical(pkg_failed[sort(names(pkg_failed))],
                  mine_failed[sort(names(mine_failed))])
cat("\n-- R-hat screen vs bayesnec::rhat() on manec_example --\n")
print(data.frame(model = mine$model, max_rhat = round(mine$max_rhat, 4),
                 mine = unname(mine_failed),
                 bayesnec = unname(pkg_failed[mine$model])), row.names = FALSE)
gate_rhat <- same
cat("agree:", same, "\n")

if (all(g$matches_plan) && gate_rhat) {
  cat("\nGATE 1 PASSED -- candidate sets are as planned and the R-hat screen",
      "matches the package.\n")
} else {
  cat("\nGATE 1 FAILED --",
      if (!all(g$matches_plan))
        paste("candidate sets:", paste(g$arm[!g$matches_plan], collapse = ", ")),
      if (!gate_rhat) "the R-hat screen disagrees with bayesnec::rhat()",
      "\nUpdate the plan and the budget before running the sweep.\n")
  quit(status = 1)
}
