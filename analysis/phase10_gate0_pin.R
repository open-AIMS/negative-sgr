## Phase 10, Gate 0 -- the pin is deterministic where it must be and unchanged
## where it must be.
##
## This gate blocks the sweep. Both halves must pass.
##
## HALF 1, determinism. Model-averaged `ecx()`, `nsec()` and `posterior_epred()`
## must return identical values on repeated calls. Under the Phase 1-9 pin they
## did not: the component draw index was resampled with an unseeded `sample()`
## at every call site, and the variation landed almost entirely on the interval.
## Coverage is an interval property, so without this the phase would be
## measuring bayesnec's resampling noise on precisely the arms under test.
##
## HALF 2, invariance. A stored SINGLE-model fit must return bit-identical
## endpoints under the new pin. #216 touches the `bayesmanecfit` path only, and
## that is what allows Phase 10's averaged results to be compared with the
## 500-iteration single-model results already in hand. If this half fails, the
## comparison straddles two packages and the single-model arms must be re-run on
## the three cells before anything is compared.
##
## Half 2 is run in a SEPARATE R process per pin, because pkgload::load_all()
## cannot cleanly swap two versions of the same package inside one session --
## the second load leaves S3 methods from the first registered, and the failure
## is silent.
##
## Run: Rscript analysis/phase10_gate0_pin.R

source("R/setup.R"); source("R/setup_phase10.R")

OUT <- "analysis/phase10_gate0.csv"
gates <- list()

## ------------------------------------------------- half 1: determinism -------
load_bayesnec_p10()
data(manec_example, package = "bayesnec")

nd <- bayesnec::bnec_newdata(manec_example, resolution = 20)
rep_ecx  <- replicate(6, suppressMessages(bayesnec::ecx(manec_example, ecx_val = 10)))
rep_nsec <- replicate(6, suppressMessages(bayesnec::nsec(manec_example)))
rep_pe   <- replicate(4, brms::posterior_epred(manec_example, newdata = nd)[1, 1])

## Identity across calls, checked on every summary quantity rather than on the
## median alone: the median was the stable end even before the fix, and the
## lower bound -- the end a protective concentration is read off -- was the
## unstable one. A gate that checked only the point estimate would have passed
## on the broken pin.
stable <- function(m) all(apply(m, 1, function(z) length(unique(z)) == 1L))
gates[[1]] <- gate_result("P10-G0-1", "model-averaged ecx() identical over 6 calls",
                          stable(rep_ecx),
                          paste0("spread ", paste(signif(apply(rep_ecx, 1, function(z) diff(range(z))), 3),
                                                  collapse = "/")))
gates[[2]] <- gate_result("P10-G0-2", "model-averaged nsec() identical over 6 calls",
                          stable(rep_nsec),
                          paste0("spread ", paste(signif(apply(rep_nsec, 1, function(z) diff(range(z))), 3),
                                                  collapse = "/")))
gates[[3]] <- gate_result("P10-G0-3", "model-averaged posterior_epred() identical over 4 calls",
                          length(unique(rep_pe)) == 1L,
                          paste0("values ", paste(signif(unique(rep_pe), 8), collapse = ", ")))

## A seed set by the caller must NOT be disturbed. The study seeds every
## iteration and forks workers from that state; a helper that called set.seed()
## outright would silently reseed the sweep from inside an estimate.
set.seed(4242); before <- runif(1)
set.seed(4242); invisible(suppressMessages(bayesnec::ecx(manec_example, ecx_val = 10)))
after <- runif(1)
gates[[4]] <- gate_result("P10-G0-4", "ecx() leaves the caller's RNG stream untouched",
                          identical(before, after),
                          paste0("before ", signif(before, 10), " after ", signif(after, 10)))

## ------------------------------------------------- half 2: invariance --------
## Endpoints from a stored single-model fit, computed under each pin in its own
## process and compared as text so the check is exact rather than approximate.
probe <- tempfile(fileext = ".R")
writeLines(c(
  'args <- commandArgs(trailingOnly = TRUE)',
  'suppressMessages(pkgload::load_all(args[1], quiet = TRUE, export_all = FALSE))',
  'f <- readRDS("analysis/phase1_fits.rds")$b2_fit',
  'v <- c(suppressMessages(bayesnec::ecx(f, ecx_val = 10, type = "absolute")),',
  '       suppressMessages(bayesnec::ecx(f, ecx_val = 50, type = "absolute")),',
  '       suppressMessages(bayesnec::nsec(f)))',
  'cat(paste(format(v, digits = 17), collapse = "|"), "\n")'
), probe)

run_probe <- function(wt) {
  o <- system2("Rscript", c(probe, wt), stdout = TRUE, stderr = FALSE)
  trimws(utils::tail(o, 1))
}
old_pin <- run_probe(BAYESNEC_WORKTREE)
new_pin <- run_probe(BAYESNEC_WORKTREE_P10)

gates[[5]] <- gate_result(
  "P10-G0-5", "single-model ecx()/nsec() bit-identical across the two pins",
  identical(old_pin, new_pin) && nzchar(old_pin),
  paste0("old: ", substr(old_pin, 1, 60), " | new: ", substr(new_pin, 1, 60)))

## ------------------------------------------------------------- report --------
g <- do.call(rbind, gates)
utils::write.csv(g, OUT, row.names = FALSE)
print(g, row.names = FALSE)
cat("\n")
if (all(g$passed)) {
  cat("GATE 0 PASSED -- the Phase 10 pin is sound.\n")
} else {
  cat("GATE 0 FAILED on:", paste(g$id[!g$passed], collapse = ", "),
      "\nDo not run the sweep.\n")
  quit(status = 1)
}
