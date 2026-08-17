## Unattended close-out, run after the Phase 5 sweep exits.
##
## Everything here is mechanical and safe to run without supervision. The
## interpretive steps -- folding the R-axis result into the plan's conclusions
## and RESUME.md -- are deliberately NOT here, because they are judgements about
## what the numbers mean rather than transformations of them.
##
## Invoked by analysis/overnight_finish.sh, which waits for the sweep process to
## exit first. Run: Rscript analysis/overnight_finish.R

options(width = 200)
cat("=== overnight close-out started", format(Sys.time()), "===\n\n")

## ----------------------------------------------------- verify the artefact ---
## Trap 9: check the artefact, not the exit code. A sweep has previously
## reported [done] on a cell where 239 of 240 iterations had failed, so the cell
## is inspected here rather than trusted.
CELLS <- c("d4.0_t2.0_R17.0_s8.1", "d4.0_t2.0_R73.0_s8.1")
ok_all <- TRUE
for (cl in CELLS) {
  p <- file.path("analysis/phase5", paste0(cl, ".rds"))
  if (!file.exists(p)) {
    cat("MISSING:", p, "\n"); ok_all <- FALSE; next
  }
  r <- readRDS(p)
  n_it <- length(unique(r$iteration))
  n_arm <- length(unique(r$arm))
  cat(sprintf("%s  written %s\n", cl,
              format(file.info(p)$mtime, "%Y-%m-%d %H:%M:%S")))
  cat(sprintf("  rows %d | iterations %d | arms %d | ok %d | NA estimates %d | f_neg %.4f\n",
              nrow(r), n_it, n_arm, sum(r$ok), sum(is.na(r$estimate)),
              mean(r$f_neg)))
  ## 240 iterations x 6 arms x 3 endpoints. Anything short means the cell was
  ## written from a partial run and must not be reported.
  if (n_it != 240 || n_arm != 6 || nrow(r) != 4320) {
    cat("  *** INCOMPLETE -- do not report this cell ***\n"); ok_all <- FALSE
  }
}
cat("\nsweep log [done]/[abort] lines:\n")
lg <- readLines("analysis/phase5_run_rcells.log", warn = FALSE)
cat(paste(grep("^\\[done\\]|^\\[abort\\]", lg, value = TRUE), collapse = "\n"),
    "\n")

if (!ok_all) {
  cat("\nOne or more cells are incomplete. Report NOT regenerated; the existing\n",
      "phase5_metrics.csv is left untouched so nothing stale is written over a\n",
      "good file. Re-launch the sweep -- completed cells are skipped.\n")
} else {
  ## ---------------------------------------------------------- the report -----
  cat("\n=== re-running phase5_report.R over all cells ===\n")
  ## `system()` with shell redirection rather than system2(): system2's stderr=
  ## takes a filename or TRUE, and pointing both streams at one file through it
  ## is fiddly, whereas 2>&1 merges them in the order they were written.
  st <- system("Rscript analysis/phase5_report.R > analysis/phase5_report.log 2>&1")
  cat("phase5_report.R exit status:", st, "\n")
  cat("wrote:", paste(c("analysis/phase5_metrics.csv",
                        "analysis/phase5_metrics_cleanfits.csv"),
                      collapse = ", "), "\n")

  ## ------------------------------------------------- the R-axis read-out -----
  ## Pulled out separately because it is the one open question left for the
  ## morning: whether the noise axis earns a place in the paper. Written to its
  ## own file so the next session does not have to recompute it.
  m <- read.csv("analysis/phase5_metrics.csv")
  s <- m[m$delta == 4 & m$top_factor == 2 &
           m$endpoint %in% c("ErC10", "ErC50"), ]
  w <- stats::reshape(s[, c("R", "arm", "endpoint", "rel_bias")],
                      idvar = c("arm", "endpoint"), timevar = "R",
                      direction = "wide")
  w <- w[order(w$endpoint, w$arm), ]
  cov <- stats::reshape(s[, c("R", "arm", "endpoint", "coverage")],
                        idvar = c("arm", "endpoint"), timevar = "R",
                        direction = "wide")
  cov <- cov[order(cov$endpoint, cov$arm), ]

  sink("analysis/phase5_r_axis.txt")
  cat("R axis at delta = 4, top_factor = 2. R rises => control CV falls, so\n",
      "these columns are a NOISE sweep, not a test of control fold-change:\n",
      "the generating model is exactly scale-equivariant in the growth rate.\n\n")
  cat("--- relative bias ---\n")
  print(w, digits = 3, row.names = FALSE)
  cat("\n--- coverage ---\n")
  print(cov, digits = 3, row.names = FALSE)
  cat("\n--- f_neg and sigma_ratio (should be near-constant across R) ---\n")
  print(unique(s[, c("R", "f_neg", "sigma_ratio")]), digits = 4,
        row.names = FALSE)
  cat("\nWhat to look for: arms A/C/D collapsing toward zero bias as noise\n",
      "falls while B1/B2/B3 approach a NON-ZERO asymptote. A quantity that\n",
      "does not vanish as noise vanishes is misspecification, not estimation\n",
      "error -- which is a stronger argument than the delta gradient.\n")
  sink()
  cat("\nwrote: analysis/phase5_r_axis.txt\n")
  cat(readLines("analysis/phase5_r_axis.txt"), sep = "\n")
}

## ------------------------------------------------------------- renv.lock -----
## `lockfile_create()` + `lockfile_write()` rather than `renv::init()`.
##
## init() would create a project library, rewrite .Rprofile and restart R --
## changing how every future R session in this project behaves. That is a
## decision for the user, not something to do to their environment unattended
## overnight. This route records the dependency state and touches nothing else.
## Converting the project to renv proper remains an open choice.
cat("\n=== writing renv.lock (record only; not renv::init()) ===\n")
lk <- try({
  lf <- renv::lockfile_create(project = ".", type = "implicit")
  renv::lockfile_write(lf, file = "renv.lock")
}, silent = TRUE)
if (inherits(lk, "try-error")) {
  cat("renv.lock FAILED:\n", as.character(lk), "\n")
} else {
  cat("wrote renv.lock with",
      length(jsonlite::fromJSON("renv.lock")$Packages), "packages\n")
}

cat("\n=== recording sessionInfo as a belt-and-braces companion ===\n")
writeLines(capture.output(utils::sessionInfo()), "analysis/session_info.txt")
cat("wrote analysis/session_info.txt\n")

cat("\n=== overnight close-out finished", format(Sys.time()), "===\n")
