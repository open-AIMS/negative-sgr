## Phase 2 -- data preparation and diagnostics. No models are fitted here.

source("R/setup.R")
source("R/data_prep.R")
source("R/diagnostics.R")
source("R/figures.R")

dats <- lapply(dataset_names(), read_sgr)
names(dats) <- dataset_names()

## ---------------------------------------------------------------- summary ---
summary_tab <- do.call(rbind, lapply(dats, dataset_summary))
utils::write.csv(summary_tab, "analysis/dataset_summary.csv", row.names = FALSE)

cat("\n===== DATASET SUMMARY =====\n")
options(width = 200)
print(t(summary_tab))

## The two orderings the study has to distinguish.
cat("\n-- ordering by Delta (discarded effect fraction) --\n")
print(summary_tab[order(-summary_tab$delta), c("dataset", "delta", "R")])
cat("\n-- ordering by R --\n")
print(summary_tab[order(summary_tab$R), c("dataset", "R", "delta")])

## ------------------------------------------------------------ LOD check -----
## The LOD asserted in dataset_meta() is confirmed against the data rather than
## trusted: the bound it implies must equal the SGR of the rows that were
## detected exactly at the limit.
cat("\n===== LOD CONFIRMATION =====\n")
lod_check <- do.call(rbind, lapply(dats, function(d) {
  ds <- unique(d$dataset)
  lod <- dataset_meta()$lod[dataset_meta()$dataset == ds]
  if (is.na(lod)) return(NULL)
  design <- recover_design(d)
  bound <- lod_bound(lod, design$n_0, design$t)
  at_lod <- d$sgr[d$density == lod]
  data.frame(dataset = ds, lod = lod, bound = bound,
             n_at_lod = length(at_lod),
             observed_at_lod = if (length(at_lod)) unique(at_lod)[1] else NA,
             max_abs_discrepancy = if (length(at_lod))
               max(abs(at_lod - bound)) else NA,
             n_below_lod = sum(d$density == 0))
}))
print(lod_check)

## ------------------------------------------------------ preparation check ---
## What each preparation does to each dataset, so the arm definitions can be
## read off rather than inferred.
cat("\n===== PREPARATIONS =====\n")
prep_tab <- do.call(rbind, lapply(dats, function(d) {
  ds <- unique(d$dataset)
  do.call(rbind, lapply(c("raw", "bound", "supplied", "floored", "censored"), function(p) {
    a <- prepare_sgr(d, p)
    data.frame(dataset = ds, preparation = p, n = nrow(a),
               n_negative = sum(a$y < 0), n_zero = sum(a$y == 0),
               n_censored = sum(a$cens != "none"),
               min_y = min(a$y),
               identical_to_supplied = isTRUE(all.equal(a$y, a$sgr)))
  }))
}))
print(prep_tab, row.names = FALSE)

## Which preparation reproduces the delivered CSV. The delivered data turn out
## NOT to be fully floored: `supplied` (undetected rows -> 0, measured negatives
## kept) matches, `floored` does not. The 100% cap is applied downstream in the
## labs' curve fitting, not in the file, so `floored` is our application of the
## convention and the study must say so.
cat("\nwhich preparation reproduces the delivered SGR column?\n")
print(prep_tab[prep_tab$identical_to_supplied,
               c("dataset", "preparation")], row.names = FALSE)

## ------------------------------------------------------------- reversals ----
cat("\n===== REVERSALS (dominance >= 0.85 of replicate pairs) =====\n")
for (nm in names(dats)) {
  rv <- reversals(dats[[nm]])
  if (nrow(rv)) {
    cat(nm, ":\n"); print(rv, row.names = FALSE)
  } else {
    cat(nm, ": none\n")
  }
}

## --------------------------------------------------------------- figures ----
dir.create("analysis/figures", showWarnings = FALSE, recursive = TRUE)
ggplot2::ggsave("analysis/figures/all_datasets.png",
                plot_all_datasets(dats), width = 8, height = 6, dpi = 200)
ggplot2::ggsave("analysis/figures/r_salina_ordering_pathology.png",
                plot_ordering_pathology(dats$r_salina),
                width = 7, height = 5, dpi = 200)
ggplot2::ggsave("analysis/figures/r_salina2_ordering_pathology.png",
                plot_ordering_pathology(dats$r_salina2),
                width = 7, height = 5, dpi = 200)
cat("\nFigures written to analysis/figures/\n")
