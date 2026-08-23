## Every number the vignette's model-averaging section quotes, generated.
##
## Same discipline as `vignette_tables.R`: the vignette must never carry
## hand-typed results. Re-run this after any change to the Phase 10 sweep --
## in particular after the 101-200 top-up -- and update the section from its
## output rather than editing the prose numbers in place.
##
## Run: Rscript analysis/phase10_vignette_numbers.R > analysis/phase10_vignette_numbers.txt

source("R/setup.R"); source("R/setup_phase10.R"); source("R/data_prep.R")
source("R/simulate.R"); source("R/metrics.R"); source("R/model_average.R")

fs <- list.files("analysis/phase10", pattern = "[.]rds$", full.names = TRUE)
bl <- lapply(fs, readRDS)
e <- do.call(rbind, lapply(bl, `[[`, "endpoints"))
w <- do.call(rbind, lapply(bl, `[[`, "weights"))
r <- do.call(rbind, lapply(bl, `[[`, "rhat"))
ec <- e[e$stage == "converged" & e$ok & !is.na(e$estimate), ]

lab <- c("d4.0_t2.0_R73.0_s8.1" = "cell 12 (control CV 1.9%)",
         "d4.0_t2.0_R2.3_s8.1"  = "cell 8  (control CV 9.6%)",
         "d8.0_t2.0_R2.3_s8.1"  = "cell 9  (CV 9.6%, delta 8)")

cat("PHASE 10 -- numbers quoted by the vignette\n")
cat("iterations per cell:",
    paste(tapply(ec$iteration, ec$cell, function(z) length(unique(z))),
          collapse = " / "), "\n")
cat("total endpoint rows:", nrow(e), "| failures:", sum(!e$ok),
    "| arm-iterations with no surviving model:",
    sum(!e$ok & !is.na(e$n_kept) & e$n_kept == 0), "\n\n")

for (cl in names(lab)) {
  z <- ec[ec$cell == cl, ]
  if (!nrow(z)) next
  iters <- sort(unique(z$iteration))
  s <- do.call(rbind, lapply(sweep_files(cl), readRDS))
  s <- s[s$arm %in% P10_ARMS & s$ok & !is.na(s$estimate) & s$iteration %in% iters, ]
  k <- c("arm", "endpoint", "iteration")
  m <- merge(z[, c(k, "lower", "upper", "estimate", "truth")],
             s[, c(k, "lower", "upper", "estimate")], by = k,
             suffixes = c("_a", "_s"))
  out <- do.call(rbind, lapply(split(m, ~ arm + endpoint, drop = TRUE),
    function(d) {
      ma <- mcse_summary(d$estimate_a, d$truth[1], d$lower_a, d$upper_a)
      ms <- mcse_summary(d$estimate_s, d$truth[1], d$lower_s, d$upper_s)
      data.frame(arm = d$arm[1], endpoint = d$endpoint[1],
                 bias_s = round(100 * ms$bias / d$truth[1], 1),
                 bias_a = round(100 * ma$bias / d$truth[1], 1),
                 bias_a_mcse = round(100 * ma$bias_mcse / d$truth[1], 2),
                 cov_s = ms$coverage, cov_a = ma$coverage,
                 rmse_s = round(ms$rmse, 2), rmse_a = round(ma$rmse, 2),
                 width_ratio = round(median((d$upper_a - d$lower_a) /
                                            (d$upper_s - d$lower_s)), 2))
    }))
  cat("== ", lab[[cl]], " (n = ", length(iters), ") ==\n", sep = "")
  for (ep in c("ErC10", "ErC50", "NSEC")) {
    y <- out[out$endpoint == ep, ]; y <- y[match(P10_ARMS, y$arm), ]
    cat("  ", ep, "\n", sep = "")
    print(y[, c("arm", "bias_s", "bias_a", "bias_a_mcse", "cov_s", "cov_a",
                "rmse_s", "rmse_a", "width_ratio")], row.names = FALSE)
  }
  ## Mean stacking weight on the generating model, completed against the full
  ## grid so that iterations where a model was dropped count as weight zero.
  wc <- w[w$stage == "converged" & w$cell == cl, ]
  g <- do.call(rbind, lapply(split(wc, ~ arm, drop = TRUE), function(d)
    expand.grid(arm = d$arm[1], model = unique(d$model),
                iteration = unique(d$iteration), stringsAsFactors = FALSE)))
  wf <- merge(g, wc[, c("arm", "model", "iteration", "wi")], all.x = TRUE)
  wf$wi[is.na(wf$wi)] <- 0
  mw <- aggregate(wi ~ arm + model, wf, mean)
  cat("  weight on nec4param: ",
      paste(sprintf("%s %.3f", P10_ARMS, sapply(P10_ARMS, function(a) {
        v <- mw$wi[mw$arm == a & mw$model == "nec4param"]
        if (length(v)) v else 0 })), collapse = "  "), "\n")
  top <- do.call(rbind, lapply(split(mw, ~ arm, drop = TRUE), function(d)
    d[order(-d$wi), ][1, ]))
  cat("  top-weighted model:   ",
      paste(sprintf("%s %s %.2f", top$arm, top$model, top$wi), collapse = "  "),
      "\n")
  nn <- aggregate(dropped ~ arm + iteration, r[r$cell == cl, ], sum)
  cat("  models dropped/arm-iteration: ",
      paste(sprintf("%s %.2f", P10_ARMS,
                    sapply(P10_ARMS, function(a) mean(nn$dropped[nn$arm == a]))),
            collapse = "  "), "\n\n")
}
