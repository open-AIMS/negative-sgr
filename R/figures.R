## Figures. Every function returns a ggplot object; saving is the caller's job.

library(ggplot2)

#' Put a control at zero concentration onto a log axis
#'
#' Concentration-response designs include a zero control, which a log axis
#' cannot show. The convention used throughout: place the control one geometric
#' step below the lowest tested concentration and label the break "0".
log_x_with_control <- function(x) {
  pos <- sort(unique(x[x > 0]))
  step <- if (length(pos) > 1) exp(mean(diff(log(pos)))) else 10
  ctrl_at <- min(pos) / step
  out <- x
  out[x == 0] <- ctrl_at
  attr(out, "control_at") <- ctrl_at
  out
}

#' The detection-floor ordering pathology in r_salina
#'
#' At 7.5 ug/L three replicates were counted at or near the detection limit and
#' report SGR down to -1.99. At 10 and 15 ug/L nothing was detected at all, and
#' every replicate is recorded as 0.000 -- so the treatments where the
#' population was wiped out most completely are recorded as the *least*
#' affected. Any fit to the as-supplied column is fitting that artefact.
#'
#' The figure shows the supplied values, the LOD bound the undetected rows
#' actually imply, and an arrow from one to the other.
plot_ordering_pathology <- function(dat, meta = dataset_meta()) {
  ds <- unique(dat$dataset)
  bnd <- prepare_sgr(dat, "censored", meta = meta)
  bound <- attr(bnd, "lod_bound")
  lod <- meta$lod[meta$dataset == ds]

  pl <- dat
  pl$xpos <- log_x_with_control(dat$x)
  ctrl_at <- attr(pl$xpos, "control_at")
  pl$state <- ifelse(dat$density == 0, "not detected (supplied as 0.000)",
                     ifelse(dat$density <= lod, "detected at the limit",
                            "counted"))
  moved <- pl[pl$density == 0, ]
  moved$yend <- bound

  ggplot(pl, aes(x = xpos, y = sgr)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
    geom_hline(yintercept = bound, linetype = "dashed", linewidth = 0.4,
               colour = "firebrick") +
    geom_segment(data = moved,
                 aes(x = xpos, xend = xpos, y = sgr, yend = yend),
                 arrow = grid::arrow(length = grid::unit(0.12, "cm")),
                 colour = "firebrick", alpha = 0.6, linewidth = 0.3) +
    geom_point(aes(shape = state, colour = state), size = 2.2,
               position = position_jitter(width = 0.03, height = 0)) +
    annotate("text", x = ctrl_at, y = bound, vjust = -0.6, hjust = 0,
             size = 3, colour = "firebrick",
             label = paste0("left-censoring bound at ", lod,
                            " cells/mL = ", signif(bound, 4))) +
    scale_x_log10(breaks = c(ctrl_at, sort(unique(dat$x[dat$x > 0]))),
                  labels = c("0", sort(unique(dat$x[dat$x > 0])))) +
    scale_colour_manual(values = c("counted" = "grey20",
                                   "detected at the limit" = "steelblue",
                                   "not detected (supplied as 0.000)" = "firebrick")) +
    scale_shape_manual(values = c("counted" = 16,
                                  "detected at the limit" = 17,
                                  "not detected (supplied as 0.000)" = 1)) +
    labs(x = "Concentration", y = "Specific growth rate (per day)",
         colour = NULL, shape = NULL,
         title = paste0(ds, ": what the supplied SGR column does to ",
                        "undetected populations")) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom", legend.direction = "vertical",
          panel.grid.minor = element_blank())
}

#' Raw SGR against concentration for all four datasets
plot_all_datasets <- function(dats) {
  pl <- do.call(rbind, lapply(dats, function(d) {
    d$xpos <- log_x_with_control(d$x)
    d$undefined <- d$density == 0
    d
  }))
  ggplot(pl, aes(x = xpos, y = sgr, shape = undefined)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
    geom_point(size = 1.6, alpha = 0.8) +
    scale_x_log10() +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1),
                       labels = c("counted", "not detected")) +
    facet_wrap(~dataset, scales = "free") +
    labs(x = "Concentration (control plotted one step below the lowest dose)",
         y = "Specific growth rate (per day)", shape = NULL) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
}

#' Fitted curves for every arm of one dataset, overlaid
plot_arm_curves <- function(arm_fits, dat, resolution = 200) {
  x_rng <- range(dat$x)
  x_seq <- seq(x_rng[1], x_rng[2], length.out = resolution)
  curves <- do.call(rbind, lapply(arm_fits, function(a) {
    if (identical(a$arm, "SQ")) return(NULL)
    p <- stats::fitted(a$fit, newdata = data.frame(x = x_seq), re_formula = NA)
    data.frame(arm = a$arm, x = x_seq, y = p[, "Estimate"],
               lo = p[, "Q2.5"], hi = p[, "Q97.5"])
  }))
  pl <- dat
  pl$undefined <- dat$density == 0
  ggplot() +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
    geom_ribbon(data = curves, aes(x = x, ymin = lo, ymax = hi, fill = arm),
                alpha = 0.12) +
    geom_line(data = curves, aes(x = x, y = y, colour = arm), linewidth = 0.6) +
    geom_point(data = pl, aes(x = x, y = sgr, shape = undefined),
               size = 1.5, alpha = 0.7) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), guide = "none") +
    labs(x = "Concentration", y = "Specific growth rate (per day)",
         colour = "Arm", fill = "Arm",
         title = unique(dat$dataset)) +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank())
}
