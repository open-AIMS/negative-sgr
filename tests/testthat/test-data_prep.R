## The known data issues, encoded so that a change to the CSVs or to the
## preparation logic is caught rather than absorbed.

dats <- lapply(dataset_names(), read_sgr, path = test_path("../../data-raw"))
names(dats) <- dataset_names()

test_that("the test design recovers exactly from density/SGR pairs", {
  # An R^2 of 1 means the supplied SGR column was computed deterministically
  # from these densities and a single inoculum density, so t and n_0 are facts
  # about the file, not estimates.
  for (nm in names(dats)) {
    d <- recover_design(dats[[nm]])
    expect_gt(d$r_squared, 0.9999)
  }
  expect_equal(recover_design(dats$c_proliferum)$t, 7, tolerance = 1e-3)
  expect_equal(recover_design(dats$c_proliferum)$n_0, 7968, tolerance = 1e-3)
  expect_equal(recover_design(dats$r_salina)$t, 3, tolerance = 1e-3)
  expect_equal(recover_design(dats$r_salina)$n_0, 3871, tolerance = 1e-3)
})

test_that("both Rhodomonas LODs are 10, and the data are consistent with that", {
  # The LOD is a protocol fact. The data can only ever bound it from above, so
  # this asserts consistency, not recovery. The earlier version of this test
  # asserted lod == min_detected_density() and called that confirmation "without
  # needing the protocol" -- which is the same inference twice, and is how
  # r_salina2 came to be recorded at 100.
  meta <- dataset_meta()
  expect_equal(meta$lod[meta$dataset == "r_salina"], 10)
  expect_equal(meta$lod[meta$dataset == "r_salina2"], 10)

  for (nm in c("r_salina", "r_salina2")) {
    d <- dats[[nm]]
    lod <- meta$lod[meta$dataset == nm]
    # Nothing may be recorded below the limit except a non-detect.
    expect_lte(lod, min_detected_density(d))
    # The reporting grid is 10 in both tests, which is what rules out 100 for
    # r_salina2: a limit of 100 cannot explain recorded densities of 230, 290.
    pos <- d$density[d$density > 0]
    expect_true(all(pos %% 10 == 0))
    expect_false(all(pos %% 100 == 0))
  }

  # r_salina is the test where a replicate actually landed on the limit, so
  # there the bound is directly visible in the data.
  d <- dats$r_salina
  design <- recover_design(d)
  at_lod <- d$sgr[d$density == 10]
  expect_gt(length(at_lod), 0)
  expect_equal(unname(at_lod), rep(lod_bound(10, design$n_0, design$t),
                                   length(at_lod)), tolerance = 1e-3)

  # r_salina2 is the test where none did: its smallest detected density is 100,
  # well above the limit, so the bound is an extrapolation from the protocol.
  expect_equal(min_detected_density(dats$r_salina2), 100)
  expect_equal(sum(dats$r_salina2$density == 10), 0)
})

test_that("zeros in the Rhodomonas SGR columns are undefined, not measured", {
  # Every zero in the supplied SGR column corresponds to a zero cell density,
  # where SGR is log(0) -- undefined, not zero.
  for (nm in c("r_salina", "r_salina2")) {
    d <- dats[[nm]]
    expect_identical(which(d$sgr == 0), which(d$density == 0))
  }
  # The Cladocopium sets have no such rows at all.
  for (nm in c("c_proliferum", "c_proliferum2")) {
    expect_equal(sum(dats[[nm]]$density == 0), 0)
    expect_equal(sum(dats[[nm]]$sgr == 0), 0)
  }
})

test_that("the supplied CSVs are NOT fully floored", {
  # The labs substituted 0 only for undetected rows; measured negatives survive
  # in the delivered data. The 100% cap is applied downstream in their curve
  # fitting, not in the file. `supplied` reproduces the column exactly;
  # `floored` is our application of the convention and differs from it.
  for (nm in c("r_salina", "r_salina2")) {
    sup <- prepare_sgr(dats[[nm]], "supplied")
    expect_equal(sup$y, sup$sgr)
    expect_gt(sum(sup$y < 0), 0)
    fl <- prepare_sgr(dats[[nm]], "floored")
    expect_false(isTRUE(all.equal(fl$y, fl$sgr)))
  }
  # On the Cladocopium sets `supplied` is the untouched column (no zero rows),
  # and flooring is a real change.
  for (nm in c("c_proliferum", "c_proliferum2")) {
    expect_equal(prepare_sgr(dats[[nm]], "supplied")$y, dats[[nm]]$sgr)
    fl <- prepare_sgr(dats[[nm]], "floored")
    expect_false(isTRUE(all.equal(fl$y, fl$sgr)))
    expect_equal(min(fl$y), 0)
  }
})

test_that("only undetected rows are censored, not rows at the limit", {
  # A density of exactly the LOD was detected. Its response value equals the
  # bound, but it is an observation. Several rows therefore share one response
  # value with different censoring codes, which is correct.
  d <- dats$r_salina
  cn <- prepare_sgr(d, "censored")
  bound <- attr(cn, "lod_bound")
  expect_equal(sum(cn$cens == "left"), sum(d$density == 0))
  at_bound <- which(abs(cn$y - bound) < 1e-8)
  expect_true(any(cn$cens[at_bound] == "left"))
  expect_true(any(cn$cens[at_bound] == "none"))
})

test_that("the raw preparation drops undefined rows rather than inventing them", {
  for (nm in c("r_salina", "r_salina2")) {
    rw <- prepare_sgr(dats[[nm]], "raw")
    expect_equal(nrow(rw), nrow(dats[[nm]]) - sum(dats[[nm]]$density == 0))
    expect_true(all(rw$density > 0))
  }
  for (nm in c("c_proliferum", "c_proliferum2")) {
    expect_equal(nrow(prepare_sgr(dats[[nm]], "raw")), nrow(dats[[nm]]))
  }
})

test_that("c_proliferum reverses in the declining limb at 15 -> 20", {
  # 22 of 25 replicate pairs at 20 ug/L sit above those at 15 ug/L, and mean
  # final density more than doubles (ratio 2.28). Not noise; flagged, not
  # modelled around. Dominance is 0.88 rather than 1 because the treatments
  # overlap at their edges, which is why the diagnostic reports the statistic
  # instead of a bare TRUE/FALSE.
  rv <- reversals(dats$c_proliferum)
  in_dec <- rv[rv$in_decline, ]
  expect_equal(nrow(in_dec), 1)
  expect_equal(in_dec$from, 15)
  expect_equal(in_dec$to, 20)
  expect_equal(in_dec$dominance, 0.88)
  expect_gt(in_dec$density_ratio, 2)
})

test_that("r_salina's detection floor inverts the dose ordering", {
  # At 7.5 ug/L replicates report SGR down to -1.99; at 10 and 15 every
  # replicate is recorded as 0.000. The most affected treatments are recorded
  # as the least affected.
  d <- dats$r_salina
  expect_lt(min(d$sgr[d$x == 7.5]), -1.9)
  expect_true(all(d$sgr[d$x %in% c(10, 15)] == 0))
  expect_true(all(d$density[d$x %in% c(10, 15)] == 0))
})

test_that("the discarded effect fraction orders differently from R", {
  # The study's discriminating claim: Delta and R do not order these datasets
  # the same way, so the four datasets can (weakly) tell them apart.
  s <- do.call(rbind, lapply(dats, dataset_summary))
  by_delta <- s$dataset[order(-s$delta)]
  by_r <- s$dataset[order(s$R)]
  expect_equal(by_delta[1], "c_proliferum")
  expect_false(identical(by_delta, by_r))
})

test_that("percent_inhibition caps at the control and floors at zero", {
  pi_dat <- percent_inhibition(dats$c_proliferum)
  expect_true(all(pi_dat$y >= 0))
  expect_equal(sum(pi_dat$y == 0), sum(dats$c_proliferum$sgr < 0))
})

test_that("truncate_at_crossing keeps only concentrations below the crossing", {
  d <- prepare_sgr(dats$c_proliferum, "raw")
  tr <- suppressWarnings(truncate_at_crossing(d, 12))
  expect_true(all(tr$x < 12))
  expect_equal(attr(tr, "crossing"), 12)
  expect_warning(truncate_at_crossing(d, 0.05))
})
