library(testthat)
source("../R/data_prep.R")
source("../R/diagnostics.R")
source("../R/simulate.R")
source("../R/metrics.R")
# Phase 10 helpers. No top-level side effects and no bayesnec dependency at
# source time, so they load here like the rest; the parts that need a fitted
# object are checked by the Phase 10 gates instead.
source("../R/setup_phase10.R")
source("../R/model_average.R")
test_dir("testthat")
