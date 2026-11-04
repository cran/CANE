# tests/testthat/test-CANE.R
# ---------------------------------------------------------------------------
# Tests for the CANE package (Comprehensive Groups of Experiments Analysis
# for Numerous Environments), version 0.2.0
#
# Scope of these tests:
#   * every exported pooled-analysis function runs on its documented example
#     data without error;
#   * each returns the components promised in the reference manual;
#   * each accepts all three documented multiple comparison test codes
#     (1 = Tukey's HSD, 2 = DMRT, 3 = LSD);
#   * results are deterministic across repeated calls.
#
# No ANOVA sums of squares, F values or p values are hardcoded here. Those
# should be locked in with expect_snapshot() only after you have inspected
# and verified the output yourself (see the template at the bottom of this
# file).
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Run an expression while swallowing printed output, messages and warnings,
# so that the test report stays readable. Errors still propagate and fail
# the test, which is what we want.
run_quiet <- function(expr) {
  res <- NULL
  invisible(utils::capture.output(
    res <- suppressWarnings(suppressMessages(expr))
  ))
  res
}

# Components returned by every pooled function (see ?PooledCRD, Value)
core_components <- c(
  "Design",
  "Individual_ANOVA",
  "Location_wise",
  "Bartlett_Test",
  "Levene_Test",
  "Normality_Test",
  "Transformation",
  "Pooled_ANOVA",
  "Treatments_Comparison",
  "Environment_Key",
  "Settings"
)

# PooledSPD and PooledStPD return the same components; their
# Treatments_Comparison is a list with one comparison per error stratum
# (see ?PooledSPD and ?PooledStPD, Value)
strip_components <- core_components
spd_strata  <- c("MainPlot", "SubPlot")
stpd_strata <- c("RowFactor", "ColumnFactor", "Interaction")

# Check that a stratified Treatments_Comparison has the documented parts
expect_strata <- function(out, strata) {
  tc <- out$Treatments_Comparison
  expect_identical(names(tc), strata)
  for (nm in strata) expect_s3_class(tc[[nm]], "cane_comparison")
  invisible(out)
}

# Shared structural check applied to the output of every pooled function
expect_pooled_output <- function(out, expected_names) {
  expect_s3_class(out, "cane")
  expect_true(is.list(out))
  # setdiff() reports exactly which component is absent if this fails
  expect_identical(setdiff(expected_names, names(out)), character(0))
  for (nm in expected_names) {
    expect_false(is.null(out[[nm]]))
  }
  expect_gt(NROW(out$Pooled_ANOVA), 0)
  invisible(out)
}


# ---------------------------------------------------------------------------
# Example datasets (taken verbatim from the CANE reference manual)
# ---------------------------------------------------------------------------

crd_data <- function() {
  data.frame(
    Loc       = factor(rep(c("L1", "L2"), each = 9)),
    Treatment = factor(rep(c("T1", "T2", "T3"), times = 6)),
    Yield     = c(18, 7, 11, 10, 19, 12, 15, 8, 13,
                  18, 5, 11,  7, 19, 21, 22, 9, 14)
  )
}

rbd_data <- function() {
  data.frame(
    Loc       = factor(rep(c("L1", "L2"), each = 9)),
    Rep       = factor(rep(c("R1", "R2", "R3"), each = 3, times = 2)),
    Treatment = factor(rep(c("T1", "T2", "T3"), times = 6)),
    Yield     = c(18, 7, 11, 10, 19, 12, 15, 8, 13,
                  18, 5, 11,  7, 19, 21, 22, 9, 14)
  )
}

lsd_data <- function() {
  data.frame(
    Location  = factor(c(rep("Agra", 16), rep("Bihar", 16))),
    Treatment = factor(c(4, 2, 3, 1, 3, 1, 4, 2, 1, 4, 2, 3, 2, 3, 1, 4,
                         4, 2, 3, 1, 3, 1, 4, 2, 1, 4, 2, 3, 2, 3, 1, 4)),
    Row       = factor(rep(1:4, each = 4, times = 2)),
    Column    = factor(rep(1:4, times = 8)),
    Yield     = c(29.1, 18.9, 29.4,  5.7, 16.4, 10.2, 21.2, 19.1,
                  5.4, 38.8, 24.0, 37.0, 24.9, 41.7,  9.5, 28.9,
                  13.0, 27.0,  7.0, 24.0, 10.0, 13.0, 41.0, 39.0,
                  42.0, 26.0, 19.0, 20.0, 20.0, 35.0,  8.0, 38.0)
  )
}

# The split plot and strip plot examples in the manual share the same layout
# and the same response values, so one generator serves both after renaming.
spd_data <- function() {
  data.frame(
    Location    = factor(rep(c("Londan", "Agumbe"), each = 12)),
    Replication = factor(rep(c(1, 2), each = 6, times = 2)),
    MainPlot    = factor(rep(c(1, 2), each = 3, times = 4)),
    SubPlot     = factor(rep(1:3, times = 8)),
    Yield       = c(4940, 4810, 5150, 4900, 4920, 5070,
                    4830, 5110, 4920, 5020, 5110, 5230,
                    4964, 4997, 5011, 5102, 4858, 4888,
                    5100, 5165, 4965, 5113, 5086, 5176)
  )
}

stpd_data <- function() {
  df <- spd_data()
  names(df)[names(df) == "MainPlot"] <- "RowFactor"
  names(df)[names(df) == "SubPlot"]  <- "ColumnFactor"
  df
}


# ---------------------------------------------------------------------------
# Exported API
# ---------------------------------------------------------------------------

test_that("all five pooled analysis functions are exported", {
  exported <- getNamespaceExports("CANE")
  expect_true(all(c("PooledCRD", "PooledLSD", "PooledRBD",
                    "PooledSPD", "PooledStPD") %in% exported))
})

test_that("exported functions take the documented arguments", {
  expect_identical(
    names(formals(PooledCRD)),
    c("data", "Response", "Location", "Treatment", "alpha", "Mult_Comp_Test",
      "transform", "env_effect")
  )
  expect_identical(
    names(formals(PooledRBD)),
    c("data", "Response", "Location", "Treatment", "Replication",
      "alpha", "Mult_Comp_Test", "transform", "env_effect")
  )
  expect_identical(
    names(formals(PooledLSD)),
    c("data", "Response", "Location", "Treatment", "Row", "Column",
      "alpha", "Mult_Comp_Test", "transform", "env_effect")
  )
  expect_identical(
    names(formals(PooledSPD)),
    c("data", "Response", "Location", "Replication", "MainPlot", "SubPlot",
      "alpha", "Mult_Comp_Test", "transform", "env_effect")
  )
  expect_identical(
    names(formals(PooledStPD)),
    c("data", "Response", "Location", "Replication", "RowFactor",
      "ColumnFactor", "alpha", "Mult_Comp_Test", "transform", "env_effect")
  )
})


# ---------------------------------------------------------------------------
# PooledCRD
# ---------------------------------------------------------------------------

test_that("PooledCRD returns the documented components", {
  df  <- crd_data()
  out <- run_quiet(PooledCRD(df, "Yield", "Loc", "Treatment", 0.05, 1))
  expect_pooled_output(out, core_components)
})

test_that("PooledCRD accepts Tukey's HSD, DMRT and LSD", {
  df <- crd_data()
  for (mct in 1:3) {
    out <- run_quiet(PooledCRD(df, "Yield", "Loc", "Treatment", 0.05, mct))
    expect_pooled_output(out, core_components)
  }
})


# ---------------------------------------------------------------------------
# PooledRBD
# ---------------------------------------------------------------------------

test_that("PooledRBD returns the documented components", {
  df  <- rbd_data()
  out <- run_quiet(PooledRBD(df, "Yield", "Loc", "Treatment", "Rep", 0.05, 1))
  expect_pooled_output(out, core_components)
})

test_that("PooledRBD accepts Tukey's HSD, DMRT and LSD", {
  df <- rbd_data()
  for (mct in 1:3) {
    out <- run_quiet(
      PooledRBD(df, "Yield", "Loc", "Treatment", "Rep", 0.05, mct)
    )
    expect_pooled_output(out, core_components)
  }
})


# ---------------------------------------------------------------------------
# PooledLSD
# ---------------------------------------------------------------------------

test_that("PooledLSD returns the documented components", {
  df  <- lsd_data()
  out <- run_quiet(
    PooledLSD(df, "Yield", "Location", "Treatment", "Row", "Column", 0.05, 1)
  )
  expect_pooled_output(out, core_components)
})

test_that("PooledLSD accepts Tukey's HSD, DMRT and LSD", {
  df <- lsd_data()
  for (mct in 1:3) {
    out <- run_quiet(
      PooledLSD(df, "Yield", "Location", "Treatment", "Row", "Column",
                0.05, mct)
    )
    expect_pooled_output(out, core_components)
  }
})


# ---------------------------------------------------------------------------
# PooledSPD
# ---------------------------------------------------------------------------

test_that("PooledSPD returns the documented components", {
  df  <- spd_data()
  out <- run_quiet(
    PooledSPD(df, "Yield", "Location", "Replication", "MainPlot", "SubPlot",
              0.05, 1)
  )
  expect_pooled_output(out, core_components)
  expect_strata(out, spd_strata)
})

test_that("PooledSPD accepts Tukey's HSD, DMRT and LSD", {
  df <- spd_data()
  for (mct in 1:3) {
    out <- run_quiet(
      PooledSPD(df, "Yield", "Location", "Replication", "MainPlot", "SubPlot",
                0.05, mct)
    )
    expect_pooled_output(out, core_components)
    expect_strata(out, spd_strata)
  }
})


# ---------------------------------------------------------------------------
# PooledStPD
# ---------------------------------------------------------------------------

test_that("PooledStPD returns the documented components", {
  df  <- stpd_data()
  out <- run_quiet(
    PooledStPD(df, "Yield", "Location", "Replication", "RowFactor",
               "ColumnFactor", 0.05, 1)
  )
  expect_pooled_output(out, strip_components)
  expect_strata(out, stpd_strata)
})

test_that("PooledStPD accepts Tukey's HSD, DMRT and LSD", {
  df <- stpd_data()
  for (mct in 1:3) {
    out <- run_quiet(
      PooledStPD(df, "Yield", "Location", "Replication", "RowFactor",
                 "ColumnFactor", 0.05, mct)
    )
    expect_pooled_output(out, strip_components)
    expect_strata(out, stpd_strata)
  }
})


# ---------------------------------------------------------------------------
# Determinism
# ---------------------------------------------------------------------------

test_that("repeated calls give identical pooled ANOVA tables", {
  df   <- rbd_data()
  one  <- run_quiet(PooledRBD(df, "Yield", "Loc", "Treatment", "Rep", 0.05, 1))
  two  <- run_quiet(PooledRBD(df, "Yield", "Loc", "Treatment", "Rep", 0.05, 1))
  expect_equal(one$Pooled_ANOVA, two$Pooled_ANOVA)
  expect_equal(one$Bartlett_Test, two$Bartlett_Test)
})

test_that("the choice of multiple comparison test does not alter the ANOVA", {
  df    <- crd_data()
  tukey <- run_quiet(PooledCRD(df, "Yield", "Loc", "Treatment", 0.05, 1))
  lsd   <- run_quiet(PooledCRD(df, "Yield", "Loc", "Treatment", 0.05, 3))
  expect_equal(tukey$Pooled_ANOVA, lsd$Pooled_ANOVA)
})


# ---------------------------------------------------------------------------
# Optional: snapshot tests
# ---------------------------------------------------------------------------
# Once you have run the examples yourself and confirmed the numerical output
# is correct, uncomment the block below. On the first run testthat records the
# printed output under tests/testthat/_snaps/CANE.md; on later runs it flags
# any change. This locks in the numbers without anyone having to type them.
#
# test_that("PooledRBD output is stable", {
#   df <- rbd_data()
#   expect_snapshot(
#     print(PooledRBD(df, "Yield", "Loc", "Treatment", "Rep", 0.05, 1))
#   )
# })
