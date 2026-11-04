## Regression tests for the defects corrected in version 0.2.0 --------------

test_that("each observation is divided by the error mean square of its OWN environment", {
  ## In 0.1.1 PooledRBD assigned the error mean square with
  ##   ave(y, Loc, FUN = function(x) MSE[as.character(unique(Loc[y == x]))])
  ## The comparison y == x recycled a group-length vector against the full
  ## response, so rows received another environment's error mean square
  ## whenever response values coincided across environments. The help example
  ## of 0.1.1 had 9 of 18 rows wrong. This test pins the corrected behaviour.
  d <- data.frame(
    Loc = factor(rep(c("L1", "L2"), each = 9)),
    Rep = factor(rep(c("R1", "R2", "R3"), each = 3, times = 2)),
    Treatment = factor(rep(c("T1", "T2", "T3"), times = 6)),
    ## deliberately repeated values across the two locations
    Yield = c(18, 7, 11, 10, 19, 12, 15, 8, 13,
              18, 5, 11, 7, 19, 21, 22, 9, 14))
  o <- PooledRBD(d, "Yield", "Loc", "Treatment", "Rep", 0.05, 1)
  ems <- o$Transformation$error_mean_squares
  ref <- vapply(levels(d$Loc), function(l) {
    s <- droplevels(d[d$Loc == l, ])
    f <- stats::lm(Yield ~ Rep + Treatment, data = s)
    stats::deviance(f) / stats::df.residual(f)
  }, numeric(1))
  expect_equal(ems[levels(d$Loc)], ref[levels(d$Loc)])
  ## and the transformation itself must be a pure per-row lookup
  o2 <- PooledRBD(d, "Yield", "Loc", "Treatment", "Rep", 0.05, 1, transform = "always")
  expect_true(o2$Transformation$applied)
  expect_equal(o2$Pooled_ANOVA$`Mean Sq`[o2$Pooled_ANOVA$Source == "Pooled error"], 1)
})

test_that("results are invariant to the row order of the input data", {
  ## In 0.1.1 PooledStPD built the transformed response with
  ##   unlist(tapply(y, Location, function(Y) Y / sqrt(var(Y))))
  ## which reassembles values in factor-level order, not input order, so every
  ## row was misaligned unless the data happened to be sorted by location.
  data(stpd_example); data(spd_example); data(sugarcane)
  set.seed(99)
  chk <- function(f, dat) {
    a <- f(dat)
    b <- f(dat[sample(nrow(dat)), , drop = FALSE])
    expect_equal(a$Pooled_ANOVA$`Sum Sq`, b$Pooled_ANOVA$`Sum Sq`)
    expect_equal(a$Pooled_ANOVA$`Pr(>F)`, b$Pooled_ANOVA$`Pr(>F)`)
    expect_equal(a$Transformation$error_mean_squares, b$Transformation$error_mean_squares)
  }
  chk(function(x) PooledStPD(x, "Yield", "Location", "Replication", "RowFactor",
                             "ColumnFactor", 0.05, 1, transform = "always"), stpd_example)
  chk(function(x) PooledSPD(x, "Yield", "Location", "Replication", "MainPlot",
                            "SubPlot", 0.05, 1, transform = "always"), spd_example)
  chk(function(x) PooledRBD(x, "ccsy", "env", "gen", "rep", 0.05, 1), sugarcane)
})

test_that("the StPD transformation divides by sqrt(error MS), not sqrt(var(Y))", {
  ## 0.1.1 used the total within-location variance of the raw response, which
  ## is not the quantity defined by Aitken's transformation.
  data(stpd_example)
  o <- PooledStPD(stpd_example, "Yield", "Location", "Replication", "RowFactor",
                  "ColumnFactor", 0.05, 1, transform = "always")
  ref <- vapply(levels(stpd_example$Location), function(l) {
    s <- droplevels(stpd_example[stpd_example$Location == l, ])
    f <- stats::lm(Yield ~ Replication + RowFactor + Replication:RowFactor +
                   ColumnFactor + Replication:ColumnFactor + RowFactor:ColumnFactor,
                   data = s)
    stats::deviance(f) / stats::df.residual(f)
  }, numeric(1))
  expect_equal(o$Transformation$error_mean_squares[levels(stpd_example$Location)],
               ref[levels(stpd_example$Location)])
  vy <- tapply(stpd_example$Yield, stpd_example$Location, stats::var)
  expect_false(isTRUE(all.equal(unname(o$Transformation$error_mean_squares), unname(vy))))
})

test_that("Latin square rows and columns are nested within environment", {
  data(lsd_example)
  o <- PooledLSD(lsd_example, "Yield", "Location", "Treatment", "Row", "Column",
                 0.05, 1)
  t <- o$Pooled_ANOVA
  expect_true("Row within Env" %in% t$Source)
  expect_true("Col within Env" %in% t$Source)
  ## crossing rows with environment, as 0.1.1 did, would leave 3 df for Row
  expect_gt(t$Df[t$Source == "Row within Env"], nlevels(lsd_example$Row) - 1L)
})

test_that("alpha reaches the multiple comparison procedure", {
  data(crd_example)
  a5 <- PooledCRD(crd_example, "Yield", "Loc", "Treatment", 0.05, 1)
  a1 <- PooledCRD(crd_example, "Yield", "Loc", "Treatment", 0.01, 1)
  expect_equal(a5$Treatments_Comparison$parameters$alpha, 0.05)
  expect_equal(a1$Treatments_Comparison$parameters$alpha, 0.01)
  expect_gt(a1$Treatments_Comparison$statistics$Critical_difference,
            a5$Treatments_Comparison$statistics$Critical_difference)
})

test_that("an out-of-range Mult_Comp_Test is an error, never a silent NULL", {
  data(crd_example)
  expect_error(PooledCRD(crd_example, "Yield", "Loc", "Treatment", 0.05, 5),
               "Mult_Comp_Test")
  expect_error(PooledCRD(crd_example, "Yield", "Loc", "Treatment", 0.05, 0),
               "Mult_Comp_Test")
  expect_error(PooledCRD(crd_example, "Yield", "Loc", "Treatment", 0.05, NA),
               "Mult_Comp_Test")
})

test_that("LSD with many treatments warns about the family-wise error rate", {
  data(sugarcane)
  w <- testthat::capture_warnings(PooledRBD(sugarcane, "ccsy", "env", "gen", "rep", 0.05, 3))
  expect_true(any(grepl("family-wise", w)))
})

test_that("a significant interaction triggers the interpretation note", {
  data(sugarcane)
  o <- PooledRBD(sugarcane, "ccsy", "env", "gen", "rep", 0.05, 1)
  expect_s3_class(o$GEI_Note, "cane_note")
  expect_match(o$GEI_Note$message, "Location_wise")
})

test_that("input validation rejects malformed calls with informative messages", {
  data(crd_example)
  expect_error(PooledCRD(crd_example, "NoSuchCol", "Loc", "Treatment", 0.05, 1),
               "was not found in 'data'")
  expect_error(PooledCRD(crd_example, "Yield", "Loc", "Loc", 0.05, 1),
               "more than one argument")
  expect_error(PooledCRD(as.list(crd_example), "Yield", "Loc", "Treatment", 0.05, 1),
               "must be a data frame")
  expect_error(PooledCRD(crd_example, "Yield", "Loc", "Treatment", 1.5, 1), "alpha")
  dchr <- crd_example; dchr$Label <- as.character(dchr$Treatment)
  expect_error(PooledCRD(dchr, "Label", "Loc", "Treatment", 0.05, 1), "must be numeric")
  one <- crd_example[crd_example$Loc == levels(crd_example$Loc)[1], ]
  expect_error(PooledCRD(one, "Yield", "Loc", "Treatment", 0.05, 1),
               "at least two environments")
})

test_that("split plot and strip plot refuse unreplicated series", {
  data(spd_example)
  d1 <- spd_example[spd_example$Replication == "R1", ]
  expect_error(suppressWarnings(PooledSPD(d1, "Yield", "Location", "Replication",
                         "MainPlot", "SubPlot", 0.05, 1)), "at least two replications")
})

test_that("transform argument overrides the automatic decision", {
  data(sugarcane)
  expect_false(PooledRBD(sugarcane, "ccsy", "env", "gen", "rep", 0.05, 1,
                         transform = "never")$Transformation$applied)
  expect_true(PooledRBD(sugarcane, "ccsy", "env", "gen", "rep", 0.05, 1,
                        transform = "always")$Transformation$applied)
  data(crd_example)
  expect_true(PooledCRD(crd_example, "Yield", "Loc", "Treatment", 0.05, 1,
                        transform = "always")$Transformation$applied)
})

test_that("env_effect = 'random' tests the treatment effect against the interaction", {
  data(sugarcane)
  rf <- PooledRBD(sugarcane, "ccsy", "env", "gen", "rep", 0.05, 1,
                  env_effect = "fixed")$Pooled_ANOVA
  rr <- PooledRBD(sugarcane, "ccsy", "env", "gen", "rep", 0.05, 1,
                  env_effect = "random")$Pooled_ANOVA
  expect_identical(rf$`Error term`[rf$Source == "Treatment"], "Pooled error")
  expect_identical(rr$`Error term`[rr$Source == "Treatment"], "Treatment x Env")
  expect_equal(rr$`F value`[rr$Source == "Treatment"],
               rr$`Mean Sq`[rr$Source == "Treatment"] /
               rr$`Mean Sq`[rr$Source == "Treatment x Env"])
})

test_that("split plot comparisons use the stratum of the factor compared", {
  data(spd_example)
  o <- PooledSPD(spd_example, "Yield", "Location", "Replication", "MainPlot",
                 "SubPlot", 0.05, 1)
  t <- o$Pooled_ANOVA
  mp <- o$Treatments_Comparison$MainPlot
  sp <- o$Treatments_Comparison$SubPlot
  expect_identical(mp$parameters$error_stratum, "Main plot error")
  expect_identical(sp$parameters$error_stratum, "Subplot error")
  expect_equal(mp$statistics$MSerror, t$`Mean Sq`[t$Source == "Main plot error"])
  expect_equal(sp$statistics$MSerror, t$`Mean Sq`[t$Source == "Subplot error"])
  expect_false(isTRUE(all.equal(mp$statistics$MSerror, sp$statistics$MSerror)))
})

test_that("strip plot comparisons use all three strata", {
  data(stpd_example)
  o <- PooledStPD(stpd_example, "Yield", "Location", "Replication", "RowFactor",
                  "ColumnFactor", 0.05, 1)
  tc <- o$Treatments_Comparison
  expect_identical(tc$RowFactor$parameters$error_stratum,    "Row strip error")
  expect_identical(tc$ColumnFactor$parameters$error_stratum, "Column strip error")
  expect_identical(tc$Interaction$parameters$error_stratum,  "Intersection error")
})

test_that("transformed-scale output also reports means on the original scale", {
  data(sugarcane)
  o <- PooledRBD(sugarcane, "ccsy", "env", "gen", "rep", 0.05, 1)
  g <- o$Treatments_Comparison$groups
  expect_true("Mean_original_scale" %in% names(g))
  expect_identical(o$Treatments_Comparison$parameters$scale, "Aitken-transformed")
  raw <- tapply(sugarcane$ccsy, sugarcane$gen, mean)
  expect_equal(g$Mean_original_scale, unname(raw[g$Treatment]))
})

test_that("missing values are dropped with a warning and the analysis still runs", {
  data(crd_example)
  d <- crd_example; d$Yield[3] <- NA
  w <- testthat::capture_warnings(o <- PooledCRD(d, "Yield", "Loc", "Treatment", 0.05, 1))
  expect_true(any(grepl("missing values", w)))
  expect_s3_class(o, "cane")
})

test_that("print and summary methods run without error", {
  data(sugarcane); data(stpd_example)
  o  <- PooledRBD(sugarcane, "ccsy", "env", "gen", "rep", 0.05, 1)
  o2 <- PooledStPD(stpd_example, "Yield", "Location", "Replication", "RowFactor",
                   "ColumnFactor", 0.05, 1)
  expect_output(print(o))
  expect_output(print(o2))
  expect_invisible(s <- summary(o))
  expect_named(s, c("Pooled_ANOVA", "Bartlett_Test", "Levene_Test", "Groups"))
  expect_output(print(o$Pooled_ANOVA))
  expect_output(print(o$Treatments_Comparison))
  expect_output(print(o$Bartlett_Test))
})
