## Regression tests for the corrections made in response to the review of
## the CANE software note (APPS-D-26-00091) ---------------------------------

## Bartlett's statistic computed by hand from mean squares and their df
hand_bartlett <- function(ms, df) {
  N <- sum(df); k <- length(ms); sp2 <- sum(df * ms) / N
  (N * log(sp2) - sum(df * log(ms))) / (1 + (sum(1 / df) - 1 / N) / (3 * (k - 1)))
}

test_that("GEI_Note checks every treatment-by-environment term (split plot)", {
  data(spd_example)
  d <- spd_example
  ## build a strong main-plot-by-environment interaction into the data
  shift <- ifelse(d$Location == levels(d$Location)[1] &
                    d$MainPlot == levels(d$MainPlot)[1], 400, -400)
  d$Yield <- d$Yield + shift
  o <- suppressWarnings(PooledSPD(d, "Yield", "Location", "Replication",
                                  "MainPlot", "SubPlot", 0.05, 1))
  t <- o$Pooled_ANOVA
  expect_lt(t$`Pr(>F)`[t$Source == "Main plot x Env"], 0.05)
  expect_gt(t$`Pr(>F)`[t$Source == "Subplot x Env"], 0.05)
  expect_s3_class(o$GEI_Note, "cane_note")
  expect_true("Main plot x Env" %in% o$GEI_Note$terms)
  expect_match(o$GEI_Note$message, "Main plot x Env")
  ## and no note when no interaction is significant
  o0 <- suppressWarnings(PooledSPD(spd_example, "Yield", "Location", "Replication",
                                   "MainPlot", "SubPlot", 0.05, 1))
  expect_null(o0$GEI_Note)
})

test_that("GEI_Note checks every treatment-by-environment term (strip plot)", {
  data(stpd_example)
  d <- stpd_example
  shift <- ifelse(d$Location == levels(d$Location)[1] &
                    d$RowFactor == levels(d$RowFactor)[1], 500, -500)
  d$Yield <- d$Yield + shift
  o <- suppressWarnings(PooledStPD(d, "Yield", "Location", "Replication",
                                   "RowFactor", "ColumnFactor", 0.05, 1))
  t <- o$Pooled_ANOVA
  expect_lt(t$`Pr(>F)`[t$Source == "Row x Env"], 0.05)
  expect_true("Row x Env" %in% o$GEI_Note$terms)
})

test_that("random environments: Tukey uses the same denominator as the F-test", {
  data(sugarcane)
  o <- PooledRBD(sugarcane, "ccsy", "env", "gen", "rep", 0.05, 1,
                 env_effect = "random")
  t  <- o$Pooled_ANOVA
  tc <- o$Treatments_Comparison
  expect_identical(t$`Error term`[t$Source == "Treatment"], "Treatment x Env")
  expect_identical(tc$parameters$error_stratum, "Treatment x Env")
  expect_equal(tc$statistics$MSerror, t$`Mean Sq`[t$Source == "Treatment x Env"])
  expect_equal(tc$statistics$Df, t$Df[t$Source == "Treatment x Env"])
  ## critical difference from first principles: 6 environments x 2 reps = 12
  expect_equal(tc$statistics$Critical_difference,
               stats::qtukey(0.95, 18, 85) *
                 sqrt(t$`Mean Sq`[t$Source == "Treatment x Env"] / 12))
  ## fixed environments are unchanged: pooled error
  f <- PooledRBD(sugarcane, "ccsy", "env", "gen", "rep", 0.05, 1)
  expect_identical(f$Treatments_Comparison$parameters$error_stratum, "Pooled error")
})

test_that("random environments in split plot: interaction and comparisons", {
  data(spd_example)
  w <- testthat::capture_warnings(
    o <- PooledSPD(spd_example, "Yield", "Location", "Replication", "MainPlot",
                   "SubPlot", 0.05, 1, env_effect = "random"))
  t <- o$Pooled_ANOVA
  expect_identical(t$`Error term`[t$Source == "Main plot"], "Main plot x Env")
  expect_identical(t$`Error term`[t$Source == "Subplot"], "Subplot x Env")
  expect_identical(t$`Error term`[t$Source == "Main plot x SP"], "MP x SP x Env")
  expect_equal(t$`F value`[t$Source == "Main plot x SP"],
               t$`Mean Sq`[t$Source == "Main plot x SP"] /
                 t$`Mean Sq`[t$Source == "MP x SP x Env"])
  ## Main plot x Env has 1 df with two environments: Tukey cannot be computed,
  ## so the comparison is skipped with an informative warning, not an error
  expect_true(any(grepl("No pooled comparison for MainPlot", w)))
  tc <- o$Treatments_Comparison
  expect_s3_class(tc, "cane_comparison")
  expect_identical(tc$parameters$error_stratum, "Subplot x Env")
  expect_equal(tc$statistics$MSerror, t$`Mean Sq`[t$Source == "Subplot x Env"])
})

test_that("random environments in strip plot: row x column against three-way term", {
  data(stpd_example)
  o <- suppressWarnings(PooledStPD(stpd_example, "Yield", "Location", "Replication",
                                   "RowFactor", "ColumnFactor", 0.05, 1,
                                   env_effect = "random"))
  t <- o$Pooled_ANOVA
  expect_identical(t$`Error term`[t$Source == "Row x Column"], "Row x Col x Env")
  expect_identical(o$Treatments_Comparison$Interaction$parameters$error_stratum,
                   "Row x Col x Env")
  expect_identical(o$Treatments_Comparison$ColumnFactor$parameters$error_stratum,
                   "Column x Env")
})

test_that("Levene's test uses the total error df as denominator df", {
  data(sugarcane)
  o <- PooledRBD(sugarcane, "ccsy", "env", "gen", "rep", 0.05, 1)
  expect_equal(unname(o$Levene_Test$parameter), c(5, 102))
  ## independent reference: with two replications the residuals occur in
  ## exact +/- pairs, so each environment's residual median is exactly zero
  ## and the classical Brown-Forsythe test applied to the absolute residuals
  ## of ONE replicate (18 values per environment) is the correct test
  e <- do.call(rbind, lapply(split(sugarcane, sugarcane$env), function(s) {
    s <- droplevels(s)
    data.frame(env = s$env, rep = s$rep,
               r = unname(stats::residuals(stats::lm(ccsy ~ rep + gen, data = s))))
  }))
  expect_lt(max(abs(tapply(e$r, e$env, stats::median))), 1e-10)
  e <- e[e$rep == levels(sugarcane$rep)[1], ]
  z <- abs(e$r)
  ref <- stats::anova(stats::lm(z ~ factor(e$env)))
  expect_equal(unname(o$Levene_Test$statistic), ref[1, "F value"])
  expect_equal(o$Levene_Test$p.value, ref[1, "Pr(>F)"])
  ## Latin square: 3 squares of 4 x 4, error df 6 each
  data(lsd_example)
  l <- PooledLSD(lsd_example, "Yield", "Location", "Treatment", "Row", "Column", 0.05, 1)
  expect_equal(unname(l$Levene_Test$parameter), c(2, 18))
})

test_that("Bartlett's test is reported for every error stratum", {
  data(spd_example); data(stpd_example)
  o <- suppressWarnings(PooledSPD(spd_example, "Yield", "Location", "Replication",
                                  "MainPlot", "SubPlot", 0.05, 1))
  os <- o$Bartlett_Test$other_strata
  expect_identical(os$Stratum, "Main plot error")
  ms <- vapply(o$Individual_ANOVA, function(a) a$`Mean Sq`[a$Source == "Main plot error"],
               numeric(1))
  df <- vapply(o$Individual_ANOVA, function(a) a$Df[a$Source == "Main plot error"],
               numeric(1))
  expect_equal(os$statistic, hand_bartlett(ms, df))
  expect_equal(os$df, length(ms) - 1)
  expect_match(o$Transformation$reason, "subplot error mean square")
  s <- suppressWarnings(PooledStPD(stpd_example, "Yield", "Location", "Replication",
                                   "RowFactor", "ColumnFactor", 0.05, 1))
  expect_identical(s$Bartlett_Test$other_strata$Stratum,
                   c("Row strip error", "Column strip error"))
  ## single-stratum designs have none
  data(sugarcane)
  expect_null(PooledRBD(sugarcane, "ccsy", "env", "gen", "rep")$Bartlett_Test$other_strata)
  expect_output(print(o$Bartlett_Test), "Main plot error")
})

test_that("a column subset of an ANOVA table prints without error", {
  data(sugarcane)
  o <- PooledRBD(sugarcane, "ccsy", "env", "gen", "rep", 0.05, 1)
  expect_output(print(o$Pooled_ANOVA[, c("Source", "Error term")]), "Rep within Env")
  expect_output(print(o$Pooled_ANOVA[1:2, ]), "Environment")
})
