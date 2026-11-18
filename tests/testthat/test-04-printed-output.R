## Tests for the changes of version 0.2.1 ----------------------------------

test_that("print() shows the environment-wise ANOVA and comparison for every environment", {
  data(sugarcane, package = "CANE", envir = environment())
  o <- PooledRBD(sugarcane, "ccsy", "env", "gen", "rep", 0.05, 1)
  out <- utils::capture.output(print(o))
  for (e in levels(sugarcane$env))
    expect_true(any(out == paste0("Environment: ", e)))
  expect_identical(sum(grepl("^Multiple comparison of means", out)),
                   nlevels(sugarcane$env) + 1L)        # six environments + pooled
  expect_identical(sum(grepl("CV = ", out, fixed = TRUE)), nlevels(sugarcane$env))
  expect_true(any(grepl("HSD = ", out, fixed = TRUE)))
  expect_false(any(grepl("INTERPRETATION", out, fixed = TRUE)))
  ## the GEI note is still returned in the object
  expect_s3_class(o$GEI_Note, "cane_note")
})

test_that("CV and critical values of the environment-wise comparisons are correct", {
  data(sugarcane, package = "CANE", envir = environment())
  for (k in 1:3) {
    o <- suppressWarnings(PooledRBD(sugarcane, "ccsy", "env", "gen", "rep", 0.05, k))
    for (e in levels(sugarcane$env)) {
      s <- droplevels(sugarcane[sugarcane$env == e, ])
      a <- anova(lm(ccsy ~ rep + gen, data = s))
      mse <- a["Residuals", "Mean Sq"]; dfe <- a["Residuals", "Df"]
      lw <- o$Location_wise[[e]]
      expect_equal(lw$statistics$CV, 100 * sqrt(mse) / mean(s$ccsy))
      cd <- switch(k,
                   qtukey(0.95, 18, dfe) * sqrt(mse / 2),
                   qtukey(1 - (1 - (1 - 0.05)^17), 18, dfe) * sqrt(mse / 2),
                   qt(0.975, dfe) * sqrt(2 * mse / 2))
      expect_equal(lw$statistics$Critical_difference, cd)
    }
  }
})

test_that("Duncan's test reports a critical range for every span", {
  data(lsd_example, package = "CANE", envir = environment())
  o <- PooledLSD(lsd_example, "Yield", "Location", "Treatment", "Row", "Column", 0.05, 2)
  cr <- o$Location_wise[[1]]$critical_ranges
  st <- o$Location_wise[[1]]$statistics
  expect_identical(cr$p, 2:4)
  ap <- 1 - (1 - 0.05)^(cr$p - 1)
  expect_equal(cr$Studentized_range, qtukey(1 - ap, cr$p, st$Df))
  expect_equal(cr$Critical_range, cr$Studentized_range * sqrt(st$MSerror / 4))
  expect_null(PooledLSD(lsd_example, "Yield", "Location", "Treatment", "Row", "Column",
                        0.05, 1)$Location_wise[[1]]$critical_ranges)
})

test_that("sugarcane genotypes are coded G1 to G18 in that order", {
  data(sugarcane, package = "CANE", envir = environment())
  expect_identical(levels(sugarcane$gen), paste0("G", 1:18))
  expect_identical(as.vector(table(sugarcane$gen)), rep(12L, 18))
})
