test_that("pooled SPD agrees with aov() Error() strata", {
  data(spd_example)
  o <- PooledSPD(spd_example, "Yield", "Location", "Replication", "MainPlot",
                 "SubPlot", 0.05, 1)
  t <- o$Pooled_ANOVA
  g <- function(s, c) t[[c]][t$Source == s]
  ref <- suppressWarnings(summary(stats::aov(
    Yield ~ Location * MainPlot * SubPlot + Error(Location:Replication / MainPlot),
    data = spd_example)))
  s1 <- ref[["Error: Location:Replication"]][[1]]
  s2 <- ref[["Error: Location:Replication:MainPlot"]][[1]]
  s3 <- ref[["Error: Within"]][[1]]
  rownames(s1) <- trimws(rownames(s1)); rownames(s2) <- trimws(rownames(s2))
  rownames(s3) <- trimws(rownames(s3))
  expect_equal(g("Environment", "Sum Sq"),     s1["Location", "Sum Sq"])
  expect_equal(g("Main plot", "Sum Sq"),       s2["MainPlot", "Sum Sq"])
  expect_equal(g("Main plot error", "Sum Sq"), s2["Residuals", "Sum Sq"])
  ## the main plot effect must be tested against the MAIN PLOT error
  expect_equal(g("Main plot", "F value"),      s2["MainPlot", "F value"])
  expect_equal(g("Main plot", "Pr(>F)"),       s2["MainPlot", "Pr(>F)"])
  expect_equal(g("Subplot", "F value"),        s3["SubPlot", "F value"])
  expect_equal(g("Subplot error", "Sum Sq"),   s3["Residuals", "Sum Sq"])
  expect_identical(g("Main plot", "Error term"), "Main plot error")
  expect_identical(g("Subplot", "Error term"),   "Subplot error")
})

test_that("pooled StPD agrees with aov() three-stratum Error()", {
  data(stpd_example)
  o <- PooledStPD(stpd_example, "Yield", "Location", "Replication", "RowFactor",
                  "ColumnFactor", 0.05, 1)
  t <- o$Pooled_ANOVA
  g <- function(s, c) t[[c]][t$Source == s]
  ref <- suppressWarnings(summary(stats::aov(
    Yield ~ Location * RowFactor * ColumnFactor +
      Error(Location:Replication / (RowFactor + ColumnFactor)),
    data = stpd_example)))
  nm <- names(ref)
  B <- ref[[grep("Replication:RowFactor", nm)]][[1]]
  C <- ref[[grep("Replication:ColumnFactor", nm)]][[1]]
  D <- ref[["Error: Within"]][[1]]
  rownames(B) <- trimws(rownames(B)); rownames(C) <- trimws(rownames(C))
  rownames(D) <- trimws(rownames(D))
  expect_equal(g("Row factor", "F value"),         B["RowFactor", "F value"])
  expect_equal(g("Row strip error", "Sum Sq"),     B["Residuals", "Sum Sq"])
  expect_equal(g("Column factor", "F value"),      C["ColumnFactor", "F value"])
  expect_equal(g("Column strip error", "Sum Sq"),  C["Residuals", "Sum Sq"])
  expect_equal(g("Row x Column", "Sum Sq"),        D["RowFactor:ColumnFactor", "Sum Sq"])
  expect_equal(g("Intersection error", "Sum Sq"),  D["Residuals", "Sum Sq"])
  expect_identical(g("Row factor", "Error term"),    "Row strip error")
  expect_identical(g("Column factor", "Error term"), "Column strip error")
})

test_that("degrees of freedom follow design theory and partition the total", {
  data(spd_example); data(stpd_example); data(lsd_example); data(sugarcane)
  sp <- PooledSPD(spd_example, "Yield", "Location", "Replication", "MainPlot",
                  "SubPlot", 0.05, 1)$Pooled_ANOVA
  l <- 2; r <- 4; a <- 2; b <- 3
  expect_equal(sp$Df[sp$Source == "Rep within Env"],  l * (r - 1))
  expect_equal(sp$Df[sp$Source == "Main plot error"], l * (r - 1) * (a - 1))
  expect_equal(sp$Df[sp$Source == "Subplot error"],   l * a * (r - 1) * (b - 1))
  expect_equal(sum(sp$Df[sp$Source != "Total"]), nrow(spd_example) - 1L)

  st <- PooledStPD(stpd_example, "Yield", "Location", "Replication", "RowFactor",
                   "ColumnFactor", 0.05, 1)$Pooled_ANOVA
  expect_equal(st$Df[st$Source == "Row strip error"],    l * (r - 1) * (a - 1))
  expect_equal(st$Df[st$Source == "Column strip error"], l * (r - 1) * (b - 1))
  expect_equal(st$Df[st$Source == "Intersection error"], l * (r - 1) * (a - 1) * (b - 1))
  expect_equal(sum(st$Df[st$Source != "Total"]), nrow(stpd_example) - 1L)

  ls <- PooledLSD(lsd_example, "Yield", "Location", "Treatment", "Row", "Column",
                  0.05, 1)$Pooled_ANOVA
  L <- 3; tt <- 4
  expect_equal(ls$Df[ls$Source == "Row within Env"], L * (tt - 1))
  expect_equal(ls$Df[ls$Source == "Col within Env"], L * (tt - 1))
  expect_equal(ls$Df[ls$Source == "Pooled error"],   L * (tt - 1) * (tt - 2))

  rb <- PooledRBD(sugarcane, "ccsy", "env", "gen", "rep", 0.05, 1)$Pooled_ANOVA
  expect_equal(rb$Df[rb$Source == "Rep within Env"],  6 * (2 - 1))
  expect_equal(rb$Df[rb$Source == "Treatment x Env"], (6 - 1) * (18 - 1))
  expect_equal(rb$Df[rb$Source == "Pooled error"],    6 * (2 - 1) * (18 - 1))
  expect_equal(sum(rb$Df[rb$Source != "Total"]), nrow(sugarcane) - 1L)
})

test_that("Bartlett uses the error mean squares and their error df", {
  data(sugarcane)
  o <- PooledRBD(sugarcane, "ccsy", "env", "gen", "rep", 0.05, 1)
  ems <- o$Transformation$error_mean_squares; dfe <- o$Transformation$error_df
  k <- length(ems); N <- sum(dfe); sp2 <- sum(dfe * ems) / N
  chi <- N * log(sp2) - sum(dfe * log(ems))
  Cc <- 1 + (sum(1 / dfe) - 1 / N) / (3 * (k - 1))
  expect_equal(unname(o$Bartlett_Test$statistic), chi / Cc)
  expect_equal(unname(o$Bartlett_Test$parameter), k - 1L)
  expect_true(all(dfe == 17L))   # error df, not n_j - 1 = 35
})

test_that("Aitken transformation satisfies its defining identity", {
  data(sugarcane)
  o <- PooledRBD(sugarcane, "ccsy", "env", "gen", "rep", 0.05, 1)
  expect_true(o$Transformation$applied)
  t <- o$Pooled_ANOVA
  expect_equal(t$`Mean Sq`[t$Source == "Pooled error"], 1)
  expect_equal(t$`Sum Sq`[t$Source == "Pooled error"],
               as.numeric(t$Df[t$Source == "Pooled error"]))
})

test_that("Tukey and LSD reproduce the base-R reference distributions", {
  set.seed(7)
  w <- data.frame(E = factor(rep(c("A", "B"), each = 40)),
                  T = factor(rep(paste0("T", 1:4), 20)), y = rnorm(80, 10, 1))
  o <- PooledCRD(w, "y", "E", "T", 0.05, 1, transform = "never")
  tc <- o$Treatments_Comparison
  MSE <- tc$statistics$MSerror; DFE <- tc$statistics$Df
  expect_equal(tc$statistics$Critical_difference,
               stats::qtukey(0.95, 4, DFE) * sqrt(MSE / 20))
  fit <- stats::aov(y ~ E + T + E:T, data = w)
  th <- stats::TukeyHSD(fit, "T")$T
  expect_equal(sort(unname(tc$comparison$p.value)), sort(unname(th[, "p adj"])))

  ol <- suppressWarnings(PooledCRD(w, "y", "E", "T", 0.05, 3, transform = "never"))
  expect_equal(ol$Treatments_Comparison$statistics$Critical_difference,
               stats::qt(0.975, DFE) * sqrt(2 * MSE / 20))
})

test_that("Duncan's two-mean critical range equals Fisher's LSD", {
  set.seed(11)
  w <- data.frame(E = factor(rep(c("A", "B"), each = 30)),
                  T = factor(rep(c("T1", "T2", "T3"), 20)), y = rnorm(60))
  od <- PooledCRD(w, "y", "E", "T", 0.05, 2, transform = "never")
  MSE <- od$Treatments_Comparison$statistics$MSerror
  DFE <- od$Treatments_Comparison$statistics$Df
  R2  <- stats::qtukey(0.95, 2, DFE) * sqrt(MSE / 20)
  expect_equal(R2, stats::qt(0.975, DFE) * sqrt(2 * MSE / 20))
})
