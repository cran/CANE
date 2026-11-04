## ----setup, include = FALSE------------------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
options(width = 95)

## --------------------------------------------------------------------------------------------
library(CANE)
chk <- function(label, a, b, tol = 1e-6) {
  ok <- isTRUE(all.equal(as.vector(a), as.vector(b), tolerance = tol))
  cat(sprintf("%-58s %s\n", label, if (ok) "PASS" else "FAIL"))
  invisible(ok)
}
results <- logical(0)

## --------------------------------------------------------------------------------------------
data(spd_example)
o <- PooledSPD(spd_example, "Yield", "Location", "Replication", "MainPlot",
               "SubPlot", 0.05, 1)
tab <- o$Pooled_ANOVA
g <- function(s, c) tab[[c]][tab$Source == s]

ref <- suppressWarnings(summary(aov(
  Yield ~ Location * MainPlot * SubPlot + Error(Location:Replication / MainPlot),
  data = spd_example)))
s1 <- ref[["Error: Location:Replication"]][[1]]
s2 <- ref[["Error: Location:Replication:MainPlot"]][[1]]
s3 <- ref[["Error: Within"]][[1]]
for (z in c("s1", "s2", "s3")) {
  tmp <- get(z); rownames(tmp) <- trimws(rownames(tmp)); assign(z, tmp)
}

results <- c(results,
  chk("Environment SS",            g("Environment", "Sum Sq"),     s1["Location", "Sum Sq"]),
  chk("Rep within Env SS",         g("Rep within Env", "Sum Sq"),  s1["Residuals", "Sum Sq"]),
  chk("Environment F",             g("Environment", "F value"),    s1["Location", "F value"]),
  chk("Main plot SS",              g("Main plot", "Sum Sq"),       s2["MainPlot", "Sum Sq"]),
  chk("Main plot x Env SS",        g("Main plot x Env", "Sum Sq"), s2["Location:MainPlot", "Sum Sq"]),
  chk("Main plot error SS",        g("Main plot error", "Sum Sq"), s2["Residuals", "Sum Sq"]),
  chk("Main plot F vs MP error",   g("Main plot", "F value"),      s2["MainPlot", "F value"]),
  chk("Main plot P",               g("Main plot", "Pr(>F)"),       s2["MainPlot", "Pr(>F)"]),
  chk("Subplot SS",                g("Subplot", "Sum Sq"),         s3["SubPlot", "Sum Sq"]),
  chk("Subplot F vs SP error",     g("Subplot", "F value"),        s3["SubPlot", "F value"]),
  chk("Subplot x Env SS",          g("Subplot x Env", "Sum Sq"),   s3["Location:SubPlot", "Sum Sq"]),
  chk("Main plot x SP SS",         g("Main plot x SP", "Sum Sq"),  s3["MainPlot:SubPlot", "Sum Sq"]),
  chk("MP x SP x Env SS",          g("MP x SP x Env", "Sum Sq"),   s3["Location:MainPlot:SubPlot", "Sum Sq"]),
  chk("Subplot error SS",          g("Subplot error", "Sum Sq"),   s3["Residuals", "Sum Sq"]))

## --------------------------------------------------------------------------------------------
F_right <- g("Main plot", "Mean Sq") / g("Main plot error", "Mean Sq")
F_wrong <- g("Main plot", "Mean Sq") / g("Subplot error", "Mean Sq")
data.frame(
  denominator = c("main plot error (correct)", "subplot error (incorrect)"),
  F = c(F_right, F_wrong),
  df2 = c(g("Main plot error", "Df"), g("Subplot error", "Df")),
  P = c(pf(F_right, g("Main plot", "Df"), g("Main plot error", "Df"), lower.tail = FALSE),
        pf(F_wrong, g("Main plot", "Df"), g("Subplot error", "Df"), lower.tail = FALSE)))

## --------------------------------------------------------------------------------------------
data(stpd_example)
o2 <- PooledStPD(stpd_example, "Yield", "Location", "Replication", "RowFactor",
                 "ColumnFactor", 0.05, 1)
t2 <- o2$Pooled_ANOVA
gg <- function(s, c) t2[[c]][t2$Source == s]
r2 <- suppressWarnings(summary(aov(
  Yield ~ Location * RowFactor * ColumnFactor +
    Error(Location:Replication / (RowFactor + ColumnFactor)), data = stpd_example)))
nm <- names(r2)
B <- r2[[grep("Replication:RowFactor", nm)]][[1]]
C <- r2[[grep("Replication:ColumnFactor", nm)]][[1]]
D <- r2[["Error: Within"]][[1]]
for (z in c("B", "C", "D")) { tmp <- get(z); rownames(tmp) <- trimws(rownames(tmp)); assign(z, tmp) }

results <- c(results,
  chk("Row factor SS",         gg("Row factor", "Sum Sq"),         B["RowFactor", "Sum Sq"]),
  chk("Row factor F",          gg("Row factor", "F value"),        B["RowFactor", "F value"]),
  chk("Row strip error SS",    gg("Row strip error", "Sum Sq"),    B["Residuals", "Sum Sq"]),
  chk("Column factor SS",      gg("Column factor", "Sum Sq"),      C["ColumnFactor", "Sum Sq"]),
  chk("Column factor F",       gg("Column factor", "F value"),     C["ColumnFactor", "F value"]),
  chk("Column strip error SS", gg("Column strip error", "Sum Sq"), C["Residuals", "Sum Sq"]),
  chk("Row x Column SS",       gg("Row x Column", "Sum Sq"),       D["RowFactor:ColumnFactor", "Sum Sq"]),
  chk("Row x Col x Env SS",    gg("Row x Col x Env", "Sum Sq"),    D["Location:RowFactor:ColumnFactor", "Sum Sq"]),
  chk("Intersection error SS", gg("Intersection error", "Sum Sq"), D["Residuals", "Sum Sq"]))

## --------------------------------------------------------------------------------------------
d <- spd_example
y <- d$Yield; E <- d$Location; P <- d$Replication; A <- d$MainPlot; B2 <- d$SubPlot
l <- nlevels(E); r <- nlevels(P); a <- nlevels(A); b <- nlevels(B2); G <- mean(y)
m <- function(...) tapply(y, list(...), mean)
mE <- m(E); mA <- m(A); mB <- m(B2); mEP <- m(E, P); mEA <- m(E, A)
mEB <- m(E, B2); mAB <- m(A, B2); mEPA <- m(E, P, A); mEAB <- m(E, A, B2)
acc <- function(n, f, ...) { gr <- expand.grid(...); n * sum(apply(gr, 1, function(z) f(z)^2)) }

hand <- c(
  "Environment"     = acc(r*a*b, function(z) mE[z[1]] - G, j = 1:l),
  "Rep within Env"  = acc(a*b,   function(z) mEP[z[1], z[2]] - mE[z[1]], j = 1:l, k = 1:r),
  "Main plot"       = acc(l*r*b, function(z) mA[z[1]] - G, i = 1:a),
  "Main plot x Env" = acc(r*b,   function(z) mEA[z[1], z[2]] - mE[z[1]] - mA[z[2]] + G,
                          j = 1:l, i = 1:a),
  "Main plot error" = acc(b, function(z) mEPA[z[1], z[2], z[3]] - mEP[z[1], z[2]] -
                            mEA[z[1], z[3]] + mE[z[1]], j = 1:l, k = 1:r, i = 1:a),
  "Subplot"         = acc(l*r*a, function(z) mB[z[1]] - G, q = 1:b),
  "Subplot x Env"   = acc(r*a,   function(z) mEB[z[1], z[2]] - mE[z[1]] - mB[z[2]] + G,
                          j = 1:l, q = 1:b),
  "Main plot x SP"  = acc(l*r,   function(z) mAB[z[1], z[2]] - mA[z[1]] - mB[z[2]] + G,
                          i = 1:a, q = 1:b),
  "MP x SP x Env"   = acc(r, function(z) mEAB[z[1], z[2], z[3]] - mEA[z[1], z[2]] -
                            mEB[z[1], z[3]] - mAB[z[2], z[3]] + mE[z[1]] + mA[z[2]] +
                            mB[z[3]] - G, j = 1:l, i = 1:a, q = 1:b))
for (s in names(hand)) results <- c(results, chk(paste0(s, " SS (by hand)"), g(s, "Sum Sq"), hand[[s]]))
results <- c(results, chk("Subplot error SS (by subtraction)",
  g("Subplot error", "Sum Sq"), sum((y - G)^2) - sum(hand)))

## --------------------------------------------------------------------------------------------
results <- c(results,
  chk("SPD: rep within env df",      g("Rep within Env", "Df"),  l*(r-1)),
  chk("SPD: main plot error df",     g("Main plot error", "Df"), l*(r-1)*(a-1)),
  chk("SPD: subplot error df",       g("Subplot error", "Df"),   l*a*(r-1)*(b-1)),
  chk("SPD: df partition complete",  sum(tab$Df[tab$Source != "Total"]), nrow(spd_example)-1),
  chk("StPD: row strip error df",    gg("Row strip error", "Df"),    2*(4-1)*(2-1)),
  chk("StPD: column strip error df", gg("Column strip error", "Df"), 2*(4-1)*(3-1)),
  chk("StPD: intersection error df", gg("Intersection error", "Df"), 2*(4-1)*(2-1)*(3-1)),
  chk("StPD: df partition complete", sum(t2$Df[t2$Source != "Total"]), nrow(stpd_example)-1))

data(lsd_example)
tl <- PooledLSD(lsd_example, "Yield", "Location", "Treatment", "Row", "Column",
                0.05, 1)$Pooled_ANOVA
results <- c(results,
  chk("LSD: row within env df",      tl$Df[tl$Source == "Row within Env"], 3*(4-1)),
  chk("LSD: column within env df",   tl$Df[tl$Source == "Col within Env"], 3*(4-1)),
  chk("LSD: error df",               tl$Df[tl$Source == "Pooled error"],   3*(4-1)*(4-2)),
  chk("LSD: df partition complete",  sum(tl$Df[tl$Source != "Total"]), nrow(lsd_example)-1))

data(sugarcane)
tr <- PooledRBD(sugarcane, "ccsy", "env", "gen", "rep", 0.05, 1)$Pooled_ANOVA
results <- c(results,
  chk("RBD: rep within env df",      tr$Df[tr$Source == "Rep within Env"],  6*(2-1)),
  chk("RBD: treatment x env df",     tr$Df[tr$Source == "Treatment x Env"], 5*17),
  chk("RBD: error df",               tr$Df[tr$Source == "Pooled error"],    6*(2-1)*17),
  chk("RBD: df partition complete",  sum(tr$Df[tr$Source != "Total"]), nrow(sugarcane)-1))

## ----eval = requireNamespace("nlme", quietly = TRUE)-----------------------------------------
s <- spd_example
s$EP  <- interaction(s$Location, s$Replication, drop = TRUE)
s$EPA <- interaction(s$Location, s$Replication, s$MainPlot, drop = TRUE)
fm <- nlme::lme(Yield ~ Location * MainPlot * SubPlot, random = ~ 1 | EP / EPA,
                data = s, method = "REML")
vc <- suppressWarnings(as.numeric(nlme::VarCorr(fm)[, "Variance"]))
vc <- vc[!is.na(vc)]
sig_rep <- vc[1]; sig_mp <- vc[2]; sig_e <- vc[length(vc)]
cat(sprintf("REML: sigma^2_rep = %.4f   sigma^2_mainplot = %.4f   sigma^2_e = %.4f\n",
            sig_rep, sig_mp, sig_e))
results <- c(results,
  chk("subplot error MS == REML residual variance", g("Subplot error", "Mean Sq"),
      sig_e, tol = 1e-4),
  chk("main plot error MS == sigma^2_e + b sigma^2_mp", g("Main plot error", "Mean Sq"),
      sig_e + b * sig_mp, tol = 1e-3))

## --------------------------------------------------------------------------------------------
ob <- PooledRBD(sugarcane, "ccsy", "env", "gen", "rep", 0.05, 1)
ems <- ob$Transformation$error_mean_squares
dfe <- ob$Transformation$error_df
k <- length(ems); N <- sum(dfe); sp2 <- sum(dfe * ems) / N
chi <- N * log(sp2) - sum(dfe * log(ems))
Cc  <- 1 + (sum(1/dfe) - 1/N) / (3 * (k - 1))
results <- c(results,
  chk("Bartlett K-squared (by hand)", ob$Bartlett_Test$statistic, chi / Cc),
  chk("Bartlett p-value (by hand)",   ob$Bartlett_Test$p.value,
      pchisq(chi / Cc, k - 1, lower.tail = FALSE)),
  chk("pooled error variance",        ob$Bartlett_Test$pooled_error_variance, sp2))

res_list <- lapply(levels(sugarcane$env), function(e) {
  s <- droplevels(sugarcane[sugarcane$env == e, ])
  residuals(lm(ccsy ~ rep + gen, data = s))
})
rv <- unlist(res_list)
gl <- rep(levels(sugarcane$env), times = vapply(res_list, length, 1L))
data.frame(
  method = c("CANE 0.2.0: on error mean squares, error df",
             "CANE 0.1.1: bartlett.test() on stacked residuals"),
  statistic = c(unname(ob$Bartlett_Test$statistic),
                unname(bartlett.test(rv, gl)$statistic)),
  df_per_env = c(unique(dfe), length(res_list[[1]]) - 1))

## --------------------------------------------------------------------------------------------
pe <- ob$Pooled_ANOVA
results <- c(results,
  chk("pooled error MS == 1",  pe$`Mean Sq`[pe$Source == "Pooled error"], 1),
  chk("pooled error SS == df", pe$`Sum Sq`[pe$Source == "Pooled error"],
                               pe$Df[pe$Source == "Pooled error"]))

## --------------------------------------------------------------------------------------------
set.seed(7)
w <- data.frame(E = factor(rep(c("A", "B"), each = 40)),
                T = factor(rep(paste0("T", 1:4), 20)), y = rnorm(80, 10, 1))
oc <- PooledCRD(w, "y", "E", "T", 0.05, 1, transform = "never")
tc <- oc$Treatments_Comparison
MSE <- tc$statistics$MSerror; DFE <- tc$statistics$Df
th <- TukeyHSD(aov(y ~ E + T + E:T, data = w), "T")$T
results <- c(results,
  chk("Tukey critical difference", tc$statistics$Critical_difference,
      qtukey(0.95, 4, DFE) * sqrt(MSE / 20)),
  chk("Tukey p-values == stats::TukeyHSD",
      sort(as.vector(tc$comparison$p.value)), sort(as.vector(th[, "p adj"]))))

ol <- suppressWarnings(PooledCRD(w, "y", "E", "T", 0.05, 3, transform = "never"))
mu <- tapply(w$y, w$T, mean); pr <- combn(names(mu), 2)
tstat <- abs(mu[pr[1, ]] - mu[pr[2, ]]) / sqrt(2 * MSE / 20)
ord <- match(paste(pr[1, ], "-", pr[2, ]), ol$Treatments_Comparison$comparison$Comparison)
results <- c(results,
  chk("LSD critical difference", ol$Treatments_Comparison$statistics$Critical_difference,
      qt(0.975, DFE) * sqrt(2 * MSE / 20)),
  chk("LSD p-values == unadjusted t test",
      ol$Treatments_Comparison$comparison$p.value[ord],
      2 * pt(tstat, DFE, lower.tail = FALSE)),
  chk("Duncan R_2 == Fisher LSD", qtukey(0.95, 2, DFE) * sqrt(MSE / 20),
      qt(0.975, DFE) * sqrt(2 * MSE / 20)))

## ----eval = requireNamespace("multcompView", quietly = TRUE)---------------------------------
tcs <- ob$Treatments_Comparison
muS <- setNames(tcs$means$Mean, tcs$means$Treatment)
prS <- combn(names(muS), 2)
pvS <- ptukey(abs(muS[prS[1, ]] - muS[prS[2, ]]) /
              sqrt(tcs$statistics$MSerror / tcs$parameters$r),
              tcs$parameters$ntr, tcs$statistics$Df, lower.tail = FALSE)
names(pvS) <- paste(prS[1, ], prS[2, ], sep = "-")
refL <- multcompView::multcompLetters(pvS, threshold = 0.05)$Letters
mine <- setNames(tcs$groups$Groups, tcs$groups$Treatment)
rel <- function(lab) {
  n <- names(lab)
  outer(n, n, Vectorize(function(i, j)
    any(strsplit(lab[[i]], "")[[1]] %in% strsplit(lab[[j]], "")[[1]])))
}
results <- c(results,
  chk("grouping relation == multcompView", sum(rel(as.list(mine)) !=
      rel(as.list(refL[names(mine)]))), 0))

## --------------------------------------------------------------------------------------------
cat(sprintf("%d of %d checks passed.\n", sum(results), length(results)))
if (!all(results)) {
  cat("FAILED:\n"); print(which(!results))
}
stopifnot(all(results))

