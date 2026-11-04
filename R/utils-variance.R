## Homogeneity of error variances -------------------------------------------

#' @keywords internal
#' @noRd
.cane_bartlett <- function(ems, dfe, alpha) {
  ## Bartlett's test computed on the environment-wise ERROR mean squares,
  ## each carrying its own error degrees of freedom. This is the statistic
  ## defined by Bartlett (1937) and reproduced as Equation 3 of the CANE
  ## software note. It is NOT the same as applying stats::bartlett.test() to
  ## the stacked residuals, which uses n_j - 1 in place of the error degrees
  ## of freedom and therefore both under-estimates each error variance and
  ## over-states its precision.
  ok <- is.finite(ems) & ems > 0 & dfe > 0
  if (sum(ok) < 2L)
    .cane_stop("Bartlett's test needs at least two environments with a positive error ",
               "mean square and non-zero error degrees of freedom.")
  if (any(!ok))
    warning("Environment(s) ", paste(names(ems)[!ok], collapse = ", "),
            " have a zero or non-finite error mean square and were excluded from ",
            "Bartlett's test. This usually means the design is saturated in that ",
            "environment (no residual degrees of freedom).", call. = FALSE)
  ems <- ems[ok]; dfe <- dfe[ok]
  k  <- length(ems)
  N  <- sum(dfe)
  sp2 <- sum(dfe * ems) / N
  chi <- N * log(sp2) - sum(dfe * log(ems))
  C   <- 1 + (sum(1 / dfe) - 1 / N) / (3 * (k - 1))
  stat <- chi / C
  p <- stats::pchisq(stat, df = k - 1L, lower.tail = FALSE)
  out <- list(statistic = c("Bartlett's K-squared" = stat),
              parameter = c(df = k - 1L),
              p.value = p,
              method = "Bartlett's test for homogeneity of error variances (on error mean squares)",
              data.name = "environment-wise error mean squares",
              error_mean_squares = ems,
              error_df = dfe,
              pooled_error_variance = sp2,
              correction_factor = C,
              variance_ratio = max(ems) / min(ems),
              alpha = alpha,
              heterogeneous = isTRUE(p < alpha))
  class(out) <- c("cane_bartlett", "htest")
  out
}

#' @keywords internal
#' @noRd
.cane_levene <- function(resid_list, dfe = NULL) {
  ## Median-centred (Brown-Forsythe) form of Levene's test, applied to the
  ## environment-wise residuals. Robust to departures from normality and
  ## therefore a useful companion to Bartlett's test, which is not.
  ##
  ## Degrees of freedom. The residuals of environment j are n_j numbers that
  ## carry only nu_j (the error degrees of freedom of that environment's
  ## ANOVA) independent pieces of information. Treating them as n_j
  ## independent values, as the classical one-way ANOVA on |residual - median|
  ## does, overstates the denominator degrees of freedom (N - l instead of
  ## sum(nu_j)) and inflates F; with two replications the residuals occur in
  ## exact +/- pairs, so every absolute deviation is counted twice. The
  ## within-environment mean square is therefore computed on sum(nu_j)
  ## degrees of freedom. This reduces to the classical Brown-Forsythe test
  ## when the within-environment model contains only a mean (nu_j = n_j - 1).
  ## With two replications each environment's residual median is exactly zero
  ## and the test equals the classical one applied to the absolute residuals
  ## of a single replicate.
  g <- rep(names(resid_list), times = vapply(resid_list, length, integer(1)))
  r <- unlist(resid_list, use.names = FALSE)
  g <- factor(g, levels = names(resid_list))
  if (nlevels(g) < 2L) return(NULL)
  z <- abs(r - stats::ave(r, g, FUN = stats::median))
  fit <- stats::anova(stats::lm(z ~ g))
  d1 <- fit[1L, "Df"]
  d2 <- if (is.null(dfe)) fit[2L, "Df"] else sum(dfe[names(resid_list)])
  stat <- (fit[1L, "Sum Sq"] / d1) / (fit[2L, "Sum Sq"] / d2)
  out <- list(statistic = c(F = stat),
              parameter = c(num.df = d1, denom.df = d2),
              p.value = stats::pf(stat, d1, d2, lower.tail = FALSE),
              method = paste0("Levene's test (median-centred, Brown-Forsythe) on environment-wise ",
                              "residuals, denominator df = total error df"),
              data.name = "environment-wise residuals")
  class(out) <- "htest"
  out
}

#' @keywords internal
#' @noRd
.cane_bartlett_strata <- function(ind, strata, alpha) {
  ## Bartlett's test for every error stratum other than the residual one, in
  ## designs that have more than one (split plot, strip plot). The residual
  ## stratum is tested by the main Bartlett_Test and drives the Aitken
  ## transformation; the other strata are reported here so that heterogeneity
  ## of the main plot or strip errors across environments is not overlooked.
  if (!length(strata)) return(NULL)
  rows <- lapply(names(strata), function(lab) {
    ms <- vapply(ind, function(tab) {
      st <- .cane_stratum(tab, strata[[lab]]); if (is.null(st)) NA_real_ else st$MS
    }, numeric(1))
    df <- vapply(ind, function(tab) {
      st <- .cane_stratum(tab, strata[[lab]]); if (is.null(st)) NA_real_ else st$Df
    }, numeric(1))
    ok <- is.finite(ms) & ms > 0 & is.finite(df) & df > 0
    if (sum(ok) < 2L)
      return(data.frame(Stratum = lab, statistic = NA_real_, df = NA_real_,
                        p.value = NA_real_, variance_ratio = NA_real_,
                        heterogeneous = NA, stringsAsFactors = FALSE))
    b <- suppressWarnings(.cane_bartlett(ms[ok], df[ok], alpha))
    data.frame(Stratum = lab, statistic = unname(b$statistic),
               df = unname(b$parameter), p.value = b$p.value,
               variance_ratio = b$variance_ratio, heterogeneous = b$heterogeneous,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' @keywords internal
#' @noRd
.cane_normality <- function(resid_list) {
  ## Shapiro-Wilk test of normality of the residuals of each environment.
  ## Reported for information: Bartlett's test assumes normality, so a
  ## rejection here is a reason to prefer the Levene result.
  do.call(rbind, lapply(names(resid_list), function(nm) {
    r <- resid_list[[nm]]
    r <- r[is.finite(r)]
    if (length(r) < 3L || length(unique(r)) < 3L)
      return(data.frame(Environment = nm, n = length(r), W = NA_real_, p.value = NA_real_))
    st <- stats::shapiro.test(r)
    data.frame(Environment = nm, n = length(r),
               W = unname(st$statistic), p.value = st$p.value)
  }))
}

#' @keywords internal
#' @noRd
.cane_decide_transform <- function(bart, lev, transform, alpha) {
  switch(transform,
    always = list(apply = TRUE,
                  reason = "transform = \"always\": Aitken's transformation applied regardless of the test outcome."),
    never  = list(apply = FALSE,
                  reason = "transform = \"never\": pooled analysis carried out on the original scale regardless of the test outcome."),
    auto   = {
      ap <- isTRUE(bart$p.value < alpha)
      list(apply = ap,
           reason = paste0("transform = \"auto\": Bartlett's test gave P = ",
                           format.pval(bart$p.value, digits = 4), " against alpha = ", alpha,
                           ", so the error variances were judged ",
                           if (ap) "heterogeneous and Aitken's transformation was applied."
                           else "homogeneous and the original scale was retained.",
                           if (!is.null(lev)) paste0(" Levene's test (robust to non-normality) gave P = ",
                                                     format.pval(lev$p.value, digits = 4), ".") else ""))
    })
}
