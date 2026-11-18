#' @title Print and Summary Methods for CANE Objects
#'
#' @description Compact display of the objects returned by the pooled analysis
#'   functions of the package.
#'
#' @param x,object An object returned by \code{\link{PooledCRD}},
#'   \code{\link{PooledRBD}}, \code{\link{PooledLSD}}, \code{\link{PooledSPD}}
#'   or \code{\link{PooledStPD}}.
#' @param digits Number of significant digits used when printing.
#' @param ... Further arguments, currently ignored.
#'
#' @return \code{print} methods return their argument invisibly.
#'   \code{summary.cane} returns, invisibly, a list with the pooled ANOVA
#'   table, the variance-homogeneity tests and the pooled treatment groups.
#'
#' @name cane-methods
NULL

#' @rdname cane-methods
#' @export
print.cane_anova <- function(x, digits = 4, ...) {
  need <- c("Source", "Df", "Sum Sq", "Mean Sq", "F value", "Pr(>F)", "Error term")
  if (!all(need %in% names(x))) {
    ## a subset of the columns (or rows dropped by the user): print as a plain data frame
    print(as.data.frame(unclass(x), check.names = FALSE, stringsAsFactors = FALSE), ...)
    return(invisible(x))
  }
  fmt <- function(v, d) ifelse(is.na(v), "", formatC(signif(v, d), format = "fg",
                                                     digits = d, flag = "#"))
  y <- data.frame(
    Source    = x$Source,
    Df        = x$Df,
    `Sum Sq`  = fmt(x$`Sum Sq`, max(5L, digits)),
    `Mean Sq` = fmt(x$`Mean Sq`, digits),
    `F value` = fmt(x$`F value`, digits),
    `Pr(>F)`  = ifelse(is.na(x$`Pr(>F)`), "",
                       format.pval(x$`Pr(>F)`, digits = 3, eps = 1e-16)),
    Sig = ifelse(is.na(x$`Pr(>F)`), "",
           ifelse(x$`Pr(>F)` < 0.001, "***",
           ifelse(x$`Pr(>F)` < 0.01,  "**",
           ifelse(x$`Pr(>F)` < 0.05,  "*",
           ifelse(x$`Pr(>F)` < 0.1,   ".", " "))))),
    Error = ifelse(is.na(x$`Error term`), "", x$`Error term`),
    check.names = FALSE, stringsAsFactors = FALSE)
  print(y, row.names = FALSE, right = FALSE)
  cat("---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  cat("'Error' names the mean square used as the denominator of each F-test.\n")
  invisible(x)
}

#' @rdname cane-methods
#' @export
print.cane_comparison <- function(x, digits = 4, ...) {
  p <- x$parameters
  s <- x$statistics
  cat("Multiple comparison of means\n")
  cat("  Procedure     : ", p$test, "\n", sep = "")
  cat("  Error stratum : ", p$error_stratum, "  (MSerror = ",
      signif(s$MSerror, digits), " on ", s$Df, " df)\n", sep = "")
  cat("  Scale         : ", p$scale, "\n", sep = "")
  cat("  alpha = ", p$alpha, " | treatments = ", p$ntr,
      " | replication = ", signif(p$r, digits), "\n", sep = "")
  cat("  Grand mean = ", signif(s$Mean, digits), sep = "")
  if (identical(p$scale, "original"))
    cat(" | CV = ", formatC(s$CV, format = "f", digits = 2), " %", sep = "")
  cat("\n")
  .cane_print_critical(x, digits)
  cat("\n")
  g <- x$groups
  num <- vapply(g, is.numeric, logical(1))
  g[num] <- lapply(g[num], function(z) round(z, digits))
  print(g, row.names = FALSE)
  cat("\nMeans sharing a letter do not differ significantly at alpha = ", p$alpha, ".\n", sep = "")
  invisible(x)
}

#' @keywords internal
#' @noRd
.cane_print_critical <- function(x, digits = 4) {
  p <- x$parameters
  s <- x$statistics
  if (grepl("Tukey", p$test)) {
    cat("  Studentized range q(", p$alpha, "; ", p$ntr, ", ", s$Df, ") = ",
        signif(p$critical_value, digits), " | HSD = ",
        signif(s$Critical_difference, digits), "\n", sep = "")
  } else if (grepl("LSD", p$test)) {
    cat("  t(", p$alpha / 2, "; ", s$Df, ") = ", signif(p$critical_value, digits),
        " | LSD (CD) = ", signif(s$Critical_difference, digits), "\n", sep = "")
  } else {
    cr <- x$critical_ranges
    if (is.null(cr)) {
      cat("  Largest critical range = ", signif(s$Critical_difference, digits), "\n", sep = "")
    } else {
      cat("  DMRT shortest significant ranges R_p (p = number of means spanned):\n")
      txt <- paste0("p = ", cr$p, ": ", signif(cr$Critical_range, digits))
      cat(strwrap(paste(txt, collapse = "; "), width = 76, prefix = "    "), sep = "\n")
    }
  }
  invisible(NULL)
}

#' @rdname cane-methods
#' @export
print.cane_bartlett <- function(x, digits = 4, ...) {
  cat("\n\t", x$method, "\n\n", sep = "")
  cat("data:  ", x$data.name, "\n", sep = "")
  cat("Bartlett's K-squared = ", signif(x$statistic, digits),
      ", df = ", x$parameter,
      ", p-value ", if (x$p.value < 2.2e-16) "< 2.2e-16" else
        paste("=", format.pval(x$p.value, digits = digits)), "\n", sep = "")
  cat("ratio of largest to smallest error mean square = ",
      signif(x$variance_ratio, digits), "\n", sep = "")
  os <- x$other_strata
  if (!is.null(os) && nrow(os)) {
    cat("Other error strata (Bartlett's test on each stratum's error mean squares):\n")
    for (k in seq_len(nrow(os)))
      cat("  ", os$Stratum[k], ": K-squared = ", signif(os$statistic[k], digits),
          ", df = ", os$df[k], ", p-value = ",
          format.pval(os$p.value[k], digits = digits),
          ", ratio = ", signif(os$variance_ratio[k], digits), "\n", sep = "")
  }
  invisible(x)
}

#' @rdname cane-methods
#' @export
print.cane_note <- function(x, ...) {
  cat(strwrap(x$message, width = 78, prefix = "  "), sep = "\n")
  invisible(x)
}

#' @rdname cane-methods
#' @export
print.cane <- function(x, digits = 4, ...) {
  line <- strrep("=", 78)
  cat(line, "\n", x$Design, "\n", line, "\n", sep = "")
  cat("CANE ", x$Settings$CANE_version, " | R ", x$Settings$R_version,
      " | alpha = ", x$Settings$alpha,
      " | environments treated as ", x$Settings$env_effect, "\n\n", sep = "")

  cat("ENVIRONMENT KEY AND ERROR MEAN SQUARES\n")
  k <- x$Environment_Key
  k$Error_MS <- signif(k$Error_MS, digits)
  print(k, row.names = FALSE)

  cat("\nENVIRONMENT-WISE ANALYSIS OF VARIANCE AND COMPARISON OF MEANS",
      "  (original scale)\n", sep = "")
  for (e in names(x$Individual_ANOVA)) {
    cat("\n", strrep("-", 78), "\nEnvironment: ", e, "\n", strrep("-", 78), "\n", sep = "")
    print(x$Individual_ANOVA[[e]], digits = digits)
    lw <- x$Location_wise[[e]]
    if (inherits(lw, "cane_comparison")) {
      cat("\n"); print(lw, digits = digits)
    } else if (length(lw)) {
      for (nm in names(lw)) {
        cat("\n-- ", nm, " --\n", sep = ""); print(lw[[nm]], digits = digits)
      }
    } else {
      cat("\nNo comparison of means is available for this environment.\n")
    }
  }

  cat("\nHOMOGENEITY OF ERROR VARIANCES\n")
  print(x$Bartlett_Test, digits = digits)
  if (!is.null(x$Levene_Test))
    cat("Levene (median-centred, robust to non-normality): F = ",
        signif(x$Levene_Test$statistic, digits), " on ",
        x$Levene_Test$parameter[1], " and ", x$Levene_Test$parameter[2],
        " df, p-value = ", format.pval(x$Levene_Test$p.value, digits = digits),
        "\n", sep = "")
  if (!is.null(x$Normality_Test)) {
    nt <- x$Normality_Test
    bad <- nt$Environment[!is.na(nt$p.value) & nt$p.value < 0.05]
    if (length(bad))
      cat("Shapiro-Wilk rejects normality of the residuals in environment(s) ",
          paste(bad, collapse = ", "),
          ". Bartlett's test is sensitive to non-normality; prefer the Levene result there.\n",
          sep = "")
  }

  cat("\nTRANSFORMATION\n")
  cat(strwrap(x$Transformation$reason, width = 78, prefix = "  "), sep = "\n")
  if (!is.null(x$Transformation$diagnostic))
    cat(strwrap(x$Transformation$diagnostic$note, width = 78, prefix = "  "), sep = "\n")

  cat("\nPOOLED ANALYSIS OF VARIANCE  (response on the ",
      x$Transformation$scale, " scale)\n", sep = "")
  print(x$Pooled_ANOVA, digits = digits)

  cat("\nPOOLED COMPARISON OF MEANS\n")
  tc <- x$Treatments_Comparison
  if (inherits(tc, "cane_comparison")) print(tc, digits = digits) else
    for (nm in names(tc)) { cat("\n-- ", nm, " --\n", sep = ""); print(tc[[nm]], digits = digits) }

  cat("\nComponents available: ", paste(setdiff(names(x), "Call"), collapse = ", "),
      "\n", sep = "")
  invisible(x)
}

#' @rdname cane-methods
#' @export
summary.cane <- function(object, digits = 4, ...) {
  print(object, digits = digits)
  tc <- object$Treatments_Comparison
  invisible(list(Pooled_ANOVA = object$Pooled_ANOVA,
                 Bartlett_Test = object$Bartlett_Test,
                 Levene_Test = object$Levene_Test,
                 Groups = if (inherits(tc, "cane_comparison")) tc$groups
                          else lapply(tc, `[[`, "groups")))
}
