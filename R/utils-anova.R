## Pooled ANOVA with explicit error strata ----------------------------------

#' @keywords internal
#' @noRd
.cane_key1 <- function(v) paste(sort(v), collapse = "|")

#' @keywords internal
#' @noRd
.cane_key <- function(x) vapply(strsplit(x, ":", fixed = TRUE), .cane_key1, character(1))

#' @keywords internal
#' @noRd
.cane_find <- function(rn, vars) {
  idx <- which(.cane_key(rn) == .cane_key1(vars))
  if (length(idx) == 1L) idx else NA_integer_
}

#' @keywords internal
#' @noRd
.cane_fit <- function(response, terms, data) {
  f <- stats::as.formula(paste0("`", response, "` ~ ", paste(terms, collapse = " + ")))
  stats::lm(f, data = data)
}

#' @keywords internal
#' @noRd
.cane_term_label <- function(vars) paste0("`", vars, "`", collapse = ":")

#' @keywords internal
#' @noRd
.cane_strata_anova <- function(fit, spec, residual_label = "Residual") {
  ## spec: named list; each element has
  ##   vars  : character vector of variables forming the effect
  ##   denom : character vector of variables forming the denominator, or "Residuals"
  ##   label : row label to print
  a  <- stats::anova(fit)
  rn <- gsub("`", "", trimws(rownames(a)), fixed = TRUE)
  resid_i <- which(rn == "Residuals")

  get <- function(vars) if (identical(vars, "Residuals")) resid_i else .cane_find(rn, vars)

  rows <- lapply(spec, function(s) {
    i <- get(s$vars)
    if (is.na(i))
      .cane_stop("Internal error: effect '", paste(s$vars, collapse = ":"),
                 "' is not estimable in the fitted model. This happens when the data are ",
                 "unbalanced or a factor level is absent from some environment. Please check ",
                 "the structure of 'data', or report a reproducible example at ",
                 "https://github.com/Vinayaka11442/CANE/issues")
    j <- get(s$denom)
    Df <- a[i, "Df"]; SS <- a[i, "Sum Sq"]; MS <- SS / Df
    if (is.na(j) || isTRUE(i == j) || a[j, "Df"] < 1L || a[j, "Sum Sq"] <= 0) {
      Fv <- NA_real_; Pv <- NA_real_; elab <- NA_character_
    } else {
      dMS <- a[j, "Sum Sq"] / a[j, "Df"]
      Fv  <- MS / dMS
      Pv  <- stats::pf(Fv, Df, a[j, "Df"], lower.tail = FALSE)
      elab <- if (identical(s$denom, "Residuals")) residual_label else {
        k <- which(vapply(spec, function(z) identical(.cane_key1(z$vars),
                                                      .cane_key1(s$denom)), logical(1)))
        if (length(k)) spec[[k[1L]]]$label else paste(s$denom, collapse = ":")
      }
    }
    data.frame(Source = s$label, Df = Df, `Sum Sq` = SS, `Mean Sq` = MS,
               `F value` = Fv, `Pr(>F)` = Pv, `Error term` = elab,
               Key = .cane_key1(s$vars),
               check.names = FALSE, stringsAsFactors = FALSE)
  })

  tab <- do.call(rbind, rows)
  tab <- rbind(tab, data.frame(
    Source = residual_label, Df = a[resid_i, "Df"], `Sum Sq` = a[resid_i, "Sum Sq"],
    `Mean Sq` = a[resid_i, "Sum Sq"] / a[resid_i, "Df"], `F value` = NA_real_,
    `Pr(>F)` = NA_real_, `Error term` = NA_character_, Key = "Residuals",
    check.names = FALSE, stringsAsFactors = FALSE))
  tot_df <- stats::nobs(fit) - 1L
  eff_df <- sum(tab$Df)
  tab <- rbind(tab, data.frame(
    Source = "Total", Df = tot_df, `Sum Sq` = sum(a[, "Sum Sq"]), `Mean Sq` = NA_real_,
    `F value` = NA_real_, `Pr(>F)` = NA_real_, `Error term` = NA_character_,
    Key = "Total", check.names = FALSE, stringsAsFactors = FALSE))
  if (abs(eff_df - tot_df) > 1e-8)
    warning("The degrees of freedom of the reported effects sum to ", eff_df,
            " but the total is ", tot_df,
            ". A term is aliased or non-estimable, which happens when the data are ",
            "unbalanced or a factor level is empty in some environment.", call. = FALSE)
  rownames(tab) <- NULL
  attr(tab, "anova_raw") <- a
  class(tab) <- c("cane_anova", "data.frame")
  tab
}

#' @keywords internal
#' @noRd
.cane_stratum <- function(tab, vars) {
  k <- if (identical(vars, "Residuals")) "Residuals" else .cane_key1(vars)
  i <- which(tab$Key == k)[1L]
  if (is.na(i)) return(NULL)
  list(MS = tab$`Mean Sq`[i], Df = tab$Df[i], label = tab$Source[i])
}

#' @keywords internal
#' @noRd
.cane_gei_note <- function(tab, gei_label, alpha) {
  ## gei_label: labels of every term involving the environment and a treatment
  ## factor (for example, in a split plot design, main plot x environment,
  ## subplot x environment and the three-way term). The note is issued if any
  ## of them is significant, because any one of them means that the pooled
  ## means of the factor(s) concerned average over inconsistent rankings.
  p <- vapply(gei_label, function(lab) {
    i <- which(tab$Source == lab)
    if (length(i)) tab$`Pr(>F)`[i[1L]] else NA_real_
  }, numeric(1))
  sig <- !is.na(p) & p < alpha
  if (!any(sig)) return(NULL)
  ps <- p[sig]
  fp <- vapply(ps, format.pval, character(1), digits = 4)
  fp <- ifelse(startsWith(fp, "<"), paste0("P ", fp), paste0("P = ", fp))
  terms <- paste0(names(ps), " (", fp, ")")
  msg <- paste0(
    if (length(ps) == 1L) "The term " else "The terms ",
    paste(terms, collapse = ", "),
    if (length(ps) == 1L) " is" else " are",
    " significant at alpha = ", alpha,
    ". Treatment rankings are therefore not consistent across environments, and the pooled ",
    "treatment means in 'Treatments_Comparison' average over that inconsistency. Base any ",
    "recommendation on the environment-wise comparisons in 'Location_wise' as well, and ",
    "consider a stability or multiplicative-model analysis (for example AMMI or GGE) before ",
    "advancing a treatment on its pooled mean alone.")
  structure(list(significant = TRUE, terms = names(ps), p.value = ps, message = msg),
            class = "cane_note")
}
