## Internal input validation -------------------------------------------------

.cane_stop <- function(...) stop(paste0(...), call. = FALSE)

#' @keywords internal
#' @noRd
.cane_check <- function(data, Response, factors, alpha, Mult_Comp_Test,
                        transform = c("auto", "always", "never")) {

  if (!is.data.frame(data)) .cane_stop("'data' must be a data frame.")
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  if (nrow(data) == 0L) .cane_stop("'data' has no rows.")

  nms <- c(Response = Response, factors)
  for (i in seq_along(nms)) {
    v <- nms[[i]]
    lab <- names(nms)[i]
    if (!is.character(v) || length(v) != 1L || is.na(v))
      .cane_stop("Argument '", lab, "' must be a single column name given as a character string.")
    if (!v %in% names(data))
      .cane_stop("Column '", v, "' (argument '", lab, "') was not found in 'data'. ",
                 "Available columns: ", paste(names(data), collapse = ", "), ".")
  }
  if (anyDuplicated(unlist(nms)))
    .cane_stop("The same column has been supplied for more than one argument: ",
               paste(unique(unlist(nms)[duplicated(unlist(nms))]), collapse = ", "), ".")

  y <- data[[Response]]
  if (!is.numeric(y))
    .cane_stop("The response column '", Response, "' must be numeric, but it is of class '",
               paste(class(y), collapse = "/"), "'.")
  if (all(is.na(y))) .cane_stop("The response column '", Response, "' is entirely missing.")

  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) || alpha <= 0 || alpha >= 1)
    .cane_stop("'alpha' must be a single number strictly between 0 and 1.")

  if (length(Mult_Comp_Test) != 1L || is.na(Mult_Comp_Test) ||
      !(as.character(Mult_Comp_Test) %in% c("1", "2", "3")))
    .cane_stop("'Mult_Comp_Test' must be 1 (Tukey's HSD), 2 (DMRT) or 3 (LSD). ",
               "Received: ", paste(deparse(Mult_Comp_Test), collapse = ""), ".")
  Mult_Comp_Test <- as.integer(as.character(Mult_Comp_Test))

  transform <- match.arg(transform)

  ## drop incomplete cases on the variables actually used
  used <- unique(c(Response, unlist(nms)))
  keep <- stats::complete.cases(data[, used, drop = FALSE])
  ndrop <- sum(!keep)
  if (ndrop > 0L) {
    if (sum(keep) == 0L) .cane_stop("Every row has a missing value in one of the variables used.")
    warning("Dropped ", ndrop, " row(s) with missing values in the analysis variables. ",
            "CANE requires complete, balanced data; please verify that the remaining structure ",
            "is still balanced.", call. = FALSE)
    data <- data[keep, , drop = FALSE]
  }

  ## coerce every classification variable to a factor, dropping unused levels
  for (v in unlist(nms)[-1L]) data[[v]] <- droplevels(factor(data[[v]]))

  Location <- factor(data[[factors[["Location"]]]])
  if (nlevels(Location) < 2L)
    .cane_stop("Pooled analysis requires at least two environments, but column '",
               factors[["Location"]], "' has ", nlevels(Location), " level(s).")

  list(data = data, alpha = alpha, Mult_Comp_Test = Mult_Comp_Test, transform = transform)
}

#' @keywords internal
#' @noRd
.cane_check_balance <- function(data, Location, cross, what) {
  tab <- table(data[[Location]], interaction(data[, cross, drop = FALSE], drop = TRUE))
  if (length(unique(as.vector(tab))) != 1L) {
    warning("The data are not balanced across environments: the number of observations per ",
            what, " combination ranges from ", min(tab), " to ", max(tab),
            ". CANE implements the classical fixed-effects pooled analysis, which assumes a ",
            "balanced structure; sums of squares and F-tests for unbalanced data are ",
            "sequential and should be interpreted with care. Consider a mixed-model ",
            "approach (for example nlme or lme4) for genuinely unbalanced series.",
            call. = FALSE)
    return(FALSE)
  }
  TRUE
}

#' @keywords internal
#' @noRd
.cane_lsd_warning <- function(ntr, Mult_Comp_Test) {
  if (Mult_Comp_Test == 3L && ntr > 3L)
    warning("The least significant difference (LSD) test was requested for ", ntr,
            " treatments, which implies ", choose(ntr, 2),
            " unadjusted pairwise comparisons. Fisher's LSD controls the comparison-wise, ",
            "not the family-wise, error rate: with ", ntr,
            " means the probability of at least one false declaration of significance is ",
            "far above alpha. The LSD option is intended for a small number of pre-planned ",
            "comparisons. For all-pairwise inference use Mult_Comp_Test = 1 (Tukey's HSD).",
            call. = FALSE)
  invisible(NULL)
}
