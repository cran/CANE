#' @title Pooled Analysis of a Series of Strip Plot Designs
#'
#' @description
#' Performs the classical pooled analysis of variance for a strip plot (split
#' block) design repeated over several environments, with three error strata
#' handled explicitly.
#'
#' @details
#' In a strip plot design the two factors are applied in perpendicular strips,
#' which creates \strong{three} error terms: a row strip error, a column strip
#' error and an intersection error. The row factor must be tested against the
#' row strip error \eqn{Environment \times Replication \times Row}, the column
#' factor against the column strip error
#' \eqn{Environment \times Replication \times Column}, and their interaction
#' against the intersection (residual) error. The interaction is therefore
#' estimated with the greatest precision, which is the characteristic merit of
#' the design.
#'
#' With \eqn{l} environments, \eqn{r} replications, \eqn{a} row levels and
#' \eqn{b} column levels the three errors carry \eqn{l(r-1)(a-1)},
#' \eqn{l(r-1)(b-1)} and \eqn{l(r-1)(a-1)(b-1)} degrees of freedom.
#'
#' Homogeneity of error variances is assessed for all three strata: the main
#' Bartlett test uses the intersection error mean squares, and the row strip
#' and column strip error mean squares are tested separately and reported in
#' \code{Bartlett_Test$other_strata}.
#'
#' Aitken's transformation, when applied, divides each observation by the
#' square root of the \emph{intersection error mean square} of its own
#' environment, matching the definition used by the other functions in the
#' package. The lookup is by environment label, so it is unaffected by the row
#' order of \code{data}.
#'
#' @inheritParams PooledCRD
#' @param Replication Character string naming the replication factor.
#' @param RowFactor Character string naming the factor applied in horizontal
#'   strips.
#' @param ColumnFactor Character string naming the factor applied in vertical
#'   strips.
#'
#' @return An object of class \code{"cane"}. The structure is as described in
#'   \code{\link{PooledCRD}}, except that \code{Location_wise} and
#'   \code{Treatments_Comparison} are each a list with components
#'   \code{RowFactor}, \code{ColumnFactor} and \code{Interaction}, computed
#'   from the row strip, column strip and intersection error strata
#'   respectively.
#'
#' @inherit PooledCRD references
#' @seealso \code{\link{PooledCRD}}, \code{\link{PooledRBD}},
#'   \code{\link{PooledLSD}}, \code{\link{PooledSPD}}
#'
#' @examples
#' data(stpd_example)
#' out <- PooledStPD(stpd_example, Response = "Yield", Location = "Location",
#'                   Replication = "Replication", RowFactor = "RowFactor",
#'                   ColumnFactor = "ColumnFactor", alpha = 0.05,
#'                   Mult_Comp_Test = 1)
#' out
#' out$Pooled_ANOVA
#'
#' @export
PooledStPD <- function(data, Response, Location, Replication, RowFactor,
                       ColumnFactor, alpha = 0.05, Mult_Comp_Test = 1,
                       transform = c("auto", "always", "never"),
                       env_effect = c("fixed", "random")) {

  env_effect <- match.arg(env_effect)
  chk <- .cane_check(data, Response,
                     c(Location = Location, Replication = Replication,
                       RowFactor = RowFactor, ColumnFactor = ColumnFactor),
                     alpha, Mult_Comp_Test, transform)
  dat <- chk$data
  .cane_check_balance(dat, Location, c(Replication, RowFactor, ColumnFactor),
                      "replication-by-row-by-column")

  nrep <- min(tapply(dat[[Replication]], dat[[Location]],
                     function(z) nlevels(droplevels(z))))
  if (nrep < 2L)
    .cane_stop("A strip plot design needs at least two replications per environment to ",
               "provide row strip and column strip errors, but the smallest number of ",
               "replications in an environment is ", nrep, ".")
  if (nrep == 2L)
    warning("Only two replications per environment: the row strip and column strip errors ",
            "have very few degrees of freedom and the corresponding F-tests will have ",
            "little power.", call. = FALSE)

  design <- list(
    name = "Strip plot design (StPD) pooled over environments",
    within_resid = "Intersection error", pooled_resid = "Intersection error",
    other_strata = list(`Row strip error` = c(Replication, RowFactor),
                        `Column strip error` = c(Replication, ColumnFactor)),
    randomise = function(sp) {
      sp$RowFactor$denom <- c(Location, RowFactor)
      sp$ColFactor$denom <- c(Location, ColumnFactor)
      sp$RxC$denom       <- c(Location, RowFactor, ColumnFactor)
      sp })

  within_terms <- c(Replication, RowFactor,
                    .cane_term_label(c(Replication, RowFactor)),
                    ColumnFactor,
                    .cane_term_label(c(Replication, ColumnFactor)),
                    .cane_term_label(c(RowFactor, ColumnFactor)))
  within_spec <- list(
    Replication = list(vars = Replication, denom = "Residuals", label = "Replication"),
    RowFactor   = list(vars = RowFactor, denom = c(Replication, RowFactor),
                       label = "Row factor"),
    RowErr      = list(vars = c(Replication, RowFactor), denom = "Residuals",
                       label = "Row strip error"),
    ColFactor   = list(vars = ColumnFactor, denom = c(Replication, ColumnFactor),
                       label = "Column factor"),
    ColErr      = list(vars = c(Replication, ColumnFactor), denom = "Residuals",
                       label = "Column strip error"),
    RxC         = list(vars = c(RowFactor, ColumnFactor), denom = "Residuals",
                       label = "Row x Column"))

  pooled_terms <- c(Location,
                    .cane_term_label(c(Location, Replication)),
                    RowFactor,
                    .cane_term_label(c(Location, RowFactor)),
                    .cane_term_label(c(Location, Replication, RowFactor)),
                    ColumnFactor,
                    .cane_term_label(c(Location, ColumnFactor)),
                    .cane_term_label(c(Location, Replication, ColumnFactor)),
                    .cane_term_label(c(RowFactor, ColumnFactor)),
                    .cane_term_label(c(Location, RowFactor, ColumnFactor)))
  pooled_spec <- list(
    Environment = list(vars = Location, denom = c(Location, Replication),
                       label = "Environment"),
    RepInEnv    = list(vars = c(Location, Replication), denom = "Residuals",
                       label = "Rep within Env"),
    RowFactor   = list(vars = RowFactor,
                       denom = c(Location, Replication, RowFactor),
                       label = "Row factor"),
    RxE         = list(vars = c(Location, RowFactor),
                       denom = c(Location, Replication, RowFactor),
                       label = "Row x Env"),
    RowErr      = list(vars = c(Location, Replication, RowFactor),
                       denom = "Residuals", label = "Row strip error"),
    ColFactor   = list(vars = ColumnFactor,
                       denom = c(Location, Replication, ColumnFactor),
                       label = "Column factor"),
    CxE         = list(vars = c(Location, ColumnFactor),
                       denom = c(Location, Replication, ColumnFactor),
                       label = "Column x Env"),
    ColErr      = list(vars = c(Location, Replication, ColumnFactor),
                       denom = "Residuals", label = "Column strip error"),
    RxC         = list(vars = c(RowFactor, ColumnFactor), denom = "Residuals",
                       label = "Row x Column"),
    RxCxE       = list(vars = c(Location, RowFactor, ColumnFactor),
                       denom = "Residuals", label = "Row x Col x Env"))

  dat[[".cane_interaction"]] <- interaction(dat[[RowFactor]], dat[[ColumnFactor]],
                                            sep = ":", drop = TRUE)

  res <- .cane_engine(dat, Response, Location, chk$alpha, chk$Mult_Comp_Test,
                      chk$transform, env_effect,
                      within_terms, within_spec, pooled_terms, pooled_spec,
                      compare = list(
                        list(name = "RowFactor", factor = RowFactor, spec = "RowFactor",
                             within_error = c(Replication, RowFactor),
                             pooled_error = c(Location, Replication, RowFactor)),
                        list(name = "ColumnFactor", factor = ColumnFactor, spec = "ColFactor",
                             within_error = c(Replication, ColumnFactor),
                             pooled_error = c(Location, Replication, ColumnFactor)),
                        list(name = "Interaction", factor = ".cane_interaction", spec = "RxC",
                             within_error = "Residuals",
                             pooled_error = "Residuals")),
                      gei_label = c("Row x Env", "Column x Env", "Row x Col x Env"),
                      design = design)
  res$Call <- match.call()
  res
}
