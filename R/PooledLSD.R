#' @title Pooled Analysis of a Series of Latin Square Designs
#'
#' @description
#' Performs the classical pooled analysis of variance for a Latin square
#' design repeated over several environments, with rows and columns correctly
#' nested within environment.
#'
#' @details
#' The model fitted within environment \eqn{j} is
#' \deqn{Y_{ijkm} = \mu_j + \rho_{jk} + \gamma_{jm} + \tau_{ij} + \epsilon_{ijkm}}
#' with \eqn{\rho_{jk}} the row effect and \eqn{\gamma_{jm}} the column effect.
#' In the pooled model rows and columns are \strong{nested within environment},
#' because each environment carries its own Latin square: row 1 at one site is
#' a different strip of land from row 1 at another. Nesting gives
#' \eqn{l(t-1)} degrees of freedom each for rows and columns. Treating rows and
#' columns as crossed with environment instead would absorb
#' \eqn{(l-1)(t-1)} genuine blocking degrees of freedom into the error term
#' and bias every F-test.
#'
#' @inheritParams PooledCRD
#' @param Row Character string naming the row factor of the square.
#' @param Column Character string naming the column factor of the square.
#'
#' @inherit PooledCRD return references
#' @seealso \code{\link{PooledCRD}}, \code{\link{PooledRBD}},
#'   \code{\link{PooledSPD}}, \code{\link{PooledStPD}}
#'
#' @examples
#' data(lsd_example)
#' out <- PooledLSD(lsd_example, Response = "Yield", Location = "Location",
#'                  Treatment = "Treatment", Row = "Row", Column = "Column",
#'                  alpha = 0.05, Mult_Comp_Test = 1)
#' out
#' out$Pooled_ANOVA
#'
#' @export
PooledLSD <- function(data, Response, Location, Treatment, Row, Column,
                      alpha = 0.05, Mult_Comp_Test = 1,
                      transform = c("auto", "always", "never"),
                      env_effect = c("fixed", "random")) {

  env_effect <- match.arg(env_effect)
  chk <- .cane_check(data, Response,
                     c(Location = Location, Treatment = Treatment,
                       Row = Row, Column = Column),
                     alpha, Mult_Comp_Test, transform)
  dat <- chk$data
  .cane_check_balance(dat, Location, c(Row, Column), "row-by-column")

  nt <- nlevels(dat[[Treatment]])
  if (nlevels(dat[[Row]]) != nt || nlevels(dat[[Column]]) != nt)
    warning("A Latin square requires as many rows and columns as treatments, but there are ",
            nt, " treatments, ", nlevels(dat[[Row]]), " rows and ",
            nlevels(dat[[Column]]), " columns. The analysis will proceed, but check that ",
            "the design really is a Latin square.", call. = FALSE)

  design <- list(
    name = "Latin square design pooled over environments",
    within_resid = "Error", pooled_resid = "Pooled error",
    randomise = function(sp) { sp$Treatment$denom <- c(Location, Treatment); sp })

  within_terms <- c(Row, Column, Treatment)
  within_spec <- list(
    Row       = list(vars = Row,       denom = "Residuals", label = "Row"),
    Column    = list(vars = Column,    denom = "Residuals", label = "Column"),
    Treatment = list(vars = Treatment, denom = "Residuals", label = "Treatment"))

  pooled_terms <- c(Location,
                    .cane_term_label(c(Location, Row)),
                    .cane_term_label(c(Location, Column)),
                    Treatment,
                    .cane_term_label(c(Location, Treatment)))
  pooled_spec <- list(
    Environment = list(vars = Location, denom = "Residuals", label = "Environment"),
    RowInEnv    = list(vars = c(Location, Row), denom = "Residuals",
                       label = "Row within Env"),
    ColInEnv    = list(vars = c(Location, Column), denom = "Residuals",
                       label = "Col within Env"),
    Treatment   = list(vars = Treatment, denom = "Residuals", label = "Treatment"),
    TxE         = list(vars = c(Location, Treatment), denom = "Residuals",
                       label = "Treatment x Env"))

  res <- .cane_engine(dat, Response, Location, chk$alpha, chk$Mult_Comp_Test,
                      chk$transform, env_effect,
                      within_terms, within_spec, pooled_terms, pooled_spec,
                      compare = list(list(name = "Treatment", factor = Treatment, spec = "Treatment",
                                          within_error = "Residuals",
                                          pooled_error = "Residuals")),
                      gei_label = "Treatment x Env", design = design)
  res$Call <- match.call()
  res
}
