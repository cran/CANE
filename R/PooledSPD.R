#' @title Pooled Analysis of a Series of Split Plot Designs
#'
#' @description
#' Performs the classical pooled analysis of variance for a split plot design
#' (SPD) repeated over several environments, with every effect tested against
#' its proper error term.
#'
#' @details
#' A split plot design has two error strata. Main plot effects belong to the
#' main plot stratum and must be tested against the main plot error,
#' \eqn{Environment \times Replication \times MainPlot}; subplot effects and
#' all interactions involving the subplot factor belong to the subplot stratum
#' and are tested against the subplot (residual) error. Testing a main plot
#' effect against the subplot error, as happens if a single-stratum model is
#' fitted and its residual mean square is used for everything, gives an
#' F-ratio with the wrong denominator and the wrong denominator degrees of
#' freedom.
#'
#' The pooled partition implemented here is
#' \deqn{Environment,\; Rep(Env),\; A,\; A \times Env,\; \mathrm{main\ plot\ error},}
#' \deqn{B,\; B \times Env,\; A \times B,\; A \times B \times Env,\;
#'   \mathrm{subplot\ error},}
#' where \eqn{A} is the main plot factor and \eqn{B} the subplot factor. With
#' \eqn{l} environments, \eqn{r} replications, \eqn{a} main plot levels and
#' \eqn{b} subplot levels the main plot error carries \eqn{l(r-1)(a-1)}
#' degrees of freedom and the subplot error \eqn{la(r-1)(b-1)}.
#'
#' Mean separation is carried out in the stratum appropriate to each factor:
#' main plot means are compared using the main plot error mean square and its
#' degrees of freedom, subplot means using the subplot error. The stratum used
#' is recorded in the \code{parameters} component of each comparison.
#'
#' Homogeneity of error variances is assessed for both strata. The main
#' Bartlett test, and the Aitken transformation when it is applied, use the
#' subplot error mean square of each environment; the main plot error mean
#' squares are tested separately and reported in
#' \code{Bartlett_Test$other_strata}, with a note in the transformation report
#' if they are heterogeneous.
#'
#' At least two replications per environment are needed for a main plot error
#' to exist, and three or more are strongly advised: with \eqn{r = 2} and
#' \eqn{a = 2} the main plot error has only \eqn{l} degrees of freedom and the
#' main plot test has almost no power.
#'
#' @inheritParams PooledCRD
#' @param Replication Character string naming the replication factor.
#' @param MainPlot Character string naming the main plot (whole plot) factor,
#'   that is the factor requiring the larger experimental units.
#' @param SubPlot Character string naming the subplot factor.
#'
#' @return An object of class \code{"cane"}. The structure is as described in
#'   \code{\link{PooledCRD}}, except that \code{Location_wise} and
#'   \code{Treatments_Comparison} are each a list with components
#'   \code{MainPlot} and \code{SubPlot}, computed from the main plot and
#'   subplot error strata respectively.
#'
#' @inherit PooledCRD references
#' @seealso \code{\link{PooledCRD}}, \code{\link{PooledRBD}},
#'   \code{\link{PooledLSD}}, \code{\link{PooledStPD}}
#'
#' @examples
#' data(spd_example)
#' out <- PooledSPD(spd_example, Response = "Yield", Location = "Location",
#'                  Replication = "Replication", MainPlot = "MainPlot",
#'                  SubPlot = "SubPlot", alpha = 0.05, Mult_Comp_Test = 1)
#' out
#' out$Pooled_ANOVA
#' out$Treatments_Comparison$MainPlot$parameters
#'
#' @export
PooledSPD <- function(data, Response, Location, Replication, MainPlot, SubPlot,
                      alpha = 0.05, Mult_Comp_Test = 1,
                      transform = c("auto", "always", "never"),
                      env_effect = c("fixed", "random")) {

  env_effect <- match.arg(env_effect)
  chk <- .cane_check(data, Response,
                     c(Location = Location, Replication = Replication,
                       MainPlot = MainPlot, SubPlot = SubPlot),
                     alpha, Mult_Comp_Test, transform)
  dat <- chk$data
  .cane_check_balance(dat, Location, c(Replication, MainPlot, SubPlot),
                      "replication-by-main-plot-by-subplot")

  nrep <- min(tapply(dat[[Replication]], dat[[Location]],
                     function(z) nlevels(droplevels(z))))
  if (nrep < 2L)
    .cane_stop("A split plot design needs at least two replications per environment to ",
               "provide a main plot error, but the smallest number of replications in an ",
               "environment is ", nrep, ".")
  if (nrep == 2L)
    warning("Only two replications per environment: the main plot error has very few ",
            "degrees of freedom and the main plot F-test will have little power. ",
            "Interpret a non-significant main plot effect with caution.", call. = FALSE)

  design <- list(
    name = "Split plot design (SPD) pooled over environments",
    within_resid = "Subplot error", pooled_resid = "Subplot error",
    other_strata = list(`Main plot error` = c(Replication, MainPlot)),
    randomise = function(sp) {
      sp$MainPlot$denom <- c(Location, MainPlot)
      sp$SubPlot$denom  <- c(Location, SubPlot)
      sp$MxS$denom      <- c(Location, MainPlot, SubPlot)
      sp })

  ## ---- within environment: two strata ----
  within_terms <- c(Replication, MainPlot,
                    .cane_term_label(c(Replication, MainPlot)),
                    SubPlot, .cane_term_label(c(MainPlot, SubPlot)))
  within_spec <- list(
    Replication = list(vars = Replication, denom = c(Replication, MainPlot),
                       label = "Replication"),
    MainPlot    = list(vars = MainPlot, denom = c(Replication, MainPlot),
                       label = "Main plot"),
    MainErr     = list(vars = c(Replication, MainPlot), denom = "Residuals",
                       label = "Main plot error"),
    SubPlot     = list(vars = SubPlot, denom = "Residuals", label = "Subplot"),
    MxS         = list(vars = c(MainPlot, SubPlot), denom = "Residuals",
                       label = "Main plot x SP"))

  ## ---- pooled: two strata ----
  pooled_terms <- c(Location,
                    .cane_term_label(c(Location, Replication)),
                    MainPlot,
                    .cane_term_label(c(Location, MainPlot)),
                    .cane_term_label(c(Location, Replication, MainPlot)),
                    SubPlot,
                    .cane_term_label(c(Location, SubPlot)),
                    .cane_term_label(c(MainPlot, SubPlot)),
                    .cane_term_label(c(Location, MainPlot, SubPlot)))
  pooled_spec <- list(
    Environment = list(vars = Location, denom = c(Location, Replication),
                       label = "Environment"),
    RepInEnv    = list(vars = c(Location, Replication),
                       denom = c(Location, Replication, MainPlot),
                       label = "Rep within Env"),
    MainPlot    = list(vars = MainPlot, denom = c(Location, Replication, MainPlot),
                       label = "Main plot"),
    MxE         = list(vars = c(Location, MainPlot),
                       denom = c(Location, Replication, MainPlot),
                       label = "Main plot x Env"),
    MainErr     = list(vars = c(Location, Replication, MainPlot), denom = "Residuals",
                       label = "Main plot error"),
    SubPlot     = list(vars = SubPlot, denom = "Residuals", label = "Subplot"),
    SxE         = list(vars = c(Location, SubPlot), denom = "Residuals",
                       label = "Subplot x Env"),
    MxS         = list(vars = c(MainPlot, SubPlot), denom = "Residuals",
                       label = "Main plot x SP"),
    MxSxE       = list(vars = c(Location, MainPlot, SubPlot), denom = "Residuals",
                       label = "MP x SP x Env"))

  res <- .cane_engine(dat, Response, Location, chk$alpha, chk$Mult_Comp_Test,
                      chk$transform, env_effect,
                      within_terms, within_spec, pooled_terms, pooled_spec,
                      compare = list(
                        list(name = "MainPlot", factor = MainPlot, spec = "MainPlot",
                             within_error = c(Replication, MainPlot),
                             pooled_error = c(Location, Replication, MainPlot)),
                        list(name = "SubPlot", factor = SubPlot, spec = "SubPlot",
                             within_error = "Residuals",
                             pooled_error = "Residuals")),
                      gei_label = c("Main plot x Env", "Subplot x Env", "MP x SP x Env"),
                      design = design)
  res$Call <- match.call()
  res
}
