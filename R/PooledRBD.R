#' @title Pooled Analysis of a Series of Randomised Block Designs
#'
#' @description
#' Performs the classical pooled analysis of variance for a randomised block
#' design (RBD) repeated over several environments. This is the workhorse
#' analysis of multi-environment variety and agronomy trials.
#'
#' @details
#' The model fitted within environment \eqn{j} is
#' \deqn{Y_{ijk} = \mu_j + \rho_{jk} + \tau_{ij} + \epsilon_{ijk}}
#' with \eqn{\rho_{jk}} the effect of replication \eqn{k} in environment
#' \eqn{j}. The pooled model is
#' \deqn{Y_{ijk} = \mu + \beta_j + \rho_{k(j)} + \tau_i + (\tau\beta)_{ij}
#'   + \epsilon_{ijk},}
#' replications being nested within environments, which gives
#' \eqn{l(r-1)} degrees of freedom for replication rather than \eqn{r-1}.
#'
#' Following the standard treatment of a series of block experiments, the
#' environment main effect is tested against the replication-within-environment
#' mean square, because replications are the units within which environments
#' are compared. Treatment and the treatment-by-environment interaction are
#' tested against the pooled error when \code{env_effect = "fixed"}. When
#' \code{env_effect = "random"} the treatment main effect is tested against the
#' treatment-by-environment interaction instead.
#'
#' @inheritParams PooledCRD
#' @param Replication Character string naming the replication (block) factor.
#'
#' @inherit PooledCRD return references
#' @seealso \code{\link{PooledCRD}}, \code{\link{PooledLSD}},
#'   \code{\link{PooledSPD}}, \code{\link{PooledStPD}}
#'
#' @examples
#' data(sugarcane)
#' out <- PooledRBD(sugarcane, Response = "ccsy", Location = "env",
#'                  Treatment = "gen", Replication = "rep",
#'                  alpha = 0.05, Mult_Comp_Test = 1)
#' out
#' out$Pooled_ANOVA
#' head(out$Treatments_Comparison$groups)
#'
#' @export
PooledRBD <- function(data, Response, Location, Treatment, Replication,
                      alpha = 0.05, Mult_Comp_Test = 1,
                      transform = c("auto", "always", "never"),
                      env_effect = c("fixed", "random")) {

  env_effect <- match.arg(env_effect)
  chk <- .cane_check(data, Response,
                     c(Location = Location, Treatment = Treatment,
                       Replication = Replication),
                     alpha, Mult_Comp_Test, transform)
  dat <- chk$data
  .cane_check_balance(dat, Location, c(Treatment, Replication),
                      "treatment-by-replication")

  design <- list(
    name = "Randomised block design (RBD) pooled over environments",
    within_resid = "Error", pooled_resid = "Pooled error",
    randomise = function(sp) { sp$Treatment$denom <- c(Location, Treatment); sp })

  within_terms <- c(Replication, Treatment)
  within_spec <- list(
    Replication = list(vars = Replication, denom = "Residuals", label = "Replication"),
    Treatment   = list(vars = Treatment,   denom = "Residuals", label = "Treatment"))

  pooled_terms <- c(Location,
                    .cane_term_label(c(Location, Replication)),
                    Treatment,
                    .cane_term_label(c(Location, Treatment)))
  pooled_spec <- list(
    Environment = list(vars = Location, denom = c(Location, Replication),
                       label = "Environment"),
    RepInEnv    = list(vars = c(Location, Replication), denom = "Residuals",
                       label = "Rep within Env"),
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
