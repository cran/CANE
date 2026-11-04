#' @title Pooled Analysis of a Series of Completely Randomised Designs
#'
#' @description
#' Performs the classical pooled analysis of variance for a completely
#' randomised design (CRD) repeated over several environments (locations,
#' seasons or years). The function carries out the environment-wise analysis,
#' tests homogeneity of the error variances, optionally applies Aitken's
#' variance-stabilising transformation, fits the pooled model including the
#' treatment-by-environment interaction, and separates treatment means.
#'
#' @details
#' The model fitted within environment \eqn{j} is
#' \deqn{Y_{ijk} = \mu_j + \tau_{ij} + \epsilon_{ijk}}
#' and the pooled model is
#' \deqn{Y_{ijk} = \mu + \beta_j + \tau_i + (\tau\beta)_{ij} + \epsilon_{ijk},}
#' where \eqn{\beta_j} is the environment effect, \eqn{\tau_i} the treatment
#' effect and \eqn{(\tau\beta)_{ij}} their interaction. A completely randomised
#' design has a single error stratum, so every effect is tested against the
#' pooled error mean square when \code{env_effect = "fixed"}. With
#' \code{env_effect = "random"} the treatment main effect is tested against the
#' treatment-by-environment interaction, which is the appropriate denominator
#' when the environments are regarded as a random sample from a target
#' population.
#'
#' Homogeneity of error variances is tested by Bartlett's statistic computed
#' on the environment-wise error mean squares, each carrying its own error
#' degrees of freedom. Because Bartlett's test is sensitive to non-normality,
#' the median-centred (Brown-Forsythe) form of Levene's test and
#' environment-wise Shapiro-Wilk tests are reported alongside it. Levene's
#' test is computed on the environment-wise residuals; because the residuals of
#' environment \eqn{j} carry only its error degrees of freedom
#' \eqn{\nu_j}{nu_j} (with two replications they occur in exact plus-minus
#' pairs), its denominator degrees of freedom are \eqn{\sum_j \nu_j}{sum nu_j}
#' rather than \eqn{N - l}. The \code{transform} argument lets the user
#' override the automatic decision.
#'
#' @param data A data frame containing the experimental data. Rows with
#'   missing values in the variables used are dropped with a warning.
#' @param Response Character string naming the numeric response variable.
#' @param Location Character string naming the factor that identifies
#'   environments (locations, seasons or years).
#' @param Treatment Character string naming the treatment factor.
#' @param alpha Significance level used for the variance-homogeneity test and
#'   for the multiple comparison procedure. Default \code{0.05}.
#' @param Mult_Comp_Test Integer selecting the multiple comparison procedure:
#'   \code{1} = Tukey's honestly significant difference (HSD) test,
#'   \code{2} = Duncan's multiple range test (DMRT),
#'   \code{3} = Fisher's least significant difference (LSD) test. Option
#'   \code{3} makes no adjustment for multiplicity and is intended for a small
#'   number of pre-planned comparisons; a warning is issued if it is requested
#'   for more than three treatments.
#' @param transform One of \code{"auto"} (default; transform only when the
#'   variances are judged heterogeneous at level \code{alpha}),
#'   \code{"always"} or \code{"never"}.
#' @param env_effect One of \code{"fixed"} (default) or \code{"random"}.
#'   With \code{"random"} each treatment main effect is tested against its
#'   interaction with environment, and in split plot and strip plot designs
#'   the interaction of the two treatment factors is tested against the
#'   three-way interaction with environment. The pooled comparison of means
#'   then uses the same mean square and degrees of freedom as the F-test of
#'   the factor compared. With few environments that interaction has few
#'   degrees of freedom; if it has fewer than two, the Tukey and Duncan
#'   comparisons for that factor are not computed and a warning is issued.
#'
#' @return An object of class \code{"cane"}, a list with components
#'   \describe{
#'     \item{\code{Design}}{Name of the design analysed.}
#'     \item{\code{Individual_ANOVA}}{List of environment-wise ANOVA tables.}
#'     \item{\code{Location_wise}}{Environment-wise multiple comparisons, each
#'       computed from the error mean square of that environment.}
#'     \item{\code{Bartlett_Test}}{Bartlett's test on the error mean squares.
#'       For split plot and strip plot designs, element \code{other_strata}
#'       gives the same test for the main plot or strip error strata.}
#'     \item{\code{Levene_Test}}{Median-centred Levene (Brown-Forsythe) test.}
#'     \item{\code{Normality_Test}}{Environment-wise Shapiro-Wilk tests.}
#'     \item{\code{Transformation}}{Whether Aitken's transformation was applied,
#'       why, the environment-wise error mean squares and the ratio of the
#'       largest to the smallest.}
#'     \item{\code{Pooled_ANOVA}}{Pooled ANOVA table, including the error term
#'       used as the denominator of each F-test.}
#'     \item{\code{Treatments_Comparison}}{Pooled treatment means, groups and
#'       all pairwise differences. When a transformation was applied, means are
#'       given on the transformed scale (on which the test was performed) and
#'       the arithmetic means on the original scale are reported alongside.}
#'     \item{\code{GEI_Note}}{Advice issued when any interaction of a treatment
#'       factor with environment is significant, naming the terms concerned.}
#'     \item{\code{Environment_Key}}{Mapping of environment levels to their
#'       error mean square and error degrees of freedom.}
#'     \item{\code{Settings}}{Arguments and software versions used.}
#'   }
#'
#' @references
#' Aitken, A. C. (1935). On least squares and linear combination of
#' observations. \emph{Proceedings of the Royal Society of Edinburgh},
#' 55, 42--48. \doi{10.1017/S0370164600014346}
#'
#' Bartlett, M. S. (1937). Properties of sufficiency and statistical tests.
#' \emph{Proceedings of the Royal Society of London A}, 160, 268--282.
#' \doi{10.1098/rspa.1937.0109}
#'
#' Cochran, W. G. and Cox, G. M. (1957). \emph{Experimental Designs},
#' 2nd ed. John Wiley and Sons, New York.
#'
#' Dean, A., Voss, D. and Draguljic, D. (2017). \emph{Design and Analysis of
#' Experiments}, 2nd ed. Springer, Cham. \doi{10.1007/978-3-319-52250-0}
#'
#' Gomez, K. A. and Gomez, A. A. (1984). \emph{Statistical Procedures for
#' Agricultural Research}, 2nd ed. John Wiley and Sons, New York.
#'
#' @seealso \code{\link{PooledRBD}}, \code{\link{PooledLSD}},
#'   \code{\link{PooledSPD}}, \code{\link{PooledStPD}}
#'
#' @examples
#' data(crd_example)
#' out <- PooledCRD(crd_example, Response = "Yield", Location = "Loc",
#'                  Treatment = "Treatment", alpha = 0.05, Mult_Comp_Test = 1)
#' out
#' out$Pooled_ANOVA
#'
#' @importFrom stats anova aov as.formula ave complete.cases deviance df.residual
#' @importFrom stats lm median nobs pchisq pf pt ptukey qt qtukey quantile
#' @importFrom stats residuals sd setNames shapiro.test terms var
#' @importFrom utils combn packageVersion
#' @export
PooledCRD <- function(data, Response, Location, Treatment, alpha = 0.05,
                      Mult_Comp_Test = 1, transform = c("auto", "always", "never"),
                      env_effect = c("fixed", "random")) {

  env_effect <- match.arg(env_effect)
  chk <- .cane_check(data, Response,
                     c(Location = Location, Treatment = Treatment),
                     alpha, Mult_Comp_Test, transform)
  dat <- chk$data
  .cane_check_balance(dat, Location, Treatment, "treatment")

  design <- list(
    name = "Completely randomised design (CRD) pooled over environments",
    within_resid = "Error", pooled_resid = "Pooled error",
    randomise = function(sp) { sp$Treatment$denom <- c(Location, Treatment)
                               sp$Treatment$label <- "Treatment"; sp })

  within_terms <- Treatment
  within_spec <- list(Treatment = list(vars = Treatment, denom = "Residuals",
                                       label = "Treatment"))
  pooled_terms <- c(Location, Treatment, .cane_term_label(c(Location, Treatment)))
  pooled_spec <- list(
    Environment = list(vars = Location, denom = "Residuals", label = "Environment"),
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
