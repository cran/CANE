#' @title CANE: Comprehensive Groups of Experiments Analysis for Numerous
#'   Environments
#'
#' @description
#' Pooled analysis of a series of designed experiments repeated over locations,
#' seasons or years, where the treatment-by-environment interaction is of
#' primary interest.
#'
#' @details
#' One design-specific function call performs the whole classical workflow:
#' \enumerate{
#'   \item analysis of variance within each environment, with every effect
#'     tested against its proper error stratum;
#'   \item formal assessment of homogeneity of the error variances, by
#'     Bartlett's test computed on the environment-wise error mean squares and
#'     by the median-centred Brown-Forsythe form of Levene's test, which is
#'     robust to non-normality, together with environment-wise Shapiro-Wilk
#'     tests of the normality that Bartlett's test assumes;
#'   \item an optional Aitken transformation, dividing each observation by the
#'     square root of the error mean square of its own environment, applied
#'     automatically, always, or never at the user's choice;
#'   \item pooled analysis of variance partitioning the treatment-by-environment
#'     interaction, with each effect tested against the mean square appropriate
#'     to its stratum, which is reported in the table; and
#'   \item multiple comparison of means by Tukey's HSD, Duncan's multiple range
#'     test or Fisher's LSD, using the error mean square and degrees of freedom
#'     of the stratum in which the factor being compared belongs.
#' }
#'
#' Five designs are supported, one function each:
#' \describe{
#'   \item{\code{\link{PooledCRD}}}{completely randomised design}
#'   \item{\code{\link{PooledRBD}}}{randomised block design}
#'   \item{\code{\link{PooledLSD}}}{Latin square design}
#'   \item{\code{\link{PooledSPD}}}{split plot design}
#'   \item{\code{\link{PooledStPD}}}{strip plot design}
#' }
#'
#' The package implements the fixed-effects pooled analysis of a balanced
#' series of experiments. It is not a substitute for a mixed-model analysis:
#' for unbalanced or partially overlapping treatment sets, for explicit
#' modelling of variance heterogeneity, or for best linear unbiased prediction,
#' use a mixed-model package such as \pkg{nlme} or \pkg{lme4}; for augmented
#' designs see \pkg{augmentedRCBD}; for multiplicative models and stability
#' analysis see \pkg{metan}.
#'
#' @section Reporting the scale of the analysis:
#' When Aitken's transformation is applied, the pooled analysis of variance and
#' the multiple comparison tests are performed on the transformed scale, and
#' the returned objects say so. The environment main effect on that scale is
#' not interpretable in the units of the original response, because each
#' environment has been rescaled by its own error standard deviation. Treatment
#' means are reported on the transformed scale, which is the scale on which the
#' test was conducted and to which the grouping letters refer, alongside the
#' arithmetic means on the original scale for interpretation.
#'
#' @keywords internal
"_PACKAGE"
