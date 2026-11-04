#' @title Sugarcane Multi-Environment Trial (Commercial Cane Sugar Yield)
#'
#' @description
#' Commercial cane sugar yield of 18 sugarcane (\emph{Saccharum} spp. hybrids)
#' genotypes evaluated in a randomised block design at two locations in
#' peninsular India over three crop cycles, giving six environments.
#'
#' @details
#' The trials were conducted by the ICAR-Sugarcane Breeding Institute at
#' Sakthi Sugars, Sivaganga (SGA) and Tiruttani Cooperative Sugar Mills,
#' Thiruvalangadu (TCSM), both in Tamil Nadu, India. The three crop cycles are
#' plant I (IP), plant II (IIP) and ratoon (R). The six environments are the
#' location-by-cycle combinations, labelled \code{SGA-IP}, \code{SGA-IIP},
#' \code{SGA-R}, \code{TCSM-IP}, \code{TCSM-IIP} and \code{TCSM-R}.
#'
#' Two replications are distributed. A third column of values present in an
#' earlier working version of these data was found, on audit, to be identically
#' the arithmetic mean of the two observed replicates in the \code{SGA-R},
#' \code{TCSM-IIP} and \code{TCSM-R} trials, for every one of the 18 genotypes.
#' Such a value is a derived quantity, not an independently randomised
#' experimental unit: including it leaves the error sum of squares of those
#' environments unchanged while doubling the error degrees of freedom, which
#' halves the error mean square; because the genotype means are unchanged, the
#' genotype sum of squares is at the same time multiplied by three halves, so
#' the genotype F-ratio is inflated by exactly three. Only the two genuine
#' replicates are therefore included, giving a fully balanced series of six
#' two-replicate randomised block experiments.
#'
#' @format A data frame with 216 rows and 6 variables:
#' \describe{
#'   \item{env}{Factor with 6 levels, the environment (location-by-crop-cycle).}
#'   \item{loc}{Factor with 2 levels, \code{SGA} and \code{TCSM}.}
#'   \item{crop}{Factor with 3 levels, \code{IP}, \code{IIP} and \code{R}.}
#'   \item{rep}{Factor with 2 levels, the replication (block) within environment.}
#'   \item{gen}{Factor with 18 levels, the genotype.}
#'   \item{ccsy}{Numeric, commercial cane sugar yield in t/ha.}
#' }
#'
#' @source ICAR-Sugarcane Breeding Institute, Coimbatore, Tamil Nadu, India.
#'   Institute Plan Project P1-15/4-1(v).
#'
#' @examples
#' data(sugarcane)
#' str(sugarcane)
#' with(sugarcane, table(env, rep))
#' round(tapply(sugarcane$ccsy, sugarcane$env, mean), 3)
"sugarcane"

#' @title Illustrative Completely Randomised Design Series
#'
#' @description A balanced series of completely randomised experiments: four
#'   treatments, three environments and five experimental units per treatment
#'   within each environment.
#'
#' @format A data frame with 60 rows and 3 variables:
#' \describe{
#'   \item{Loc}{Factor with 3 levels, the environment.}
#'   \item{Treatment}{Factor with 4 levels.}
#'   \item{Yield}{Numeric response.}
#' }
#'
#' @details Simulated from a known additive model with a treatment main effect,
#'   an environment main effect, no treatment-by-environment interaction and
#'   homogeneous error variances, so that the expected outcome of the analysis
#'   is known. Provided for documentation, examples and regression tests.
#'
#' @examples
#' data(crd_example)
#' str(crd_example)
"crd_example"

#' @title Illustrative Latin Square Design Series
#'
#' @description Three 4 x 4 Latin squares, one per environment, with a
#'   different cyclic arrangement in each environment.
#'
#' @format A data frame with 48 rows and 5 variables:
#' \describe{
#'   \item{Location}{Factor with 3 levels, the environment.}
#'   \item{Row}{Factor with 4 levels, the row of the square.}
#'   \item{Column}{Factor with 4 levels, the column of the square.}
#'   \item{Treatment}{Factor with 4 levels.}
#'   \item{Yield}{Numeric response.}
#' }
#'
#' @details Simulated from a known model containing treatment, environment, row
#'   and column effects and no treatment-by-environment interaction. Rows and
#'   columns are specific to each environment, which is why the pooled model
#'   nests them within environment.
#'
#' @examples
#' data(lsd_example)
#' with(lsd_example, table(Treatment, Row, Location))
"lsd_example"

#' @title Illustrative Split Plot Design Series
#'
#' @description A balanced series of split plot experiments: two environments,
#'   four replications, two main plot levels and three subplot levels.
#'
#' @format A data frame with 48 rows and 5 variables:
#' \describe{
#'   \item{Location}{Factor with 2 levels, the environment.}
#'   \item{Replication}{Factor with 4 levels.}
#'   \item{MainPlot}{Factor with 2 levels, the main plot (whole plot) factor.}
#'   \item{SubPlot}{Factor with 3 levels, the subplot factor.}
#'   \item{Yield}{Numeric response.}
#' }
#'
#' @details Simulated from a known model with a genuine main plot error
#'   component in addition to the subplot error, so that the two error strata
#'   of the design are both identifiable. Four replications are used so that
#'   the main plot error has adequate degrees of freedom.
#'
#' @examples
#' data(spd_example)
#' with(spd_example, table(MainPlot, SubPlot, Location))
"spd_example"

#' @title Illustrative Strip Plot Design Series
#'
#' @description A balanced series of strip plot (split block) experiments: two
#'   environments, four replications, two row strip levels and three column
#'   strip levels.
#'
#' @format A data frame with 48 rows and 5 variables:
#' \describe{
#'   \item{Location}{Factor with 2 levels, the environment.}
#'   \item{Replication}{Factor with 4 levels.}
#'   \item{RowFactor}{Factor with 2 levels, applied in horizontal strips.}
#'   \item{ColumnFactor}{Factor with 3 levels, applied in vertical strips.}
#'   \item{Yield}{Numeric response.}
#' }
#'
#' @details Simulated from a known model containing separate row strip, column
#'   strip and intersection error components, so that all three error strata of
#'   the design are identifiable.
#'
#' @examples
#' data(stpd_example)
#' with(stpd_example, table(RowFactor, ColumnFactor, Location))
"stpd_example"
