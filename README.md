
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CANE

<!-- badges: start -->
<!-- badges: end -->

**CANE** (Comprehensive Analysis of groups of experiments for Numerous
Environments) performs the classical pooled analysis of variance of a
balanced series of designed experiments repeated over locations, seasons
or years.

A single **design-specific** function call runs the whole workflow for
one design. Choose the function that matches the design actually
conducted; there is deliberately no omnibus function covering all five.

| Design                       | Function       |
|------------------------------|----------------|
| Completely randomised design | `PooledCRD()`  |
| Randomised block design      | `PooledRBD()`  |
| Latin square design          | `PooledLSD()`  |
| Split plot design            | `PooledSPD()`  |
| Strip plot design            | `PooledStPD()` |

Each call performs:

1.  analysis of variance within each environment, with every effect
    tested against the error term of its own stratum;
2.  assessment of the homogeneity of the error variances: Bartlett’s
    test on the environment-wise error mean squares (for every error
    stratum in split plot and strip plot designs), the median-centred
    Brown-Forsythe form of Levene’s test, environment-wise Shapiro-Wilk
    tests and the ratio of the largest to the smallest error mean
    square;
3.  Aitken’s transformation, applied automatically, always or never
    (`transform = "auto"`, `"always"` or `"never"`);
4.  pooled analysis of variance in which the denominator of every F-test
    is printed in an `Error` column; and
5.  multiple comparison of means (Tukey’s HSD, Duncan’s multiple range
    test or Fisher’s LSD) using the error mean square and degrees of
    freedom of the stratum in which the compared factor belongs.

The environment-wise analysis of variance and comparison of means are
printed for every environment alongside the pooled analysis, so that
environment-specific response is not hidden by the pooled grouping.

## Installation

``` r
# from CRAN
install.packages("CANE")

# development version from GitHub
# install.packages("remotes")
remotes::install_github("Vinayaka11442/CANE")
```

## Example

Commercial cane sugar yield of 18 sugarcane genotypes (coded G1 to G18) in a randomised
block design with two replications in six environments:

``` r
library(CANE)
data(sugarcane)

fit <- PooledRBD(data = sugarcane, Response = "ccsy", Location = "env",
                 Treatment = "gen", Replication = "rep",
                 alpha = 0.05, Mult_Comp_Test = 1)

fit$Environment_Key
#>      Level  n   Error_MS Error_df
#> 1   SGA-IP 36 10.7245624       17
#> 2  SGA-IIP 36 20.5435136       17
#> 3    SGA-R 36 10.6075908       17
#> 4  TCSM-IP 36  0.1605052       17
#> 5 TCSM-IIP 36  1.3576430       17
#> 6   TCSM-R 36  3.3007993       17
fit$Bartlett_Test
#> 
#>  Bartlett's test for homogeneity of error variances (on error mean squares)
#> 
#> data:  environment-wise error mean squares
#> Bartlett's K-squared = 81.17, df = 5, p-value = 4.77e-16
#> ratio of largest to smallest error mean square = 128
fit$Pooled_ANOVA
#>  Source          Df  Sum Sq Mean Sq F value Pr(>F)  Sig Error         
#>  Environment       5 30753. 6151.   1654.   2.5e-09 *** Rep within Env
#>  Rep within Env    6 22.312 3.719   3.719   0.0022  **  Pooled error  
#>  Treatment        17 503.72 29.63   29.63   < 1e-16 *** Pooled error  
#>  Treatment x Env  85 1233.5 14.51   14.51   < 1e-16 *** Pooled error  
#>  Pooled error    102 102.00 1.000                                     
#>  Total           215 32614.                                           
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 'Error' names the mean square used as the denominator of each F-test.
head(fit$Treatments_Comparison$groups)
#>   Treatment     Mean Groups Mean_original_scale
#> 1        G2 13.93197      a            15.89367
#> 2       G10 13.66819      a            15.59433
#> 3       G14 12.84143     ab            15.19867
#> 4        G5 12.63899     ab            15.27950
#> 5        G1 11.83662     bc            12.86358
#> 6       G18 11.80427    bcd            14.23675
```

`print(fit)` shows the complete analysis, including the analysis of
variance and the comparison of means (with the coefficient of variation
and the critical value of the test) within each environment. The error variances are
strongly heterogeneous here, so the pooled analysis is carried out on
the Aitken-transformed scale; the grouping letters refer to that scale,
and the arithmetic means on the original scale are reported alongside
them.

## Validation

Every design is checked against reference analyses specified
independently of the package (stratified `aov()` fits with `Error()`
strata, sums of squares accumulated from cell means, restricted maximum
likelihood in **nlme**, design theory for the degrees of freedom, and
the base-R distribution functions behind the post hoc procedures):

``` r
vignette("validation", package = "CANE")
```

## Scope

CANE implements the fixed-effects pooled analysis of a **balanced**
series. For other situations use:

- unbalanced series or partly overlapping treatment sets: mixed models
  with **lme4**, **nlme** or **sommer**, with **statgenSTA** for the
  single-trial stage;
- augmented designs: **augmentedRCBD** (archived on CRAN in April 2026
  for an administrative reason; installable from its own repository);
- AMMI, GGE biplots and stability statistics: **metan** or
  **ammistability**.

## Version 0.2.0

Version 0.2.0 corrects defects in 0.1.1 that affected numerical results,
including the error strata of the split plot and strip plot designs.
Users who have published analyses produced with 0.1.1 should re-run
them. See the NEWS file (`news(package = "CANE")`) for the full list.

## Citation

``` r
citation("CANE")
```
