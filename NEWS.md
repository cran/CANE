# CANE 0.2.1

This release changes the printed output and the coding of the case-study data.
No statistical computation has changed: every analysis of variance, test,
transformation and comparison of means gives the same numbers as in 0.2.0.

## Printed output

* `print()` now shows, for every environment, the environment-wise analysis of
  variance followed by the comparison of treatment means within that
  environment (the `Individual_ANOVA` and `Location_wise` components, which
  were already returned but were not printed). For each comparison the output
  gives the grand mean, the coefficient of variation (CV, %) and the critical
  value of the procedure: the studentized range and Tukey's HSD, the t value
  and Fisher's LSD (critical difference, CD), or Duncan's shortest significant
  ranges R_p for p = 2, ..., t means. In split plot and strip plot designs the
  CV and the critical values are given for each error stratum.

* For Duncan's multiple range test the comparison object gains a
  `critical_ranges` component with the studentized range and the critical
  range for every number of means spanned.

* The printed `INTERPRETATION` paragraph has been removed. The `GEI_Note`
  component is still returned in the object.

## Case-study data

* The 18 genotypes of the `sugarcane` dataset are now coded `G1` to `G18`
  (factor levels in that order). The observations are unchanged.

# CANE 0.2.0

This is a correctness release. Several defects in version 0.1.1 affected the
numerical results, and users who have published analyses produced with 0.1.1
are advised to re-run them. Every correction below is covered by a regression
test in `tests/testthat/`, and the pooled analyses are validated against
independently specified reference analyses in the vignette
`vignette("validation", package = "CANE")`.

## Corrected errors

* **`PooledSPD()` and `PooledStPD()` tested effects against the wrong error
  term.** Both functions fitted a single-stratum linear model and reported the
  F-tests from `anova()`, so every effect was tested against the final
  residual mean square. In a split plot design the main plot effect must be
  tested against the main plot error, and in a strip plot design the row and
  column factors must be tested against the row strip and column strip errors
  respectively. All effects are now tested against the mean square appropriate
  to their stratum, and the denominator used is printed in a new `Error term`
  column of the ANOVA table. Mean separation likewise uses the error mean
  square and degrees of freedom of the stratum in which the compared factor
  belongs.

* **`PooledRBD()` could assign the wrong environment's error mean square.**
  The error mean square was distributed over rows with
  `ave(y, Location, FUN = function(x) MSE[as.character(unique(Location[y == x]))])`.
  Inside `ave()`, `x` holds one group's observations while `y` is the complete
  response, so the comparison recycled a short vector against a long one and
  matched rows in unrelated environments whenever response values coincided.
  On the help example of 0.1.1 this produced the wrong error mean square for
  9 of the 18 rows, together with two `number of items to replace is not a
  multiple of replacement length` warnings. The failure was silent and
  data-dependent: it corrupted the Aitken transformation whenever Bartlett's
  test was significant. The assignment is now a direct lookup by environment
  label. Thanks to the user who reported this.

* **`PooledStPD()` applied the wrong transformation and could misalign every
  row.** The transformed response was computed as
  `unlist(tapply(y, Location, function(Y) Y / sqrt(var(Y))))`. This divided by
  the square root of the *total* within-environment variance of the raw
  response rather than by the square root of the *error* mean square, which is
  the quantity Aitken's transformation is defined by and which the other four
  functions used. In addition, `unlist(tapply(...))` returns values in
  factor-level order rather than in the order of the input rows, so unless the
  data were already sorted by environment the transformed response was
  assigned to the wrong rows. Both faults are fixed, and a regression test now
  checks that results are invariant to the row order of the input.

* **`PooledLSD()` crossed rows and columns with environment instead of nesting
  them.** Each environment carries its own Latin square, so row 1 at one site
  is a different strip of land from row 1 at another. The pooled model now
  nests rows and columns within environment, giving l(t-1) degrees of freedom
  to each instead of t-1. The previous model absorbed 2(l-1)(t-1) degrees of
  freedom of genuine blocking variation, (l-1)(t-1) each for rows and
  columns, into the error term and biased every F-test.

* **Bartlett's test was computed from the stacked residuals rather than from
  the error mean squares.** Version 0.1.1 called
  `bartlett.test(residuals_values, group_labels)`, which estimates each
  environment's variance as the residual sum of squares divided by n_j - 1 and
  uses n_j - 1 degrees of freedom. Both quantities are wrong: the correct
  divisor and degrees of freedom are those of the error term of that
  environment's ANOVA. The statistic is now computed directly from the
  environment-wise error mean squares and their error degrees of freedom, as
  specified by Bartlett (1937).

* **`alpha` did not reach the multiple comparison procedure.** In
  `PooledCRD()`, `PooledRBD()`, `PooledLSD()` and `PooledSPD()` the
  significance level was used for Bartlett's test but not passed to the post
  hoc test, which therefore always ran at 0.05. A single `alpha` now governs
  both.

* **An out-of-range `Mult_Comp_Test` returned `NULL` silently.**
  `PooledStPD()` dispatched with `switch()` on a numeric value, which selects
  by position and ignores the element names, so `Mult_Comp_Test = 5` returned
  `NULL` without any error. All arguments are now validated before any
  computation.

* **The interaction note checked only one interaction term in split plot and
  strip plot designs.** `GEI_Note` looked only at the subplot-by-environment
  (split plot) or row-by-column-by-environment (strip plot) term, so a
  significant main-plot-by-environment, row-by-environment or
  column-by-environment interaction produced pooled groupings with no
  warning. Every interaction of a treatment factor with environment is now
  checked, and the note names the significant terms.

* **With `env_effect = "random"`, mean separation used the pooled error.**
  The F-test of a treatment main effect used its interaction with
  environment, but the pooled comparison of means still used the pooled
  error, so the F-test and the grouping letters rested on different
  denominators. The comparison now uses the same mean square and degrees of
  freedom as the F-test of the factor compared. In split plot and strip plot
  designs the interaction of the two treatment factors is now also tested
  against its three-way interaction with environment when environments are
  random. If the denominator has fewer than two degrees of freedom the Tukey
  and Duncan comparisons cannot be computed; the comparison is then skipped
  with a warning instead of failing.

* **Levene's test overstated its denominator degrees of freedom.** The
  residuals of environment j carry only its error degrees of freedom, not
  n_j - 1; with two replications they occur in exact plus-minus pairs, so
  every absolute deviation was counted twice. The denominator degrees of
  freedom are now the total error degrees of freedom, which reduces to the
  classical Brown-Forsythe test for a mean-only model and, with two
  replications, equals the classical test applied to the absolute residuals
  of a single replicate. For the
  `sugarcane` series the statistic changes from F = 25.02 on 5 and 210 df to
  F = 12.15 on 5 and 102 df (P = 3.0e-09).

* **Printing a column subset of an ANOVA table failed.** For example,
  `out$Pooled_ANOVA[, c("Source", "Error term")]` raised an error. Such
  subsets are now printed as ordinary data frames.

* **The degrees of freedom of the pooled partition are now checked.** The
  returned table includes a `Total` row, and a warning is issued if the
  reported effects do not account for the total degrees of freedom, which
  catches aliased terms arising from unbalanced data.

## New features

* Every effect's denominator is reported in the `Error term` column of the
  ANOVA table, so the stratum used for each F-test is visible rather than
  implicit.

* `transform` argument (`"auto"`, `"always"`, `"never"`) gives the user control
  over the variance-stabilising transformation instead of relying on an
  automatic decision at a significance threshold.

* `env_effect` argument (`"fixed"`, `"random"`) selects the denominator of the
  treatment main effect. With random environments the treatment effect is
  tested against the treatment-by-environment interaction, which is the
  appropriate test when the environments sampled are to represent a wider
  target population.

* For split plot and strip plot designs, Bartlett's test is also reported for
  the main plot error, or the row strip and column strip errors, in
  `Bartlett_Test$other_strata`, and the transformation report states which
  stratum's error mean square the Aitken scaling uses and flags any other
  stratum whose variances are heterogeneous.

* Critical values of Tukey's and Duncan's procedures are computed once per
  comparison and reused, which makes the analysis of the 18-genotype
  `sugarcane` series about seven times faster; results are unchanged.

* Levene's test in its median-centred (Brown-Forsythe) form is reported
  alongside Bartlett's test. It is robust to non-normality, to which Bartlett's
  test is not. Environment-wise Shapiro-Wilk tests of the normality that
  Bartlett's test assumes are reported as well, and the printed output names
  any environment in which normality is rejected.

* The ratio of the largest to the smallest environment-wise error mean square
  is reported, since a ratio above about 3 warrants caution even when the
  formal test is not significant.

* When the treatment-by-environment interaction is significant, the returned
  object carries a `GEI_Note` component advising that the pooled treatment
  means average over inconsistent rankings and that the environment-wise
  comparisons should be consulted as well.

* When Aitken's transformation has been applied, treatment means are reported
  on the transformed scale, which is the scale on which the test was conducted
  and to which the grouping letters refer, together with the arithmetic means
  on the original scale. The scale is recorded in the `parameters` component of
  every comparison.

* Requesting Fisher's LSD for more than three treatments now produces a
  warning that it controls the comparison-wise rather than the family-wise
  error rate and is intended for pre-planned comparisons.

* Input validation: column names, response type, `alpha`, `Mult_Comp_Test`,
  duplicate arguments, number of environments, adequacy of replication for the
  strata the design requires, and balance of the series are all checked, with
  messages that name the offending argument.

* `print()` and `summary()` methods for the returned object, and `print()`
  methods for the ANOVA table, the comparison object and the Bartlett result.

* Five datasets are shipped: `sugarcane` (the case study) plus `crd_example`,
  `lsd_example`, `spd_example` and `stpd_example`. The illustrative datasets
  are simulated from known models and carry four replications, so that the
  separate error strata of the split plot and strip plot designs are
  identifiable and the documented analyses are informative.

* Two vignettes: `vignette("CANE")` for the workflow and
  `vignette("validation")` for the numerical validation of all five designs
  against independently specified reference analyses.

* A test suite with regression tests for every defect listed above.

## Changes to the case-study data

The `sugarcane` dataset now contains the two genuine replications of each
trial. An audit of the working data used for version 0.1.1 found that, in the
`SGA-R`, `TCSM-IIP` and `TCSM-R` trials, the value recorded as the third
replicate was identically the arithmetic mean of the first two replicates, for
every one of the 18 genotypes. A value computed in that way is not an
independently randomised experimental unit. Including it leaves the error sum
of squares of those environments unchanged while doubling the error degrees of
freedom, which halves the error mean square. Because the genotype means are
unchanged, the genotype sum of squares is at the same time multiplied by three
halves, so the genotype F-ratio was inflated by exactly three. Only the
observed replicates are therefore distributed.

## Dependencies

`Imports` is reduced to the base-R packages `stats` and `utils`. The
dependencies on **agricolae**, **dplyr** and **emmeans** are removed.

* **dplyr** and **emmeans** were declared and imported but never called.

* **agricolae** was used for the three multiple comparison procedures, called
  with a fitted model object. In that form the procedure reads the residual
  mean square off the model and therefore cannot be given the error stratum
  of the factor being compared, which is precisely the defect corrected above.
  Tukey's HSD, Duncan's multiple range test and Fisher's LSD are now computed
  within the package from `stats::ptukey()`, `stats::qtukey()` and
  `stats::qt()`, with the error mean square and degrees of freedom supplied
  explicitly for each stratum. The implementations are checked against
  `stats::TukeyHSD()`, against an unadjusted pooled-variance t test, and
  against the identity that Duncan's two-mean critical range equals Fisher's
  LSD; the letter-based grouping is checked against
  `multcompView::multcompLetters()`. These checks are reported in
  `vignette("validation")`.

# CANE 0.1.1

* First CRAN release.
