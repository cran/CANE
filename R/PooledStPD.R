#' @title Pooled Strip Plot Design Analysis
#' @description This function conducts a pooled analysis of variance (ANOVA) using the strip plot design (StPD) for data collected across multiple locations or years. In this design, the interaction between factors (RowFactor and ColumnFactor) is estimated with higher precision. For more details see Dean et al. (2017)<doi:10.1007/978-3-319-52250-0> and Ruíz et al. (2024)<doi:10.1007/978-3-031-65575-3>.
#' @param data A data frame containing the experimental data.
#' @param Response A numeric variable representing the dependent variable (response).
#' @param Location A factor indicating different locations or years.
#' @param Replication A factor indicating replications.
#' @param RowFactor A factor used for horizontal strips.
#' @param ColumnFactor A factor used for vertical strips.
#' @param alpha A numeric value specifying the significance level for Bartlett’s test.
#' @param Mult_Comp_Test An integer specifying the type of multiple comparison test:
#'   \itemize{
#'     \item 1 = Tukey's honestly significant difference (Tukey's HSD) test
#'     \item 2 = Duncan's multiple range test (DMRT)
#'     \item 3 = least significant difference (LSD) test
#'   }
#' @return A list containing the following components:
#' \itemize{
#'   \item \strong{Individual_ANOVA}: Summary of ANOVA results for each location or year.
#'   \item \strong{Location_wise}: Multiple comparisons of interaction of RowFactor and ColumnFactor within each location or year.
#'   \item \strong{Bartlett_Test}: Results of Bartlett's test for homogeneity of variances.
#'   \item \strong{Pooled_ANOVA}: Combined (pooled) ANOVA table across all locations or years.
#'   \item \strong{Interaction_Comparison}: Summary of pooled interaction of RowFactor and ColumnFactor comparisons using the selected multiple comparison test..
#' }
#' @references Dean A, Voss D, Draguljic D (2017)<doi:10.1007/978-3-319-52250-0>.
#'
#' Ruíz JS, López OAM, Crossa J (2024)<doi:10.1007/978-3-031-65575-3>.
#' @importFrom agricolae HSD.test duncan.test LSD.test
#' @importFrom dplyr mutate filter
#' @importFrom emmeans emmeans
#' @importFrom stats aov anova lm bartlett.test
#' @importFrom stats as.formula residuals var
#' @importFrom stats ave deviance df.residual formula resid
#' @examples
#' # Creating a sample dataset for Pooled Strip Plot Design (StPD)
#' df <- data.frame(
#'   Location = factor(rep(c("Londan", "Agumbe"), each = 12)),  # Locations
#'   Replication = factor(rep(c(1, 2), each = 6, times = 2)),  # Replications
#'   RowFactor = factor(rep(c(1, 2), each = 3, times = 4)),  # Row factor
#'   ColumnFactor = factor(rep(1:3, times = 8)),  # Column factor
#'   Yield = c(4940, 4810, 5150, 4900, 4920, 5070, 
#'             4830, 5110, 4920, 5020, 5110, 5230,
#'             4964, 4997, 5011, 5102, 4858, 4888, 
#'             5100, 5165, 4965, 5113, 5086, 5176)  # Yield values
#' )
#' 
#' # Running PooledStPD function on the dataset
#' out <- PooledStPD(df, "Yield", "Location", "Replication", "RowFactor", "ColumnFactor", 0.05, 1)
#' 
#' # Print results
#' print(out)
#' @export
PooledStPD <- function(data, Response, Location, Replication, RowFactor, ColumnFactor, alpha, Mult_Comp_Test) {
  # Ensure categorical variables are factors
  data[[Location]]        <- factor(data[[Location]])
  data[[RowFactor]]   <- factor(data[[RowFactor]])
  data[[ColumnFactor]]<- factor(data[[ColumnFactor]])
  data[[Replication]] <- factor(data[[Replication]])

  # Individual ANOVA by Location or Year
  Result1 <- by(data, data[[Location]], function(subdata) {
    # Within each location: strip-plot ANOVA
    aov(formula(paste(Response, "~", Replication, "+", RowFactor, "+",
                      Replication, ":", RowFactor, "+", ColumnFactor, "+",
                      Replication, ":", ColumnFactor, "+", RowFactor, ":", ColumnFactor)),
        data = subdata)
  })
  Result1 <- as.list(Result1)

  # Bartlett’s Test for homogeneity of variances across locations
  residuals_list  <- lapply(Result1, resid)  # residuals per location
  residuals_values<- unlist(residuals_list)
  group_labels    <- rep(names(Result1), times = sapply(residuals_list, length))
  bartlett_result <- bartlett.test(residuals_values, group_labels)

  # If variances are heterogeneous, apply Aitken’s transformation (scale by sqrt(var) of each group)
  if (bartlett_result$p.value < alpha) {
    data$Transformed_Yield <- unlist(tapply(data[[Response]], data[[Location]],
                                            function(Y) Y / sqrt(var(Y))))
    # Combined model on transformed response (with replication nested in Location)
    Im2 <- lm(formula(paste("Transformed_Yield ~", Location, "+", Location, ":", Replication, "+",
                            RowFactor, "+", Location, ":", RowFactor, "+", Location, ":", Replication, ":", RowFactor, "+",
                            ColumnFactor, "+", Location, ":", ColumnFactor, "+", Location, ":", Replication, ":", ColumnFactor, "+",
                            RowFactor, ":", ColumnFactor, "+", Location, ":", RowFactor, ":", ColumnFactor)),
              data = data)
    response_vec <- data$Transformed_Yield  # for multiple comparisons
  } else {
    # Combined model on original response
    Im2 <- lm(formula(paste(Response, "~", Location, "+", Location, ":", Replication, "+",
                            RowFactor, "+", Location, ":", RowFactor, "+", Location, ":", Replication, ":", RowFactor, "+",
                            ColumnFactor, "+", Location, ":", ColumnFactor, "+", Location, ":", Replication, ":", ColumnFactor, "+",
                            RowFactor, ":", ColumnFactor, "+", Location, ":", RowFactor, ":", ColumnFactor)),
              data = data)
    response_vec <- data[[Response]]
  }

  # Combined (pooled) ANOVA table
  anova_result <- anova(Im2)

  # Prepare interaction factor for treatments
  data$InteractionFactor <- interaction(data[[RowFactor]], data[[ColumnFactor]])

  # Multiple comparisons for RowFactor:ColumnFactor interaction at each location
  comparison_results_location <- list()
  for (loc in names(Result1)) {
    subdata <- data[data[[Location]] == loc, ]
    mod     <- Result1[[loc]]
    # Compute error DF and MS for this location
    dfe    <- df.residual(mod)
    MSe    <- deviance(mod) / dfe
    # Perform the selected multiple comparison test
    comp <- switch(Mult_Comp_Test,
                   `1` = HSD.test(subdata[[Response]], interaction(subdata[[RowFactor]], subdata[[ColumnFactor]]),
                                  DFerror = dfe, MSerror = MSe, alpha = alpha, group = TRUE),
                   `2` = duncan.test(subdata[[Response]], interaction(subdata[[RowFactor]], subdata[[ColumnFactor]]),
                                     DFerror = dfe, MSerror = MSe, alpha = alpha, group = TRUE),
                   `3` = LSD.test(subdata[[Response]], interaction(subdata[[RowFactor]], subdata[[ColumnFactor]]),
                                  DFerror = dfe, MSerror = MSe, alpha = alpha, group = TRUE),
                   stop("Invalid test selection! Use 1 for Tukey HSD, 2 for DMRT, or 3 for LSD."))
    comparison_results_location[[loc]] <- comp
  }

  # Multiple comparison for pooled interaction (across all locations)
  # Use the pooled residual DF and MS from Im2:
  dfe_pooled <- df.residual(Im2)
  MSe_pooled <- deviance(Im2) / dfe_pooled
  Compare_result <- switch(Mult_Comp_Test,
                           `1` = HSD.test(response_vec, data$InteractionFactor,
                                          DFerror = dfe_pooled, MSerror = MSe_pooled, alpha = alpha, group = TRUE),
                           `2` = duncan.test(response_vec, data$InteractionFactor,
                                             DFerror = dfe_pooled, MSerror = MSe_pooled, alpha = alpha, group = TRUE),
                           `3` = LSD.test(response_vec, data$InteractionFactor,
                                          DFerror = dfe_pooled, MSerror = MSe_pooled, alpha = alpha, group = TRUE),
                           stop("Invalid test selection! Use 1 for Tukey's HSD, 2 for DMRT, or 3 for LSD test."))

  # Return results as a list
  return(list(
    Individual_ANOVA     = lapply(Result1, summary),
    Location_wise        = comparison_results_location,
    Bartlett_Test        = bartlett_result,
    Pooled_ANOVA         = anova_result,
    Interaction_Comparison = Compare_result
  ))
}
