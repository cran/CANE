#' @title Pooled Latin Square Design Analysis
#' @description This function performs pooled analysis of variance (ANOVA) using the Latin square design for multiple locations or years. For more details see Montgomery (2017), Dean et al. (2017)<doi:10.1007/978-3-319-52250-0> and Ruíz et al. (2024)<doi:10.1007/978-3-031-65575-3>.
#' @param data A data frame containing the experimental data.
#' @param Response A numeric variable representing the dependent variable (response).
#' @param Location A factor indicating different locations or years.
#' @param Treatment A factor indicating the different treatments applied.
#' @param Row A factor indicating rows.
#' @param Column A factor indicating columns.
#' @param alpha A numeric value specifying the significance level for Bartlett’s test.
#' @param Mult_Comp_Test An integer specifying the type of multiple comparison test:
#'   \itemize{
#'     \item 1 = Tukey's honestly significant difference (Tukey's HSD) test
#'     \item 2 = Duncan's multiple range test (DMRT)
#'     \item 3 = least significant difference (LSD) test
#'   }
#' @return A list containing the following components:
#'   \itemize{
#'     \item \strong{Individual_ANOVA}: Summary of ANOVA results for each location or year.
#'     \item \strong{Location_wise}: Multiple comparisons of treatments within each location or year.
#'     \item \strong{Bartlett_Test}: Results of Bartlett's test for homogeneity of variances.
#'     \item \strong{Pooled_ANOVA}: Combined (pooled) ANOVA table across all locations or years.
#'     \item \strong{Treatments_Comparison}: Summary of pooled treatment comparisons using the selected multiple comparison test.
#'   }
#' @references Dean A, Voss D, Draguljic D (2017)<doi:10.1007/978-3-319-52250-0>.
#'
#' Montgomery DC (2017). \emph{Design and Analysis of Experiments}. John wiley & sons.
#'
#' Ruíz JS, López OAM, Crossa J (2024)<doi:10.1007/978-3-031-65575-3>.
#' @importFrom agricolae HSD.test duncan.test LSD.test
#' @importFrom dplyr mutate filter
#' @importFrom emmeans emmeans
#' @importFrom stats aov anova lm bartlett.test
#' @importFrom stats as.formula residuals var
#' @importFrom stats ave deviance df.residual formula resid
#' @examples
#' # Creating a sample dataset for Pooled Latin Square Design (LSD)
#' df <- data.frame(
#'   Location = factor(c(rep("Agra", 16), rep("Bihar", 16))),  # Locations
#'   Treatment = factor(c(4, 2, 3, 1, 3, 1, 4, 2, 1, 4, 2, 3, 2, 3, 1, 4,
#'   4, 2, 3, 1, 3, 1, 4, 2, 1, 4, 2, 3, 2, 3, 1, 4)),  # Treatments
#'   Row = factor(rep(1:4, each = 4, times = 2)),  # Row factor
#'   Column = factor(rep(1:4, times = 8)),  # Column factor
#'   Yield = c(29.1, 18.9, 29.4, 5.7, 16.4, 10.2, 21.2, 19.1, 5.4, 38.8, 24, 37, 
#'             24.9, 41.7, 9.5, 28.9, 13, 27, 7, 24, 10, 13, 41, 39, 42, 26, 19, 20, 
#'             20, 35, 8, 38)  # Yield values
#' )
#' 
#' # Running PooledLSD function on the dataset
#' out <- PooledLSD(df, "Yield", "Location", "Treatment", "Row", "Column", 0.05, 1)
#' 
#' # Print results
#' print(out)
#' @export
PooledLSD <- function(data, Response, Location, Treatment, Row, Column, alpha, Mult_Comp_Test) {

  # Convert categorical variables to factors
  data[, Location] <- factor(data[, Location])
  data[, Treatment] <- factor(data[, Treatment])
  data[, Row] <- factor(data[, Row])
  data[, Column] <- factor(data[, Column])

  # Individual Analysis for each Location
  Result1 <- by(data, data[, Location], function(x) aov(as.formula(paste(Response, "~", Row, "+", Column, "+", Treatment)), data = x))
  Result1 <- as.list(Result1)

  # Compute Mean Squared Error (MSE) for each Location
  MSE <- sapply(Result1, function(model) sum(resid(model)^2) / df.residual(model))

  # Assign MSE to data based on Location
  data$MSE <- sapply(data[, Location], function(loc) MSE[as.character(loc)])

  # Perform Bartlett’s Test for Homogeneity of Variance
  residuals_list <- lapply(Result1, residuals)
  residuals_values <- unlist(residuals_list)
  group_labels <- rep(names(Result1), times = sapply(residuals_list, length))
  bartlett_result <- bartlett.test(residuals_values, group_labels)

  # Check if transformation is needed
  if (bartlett_result$p.value < alpha) {
    data$Transformed_Yield <- data[, Response] / sqrt(as.numeric(data$MSE))
    Im2 <- lm(as.formula(paste("Transformed_Yield ~", Location, "+", Row, "+", Column, "+", Treatment, "+", Location, ":", Treatment)), data = data)
  } else {
    Im2 <- lm(as.formula(paste(Response, "~", Location, "+", Row, "+", Column, "+", Treatment, "+", Location, ":", Treatment)), data = data)
  }

  # Perform Combined ANOVA
  anova_result <- anova(Im2)

  # Store Individual ANOVA and Treatment Comparisons for each Location
  comparison_results_Location_wise <- list()
  formatted_Individual_ANOVA <- list()

  for (loc in levels(data[, Location])) {
    sub_data <- data[data[, Location] == loc, ]

    # Individual ANOVA tables
    model <- aov(as.formula(paste(Response, "~", Row, "+", Column, "+", Treatment)), data = sub_data)
    formatted_Individual_ANOVA[[loc]] <- summary(model)

    # Perform multiple comparisons based on user selection
    if (Mult_Comp_Test == 1) {
      comp <- HSD.test(model, Treatment, group = TRUE)
    } else if (Mult_Comp_Test == 2) {
      comp <- duncan.test(model, Treatment, group = TRUE)
    } else if (Mult_Comp_Test == 3) {
      comp <- LSD.test(model, Treatment, group = TRUE)
    } else {
      stop("Invalid test selection! Use 1 for Tukey's HSD, 2 for DMRT, or 3 for LSD test.")
    }

    # Store treatment comparison results
    comparison_results_Location_wise[[loc]] <- list(
      statistics = comp$statistics,
      parameters = comp$parameters,
      means = comp$means,
      comparison = comp$comparison,
      groups = comp$groups
    )
  }

  # Perform overall treatment comparison across Locations
  if (Mult_Comp_Test == 1) {
    Compare_result <- HSD.test(Im2, Treatment, group = TRUE)
  } else if (Mult_Comp_Test == 2) {
    Compare_result <- duncan.test(Im2, Treatment, group = TRUE)
  } else if (Mult_Comp_Test == 3) {
    Compare_result <- LSD.test(Im2, Treatment, group = TRUE)
  } else {
    stop("Invalid test selection! Use 1 for Tukey's HSD, 2 for DMRT, or 3 for LSD test.")
  }

  # Return results as a structured list
  return(list(
    Individual_ANOVA = formatted_Individual_ANOVA,
    Location_wise = comparison_results_Location_wise,
    Bartlett_Test = bartlett_result,
    Pooled_ANOVA = anova_result,
    Treatments_Comparison = Compare_result
  ))
}
