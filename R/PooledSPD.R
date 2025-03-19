#' @title Pooled Split Plot Design Analysis
#' @description This function conducts a pooled analysis of variance (ANOVA) using the split plot design (SPD) for data collected across multiple locations or years. In this design, subplot effects are estimated with greater precision. For more details see Dean et al. (2017)<doi:10.1007/978-3-319-52250-0> and Ruíz et al. (2024)<doi:10.1007/978-3-031-65575-3>.
#' @param data A data frame containing the experimental data.
#' @param Response A numeric variable representing the dependent variable (response).
#' @param Location A factor indicating different locations or years.
#' @param Replication A factor indicating replications.
#' @param MainPlot A factor which require larger plot sizes.
#' @param SubPlot A factor which require smaller plot sizes.
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
#' \item \strong{Location_wise}: Multiple comparisons of subplots within each location or year.
#' \item \strong{Bartlett_Test}: Results of Bartlett's test for homogeneity of variances.
#' \item \strong{Pooled_ANOVA}: Combined (pooled) ANOVA table across all locations or years.
#' \item \strong{Treatments_Comparison}: Summary of pooled subplots comparisons using the selected multiple comparison test.
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
#' # Creating a sample dataset for Pooled Split Plot Design (SPD)
#' df <- data.frame(
#'   Location = factor(rep(c("Londan", "Agumbe"), each = 12)),  # Locations
#'   Replication = factor(rep(c(1, 2), each = 6, times = 2)),  # Replications
#'   MainPlot = factor(rep(c(1, 2), each = 3, times = 4)),  # Main plot factor
#'   SubPlot = factor(rep(1:3, times = 8)),  # Sub plot factor
#'   Yield = c(4940, 4810, 5150, 4900, 4920, 5070, 
#'             4830, 5110, 4920, 5020, 5110, 5230,
#'             4964, 4997, 5011, 5102, 4858, 4888, 
#'             5100, 5165, 4965, 5113, 5086, 5176)  # Yield values
#' )
#' 
#' # Running PooledSPD function on the dataset
#' out <- PooledSPD(df, "Yield", "Location", "Replication", "MainPlot", "SubPlot", 0.05, 1)
#' 
#' # Print results
#' print(out)
#' @export
PooledSPD <- function(data, Response, Location, Replication, MainPlot, SubPlot, alpha, Mult_Comp_Test) {

  # Convert categorical variables to factors
  data[, Location] <- factor(data[, Location])
  data[, MainPlot] <- factor(data[, MainPlot])
  data[, SubPlot] <- factor(data[, SubPlot])
  data[, Replication] <- factor(data[, Replication])

  # Individual Analysis for each Location
  Result1 <- by(data, data[, Location], function(x) {
    aov(as.formula(paste(Response, "~", Replication, "+", MainPlot, "+", Replication, ":", MainPlot, "+", SubPlot, "+", MainPlot, ":", SubPlot)), data = x)
  })
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
    Im2 <- lm(as.formula(paste("Transformed_Yield ~", Location, "+", Location, ":", Replication, "+", MainPlot, "+", Location, ":", MainPlot, "+",
                               Location, ":", Replication, ":", MainPlot, "+", SubPlot, "+", Location, ":", SubPlot, "+",
                               MainPlot, ":", SubPlot, "+", Location, ":", MainPlot, ":", SubPlot)), data = data)
  } else {
    Im2 <- lm(as.formula(paste(Response, "~", Location, "+", Location, ":", Replication, "+", MainPlot, "+", Location, ":", MainPlot, "+",
                               Location, ":", Replication, ":", MainPlot, "+", SubPlot, "+", Location, ":", SubPlot, "+",
                               MainPlot, ":", SubPlot, "+", Location, ":", MainPlot, ":", SubPlot)), data = data)
  }

  # Perform Combined ANOVA
  anova_result <- anova(Im2)

  # Store Individual ANOVA and Treatment Comparisons for each Location
  comparison_results_location_wise <- list()
  formatted_Individual_ANOVA <- list()

  for (loc in levels(data[, Location])) {
    sub_data <- data[data[, Location] == loc, ]

    # Individual ANOVA tables
    model <- aov(as.formula(paste(Response, "~", Replication, "+", MainPlot, "+", Replication, "*", MainPlot, "+", SubPlot, "+", MainPlot, "*", SubPlot)), data = sub_data)
    formatted_Individual_ANOVA[[loc]] <- summary(model)

    # Perform multiple comparisons based on user selection
    if (Mult_Comp_Test == 1) {
      comp <- HSD.test(model, SubPlot, group = TRUE)
    } else if (Mult_Comp_Test == 2) {
      comp <- duncan.test(model, SubPlot, group = TRUE)
    } else if (Mult_Comp_Test == 3) {
      comp <- LSD.test(model, SubPlot, group = TRUE)
    } else {
      stop("Invalid test selection! Use 1 for Tukey's HSD, 2 for DMRT, or 3 for LSD test.")
    }

    # Store treatment comparison results
    comparison_results_location_wise[[loc]] <- list(
      statistics = comp$statistics,
      parameters = comp$parameters,
      means = comp$means,
      comparison = comp$comparison,
      groups = comp$groups
    )
  }

  # Perform overall treatment comparison across Locations
  if (Mult_Comp_Test == 1) {
    Compare_result <- HSD.test(Im2, SubPlot, group = TRUE)
  } else if (Mult_Comp_Test == 2) {
    Compare_result <- duncan.test(Im2, SubPlot, group = TRUE)
  } else if (Mult_Comp_Test == 3) {
    Compare_result <- LSD.test(Im2, SubPlot, group = TRUE)
  } else {
    stop("Invalid test selection! Use 1 for Tukey's HSD, 2 for DMRT, or 3 for LSD test.")
  }

  # Return results as a structured list
  return(list(
    Individual_ANOVA = formatted_Individual_ANOVA,
    Location_wise = comparison_results_location_wise,
    Bartlett_Test = bartlett_result,
    Pooled_ANOVA = anova_result,
    Treatments_Comparison = Compare_result
  ))
}
