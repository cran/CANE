## ----setup, include = FALSE------------------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7)
options(width = 95)

## --------------------------------------------------------------------------------------------
library(CANE)

## --------------------------------------------------------------------------------------------
data.frame(
  Design = c("Completely randomised", "Randomised block", "Latin square",
             "Split plot", "Strip plot"),
  Function = c("PooledCRD()", "PooledRBD()", "PooledLSD()", "PooledSPD()",
               "PooledStPD()"),
  `Error strata` = c(1, 1, 1, 2, 3),
  `Design arguments` = c("-", "Replication", "Row, Column",
                         "Replication, MainPlot, SubPlot",
                         "Replication, RowFactor, ColumnFactor"),
  check.names = FALSE)

## --------------------------------------------------------------------------------------------
data(sugarcane)
str(sugarcane)
with(sugarcane, table(env, rep))

## --------------------------------------------------------------------------------------------
fit <- PooledRBD(sugarcane, Response = "ccsy", Location = "env",
                 Treatment = "gen", Replication = "rep",
                 alpha = 0.05, Mult_Comp_Test = 1)

## --------------------------------------------------------------------------------------------
fit$Individual_ANOVA[["SGA-IP"]]

## --------------------------------------------------------------------------------------------
fit$Environment_Key

## --------------------------------------------------------------------------------------------
fit$Bartlett_Test
fit$Levene_Test

## --------------------------------------------------------------------------------------------
fit$Normality_Test

## --------------------------------------------------------------------------------------------
fit$Transformation$applied
cat(fit$Transformation$reason)

## --------------------------------------------------------------------------------------------
cat(fit$Transformation$diagnostic$note)

## --------------------------------------------------------------------------------------------
fit$Pooled_ANOVA

## --------------------------------------------------------------------------------------------
fit$Treatments_Comparison

## --------------------------------------------------------------------------------------------
fit$Individual_ANOVA[["TCSM-IP"]]
fit$Location_wise[["TCSM-IP"]]

## --------------------------------------------------------------------------------------------
data(spd_example)
sp <- PooledSPD(spd_example, "Yield", "Location", "Replication", "MainPlot",
                "SubPlot", 0.05, 1)
sp$Pooled_ANOVA

## --------------------------------------------------------------------------------------------
sp$Treatments_Comparison$MainPlot$statistics
sp$Treatments_Comparison$SubPlot$statistics

## --------------------------------------------------------------------------------------------
data(stpd_example)
st <- PooledStPD(stpd_example, "Yield", "Location", "Replication", "RowFactor",
                 "ColumnFactor", 0.05, 1)
st$Pooled_ANOVA
names(st$Treatments_Comparison)

## --------------------------------------------------------------------------------------------
fr <- PooledRBD(sugarcane, "ccsy", "env", "gen", "rep", 0.05, 1,
                env_effect = "random")
r <- fr$Pooled_ANOVA
r[r$Source %in% c("Treatment", "Treatment x Env"), ]

## --------------------------------------------------------------------------------------------
fr$Treatments_Comparison$parameters[, c("test", "error_stratum", "scale")]

