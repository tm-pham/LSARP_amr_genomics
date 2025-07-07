# ============================================================================ #
# FUNCTION
# Project: LSARP 
# Author: Thi Mui Pham, tmpham@hsph.harvard.edu
# Title: Function to calculate AAPC from a time trend model
# ============================================================================ #
# Load required packages
library(MASS)
library(ggplot2)
library(dplyr)

# Define a function to fit the model and calculate AAPC
calculate.aapc.function <- function(model) {
  # Extract coefficient and standard error for YEAR
  coef_year <- summary(model)$coefficients["YEAR", "Estimate"]
  se_year <- summary(model)$coefficients["YEAR", "Std. Error"]
  
  # Calculate AAPC and 95% CI
  aapc <- (exp(coef_year) - 1) * 100
  lower_bound <- (exp(coef_year - 1.96 * se_year) - 1) * 100
  upper_bound <- (exp(coef_year + 1.96 * se_year) - 1) * 100
  
  # Return the results as a list
  return(list(AAPC = aapc, lower_CI = lower_bound, upper_CI = upper_bound, model = model))
}