#!/usr/bin/env Rscript

# Model selection analysis that conclude the SecSSE pipeline.
# This script compares the fit of three SecSSE models (Constant Rates, Examined Trait-Dependent, and Concealed Trait-Dependent) 
# to determine which hypothesis best explains the evolutionary dynamics, using the Akaike Information Criterion (AIC).
# It requires three RDS objects: one from the CR model, one from the ETD model and one for the CTD model. 
# Any of these .rds files has to refer to the same taxonomic group for the comparison.
# Each RDS file must contain:
#   model
#   log_likelihood
#   number_of_parameters


# Define the input files in .rds format.
cr_file <- "CR_FIT.RDS"  # From 'SecSSE_CR_model_fit.R'
etd_file <- "ETD_FIT.RDS"  # From 'SecSSE_ETD_model_fit.R'
ctd_file <- "CTD_FIT.RDS"  # From 'SecSSE_CTD_model_fit.R'

# Define the output CSV file that will store the final comparative table (change accordingly).
output_file <- "MODEL_COMPARISON.CSV"

# Read the complete results lists exported by the previous model-fitting scripts.
cr_results <- readRDS(cr_file)
etd_results <- readRDS(etd_file)
ctd_results <- readRDS(ctd_file)

# Create the foundational dataframe by extracting the model name, Maximum Log-Likelihood, and the number of estimated parameters for each model.
comparison <- data.frame(
  model = c(
    cr_results$model,
    etd_results$model,
    ctd_results$model
  ),
  log_likelihood = c(
    cr_results$log_likelihood,
    etd_results$log_likelihood,
    ctd_results$log_likelihood
  ),
  number_of_parameters = c(
    cr_results$number_of_parameters,
    etd_results$number_of_parameters,
    ctd_results$number_of_parameters
  ),
  stringsAsFactors = FALSE)

# Calculate the Akaike Information Criterion (AIC) for each model.
# Formula: -2 * logLikelihood + 2 * number_of_parameters.
# Lower AIC indicates a better trade-off between model fit and complexity.
comparison$AIC <- (
  -2 * comparison$log_likelihood +
    2 * comparison$number_of_parameters)

# Calculate Delta AIC (the difference between each model's AIC and the lowest AIC).
# The best model will always have a delta_AIC of 0.
comparison$delta_AIC <- (
  comparison$AIC -
    min(comparison$AIC))

# Calculate the relative likelihood of each model.
# This represents the proportional evidence for each model compared to the best one.
relative_likelihood <- exp( -0.5 * comparison$delta_AIC)

# Calculate Akaike weights.
# These values normalize the relative likelihoods to sum to 1. 
# A weight of 0.85 means there is an 85% chance this is the best model among the set.
comparison$Akaike_weight <- (
  relative_likelihood /
    sum(relative_likelihood))

# Sort the rows of the dataframe from the lowest (best) to highest (worst) AIC score.
comparison <- comparison[
  order(comparison$AIC),
]

# Reset the row names sequentially (1, 2, 3...) after sorting.
row.names(comparison) <- NULL

# Add a ranking column indicating the position of each model (1 = best).
comparison$rank <- seq_len(nrow(comparison))

# Create a boolean column (TRUE/FALSE) to explicitly flag the best-fitting model.
comparison$best_model <- (comparison$rank == 1)

# Assign a qualitative evidence category based on standard Burnham & Anderson rules:
# Delta < 2: Models are statistically indistinguishable; both have substantial support.
# Delta 2-7: The model has considerably less support compared to the top model.
# Delta > 7: The model has very little to no empirical support.
comparison$support_category <- ifelse(
  comparison$delta_AIC < 2,
  "substantial support",
  ifelse(
    comparison$delta_AIC < 7,
    "considerably less support",
    "little support"))

# Save the thoroughly annotated model comparison table as a CSV file.
write.csv(comparison, output_file, row.names = FALSE)
