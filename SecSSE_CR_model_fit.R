#!/usr/bin/env Rscript

# This script fits a Constant Rates (CR) model using the SecSSE package to evaluate whether lineage diversification heterogeneity is linked to an environment parameter.
# In this case depth is considered as the default parameter as in the previous scripts, but you can use any other condition.
# It requires a Newick file, containing species names written with underscores, and the .csv file obtained after runnning the 'SecSSE_prepare_trait.R' script.
# Estimated parameters:
#   1 = one speciation rate shared across all states
#   2 = one extinction rate shared across all states
#   3 = transition rate from state 0 to state 1
#   4 = transition rate from state 1 to state 0
# The output of this script has to be used in 'SecSSE_model_comparison' script.


library(ape)
library(secsse)

# Define the input treefile in .nwk format.
tree_file <- "TREEFILE.NWK"
# Definte the .csv file containing the data for the environmental parameter selected for the analysis.
# It could be used the output file of the 'SecSSE_trait_data.R' script.
trait_file <- "TRAIT_FILE.CSV"
# Define the output files name: one in .rds format and the other in .csv format (so you can easily read it).
rds_output_file <- "OUTPUT_CR_FIT.RDS"
csv_output_file <- "OUTPUT_CR_RESULTS.CSV"

# Define the estimated total number of accepted extant species.
# To gain this number you can use online databases (e.g., WoRMS for marine species).
total_known_species <- 610   #(Change accordingly).

# Define the common starting value. 
# Running the script multiple times with different starting values (e.g., 0.1, 0.5, 1.0) 
# helps ensure the maximum likelihood optimization isn't stuck in a local optimum.
initial_rate <- 0.3

# Load the phylogenetic tree.
tree <- read.tree(tree_file)
# Check if the file contains multiple trees and just in case it uses only the first one. 
if (inherits(tree, "multiPhylo")) {
  tree <- tree[[1]]
}

# Calculate the phylogenetic sampling fraction.
sampling_fraction_value <- Ntip(tree) / total_known_species
# Since we have two observed states (0 and 1), we assign the same sampling fraction to both.
sampling_fraction <- rep(sampling_fraction_value, 2)

# Load the table containing the environmental trait data.
trait_data <- read.csv(trait_file, stringsAsFactors = FALSE, na.strings = c("NA", ""))

# Convert tree labels from underscores to spaces to match with the species name in the environment file.
# If the environment file contains species name with the underscore (e.g., Verpa_penis), skip this line.
trait_data$tree_label <- gsub(" ", "_", trait_data$species)

# Reorder the trait table to match the tree tips.
trait_data <- trait_data[
  match(tree$tip.label, trait_data$tree_label),
]

# Extract the trait column as a numeric vector and assign the tip labels as names.
# This named numeric vector is the standard input format required by secsse.
trait <- as.numeric(trait_data$trait)
names(trait) <- tree$tip.label

# Define observed and concealed states.
# Assume a binary observed trait (0 and 1).
state_names <- c(0, 1)
# Assume 2 concealed states (e.g., A and B).
num_concealed_states <- 2

# Define the speciation matrix structure.
# In a Constant Rate (CR) model, all states (observed and concealed) share the same speciation rate.
# Row 1 indicates observed state 0, Row 2 indicates observed state 1.
speciation_matrix <- rbind(
  c(0, 0, 0, 1),
  c(1, 1, 1, 1)
)

# Generate the lambda (speciation) list based on the transition matrix.
lambda_list <- secsse::create_lambda_list(
  state_names = state_names,
  num_concealed_states = num_concealed_states,
  transition_matrix = speciation_matrix,
  model = "CR")

# Create the extinction (mu) rate vector.
mu_vector <- secsse::create_mu_vector(
  state_names = state_names,
  num_concealed_states = num_concealed_states,
  model = "CR",
  lambda_list = lambda_list)

# Link the extinction rates so that a single parameter (parameter ID 2) estimates extinction across all states.
mu_vector[] <- 2

# Define directional state transitions (trait evolution rates).
# 3 represents the rate of shifting from state 0 to 1.
# 4 represents the rate of shifting from state 1 to 0.
shift_matrix <- rbind(
  c(0, 1, 3),
  c(1, 0, 4)
)

# Generate the Q matrix (transition rates) using the shift matrix.
q_matrix <- secsse::create_q_matrix(
  state_names = state_names,
  num_concealed_states = num_concealed_states,
  shift_matrix = shift_matrix,
  # diff.conceal = FALSE means transition rates between observed states do not depend on the concealed state.
  diff.conceal = FALSE)

# Combine the lambda, mu, and q structures into a single parameter list for the model.
parameter_structure <- list(
  lambda_list,
  mu_vector,
  q_matrix)

# Define which parameters are estimated (1: speciation, 2: extinction, 3 & 4: transition rates).
parameters_to_estimate <- 1:4
# Parameter ID 0 is used for parameters that should be kept strictly fixed (e.g., 0.0).
parameters_to_fix <- 0
fixed_values <- 0

# Assign the chosen initial starting value (e.g., 0.3) to all parameters being estimated.
initial_values <- rep(
  initial_rate,
  length(parameters_to_estimate))


# Fit the Constant Rates (CR) SecSSE model using Maximum Likelihood.
fit <- secsse::cla_secsse_ml(
  # The phylogenetic tree object (must be of class 'phylo').
  phy = tree,
  traits = trait,
  # The number of concealed (hidden) states used in the model.
  num_concealed_states = num_concealed_states,
  # A list containing the structures for speciation, extinction, and transition matrices.
  idparslist = parameter_structure,
  # A numeric vector indicating the IDs of the parameters to be estimated.
  idparsopt = parameters_to_estimate,
  # A numeric vector providing the initial starting values for the parameters being estimated.
  initparsopt = initial_values,
  # A numeric vector indicating the IDs of the parameters to be fixed (not estimated).
  idparsfix = parameters_to_fix,
  # A numeric vector providing the exact values for the fixed parameters
  parsfix = fixed_values,
  # Conditions the likelihood on the survival of the root to the present
  cond = "proper_cond",
  # Assigns weights to the root states based on their equilibrium frequencies
  root_state_weight = "proper_weights",
  # A numeric vector defining the proportion of extant species sampled in the tree for each state
  sampling_fraction = sampling_fraction,
  # Sets the maximum number of optimization cycles (Inf means it runs until convergence)
  num_cycles = Inf,
  # Suppresses the printing of the optimization progress to the console
  verbose = FALSE)

# Extract the maximum likelihood estimates for each parameter.
parameter_estimates <- secsse::extract_par_vals(parameter_structure, fit$MLpars)

# Assign meaningful names to the estimated parameters for clarity.
names(parameter_estimates) <- c(
  "lambda_all_states",
  "mu_all_states",
  "q_state_0_to_1",
  "q_state_1_to_0")

# Calculate model statistics: count the number of estimated parameters.
number_of_parameters <- length(parameters_to_estimate)

# Calculate the Akaike Information Criterion (AIC) for model comparison.
aic <- (
  -2 * fit$ML +
    2 * number_of_parameters
)

# Compile a comprehensive list containing the model fit and all metadata.
results <- list(
  model = "CR",
  log_likelihood = fit$ML,
  number_of_parameters = number_of_parameters,
  AIC = aic,
  parameter_estimates = parameter_estimates,
  sampling_fraction = sampling_fraction,
  total_known_species = total_known_species,
  initial_values = initial_values,
  num_cycles = Inf,
  parameter_structure = parameter_structure,
  fit = fit)

# Save the complete results object as an RDS file (useful for future plotting/analysis in R).
saveRDS(results, rds_output_file)

# Create a flattened, readable dataframe containing the key summary statistics.
results_table <- data.frame(
  model = "CR",
  log_likelihood = fit$ML,
  number_of_parameters = number_of_parameters,
  AIC = aic,
  parameter_id = parameters_to_estimate,
  parameter = names(parameter_estimates),
  estimate = as.numeric(parameter_estimates),
  initial_value = initial_values,
  sampling_fraction_state_0 = sampling_fraction[1],
  sampling_fraction_state_1 = sampling_fraction[2],
  stringsAsFactors = FALSE)

# Export the summary dataframe as a CSV file.
write.csv(results_table, csv_output_file, row.names = FALSE)
