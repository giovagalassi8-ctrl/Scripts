#!/usr/bin/env Rscript

# This script fits a Concealed Trait-Dependent (CTD) model. 
# The CTD model tests the hypothesis that diversification rate heterogeneity is driven by an unobserved (hidden) trait rather than the observed trait (e.g., bathymetric depth).
# In this case depth is considered as the default parameter as in the previous scripts, but you can use any other condition.
# It requires a Newick file, containing species names written with underscores, and the .csv file obtained after runnning the 'SecSSE_prepare_trait.R' script.
# Estimated parameters (Total = 5):
#   1 = lambda (speciation) driven by concealed state A
#   2 = lambda (speciation) driven by concealed state B
#   3 = one extinction rate shared across all states
#   4 = transition rate from observed state 0 to state 1
#   5 = transition rate from observed state 1 to state 0
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
  model = "CTD")

# Create the extinction (mu) rate vector structure.
mu_vector <- secsse::create_mu_vector(
  state_names = state_names,
  num_concealed_states = num_concealed_states,
  model = "CTD",
  lambda_list = lambda_list)

# Constrain the model to estimate a single, shared extinction rate (ID 3) across all states.
mu_vector[] <- 3

# Define directional transitions for the observed states.
# Parameter 4: transition from state 0 to 1.
# Parameter 5: transition from state 1 to 0.
shift_matrix <- rbind(
  c(0, 1, 4),
  c(1, 0, 5)
)

# Generate the Q matrix keeping transition rates independent of the hidden states.
q_matrix <- secsse::create_q_matrix(
  state_names = state_names,
  num_concealed_states = num_concealed_states,
  shift_matrix = shift_matrix,
  diff.conceal = FALSE)

# Combine speciation, extinction, and transition structures into a single list.
parameter_structure <- list(
  lambda_list,
  mu_vector,
  q_matrix)

# Define the parameters to be estimated (1 to 5) and those to be fixed (0).
parameters_to_estimate <- 1:5
parameters_to_fix <- 0
fixed_values <- 0

# Set the initial guess (0.3) for all 5 estimated parameters.
initial_values <- rep(
  initial_rate,
  length(parameters_to_estimate)
)

# Fit the Concealed Trait-Dependent (CTD) model using Maximum Likelihood.
fit <- secsse::cla_secsse_ml(
  # The phylogenetic tree object containing the branching times and topology.
  phy = tree,
  # A named numeric vector linking each tree tip to its observed state (0 = shallow, 1 = deep).
  traits = trait,
  # The number of hidden states to be considered by the model (2 in this setup).
  num_concealed_states = num_concealed_states,
  # The structured list containing the CTD rules for speciation, extinction, and transitions.
  idparslist = parameter_structure,
  # A vector listing the IDs of the 5 parameters to be estimated by the algorithm.
  idparsopt = parameters_to_estimate,
  # A vector providing the initial starting values (0.3) for the 5 parameters being estimated.
  initparsopt = initial_values,
  # A vector listing the IDs of any parameters constrained to a strictly fixed value.
  idparsfix = parameters_to_fix,
  # A vector specifying the exact numerical values for those fixed parameters (e.g., 0).
  parsfix = fixed_values,
  # Conditions the likelihood computation on the survival of the root node to the present day.
  cond = "proper_cond",
  # Computes the probabilities of the root states based on their theoretical equilibrium frequencies.
  root_state_weight = "proper_weights",
  # Applies the calculated sampling fraction to correct for missing extant species in the phylogeny.
  sampling_fraction = sampling_fraction,
  # Allows the optimization algorithm to run without cycle limits until mathematical convergence.
  num_cycles = Inf,
  # Mutes the console output during the iterative parameter estimation process.
  verbose = FALSE
)

# Extract the maximum likelihood parameter estimates.
parameter_estimates <- secsse::extract_par_vals(
  parameter_structure,
  fit$MLpars)

# Name the parameters to reflect their biological role in the CTD model.
names(parameter_estimates) <- c(
  "lambda_concealed_A",
  "lambda_concealed_B",
  "mu_all_states",
  "q_state_0_to_1",
  "q_state_1_to_0")

# Compute model statistics.
number_of_parameters <- length(parameters_to_estimate)

# Calculate the Akaike Information Criterion (AIC) for model selection.
aic <- (
  -2 * fit$ML +
    2 * number_of_parameters
)

# Store all the model information, settings, and results into a comprehensive list.
results <- list(
  model = "CTD",
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

# Save the complete results object as an RDS file.
saveRDS(results, rds_output_file)

# Create a clear dataframe to summarize the key outputs.
results_table <- data.frame(
  model = "CTD",
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

# Export the summary dataframe as a CSV table.
write.csv(results_table, csv_output_file, row.names = FALSE)
