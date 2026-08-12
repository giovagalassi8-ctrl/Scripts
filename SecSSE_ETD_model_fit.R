#!/usr/bin/env Rscript

# This script fits the parsimonious Examined Trait-Dependent (ETD) model.
# This model tests if the observed ecological trait (e.g., bathymetric depth) directly influences speciation rates, independent of hidden states.
# In this case depth is considered as the default parameter as in the previous scripts, but you can use any other condition.
# It requires a Newick file, containing species names written with underscores, and the .csv file obtained after runnning the 'SecSSE_prepare_trait.R' script.
# Estimated parameters:
#   1 = lambda (speciation) for shallow species (state 0)
#   2 = lambda (speciation) for deep species (state 1)
#   3 = one extinction rate shared across all states
#   4 = transition rate from state 0 to state 1 (shallow to deep)
#   5 = transition rate from state 1 to state 0 (deep to shallow)
# The output of this script has to be used in 'SecSSE_model_comparison' script.


library(ape)
library(secsse)

# Define the input treefile in .nwk format.
tree_file <- "TREEFILE.NWK"
# Definte the .csv file containing the data for the environmental parameter selected for the analysis.
# It could be used the output file of the 'SecSSE_trait_data.R' script.
trait_file <- "TRAIT_FILE.CSV"
# Define the output files name: one in .rds format and the other in .csv format (so you can easily read it).
rds_output_file <- "OUTPUT_ETD_FIT.RDS"
csv_output_file <- "OUTPUT_ETD_RESULTS.CSV"

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

# Define the speciation matrix for the ETD model.
# Parameter 1 is assigned to state 0 (shallow).
# Parameter 2 is assigned to state 1 (deep).
# Concealed states inherit the same rates as their observed counterparts.
speciation_matrix <- rbind(
  c(0, 0, 0, 1),
  c(1, 1, 1, 2)
)

# Generate the lambda (speciation) list based on the transition matrix.
lambda_list <- secsse::create_lambda_list(
  state_names = state_names,
  num_concealed_states = num_concealed_states,
  transition_matrix = speciation_matrix,
  model = "ETD")

# Create the extinction (mu) rate vector structure
mu_vector <- secsse::create_mu_vector(
  state_names = state_names,
  num_concealed_states = num_concealed_states,
  model = "ETD",
  lambda_list = lambda_list)

# Link the extinction rates so a single parameter (ID 3) estimates extinction for all states
mu_vector[] <- 3

# Define directional state transitions (trait evolution rates).
# Parameter 4: rate of shifting from state 0 to 1 (shallow -> deep).
# Parameter 5: rate of shifting from state 1 to 0 (deep -> shallow).
shift_matrix <- rbind(
  c(0, 1, 4),
  c(1, 0, 5)
)

# Generate the Q matrix (transition rates) using the shift matrix.

q_matrix <- secsse::create_q_matrix(
  state_names = state_names,
  num_concealed_states = num_concealed_states,
  shift_matrix = shift_matrix,
  # diff.conceal = FALSE keeps transition rates independent of the hidden states.
  diff.conceal = FALSE)

# Combine lambda, mu, and q into a single parameter list for the model fitting
parameter_structure <- list(
  lambda_list,
  mu_vector,
  q_matrix)

# Define which parameters to estimate (1 through 5) and which to fix (0).
parameters_to_estimate <- 1:5
parameters_to_fix <- 0
fixed_values <- 0

# Assign the chosen initial starting value (0.3) to all 5 estimated parameters.
initial_values <- rep(
  initial_rate,
  length(parameters_to_estimate))

# Fit the ETD model using Maximum Likelihood optimization.
fit <- secsse::cla_secsse_ml(
  # The phylogenetic tree object containing the branching times and topology.
  phy = tree,
  # A named numeric vector linking each tree tip to its observed state (0 = shallow, 1 = deep).
  traits = trait,
  # The number of hidden states to be considered by the model (2 in this specific setup).
  num_concealed_states = num_concealed_states,
  # The structured list containing the ETD rules for speciation, extinction, and transitions.
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
  verbose = FALSE)

# Extract the maximum likelihood estimates for the 5 parameters.
parameter_estimates <- secsse::extract_par_vals(parameter_structure, fit$MLpars)

# Assign biological meaning to the estimated parameter names.
names(parameter_estimates) <- c(
  "lambda_shallow",
  "lambda_deep",
  "mu_all_states",
  "q_state_0_to_1",
  "q_state_1_to_0")

# Calculate model statistics.
number_of_parameters <- length(parameters_to_estimate)

# Calculate the Akaike Information Criterion (AIC) for model comparison.
aic <- (
  -2 * fit$ML +
    2 * number_of_parameters
)

# Compile a comprehensive list containing the model fit and metadata.
results <- list(
  model = "ETD",
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

# Save the complete results object as an RDS file
saveRDS(results, rds_output_file)

# Create a clean dataframe for the summary table
results_table <- data.frame(
  model = "ETD",
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

# Export the summary dataframe as a CSV file
write.csv(results_table, csv_output_file, row.names = FALSE)
