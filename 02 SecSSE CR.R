# =============================================================================
# Script: 03b_fit_CR.R
#
# Purpose:
# Fit the parsimonious Constant Rate model.
#
# Required trait CSV columns:
#   species
#   trait
#
# Estimated parameters:
#   1 = one speciation rate shared across all states
#   2 = one extinction rate shared across all states
#   3 = transition rate from state 0 to state 1
#   4 = transition rate from state 1 to state 0
#
# Outputs:
#   analyses/SSE/Terebridae_CR_fit.rds
#   analyses/SSE/Terebridae_CR_results.csv
# =============================================================================


library(ape)
library(secsse)


# Define input and output files.
tree_file <- "00_Data/Mollusca/Conidae/Conidae.nwk"
trait_file <- "Conidae_depth_trait.csv"

rds_output_file <- "Conidae_CR_fit.rds"
csv_output_file <- "Conidae_CR_results.csv"


# Define the total number of known species.
total_known_species <- 1080


# Define the common starting value.
initial_rate <- 0.3


# Load the phylogenetic tree.
tree <- read.tree(tree_file)

if (inherits(tree, "multiPhylo")) {
  tree <- tree[[1]]
}


# Calculate the phylogenetic sampling fraction.
sampling_fraction_value <- Ntip(tree) / total_known_species

sampling_fraction <- rep(
  sampling_fraction_value,
  2
)


# Load the trait data.
trait_data <- read.csv(
  trait_file,
  stringsAsFactors = FALSE,
  na.strings = c("NA", "")
)


# Convert species names to tree-label format.
trait_data$tree_label <- gsub(
  " ",
  "_",
  trait_data$species
)


# Reorder trait values to match the tree tips.
trait_data <- trait_data[
  match(tree$tip.label, trait_data$tree_label),
]


# Extract the trait vector.
trait <- as.numeric(trait_data$trait)
names(trait) <- tree$tip.label


# Define observed and concealed states.
state_names <- c(0, 1)
num_concealed_states <- 2


# Use one shared speciation-rate indicator.
speciation_matrix <- rbind(
  c(0, 0, 0, 1),
  c(1, 1, 1, 1)
)


# Create the constant speciation-rate structure.
lambda_list <- secsse::create_lambda_list(
  state_names = state_names,
  num_concealed_states = num_concealed_states,
  transition_matrix = speciation_matrix,
  model = "CR"
)


# Create the constant extinction-rate structure.
mu_vector <- secsse::create_mu_vector(
  state_names = state_names,
  num_concealed_states = num_concealed_states,
  model = "CR",
  lambda_list = lambda_list
)


# Use one extinction rate across all states.
mu_vector[] <- 2


# Define directional state transitions.
shift_matrix <- rbind(
  c(0, 1, 3),
  c(1, 0, 4)
)


# Use the same transition rates for observed and concealed states.
q_matrix <- secsse::create_q_matrix(
  state_names = state_names,
  num_concealed_states = num_concealed_states,
  shift_matrix = shift_matrix,
  diff.conceal = FALSE
)


# Combine the model components.
parameter_structure <- list(
  lambda_list,
  mu_vector,
  q_matrix
)


# Define estimated and fixed parameters.
parameters_to_estimate <- 1:4
parameters_to_fix <- 0
fixed_values <- 0

initial_values <- rep(
  initial_rate,
  length(parameters_to_estimate)
)


# Fit the CR model.
fit <- secsse::cla_secsse_ml(
  phy = tree,
  traits = trait,
  num_concealed_states = num_concealed_states,
  idparslist = parameter_structure,
  idparsopt = parameters_to_estimate,
  initparsopt = initial_values,
  idparsfix = parameters_to_fix,
  parsfix = fixed_values,
  cond = "proper_cond",
  root_state_weight = "proper_weights",
  sampling_fraction = sampling_fraction,
  num_cycles = Inf,
  verbose = FALSE
)


# Extract parameter estimates.
parameter_estimates <- secsse::extract_par_vals(
  parameter_structure,
  fit$MLpars
)


names(parameter_estimates) <- c(
  "lambda_all_states",
  "mu_all_states",
  "q_state_0_to_1",
  "q_state_1_to_0"
)


# Calculate model statistics.
number_of_parameters <- length(parameters_to_estimate)

aic <- (
  -2 * fit$ML +
    2 * number_of_parameters
)


# Store the complete model results.
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
  fit = fit
)


# Save the complete model object.
saveRDS(
  results,
  rds_output_file
)


# Create the readable result table.
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
  stringsAsFactors = FALSE
)


# Save the readable result table.
write.csv(
  results_table,
  csv_output_file,
  row.names = FALSE
)
