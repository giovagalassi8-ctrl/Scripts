#!/usr/bin/env Rscript

# This script tests whether the summary statistics of different supermatrices differ among filtering categories (allgenes, lb, rcv),
# using a one-way ANOVA followed by a Tukey HSD post-hoc test.
# It requires a CSV file containing one row per supermatrix, with the matrix name, its filtering category and the statistics of interest
# (which can be obtained with 'AMAS.py summary', see the 'AMAS_matrix_stats.sh' script).
# The analysis is performed only on the supermatrices of one approach (by default MM, i.e. mixture model), identified by the suffix of the matrix name.
# The final table is saved in .xlsx format with the following columns:
#   variable | comparison | diff | lwr | upr | p_adj | significant
# WARNING - ANOVA assumptions:
#   1. INDEPENDENCE: the three categories share the same trimming/source dataset (block design), so the independence assumption is violated.
#   2. NORMALITY OF RESIDUALS: with the analysed data the Shapiro-Wilk test on the residuals is significant (p < 0.05) for almost all the variables,
#      so the residuals are NOT normally distributed and the assumption is violated. The script checks it and prints an automatic warning.
# The ANOVA/Tukey results should therefore be interpreted with caution: when possible, a paired non-parametric test (Friedman)
# or a Kruskal-Wallis/Dunn test is preferable.


library(tidyverse)
library(broom)
library(writexl)

# Read the CSV file containing the statistics of all the supermatrices (change with the correct file path).
data <- read.csv("02_Matrix_Stats/All_matrices_stats.csv", check.names = TRUE)

# Extract the approach (6aa, MM or PM) from the suffix of the matrix name and convert the category column into a factor.
data <- data %>%
  mutate(
    scheme   = str_extract(matrix, "(6aa|MM|PM)$"),
    category = factor(category, levels = c("allgenes", "lb", "rcv"))
  )

# Keep only the supermatrices of one approach (change accordingly).
SCHEME <- "MM"
data_filtered <- data %>% filter(scheme == SCHEME)
# Check that every category contains the same number of supermatrices (in this case 14; change accordingly).
stopifnot(all(table(data_filtered$category) == 14))

# Select the columns (statistics) to be tested (change with the names of the columns of interest).
vars_of_interest <- c("alignment.lenght", "parsimony.informative.sites",
                      "Proportion.variable.sites", "missing.percent")

# Check the ANOVA assumptions for each statistic: normality of the residuals (Shapiro-Wilk test) and homogeneity of variances (Levene test, from the car package).
# A p-value lower than 0.05 means that the assumption is violated.
assumptions <- map_dfr(vars_of_interest, function(v) {
  # Build the formula 'statistic ~ category' (the backticks handle column names containing special characters).
  formula <- as.formula(paste0("`", v, "` ~ category"))
  fit  <- aov(formula, data = data_filtered)
  sh   <- shapiro.test(residuals(fit))
  lev  <- car::leveneTest(formula, data = data_filtered)
  tibble(
    variable         = v,
    shapiro_p        = round(sh$p.value, 4),
    normal_residuals = ifelse(sh$p.value > 0.05, "OK", "VIOLATED"),
    levene_p         = round(lev$`Pr(>F)`[1], 4),
    equal_variances  = ifelse(lev$`Pr(>F)`[1] > 0.05, "OK", "VIOLATED")
  )
})
# Print the results of the assumption checks.
print(assumptions)
# Print a warning if the normality of the residuals is violated for at least one statistic.
if (any(assumptions$normal_residuals == "VIOLATED")) {
  warning("Normality of residuals violated for some variables: ",
          "the ANOVA/Tukey results should be interpreted with caution ",
          "(consider Kruskal-Wallis/Dunn or Friedman).")
}

# Calculate the descriptive statistics (number of observations, mean, standard deviation, median, minimum and maximum) of each statistic for each category.
desc_stats <- data_filtered %>%
  group_by(category) %>%
  summarise(across(all_of(vars_of_interest),
                   list(n = ~n(), mean = ~mean(.x), sd = ~sd(.x),
                        median = ~median(.x), min = ~min(.x), max = ~max(.x)),
                   .names = "{.col}__{.fn}")) %>%
  pivot_longer(-category, names_to = c("variable", "stat"), names_sep = "__") %>%
  pivot_wider(names_from = stat, values_from = value)
print(desc_stats, n = Inf)

# Perform a one-way ANOVA for each statistic, keeping only the results of the category term.
anova_results <- map_dfr(vars_of_interest, function(v) {
  formula <- as.formula(paste0("`", v, "` ~ category"))
  fit <- aov(formula, data = data_filtered)
  tidy(fit) %>%
    filter(term == "category") %>%
    mutate(variable = v, .before = 1)
})
print(anova_results)

# Perform the Tukey HSD post-hoc test for each statistic, comparing all the pairs of categories.
tukey_results <- map_dfr(vars_of_interest, function(v) {
  formula <- as.formula(paste0("`", v, "` ~ category"))
  fit <- aov(formula, data = data_filtered)
  tidy(TukeyHSD(fit)) %>%
    mutate(variable = v, .before = 1)
})

# Build the final table, renaming the columns and flagging the significant comparisons (adjusted p-value < 0.05).
results_table <- tukey_results %>%
  transmute(
    variable,
    comparison  = contrast,
    diff        = estimate,
    lwr         = conf.low,
    upr         = conf.high,
    p_adj       = adj.p.value,
    significant = ifelse(adj.p.value < 0.05, "YES", "NO")
  )

print(results_table)

# Save the final table in .xlsx format (change the sheet and file names as desired).
write_xlsx(list(Statistical_significance_result = results_table),
           path = "Statistical_significance_results_ANOVA.xlsx")

message("File saved: Statistical_significance_results_ANOVA.xlsx")
