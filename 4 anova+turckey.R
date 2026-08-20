# ============================================================
# Analisi All_matrices_stats.csv - ANOVA + Tukey HSD
#
# Riproduce il formato di Statistical_significance_results.xlsx:
# variable | comparison | diff | lwr | upr | p_adj | significant
# (i valori di riferimento nel tuo file corrispondono esattamente a
#  un'ANOVA + TukeyHSD calcolata sullo schema di ricodifica "MM", che
#  uso qui di default).
#
# ATTENZIONE - assunzioni dell'ANOVA:
# 1) INDIPENDENZA: come discusso per KW, le tre categorie condividono
#    lo stesso trimming/dataset sorgente (disegno a blocchi), quindi
#    l'assunzione di indipendenza e' violata anche qui.
# 2) NORMALITA' DEI RESIDUI: con questi dati lo Shapiro-Wilk sui
#    residui risulta significativo (p < 0.05) per quasi tutte le
#    variabili -> i residui NON sono normali, assunzione dell'ANOVA
#    violata. Lo script lo verifica e stampa un avviso automatico.
# In sintesi: l'ANOVA/Tukey qui sotto riproduce lo stesso identico
# output del tuo file di riferimento, ma se possibile e' preferibile
# un test non parametrico appaiato (Friedman) o Kruskal-Wallis/Dunn,
# gia' preparati negli script precedenti.
# ============================================================

library(tidyverse)
library(broom)
library(writexl)

# --- 1. Lettura -------------------------------------------------------
data <- read.csv("02_Matrix_Stats/All_matrices_stats.csv", check.names = TRUE)

data <- data %>%
  mutate(
    scheme   = str_extract(matrix, "(6aa|MM|PM)$"),
    category = factor(category, levels = c("allgenes", "lb", "rcv"))
  )

# --- 2. Filtro a un solo schema di ricodifica --------------------------
SCHEME <- "MM"
data_filtered <- data %>% filter(scheme == SCHEME)
stopifnot(all(table(data_filtered$category) == 14))

vars_of_interest <- c("alignment.lenght", "parsimony.informative.sites",
                      "Proportion.variable.sites", "missing.percent")

# --- 3. Verifica assunzioni ANOVA (normalita' residui, omogeneita') ----
assumptions <- map_dfr(vars_of_interest, function(v) {
  formula <- as.formula(paste0("`", v, "` ~ category"))
  fit  <- aov(formula, data = data_filtered)
  sh   <- shapiro.test(residuals(fit))
  lev  <- car::leveneTest(formula, data = data_filtered)
  tibble(
    variable        = v,
    shapiro_p       = round(sh$p.value, 4),
    normal_residui  = ifelse(sh$p.value > 0.05, "OK", "VIOLATA"),
    levene_p        = round(lev$`Pr(>F)`[1], 4),
    varianze_omog   = ifelse(lev$`Pr(>F)`[1] > 0.05, "OK", "VIOLATA")
  )
})
print(assumptions)
if (any(assumptions$normal_residui == "VIOLATA")) {
  warning("Normalita' dei residui violata per alcune variabili: ",
          "i risultati ANOVA/Tukey vanno interpretati con cautela ",
          "(considerare Kruskal-Wallis/Dunn o Friedman).")
}

# --- 4. Statistiche descrittive -----------------------------------------
desc_stats <- data_filtered %>%
  group_by(category) %>%
  summarise(across(all_of(vars_of_interest),
                   list(n = ~n(), mean = ~mean(.x), sd = ~sd(.x),
                        median = ~median(.x), min = ~min(.x), max = ~max(.x)),
                   .names = "{.col}__{.fn}")) %>%
  pivot_longer(-category, names_to = c("variable", "stat"), names_sep = "__") %>%
  pivot_wider(names_from = stat, values_from = value)
print(desc_stats, n = Inf)

# --- 5. ANOVA a una via ---------------------------------------------------
anova_results <- map_dfr(vars_of_interest, function(v) {
  formula <- as.formula(paste0("`", v, "` ~ category"))
  fit <- aov(formula, data = data_filtered)
  tidy(fit) %>%
    filter(term == "category") %>%
    mutate(variable = v, .before = 1)
})
print(anova_results)

# --- 6. Post-hoc Tukey HSD per tutte le variabili -------------------------
tukey_results <- map_dfr(vars_of_interest, function(v) {
  formula <- as.formula(paste0("`", v, "` ~ category"))
  fit <- aov(formula, data = data_filtered)
  tidy(TukeyHSD(fit)) %>%
    mutate(variable = v, .before = 1)
})

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

# --- 7. Salvataggio Excel --------------------------------------------------
write_xlsx(list(Statistical_significance_result = results_table),
           path = "Statistical_significance_results_ANOVA.xlsx")

message("File salvato: Statistical_significance_results_ANOVA.xlsx")
