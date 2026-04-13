# 0) Clear workspace
rm(list = ls())

# 1) Load libraries
library(dplyr)
library(tidyr)
library(purrr)

# 2) Read and filter data
dataset <- read.csv("C:/Users/c3371138/Dropbox/optical-imaging/AF_optical_3factor.csv") %>%
  filter(record_id >= 2001)

# 3) Create MedDiet quartiles
dataset <- dataset %>%
  mutate(
    MedDiet_Quartiles = cut(
      MedDiet.Panagiotakis,
      breaks = quantile(MedDiet.Panagiotakis, probs = seq(0, 1, 0.25), na.rm = TRUE),
      include.lowest = TRUE,
      labels = c("Q1","Q2","Q3","Q4")
    )
  )

# 4) Define continuous variables
continuous_vars <- c(
  "age", "bmi", "anu2_total_yrs_edu", "hr_average", "kJwithDF",
  "waist_hip_ratio", "dbp", "sbp", "LDL", "chol_mmoll",
  "hdl", "trigs", "glucose", "totalmvpa"
)

# 5) Function to compute summary stats (Mean ± SD)
calculate_summary <- function(data, group_var, vars) {
  cols <- c(group_var, vars)
  
  # per‐group (quartile) stats
  summary_stats <- data %>%
    dplyr::select(dplyr::all_of(cols)) %>%
    pivot_longer(
      cols      = -dplyr::all_of(group_var),
      names_to  = "variable",
      values_to = "value"
    ) %>%
    group_by(across(all_of(group_var)), variable) %>%
    summarise(
      Mean = mean(value, na.rm = TRUE),
      SD   = sd(value,   na.rm = TRUE),
      .groups = "drop"
    )
  
  # whole‐sample stats
  whole_sample <- data %>%
    dplyr::select(dplyr::all_of(vars)) %>%
    pivot_longer(
      cols      = everything(),
      names_to  = "variable",
      values_to = "value"
    ) %>%
    group_by(variable) %>%
    summarise(
      Mean = mean(value, na.rm = TRUE),
      SD   = sd(value,   na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(!!group_var := "Whole Sample")
  
  bind_rows(whole_sample, summary_stats) %>%
    arrange(variable)
}

# 6) Generate summary and round
med_diet_summary <- calculate_summary(dataset, "MedDiet_Quartiles", continuous_vars) %>%
  mutate(
    Mean = round(Mean, 2),
    SD   = round(SD,   2)
  )

# 7) Create the wide “Mean (SD)” table
wide_summary <- med_diet_summary %>%
  unite("Stat", Mean, SD, sep = " (") %>%
  mutate(Stat = paste0(Stat, ")")) %>%
  dplyr::select(MedDiet_Quartiles, variable, Stat) %>%
  pivot_wider(
    names_from  = MedDiet_Quartiles,
    values_from = Stat
  )

# 8) Run one‐way ANOVAs and extract p-values
anova_pvalues <- map_df(continuous_vars, function(var) {
  form   <- as.formula(paste(var, "~ MedDiet_Quartiles"))
  aov_mod <- aov(form, data = dataset)
  pval   <- summary(aov_mod)[[1]]["Pr(>F)"][1,1]
  tibble(variable = var, p_value = ifelse(pval < .001, "<.001", sprintf("%.3f", pval)))
})

# 9) Combine the wide table with ANOVA p-values
final_table <- wide_summary %>%
  left_join(anova_pvalues, by = "variable")

# 10) Print the final table
print(final_table)

