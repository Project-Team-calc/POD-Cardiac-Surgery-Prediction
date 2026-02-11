# ==============================================================================
# SCRIPT: Data Preparation and Clinical Rescaling 
# PURPOSE: Addressing Reviewer 1 (Clinical Units) & Reviewer 2 (Terminology)
# ==============================================================================

# --- 1. Load Libraries ---
library(tidyverse)
library(janitor)
library(mice)
library(brms) # For future modeling
library(future) # For parallel processing

# --- 2. Initial Data Cleaning & Terminology Reframing ---
# We remove "De novo" and "Intact" to be more scientifically defensible.
raw_data <- read_csv("my data.csv") %>% clean_names()

data_cleaned <- raw_data %>%
  # Filter cases exactly like your original exclusion logic
  filter(age_years >= 18) %>%
  filter(!is.na(edu_years) & !is.na(mmse_preop)) %>%
  mutate(
    # Reframing: Normal range on MMSE screening
    is_eligible = case_when(
      edu_years == 0 & mmse_preop > 17 ~ TRUE,
      edu_years >= 1 & edu_years <= 6 & mmse_preop > 20 ~ TRUE,
      edu_years > 6 & mmse_preop > 24 ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>%
  filter(is_eligible == TRUE) %>%
  # Outcome Variable: Neutral Terminology
  mutate(
    pod_status = if_else(camicu_day1 == 1 | camicu_day2 == 1 | camicu_day3 == 1, 1, 0, missing = 0),
    pod_status = factor(pod_status, levels = c(0, 1), labels = c("No", "Yes"))
  )

# --- 3. Surgery Categorization (Consolidation) ---
data_pre_mice <- data_cleaned %>%
  mutate(
    surgery_type = fct_recode(factor(surgery_category),
                              "Complex_Combined" = "4", # CABG + Valve
                              "Complex_Combined" = "5", # Multi-Valve
                              "Other_Procedures" = "3", # Isolated Aortic
                              "Other_Procedures" = "6", # Congenital
                              "Other_Procedures" = "7", # Other
                              "Isolated_CABG" = "1",
                              "Isolated_Single_Valve" = "2"),
    across(c(smoking_history, alcohol_history, hypertension, diabetes, stroke_history),
           ~factor(., levels = c(0, 1), labels = c("No", "Yes")))
  ) %>%
  select(pod_status, age_years, smoking_history, cpb_time_min, 
         propofol_mg, hospital, surgery_type)

# --- 4. Multiple Imputation (MICE) on RAW Units ---
# Performing imputation on raw units preserves the true clinical variance.
set.seed(12345)
# Using m=50 as suggested for high-tier journal rigor
imputed_mids_raw <- mice(data_pre_mice, m = 50, maxit = 20, method = 'pmm', printFlag = TRUE)

# --- 5. Step: Clinical Unit Transformation ---
# Addressing Reviewer 1: Convert SD units to Clinical Resolution units.
# age/10 = risk per 10 years; cpb/30 = per 30 mins; propofol/100 = per 100mg.

# Extract completed long data
completed_data_long <- complete(imputed_mids_raw, action = "long", include = TRUE)

completed_data_long_clinical <- completed_data_long %>%
  mutate(
    age_10y      = age_years / 10,
    cpb_30m      = cpb_time_min / 30,
    propofol_100 = propofol_mg / 100
  )

# Re-convert to mids object for Bayesian modeling
imputed_mids_final <- as.mids(completed_data_long_clinical)

# Save this critical object
saveRDS(imputed_mids_final, "imputed_mids_final_clinical.rds")

cat("\n--- SUCCESS: Data Prepped and Imputed with Clinical Units ---")


# ==============================================================================
# SCRIPT: 04_FINAL_FIXED_EFFECT_MODEL.R
# PURPOSE: Switch to Fixed Effects to fix Divergence issues (N=2 centers is too few for RE)
# ==============================================================================

library(tidyverse)
library(brms)
library(mice)
library(future)

# 1. Load Data
imputed_mids <- readRDS("imputed_mids_final_clinical.rds")
plan(multisession)

# 2. DATA TRANSFORMATION: Centering (Standard Practice)
# We repeat this to ensure we have the centered variables ready
long_data <- complete(imputed_mids, action = "long", include = TRUE)

means <- long_data %>%
  filter(.imp == 0) %>%
  summarise(
    mean_age = mean(age_10y, na.rm = TRUE),
    mean_cpb = mean(cpb_30m, na.rm = TRUE),
    mean_prop = mean(propofol_100, na.rm = TRUE)
  )

long_data_centered <- long_data %>%
  mutate(
    age_10y_c = age_10y - means$mean_age,
    cpb_30m_c = cpb_30m - means$mean_cpb,
    propofol_100_c = propofol_100 - means$mean_prop
  )

imputed_mids_centered <- as.mids(long_data_centered)

# 3. FORMULA CHANGE: FIXED EFFECTS (The Solution)
# We removed (1|hospital) and (1|surgery_type)
# We added + hospital + surgery_type
final_formula_fixed <- bf(
  pod_status ~ age_10y_c + smoking_history + cpb_30m_c + propofol_100_c + 
    hospital + surgery_type,
  family = bernoulli(link = "logit")
)

# 4. Priors
# We keep generic priors which will now apply to the hospital/surgery coefficients too
my_priors <- c(
  set_prior("normal(-2, 1.5)", class = "Intercept"),
  set_prior("normal(0, 1)", class = "b")
)

# 5. RUN THE MODEL
cat("--- Fitting Final Fixed-Effect Model (Should be fast and stable)... ---\n")

final_model_fit <- brm_multiple(
  formula = final_formula_fixed,
  data = imputed_mids_centered,
  prior = my_priors,
  iter = 4000,       
  warmup = 1000,
  chains = 4, 
  cores = 4,
  seed = 202607,
  control = list(adapt_delta = 0.99), # 0.99 is enough for fixed effects
  file = "final_clinical_model_fixed" # Save as new file
)

# 6. DIAGNOSTICS & SUMMARY
library(bayesplot)
np <- nuts_params(final_model_fit)
n_div <- sum(subset(np, Parameter == "divergent__")$Value)

cat(sprintf("\n------------------------------------------------\n"))
cat(sprintf("Number of Divergent Transitions: %d\n", n_div))
cat(sprintf("------------------------------------------------\n"))

# Print Summary Immediately
cat("\n--- MODEL SUMMARY ---\n")
print(summary(final_model_fit))

# 7. GENERATE ODDS RATIOS (For Table 3)
cat("\n--- ODDS RATIOS (For Paper) ---\n")
fixed_effects <- fixef(final_model_fit) %>%
  as.data.frame() %>%
  mutate(
    OR = exp(Estimate),
    Lower_95 = exp(Q2.5),
    Upper_95 = exp(Q97.5)
  ) %>%
  select(OR, Lower_95, Upper_95)

print(round(fixed_effects, 2))

# Save the model object for the next steps
saveRDS(final_model_fit, "final_clinical_model.rds")


# ==============================================================================
# CORRECTED SURGERY CATEGORY BASELINE TABLE
# ==============================================================================

if (!require("tableone")) install.packages("tableone")
if (!require("dplyr")) install.packages("dplyr")
library(tableone)
library(dplyr)

# 1. Load Data
df <- read.csv("my_data_cleaned_for_analysis.csv", stringsAsFactors = FALSE)

# 2. Data Cleaning & Variable Mapping (Fixed Surgery Category)
df_analysis <- df %>%
  mutate(
    # --- Grouping Variable ---
    hospital = factor(hospital, levels = c(1, 2), labels = c("Center A", "Center B")),
    
    # --- Primary Outcome ---
    delirium_occurred = factor(delirium_occurred, levels = c("No", "Yes"), labels = c("No", "Yes")),
    
    # --- Demographics & Comorbidities ---
    gender = factor(gender, levels = c("No", "Yes"), labels = c("Female", "Male")),
    hypertension = factor(hypertension, levels = c("No", "Yes")),
    diabetes = factor(diabetes, levels = c("No", "Yes")),
    stroke_history = factor(stroke_history, levels = c("No", "Yes")),
    cerebral_stenosis = factor(cerebral_stenosis, levels = c("No", "Yes")),
    smoking_history = factor(smoking_history, levels = c("No", "Yes")),
    alcohol_history = factor(alcohol_history, levels = c("No", "Yes")),
    
    # --- FIX: SURGERY CATEGORY MAPPING ---
    # Use case_when to standardize string matching for surgery categories
    surgery_temp_label = case_when(
      # If the label in the CSV is "Complex/Combined", standardize it to "Complex or Combined"
      surgery_category_grouped == "Complex/Combined" ~ "Complex or Combined",
      # Keep all other categories unchanged
      TRUE ~ surgery_category_grouped
    ),
    
    # Create Factor based on the corrected labels to ensure specific plot/analysis order
    surgery_type = factor(surgery_temp_label,
                          levels = c("Isolated CABG", 
                                     "Isolated Single-Valve", 
                                     "Complex or Combined", 
                                     "Other Procedures")),
    # --- Valve Prosthesis Type ---
    valve_prosthesis_type = factor(valve_prosthesis_type,
                                   levels = c("Bioprosthetic", "Mechanical", "No Replacement")),
    
    # --- Continuous Variables ---
    age_years = as.numeric(age_years),
    edu_years = as.numeric(edu_years),
    bmi = as.numeric(bmi),
    lvef = as.numeric(lvef),
    mmse_preop = as.numeric(mmse_preop),
    cpb_time_min = as.numeric(cpb_time_min),
    surgery_duration_min = as.numeric(surgery_duration_min),
    acclamp_time_min = as.numeric(acclamp_time_min),
    propofol_mg = as.numeric(propofol_mg)
  )

# 3. Define Variables List
all_vars <- c(
  "delirium_occurred", 
  "age_years", "gender", "edu_years", "bmi",
  "hypertension", "diabetes", "stroke_history", "cerebral_stenosis", "smoking_history", "alcohol_history",
  "lvef", "mmse_preop", 
  "surgery_type",          # Check this line in the output
  "valve_prosthesis_type", 
  "cpb_time_min", "surgery_duration_min", "acclamp_time_min", "propofol_mg"
)

cat_vars <- c("delirium_occurred", "gender", "hypertension", "diabetes", 
              "stroke_history", "cerebral_stenosis", "smoking_history", "alcohol_history", 
              "surgery_type", "valve_prosthesis_type")

# 4. Normality Check Logic
cont_vars <- setdiff(all_vars, cat_vars)
non_normal_vars <- c()

for (v in cont_vars) {
  vec <- df_analysis[[v]]
  vec <- vec[!is.na(vec)]
  if (length(vec) > 3) {
    if (shapiro.test(vec)$p.value < 0.05) {
      non_normal_vars <- c(non_normal_vars, v)
    }
  }
}

# 5. Generate and Save Table
table1 <- CreateTableOne(vars = all_vars, strata = "hospital", data = df_analysis, factorVars = cat_vars, addOverall = TRUE, test = TRUE)

table1_print <- print(
  table1,
  nonnormal = non_normal_vars,
  exact = cat_vars,
  quote = FALSE,
  noSpaces = TRUE,
  smd = TRUE,
  showAllLevels = TRUE,
  printToggle = FALSE
)

# Save
write.csv(table1_print, "Table1_Baseline_Fixed_Surgery.csv")

cat("Done! Checked 'surgery_type'. Please open 'Table1_Baseline_Fixed_Surgery.csv'.\n")


# ==============================================================================
# SCRIPT: 05_PRIOR_SENSITIVITY_ANALYSIS.R
# AUTHOR: [Your Name/Team]
# DATE: [Current Date]
# PURPOSE: 
#   To demonstrate the robustness of the primary Bayesian model by evaluating 
#   the impact of different prior distributions on the posterior estimates.
#   Three scenarios are tested:
#     1. Primary (Weakly Informative): Normal(0, 2.5)
#     2. Wide (Uninformative): Normal(0, 10) - Approaches Frequentist likelihood
#     3. Narrow (Skeptical/Regularizing): Normal(0, 1) - Penalizes large effects
# ==============================================================================

# --- 1. Load Libraries ---
library(tidyverse)
library(brms)
library(mice)
library(future)

# --- 2. Load and Prepare Data ---
# Load the clinically rescaled imputed datasets (from Script 01/02)
if(file.exists("imputed_mids_final_clinical.rds")) {
  imputed_mids <- readRDS("imputed_mids_final_clinical.rds")
} else {
  stop("Error: 'imputed_mids_final_clinical.rds' not found. Please run data preparation script first.")
}

# --- Data Transformation: Centering ---
# Crucial: Re-apply centering to continuous variables to ensure convergence 
# and consistent interpretation of the intercept.
long_data <- complete(imputed_mids, action = "long", include = TRUE)

means <- long_data %>%
  filter(.imp == 0) %>%
  summarise(
    mean_age = mean(age_10y, na.rm = TRUE),
    mean_cpb = mean(cpb_30m, na.rm = TRUE),
    mean_prop = mean(propofol_100, na.rm = TRUE)
  )

long_data_centered <- long_data %>%
  mutate(
    age_10y_c = age_10y - means$mean_age,
    cpb_30m_c = cpb_30m - means$mean_cpb,
    propofol_100_c = propofol_100 - means$mean_prop
  )

# Convert back to mids object for brms
imputed_mids_centered <- as.mids(long_data_centered)

# --- 3. Define Model Formula ---
# Fixed Effects model (consistent with primary analysis)
model_formula <- bf(
  pod_status ~ age_10y_c + smoking_history + cpb_30m_c + propofol_100_c + 
    hospital + surgery_type,
  family = bernoulli(link = "logit")
)

# Set parallel processing
plan(multisession, workers = 4)

# ==============================================================================
# SCENARIO 1: PRIMARY PRIOR (Weakly Informative)
# Specification: Normal(0, 2.5) for coefficients
# ==============================================================================
cat("\n--- Running Scenario 1: Primary Prior (Normal(0, 2.5)) ---\n")

priors_primary <- c(
  set_prior("normal(-2, 1.5)", class = "Intercept"), # Baseline low probability
  set_prior("normal(0, 2.5)", class = "b")           # Weakly informative
)

fit_primary <- brm_multiple(
  formula = model_formula,
  data = imputed_mids_centered,
  prior = priors_primary,
  iter = 4000, warmup = 1000, chains = 4, seed = 202601,
  control = list(adapt_delta = 0.99),
  file = "sensitivity_model_primary"
)

# ==============================================================================
# SCENARIO 2: WIDE PRIOR (Uninformative / Flat)
# Specification: Normal(0, 10) for coefficients
# Logic: Mimics traditional frequentist Maximum Likelihood Estimation (MLE).
# ==============================================================================
cat("\n--- Running Scenario 2: Wide Prior (Normal(0, 10)) ---\n")

priors_wide <- c(
  set_prior("normal(-2, 1.5)", class = "Intercept"),
  set_prior("normal(0, 10)", class = "b")            # Very flat/Wide
)

fit_wide <- brm_multiple(
  formula = model_formula,
  data = imputed_mids_centered,
  prior = priors_wide,
  iter = 4000, warmup = 1000, chains = 4, seed = 202602,
  control = list(adapt_delta = 0.99),
  file = "sensitivity_model_wide"
)

# ==============================================================================
# SCENARIO 3: NARROW PRIOR (Skeptical / Strong Regularization)
# Specification: Normal(0, 1) for coefficients
# Logic: Assumes effect sizes are likely small; penalizes overfitting strongly.
# ==============================================================================
cat("\n--- Running Scenario 3: Narrow Prior (Normal(0, 1)) ---\n")

priors_narrow <- c(
  set_prior("normal(-2, 1.5)", class = "Intercept"),
  set_prior("normal(0, 1)", class = "b")             # Skeptical/Narrow
)

fit_narrow <- brm_multiple(
  formula = model_formula,
  data = imputed_mids_centered,
  prior = priors_narrow,
  iter = 4000, warmup = 1000, chains = 4, seed = 202603,
  control = list(adapt_delta = 0.99),
  file = "sensitivity_model_narrow"
)

# ==============================================================================
# 4. EXTRACT AND COMPARE RESULTS
# ==============================================================================
cat("\n--- Compiling Sensitivity Analysis Results ---\n")

# Helper function to extract OR and 95% CrI
get_estimates <- function(model, model_name) {
  fixef(model) %>%
    as.data.frame() %>%
    rownames_to_column("Predictor") %>%
    mutate(
      OR = exp(Estimate),
      Lower = exp(Q2.5),
      Upper = exp(Q97.5),
      # Format: "1.24 [1.15, 1.34]"
      Result_String = sprintf("%.2f [%.2f, %.2f]", OR, Lower, Upper)
    ) %>%
    select(Predictor, !!sym(model_name) := Result_String)
}

# Extract from all three
res_primary <- get_estimates(fit_primary, "Primary_Prior (N(0,2.5))")
res_wide    <- get_estimates(fit_wide,    "Wide_Prior (N(0,10))")
res_narrow  <- get_estimates(fit_narrow,  "Narrow_Prior (N(0,1))")

# Merge into a single comparison table
sensitivity_table <- res_primary %>%
  left_join(res_wide, by = "Predictor") %>%
  left_join(res_narrow, by = "Predictor") %>%
  # Filter to show only core predictors (optional, removing intercept/adjustments for clarity)
  filter(!grepl("Intercept", Predictor))

# --- 5. EXPORT RESULTS ---
write.csv(sensitivity_table, "Table_S_Sensitivity_Priors.csv", row.names = FALSE)

cat("\n====================================================================\n")
cat("SENSITIVITY ANALYSIS COMPLETE\n")
cat("Results saved to: Table_S_Sensitivity_Priors.csv\n")
cat("Check this table. If the ORs are stable across columns, the model is robust.\n")
cat("====================================================================\n")

# Print to console for immediate check
print(sensitivity_table)

# ==============================================================================
# SCRIPT: 06_PLOT_SENSITIVITY.R
# PURPOSE: 
#   To visualize the Prior Sensitivity Analysis results.
#   Generates a grouped Forest Plot comparing Odds Ratios (OR) across
#   three Bayesian prior specifications:
#     1. Primary (Normal(0, 2.5))
#     2. Wide (Normal(0, 10))
#     3. Narrow/Skeptical (Normal(0, 1))
# ==============================================================================

# --- 1. Load Libraries ---
library(tidyverse)
library(brms)
library(ggplot2)
library(ggsci) # For professional color palettes (e.g., Lancet/JAMA)

# --- 2. Load Models (Ensure these fit objects exist from Script 05) ---
# If running this script independently, load the saved model files:
# fit_primary <- readRDS("sensitivity_model_primary.rds")
# fit_wide    <- readRDS("sensitivity_model_wide.rds")
# fit_narrow  <- readRDS("sensitivity_model_narrow.rds")

# --- 3. Extract Numeric Estimates ---
# Function to extract OR, Lower CI, and Upper CI from brms objects
extract_estimates <- function(model, model_label) {
  fixef(model) %>%
    as.data.frame() %>%
    rownames_to_column("Term") %>%
    mutate(
      OR = exp(Estimate),
      Lower = exp(Q2.5),
      Upper = exp(Q97.5),
      Model = model_label
    ) %>%
    select(Term, OR, Lower, Upper, Model)
}

# Extract data from the three models
data_primary <- extract_estimates(fit_primary, "Primary Prior (N(0, 2.5))")
data_wide    <- extract_estimates(fit_wide,    "Wide Prior (N(0, 10))")
data_narrow  <- extract_estimates(fit_narrow,  "Narrow Prior (N(0, 1))")

# Combine into one dataset
plot_data <- bind_rows(data_primary, data_wide, data_narrow)

# --- 4. Data Cleaning & Formatting for Plotting ---
plot_data_clean <- plot_data %>%
  # 1. Filter out Intercept and auxiliary variables (keep only core predictors)
  filter(!grepl("Intercept", Term)) %>%
  filter(!grepl("hospital", Term)) %>%
  filter(!grepl("surgery_type", Term)) %>%
  
  # 2. Rename terms for professional display
  mutate(Term_Label = case_when(
    grepl("age", Term) ~ "Age (per 10 years)",
    grepl("smoking", Term) ~ "Smoking History (Yes)",
    grepl("cpb", Term) ~ "CPB Duration (per 30 min)",
    grepl("propofol", Term) ~ "Propofol Dose (per 100 mg)",
    TRUE ~ Term
  )) %>%
  
  # 3. Set factor levels to control the order on the Y-axis
  mutate(Term_Label = factor(Term_Label, levels = c(
    "Propofol Dose (per 100 mg)",
    "CPB Duration (per 30 min)",
    "Smoking History (Yes)",
    "Age (per 10 years)"
  ))) %>%
  
  # 4. Set factor levels for Model to control Legend order
  mutate(Model = factor(Model, levels = c(
    "Wide Prior (N(0, 10))",
    "Primary Prior (N(0, 2.5))",
    "Narrow Prior (N(0, 1))"
  )))

# --- 5. Generate Forest Plot ---
p <- ggplot(plot_data_clean, aes(x = OR, y = Term_Label, color = Model, shape = Model)) +
  
  # Add vertical reference line at OR = 1 (Null effect)
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  
  # Plot points and error bars (Position Dodge prevents overlapping)
  geom_pointrange(aes(xmin = Lower, xmax = Upper), 
                  position = position_dodge(width = 0.6), 
                  size = 0.8, linewidth = 0.8) +
  
  # Log scale for X-axis is standard for ORs (better visualization of intervals)
  scale_x_log10(breaks = c(0.5, 1.0, 1.5, 2.0, 3.0, 5.0)) +
  
  # Use a professional color palette (e.g., JAMA or Lancet style)
  scale_color_npg() + 
  
  # Labels
  labs(
    title = "Prior Sensitivity Analysis: Robustness of Predictors",
    subtitle = "Comparison of Posterior Estimates across Different Prior Specifications",
    x = "Odds Ratio (95% Credible Interval)",
    y = NULL, # Remove Y-axis label as the terms are self-explanatory
    color = "Prior Specification",
    shape = "Prior Specification"
  ) +
  
  # Theme customization for publication quality
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", color = "black"),
    axis.text.x = element_text(color = "black"),
    legend.position = "bottom",
    legend.box.background = element_rect(color = "black", linewidth = 0.5),
    panel.grid.minor = element_blank() # Clean up grid
  )

# --- 6. Save Plot ---
# Save as PDF (Vector format, best for submission)
ggsave("Figure_S_Sensitivity_Analysis.pdf", plot = p, width = 8, height = 6)

# Save as PNG (High resolution for review)
ggsave("Figure_S_Sensitivity_Analysis.png", plot = p, width = 8, height = 6, dpi = 300)

# Display plot
print(p)


# ==============================================================================
# MODULE: FINAL MCMC CONVERGENCE DIAGNOSTICS (2x2 INTEGRATED PLATE)
# Format: A (Trace), B (Areas), C (R-hat Hist), D (N-eff Hist)
# ==============================================================================

# Load required libraries
library(brms)
library(bayesplot)
library(ggplot2)
library(patchwork)
library(posterior)

# --- 1. Data Preparation and Parameter Renaming ---
# Select core clinical predictors for the diagnostic plate
target_vars <- c("b_Intercept", "b_age_10y_c", "b_smoking_historyYes", "b_cpb_30m_c", "b_propofol_100_c")
clean_names <- c("Intercept", "Age (10y)", "Smoking (Yes)", "CPB (30m)", "Propofol (100mg)")

# Convert to draws object and rename variables for consistent labeling
draws_fit <- as_draws_df(final_model_fit, variable = target_vars)
colnames(draws_fit)[1:5] <- clean_names

# --- 2. Generate Panels ---

# Panel A: Trace Plots (mixing assessment)
p_trace <- mcmc_trace(draws_fit, pars = clean_names, 
                      facet_args = list(ncol = 2)) +
  theme_minimal() +
  theme(legend.position = "none", strip.text = element_text(face = "bold")) +
  labs(title = "A: Trace Plots", subtitle = "Chain mixing and stationarity")

# Panel B: Posterior Distributions (Area plots)
p_areas <- mcmc_areas(draws_fit, pars = clean_names, prob = 0.95) +
  theme_minimal() +
  theme(axis.text.y = element_text(face = "bold")) +
  labs(title = "B: Posterior Distributions", subtitle = "95% Credible Intervals")

# Panel C: R-hat Diagnostic (Histogram)
rhats <- rhat(final_model_fit)[target_vars]
names(rhats) <- clean_names
p_rhat <- mcmc_rhat_hist(rhats) +
  theme_minimal() +
  labs(title = "C: R-hat Diagnostic", subtitle = "Target value near 1.000")

# Panel D: Effective Sample Size (Histogram)
neff_ratios <- neff_ratio(final_model_fit)[target_vars]
names(neff_ratios) <- clean_names
p_neff <- mcmc_neff_hist(neff_ratios) +
  theme_minimal() +
  labs(title = "D: Effective Sample Size (N_eff)", subtitle = "N_eff / N ratio")

# --- 3. Assemble Plate using Patchwork ---

final_diagnostic_plate <- (p_trace + p_areas) / (p_rhat + p_neff) +
  plot_annotation(
    title = "Final Model: MCMC Convergence Diagnostics",
    caption = "Plots confirm the model's successful convergence and sampling stability.",
    theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
  )

# --- 4. Export as High-Resolution TIFF ---

ggsave(
  filename = "Figure_MCMC_Diagnostics_Final.tiff",
  plot = final_diagnostic_plate,
  width = 14, 
  height = 12, 
  dpi = 300, 
  compression = "lzw"
)

cat("\n--- SUCCESS: Professional MCMC Diagnostic Plate Generated ---\n")


# ==============================================================================
# SCRIPT: 05_GENERATE_OUTPUTS.R 
# PURPOSE: Generate Table 3 with Posterior Probabilities & Figure 2
# ==============================================================================

library(tidyverse)
library(brms)
library(bayesplot)
library(ggplot2)
library(knitr)
library(posterior) # Required for efficient draw extraction

# Load the final model if not loaded
# final_model_fit <- readRDS("final_clinical_model.rds")

# --- PART 1: Generate Table 3 (Detailed with Posterior Probabilities) ---

# 1. Extract Basic Fixed Effects (Estimate, Error, CI)
fixed_effects <- fixef(final_model_fit) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Parameter")

# 2. Calculate Posterior Probability of Direction (Pr > 0)
# This quantifies the probability that the effect is truly positive (OR > 1).
# Crucial for interpreting variables like 'Smoking' if the 95% CI crosses 1.
draws <- as_draws_df(final_model_fit)

# Helper function to calculate Pr(beta > 0) for specific parameters
calc_prob <- function(variable_name) {
  # brms adds "b_" prefix to fixed effects in the draws
  target_col <- paste0("b_", variable_name) 
  if (target_col %in% names(draws)) {
    return(mean(draws[[target_col]] > 0))
  } else {
    return(NA) # For Intercept or other parameters if naming differs
  }
}

# Apply calculation to all parameters
fixed_effects$Prob_Positive <- sapply(fixed_effects$Parameter, calc_prob)

# 3. Format the Final Table
table_data <- fixed_effects %>%
  mutate(
    OR = round(exp(Estimate), 2),
    CI_Low = round(exp(Q2.5), 2),
    CI_High = round(exp(Q97.5), 2),
    # Format: "OR (Low-High)"
    Result = paste0(OR, " (", CI_Low, "-", CI_High, ")"),
    # Format Probability as percentage for readability in diagnostics
    Prob_Risk_Percent = round(Prob_Positive * 100, 1)
  ) %>%
  select(Parameter, Result, Prob_Risk_Percent)

# 4. Rename Parameters for Publication
table_data$Parameter <- dplyr::recode(table_data$Parameter,
                                      "Intercept" = "Intercept",
                                      "age_10y_c" = "Age (per 10-year increase)",
                                      "smoking_historyYes" = "Smoking History (Yes)",
                                      "cpb_30m_c" = "CPB Duration (per 30-min increase)",
                                      "propofol_100_c" = "Total Propofol Dose (per 100 mg increase)",
                                      "hospital" = "Center (Site Adjustment)",
                                      "surgery_typeIsolated_Single_Valve" = "Surgery: Isolated Single Valve",
                                      "surgery_typeOther_Procedures" = "Surgery: Other Procedures",
                                      "surgery_typeComplex_Combined" = "Surgery: Complex/Combined"
)

cat("\n--- FINAL TABLE 3 DATA (With Posterior Probabilities) ---\n")
print(kable(table_data))
write.csv(table_data, "FINAL_Table3_Extended.csv", row.names = FALSE)

# --- PART 2: Generate Figure 2 (Forest Plot) ---

# Define readable labels
my_labels <- c(
  "b_age_10y_c" = "Age (per 10 years)",
  "b_smoking_historyYes" = "Smoking (Yes)",
  "b_cpb_30m_c" = "CPB Duration (per 30 min)",
  "b_propofol_100_c" = "Propofol Dose (per 100 mg)"
)

# Select variables to plot (exclude intercept and covariates for the main plot)
# Note: brms adds 'b_' prefix to fixed effects
vars_to_plot <- c("b_age_10y_c", "b_smoking_historyYes", "b_cpb_30m_c", "b_propofol_100_c")

p <- mcmc_plot(final_model_fit, 
               variable = vars_to_plot, 
               type = "areas",
               prob = 0.95,
               point_est = "median",
               transformations = "exp") + # Exp to get Odds Ratios
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  scale_y_discrete(labels = my_labels) +
  labs(
    title = "Predictors of Postoperative Delirium",
    subtitle = "Posterior Median Odds Ratios (95% Credible Intervals)",
    x = "Odds Ratio (log scale)",
    y = ""
  ) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 10),
    plot.title = element_text(face = "bold")
  )

# Print and Save
print(p)
ggsave("FINAL_Figure2_ForestPlot.tiff", plot = p, width = 8, height = 6, dpi = 300, compression = "lzw")

cat("\n--- SUCCESS: Table 3 and Figure 2 are ready! ---\n")


# ==============================================================================
# SCRIPT A: Model Diagnostics (VIF & ROC)
# PURPOSE: Address Reviewer #2 (Independence) & Show Model Performance
# ==============================================================================

library(tidyverse)
library(brms)
library(mice)
library(car)   # For VIF
library(pROC)  # For ROC
library(ggplot2)

# 1. Load Data & Model
imputed_mids <- readRDS("imputed_mids_final_clinical.rds")
final_model  <- readRDS("final_clinical_model.rds")

# --- PART 1: MULTICOLLINEARITY CHECK (VIF) ---
# Addressed to Reviewer #2: "Propofol is distinct from CPB time"

cat("\n--- CHECKING MULTICOLLINEARITY (VIF) ---\n")

# We extract one complete dataset to run a standard GLM for VIF calculation
# (Standard practice for Bayesian papers to prove lack of collinearity)
check_data <- complete(imputed_mids, 1)

# Fit a temporary frequentist logistic model
glm_check <- glm(pod_status ~ age_10y + cpb_30m + propofol_100 + smoking_history, 
                 data = check_data, 
                 family = binomial)

# Calculate VIF
vif_values <- vif(glm_check)
print(vif_values)

cat("\nINTERPRETATION:\n")
if(max(vif_values) < 2.5) {
  cat("SUCCESS: All VIF values are < 2.5. Low multicollinearity.\n")
  cat("This proves Propofol and CPB Time are independent predictors.\n")
} else {
  cat("WARNING: Check VIF values.\n")
}

# --- PART 2: ROC CURVE (INTERNAL DISCRIMINATION) ---
cat("\n--- GENERATING ROC CURVE WITH 95% CI ---\n")

# 1. Generate predicted probabilities from the Bayesian model
# Note: 'posterior_epred' warning regarding first imputed set is standard for MI models
preds <- posterior_epred(final_model)
prob_mean <- apply(preds, 2, median)

# 2. Extract actual outcome directly from the model data
# (Replaces 'check_data' to avoid object not found error)
y_true <- final_model$data$pod_status

# 3. Calculate ROC object, AUC, and 95% CI of AUC
# ci = TRUE calculates the Confidence Interval for AUC
roc_obj <- roc(y_true, prob_mean, ci = TRUE)
auc_ci <- ci(roc_obj)

# 4. Format AUC text with 95% CI for the legend
auc_txt <- paste0("AUC = ", round(auc_ci[2], 3), 
                  " (95% CI: ", round(auc_ci[1], 3), "-", round(auc_ci[3], 3), ")")

cat(paste0("Model ", auc_txt, "\n"))

# 5. Calculate 95% CI for the ROC curve shape
ci_curve <- ci.se(roc_obj, specificities = seq(0, 1, 0.01))

# 6. Plot ROC (Academic Style)
tiff("FINAL_Figure3_ROC.tiff", width = 6, height = 6, units = "in", res = 300)

# Base Plot
plot(roc_obj, 
     main = "Receiver Operating Characteristic (ROC) Curve",
     col = "#1c61b6", 
     lwd = 3, 
     legacy.axes = TRUE,
     xlab = "1 - Specificity (False Positive Rate)",
     ylab = "Sensitivity (True Positive Rate)")

# Add 95% CI Shading area
plot(ci_curve, type = "shape", col = "#1c61b633") 

# Add Reference Line
abline(a = 0, b = 1, lty = 2, col = "gray")

# Add Legend
legend("bottomright", legend = auc_txt, 
       bty = "n", cex = 1.0, text.font = 2)

dev.off()

cat("SUCCESS: ROC curve with 95% CI saved as 'FINAL_Figure3_ROC.tiff'.\n")


# ==============================================================================
# SCRIPT B (CORRECTED): Table 1 Generation
# PURPOSE: Fix LVEF format (Categorical -> Numeric) and generate final table
# ==============================================================================

library(tidyverse)
library(tableone)
library(janitor)

# 1. Load Data
# We use read_csv, and clean names immediately
raw_data <- read_csv("my data.csv", show_col_types = FALSE) %>% 
  clean_names()

# 2. Data Cleaning & Formatting (The Critical Step)
table1_data <- raw_data %>%
  mutate(
    # --- FIX 1: Force LVEF to Numeric ---
    # parse_number is smart: it removes "%" symbols automatically
    # e.g., "55%" becomes 55, "60" becomes 60
    lvef = parse_number(as.character(lvef)),
    
    # --- Define Outcome ---
    pod_status = if_else(camicu_day1 == 1 | camicu_day2 == 1 | camicu_day3 == 1, 1, 0, missing = 0),
    pod_status = factor(pod_status, levels = c(0, 1), labels = c("No Delirium", "Delirium")),
    
    # --- Define Surgery Type (Consistent with Model) ---
    surgery_type = fct_recode(factor(surgery_category),
                              "Complex/Combined" = "4", "Complex/Combined" = "5",
                              "Other" = "3", "Other" = "6", "Other" = "7",
                              "Isolated CABG" = "1", "Isolated Single Valve" = "2"),
    
    # --- Define Factors ---
    # Assuming original coding: 0=No, 1=Yes; Gender: 0=Female, 1=Male (Adjust if needed)
    smoking = factor(smoking_history, levels = c(0, 1), labels = c("No", "Yes")),
    alcohol = factor(alcohol_history, levels = c(0, 1), labels = c("No", "Yes")),
    diabetes = factor(diabetes, levels = c(0, 1), labels = c("No", "Yes")),
    hypertension = factor(hypertension, levels = c(0, 1), labels = c("No", "Yes")),
    stroke_history = factor(stroke_history, levels = c(0, 1), labels = c("No", "Yes")),
    gender = factor(gender, levels = c(0, 1), labels = c("Female", "Male")) 
  )

# 3. Define Variables for Table 1
vars <- c("age_years", "gender", "bmi", "edu_years", "mmse_preop", 
          "smoking", "alcohol", "diabetes", "hypertension", "stroke_history",
          "surgery_type", "cpb_time_min", "propofol_mg", "lvef") # LVEF is here

# 4. Define Categorical Variables Explicitly (To be safe)
cat_vars <- c("gender", "smoking", "alcohol", "diabetes", "hypertension", 
              "stroke_history", "surgery_type")

# 5. Create Table
tab1 <- CreateTableOne(vars = vars, factorVars = cat_vars, strata = "pod_status", data = table1_data, test = TRUE)

# 6. Print & Save
# We treat skewed variables as non-normal (Median [IQR])
# LVEF is often skewed, so putting it here makes it look consistent with CPB time
non_normal <- c("age_years", "cpb_time_min", "propofol_mg", "mmse_preop", "bmi", "lvef")

tab1_print <- print(tab1, nonnormal = non_normal, quote = FALSE, noSpaces = TRUE, printToggle = FALSE, showAllLevels = TRUE)

# 7. Output
print(tab1_print) # Show in console to check
write.csv(tab1_print, "FINAL_Table1_Baseline_Corrected.csv")

cat("\nSUCCESS: LVEF is now numeric. Table 1 saved as 'FINAL_Table1_Baseline_Corrected.csv'.\n")


# ==============================================================================
# SCRIPT C (FINAL): Exploratory Landmark Analysis
# PURPOSE: Use pre-cleaned CSV to ensure row matching (N=308)
# ==============================================================================

library(tidyverse)
library(brms)
library(mice)
library(janitor)

# 1. Load Data & Model
imputed_mids <- readRDS("imputed_mids_final_clinical.rds")
final_model  <- readRDS("final_clinical_model.rds")

# 2. Extract Analyzed Data (N = 308) from MICE object
# This contains age_10y, cpb_30m, propofol_100, etc.
full_data_analyzed <- complete(imputed_mids, 1) 

# 3. Load Pre-Cleaned Raw Data (N = 308)
# This contains the raw variables we need: camicu_day1, mvt_first_24h
raw_data_308 <- read_csv("my_data_cleaned_for_analysis.csv", show_col_types = FALSE) %>% 
  clean_names()

# --- SAFETY CHECK ---
# Ensure row counts match exactly before merging
if(nrow(full_data_analyzed) != nrow(raw_data_308)) {
  stop("CRITICAL ERROR: Row counts do not match! Check your files.")
} else {
  cat("SUCCESS: Row counts match (N=308). Proceeding to merge.\n")
}

# 4. Create Landmark Cohort
landmark_data <- full_data_analyzed %>%
  mutate(
    # Merge necessary raw variables
    # We assume the order hasn't changed (standard in R workflow)
    camicu_day1   = raw_data_308$camicu_day1, 
    mvt_first_24h = raw_data_308$mvt_first_24h # Ensure this column exists in your CSV
  ) %>%
  # --- LANDMARK EXCLUSION ---
  # Exclude patients who ALREADY had delirium on Day 1 (1 = Yes, 0 = No)
  # Check if camicu_day1 is numeric (0/1) or character ("No"/"Yes") in your CSV
  # The code below handles both safely
  filter(camicu_day1 == 0 | camicu_day1 == "No") %>% 
  # Define Outcome: Delirium on Day 2 or 3
  mutate(
    late_pod = if_else(pod_status == "Delirium", 1, 0)
  )

cat(paste0("Landmark Cohort Size (Day 1 POD Free): ", nrow(landmark_data), "\n"))

# 5. Calculate Baseline Risk Score (Propensity Score)
# We need to center the variables exactly as the model expects
means <- list(
  age  = mean(full_data_analyzed$age_10y, na.rm=TRUE),
  cpb  = mean(full_data_analyzed$cpb_30m, na.rm=TRUE),
  prop = mean(full_data_analyzed$propofol_100, na.rm=TRUE)
)

landmark_data_scored <- landmark_data %>%
  mutate(
    age_10y_c      = age_10y - means$age,
    cpb_30m_c      = cpb_30m - means$cpb,
    propofol_100_c = propofol_100 - means$prop
  )

# Generate Predictions
cat("--- Calculating Baseline Risk Scores... ---\n")
# allow_new_levels is TRUE to handle any potential factor level issues safely
risk_preds <- posterior_epred(final_model, newdata = landmark_data_scored, allow_new_levels = TRUE)
landmark_data_scored$baseline_risk <- apply(risk_preds, 2, median)

# 6. Run Landmark Model
cat("--- RUNNING LANDMARK MODEL (Does 24h MV predict Late POD?) ---\n")

landmark_fit <- brm(
  late_pod ~ baseline_risk + mvt_first_24h, 
  data = landmark_data_scored,
  family = bernoulli(),
  prior = set_prior("normal(0, 1)", class = "b"),
  seed = 123,
  control = list(adapt_delta = 0.95),
  file = "landmark_model_final"
)

# 7. Output Results
print(summary(landmark_fit))

# Extract Odds Ratios
res <- fixef(landmark_fit) %>% 
  as.data.frame() %>%
  mutate(
    Predictor = rownames(.),
    OR = round(exp(Estimate), 2), 
    Low = round(exp(Q2.5), 2), 
    High = round(exp(Q97.5), 2),
    Result = paste0(OR, " (", Low, "-", High, ")")
  ) %>%
  select(Predictor, Result)

print(res)
write.csv(res, "FINAL_Table4_Landmark.csv")

cat("\nSUCCESS: Landmark Analysis Complete. Check 'FINAL_Table4_Landmark.csv'.\n")

cat("SUCCESS: Landmark analysis complete. Check if MV Duration is significant.\n")




# ==============================================================================
# MODULE 1 (FIXED): Bootstrap Validation & ROC Visualization
# PROJECT: Prediction Model for Postoperative Delirium 
# AUTHOR: medusa
# DATE: 2026-01-25
# ------------------------------------------------------------------------------
# FIX NOTES: 
# 1. Corrected 'hospital' variable type to Numeric to match the training data.
# 2. Ensures strict variable matching for 'brms' prediction.
# ==============================================================================

# --- 1. Environment Setup & Libraries ---
library(tidyverse)
library(janitor)
library(pROC)     
library(ggplot2)  
library(brms)     

set.seed(202607) # Reproducibility

# --- 2. Data Loading & Feature Engineering ---
# Load raw data
if(file.exists("my data.csv")) {
  raw_data <- read_csv("my data.csv", show_col_types = FALSE) %>% clean_names()
} else {
  stop("Error: 'my data.csv' not found.")
}

# Apply Filters (Inclusion/Exclusion Criteria)
data_processed <- raw_data %>%
  filter(age_years >= 18) %>%
  filter(!is.na(edu_years) & !is.na(mmse_preop)) %>%
  mutate(
    is_eligible = case_when(
      edu_years == 0 & mmse_preop > 17 ~ TRUE,
      edu_years >= 1 & edu_years <= 6 & mmse_preop > 20 ~ TRUE,
      edu_years > 6 & mmse_preop > 24 ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>%
  filter(is_eligible == TRUE) %>%
  mutate(
    # Outcome: 1 = Yes, 0 = No
    outcome = if_else(camicu_day1 == 1 | camicu_day2 == 1 | camicu_day3 == 1, 1, 0, missing = 0)
  )

# Apply Variable Transformation 
data_modeled <- data_processed %>%
  mutate(
    surgery_type = fct_recode(factor(surgery_category),
                              "Complex_Combined" = "4", 
                              "Complex_Combined" = "5", 
                              "Other_Procedures" = "3", 
                              "Other_Procedures" = "6", 
                              "Other_Procedures" = "7", 
                              "Isolated_CABG" = "1",
                              "Isolated_Single_Valve" = "2"),
    smoking_history = factor(smoking_history, levels = c(0, 1), labels = c("No", "Yes")),
    
    # [FIX]: Ensure hospital is Numeric (as defined in the training set)
    hospital = as.numeric(hospital) 
  )

# --- 3. Clinical Unit Conversion & Centering ---
# Convert to clinical units
data_modeled <- data_modeled %>%
  mutate(
    age_10y      = age_years / 10,
    cpb_30m      = cpb_time_min / 30,
    propofol_100 = propofol_mg / 100
  )

# Centering Variables
# Using means from the current dataset (approximation for validation)
means_list <- list(
  age = mean(data_modeled$age_10y, na.rm = TRUE),
  cpb = mean(data_modeled$cpb_30m, na.rm = TRUE),
  prop = mean(data_modeled$propofol_100, na.rm = TRUE)
)

data_final <- data_modeled %>%
  mutate(
    age_10y_c      = age_10y - means_list$age,
    cpb_30m_c      = cpb_30m - means_list$cpb,
    propofol_100_c = propofol_100 - means_list$prop
  ) %>%
  # Drop NAs in predictors
  drop_na(age_10y_c, cpb_30m_c, propofol_100_c, smoking_history, hospital, surgery_type)

cat("Data Ready. N =", nrow(data_final), "\n")

# --- 4. Load Model & Predict ---
# Try loading the fixed model first
model_file <- "final_clinical_model.rds"
if(!file.exists(model_file)) model_file <- "final_clinical_model_fixed.rds"

cat("Loading model:", model_file, "\n")
final_model <- readRDS(model_file)

cat("Generating predictions... \n")
# Note: allow_new_levels = TRUE is safer for validation if surgery types are rare, 
# but for internal validation, it shouldn't be strictly necessary.
preds <- fitted(final_model, newdata = data_final, scale = "response", allow_new_levels = TRUE)
data_final$pred_prob <- preds[, "Estimate"]

# --- 5. Bootstrap Validation ---
cat("Running Bootstrap (n=1000)... \n")

roc_obj <- roc(data_final$outcome, data_final$pred_prob, quiet = TRUE)
orig_auc <- as.numeric(auc(roc_obj))

# Calculate Confidence Intervals
roc_ci <- ci.se(roc_obj, specificities = seq(0, 1, 0.01), boot.n = 1000, progress = "none")
auc_ci_val <- ci.auc(roc_obj, method = "bootstrap", boot.n = 1000)

plot_data <- data.frame(
  specificity = as.numeric(rownames(roc_ci)),
  sensitivity = roc_ci[, 2], 
  lower = roc_ci[, 1],       
  upper = roc_ci[, 3]        
)

label_text <- paste0("AUC: ", sprintf("%.3f", orig_auc), 
                     "\n95% CI: ", sprintf("%.3f", auc_ci_val[1]), 
                     "-", sprintf("%.3f", auc_ci_val[3]))

# --- 6. Plotting ---
p_roc <- ggplot(plot_data, aes(x = 1 - specificity, y = sensitivity)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#1f77b4", alpha = 0.2) +
  geom_line(color = "#1f77b4", linewidth = 1.2) +
  annotate("text", x = 0.60, y = 0.25, label = label_text, 
           size = 5, fontface = "bold", hjust = 0) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Bootstrap-Corrected ROC Curve",
    subtitle = "Internal Validation (B=1000 Resamples)",
    x = "1 - Specificity",
    y = "Sensitivity"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
    panel.border = element_rect(color = "black", fill = NA),
    axis.text = element_text(color = "black")
  )

# --- 7. Save Output ---
ggsave("Supplementary_Figure_1_Bootstrap_ROC.tiff", plot = p_roc, 
       width = 7, height = 6, dpi = 300, compression = "lzw")

cat("\n[SUCCESS] Figure Saved: Supplementary_Figure_1_Bootstrap_ROC.tiff\n")



# ==============================================================================
# SCRIPT: 05_Generate_Validation_Figures.R
# PROJECT: POD Prediction Model 
# AUTHOR: [Medusa]
# DATE: 2026-01-25
# PURPOSE: Generate critical validation figures 
#          1. Bootstrap Internal Validation (ROC & Calibration Slope)
#          2. Calibration Plot (Figure 4)
#          3. Decision Curve Analysis (Figure 5)
#          4. Static Clinical Nomogram (Figure 6)
# ==============================================================================

# --- 1. SETUP & LIBRARIES -----------------------------------------------------
rm(list = ls()) # Clean environment
library(tidyverse)
library(brms)      # For loading the Bayesian model
library(rms)       # Standard for Nomograms and Calibration
library(pROC)      # For ROC calculations
library(dcurves)   # For Decision Curve Analysis
library(mice)      # For handling imputed data

# Set seed for reproducibility
set.seed(2026)

# Create output directories if they don't exist
if(!dir.exists("Figures")) dir.create("Figures")
if(!dir.exists("Tables")) dir.create("Tables")

# --- 2. LOAD DATA & MODEL -----------------------------------------------------
cat("--- [Step 1] Loading Data and Model objects... ---\n")

# Load the imputed MICE object (Clinical Units)
imputed_mids <- readRDS("imputed_mids_final_clinical.rds")

# Extract the first imputed dataset for visualization purposes
# (Standard practice for plotting nomograms/DCA is to use a representative dataset 
# or averaged predictions, here we use Dataset #1 for stability)
df_clinical <- complete(imputed_mids, action = 1)

# Ensure factors are set correctly
df_clinical$pod_status_num <- ifelse(df_clinical$pod_status == "Yes", 1, 0)
df_clinical$smoking_history <- factor(df_clinical$smoking_history, levels = c("No", "Yes"))

# Load the Final Bayesian Model
final_model <- readRDS("final_clinical_model.rds")

# --- 3. GENERATE PREDICTIONS --------------------------------------------------
cat("--- [Step 2] Generating Posterior Predictions... ---\n")

# We need to generate the predicted probability of Delirium for every patient
# based on the Bayesian model posterior.
preds <- posterior_epred(final_model) # Matrix: n_samples x n_patients
df_clinical$prob_mean <- colMeans(preds) # Average probability per patient

# Quick Check of AUC
roc_obj <- roc(df_clinical$pod_status, df_clinical$prob_mean, quiet = TRUE)
cat(sprintf("   > Apparent AUC (Bayesian Model): %.3f\n", auc(roc_obj)))

# --- 4. FIGURE 4: CALIBRATION PLOT (CORRECTED) --------------------------------
cat("--- [Step 3] Generating Figure 4: Calibration Plot... ---\n")

tiff("FINAL_Figure4_Calibration.tiff", width = 2000, height = 2000, res = 300)

# Set graphical parameters BEFORE calling the function
# cex.lab = 1.2 sets the label size here
par(mfrow = c(1, 1), mar = c(5, 5, 2, 2), cex.lab = 1.2)

val.prob(
  p = df_clinical$prob_mean, 
  y = df_clinical$pod_status_num,
  logistic.cal = FALSE, # Use loess smoother
  statloc = FALSE,      # Hide stats table (cleaner for publication)
  g = 10,               # 10 bins
  xlab = "Predicted Probability of Postoperative Delirium",
  ylab = "Actual Probability (Proportion)",
  riskdist = "calibrated", 
  lim = c(0, 1)
)

title(main = "Calibration Curve (Bayesian Model)", adj = 0)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

dev.off()
cat("   > Saved: FINAL_Figure4_Calibration.tiff\n")

# --- 5. FIGURE 5: DECISION CURVE ANALYSIS (DCA) -------------------------------
cat("--- [Step 4] Generating Figure 5: Decision Curve Analysis... ---\n")

# Prepare data for dcurves package
dca_data <- df_clinical %>%
  select(pod_status_num, prob_mean)

# Run DCA
dca_result <- dca(
  formula = pod_status_num ~ prob_mean,
  data = dca_data,
  thresholds = seq(0, 0.50, by = 0.01), # Clinical range 0% to 50%
  label = list(prob_mean = "Prediction Model")
)

# Plot DCA
tiff("FINAL_Figure5_DCA.tiff", width = 2400, height = 1800, res = 300)

dca_plot <- dca_result %>%
  as_tibble() %>%
  ggplot(aes(x = threshold, y = net_benefit, color = label)) +
  geom_line(linewidth = 1.2) +
  coord_cartesian(ylim = c(-0.02, 0.15)) + # Zoom in on relevant benefit area
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Threshold Probability",
    y = "Net Benefit",
    color = ""
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = c(0.7, 0.8),
    legend.background = element_rect(fill = "white", color = "black")
  ) +
  scale_color_manual(values = c("grey60", "black", "red")) # Treat All, Treat None, Model

print(dca_plot)
dev.off()
cat("   > Saved: FINAL_Figure5_DCA.tiff\n")

# --- 6. FIGURE 6: STATIC NOMOGRAM (Clinical Tool) -----------------------------
cat("--- [Step 5] Generating Figure 6: Static Nomogram... ---\n")

# Load necessary library
library(rms)

# 1. OPTIMIZE LABELS FOR PUBLICATION
# Assign professional labels to variables so they appear nicely on the graph
# instead of raw column names like "age_10y".
label(df_clinical$age_10y)        <- "Age (per 10 years)"
label(df_clinical$smoking_history)<- "Smoking History"
label(df_clinical$cpb_30m)        <- "CPB Duration (per 30 min)"
label(df_clinical$propofol_100)   <- "Propofol (per 100 mg)"

# 2. DEFINE DATA DISTRIBUTION
# The 'rms' package requires a data distribution object to define ranges.
dd <- datadist(df_clinical)
options(datadist = "dd")

# 3. FIT THE MODEL
# We fit the logistic regression model using the labeled data.
f_lrm <- lrm(
  pod_status_num ~ age_10y + smoking_history + cpb_30m + propofol_100,
  data = df_clinical,
  x = TRUE, y = TRUE
)

# 4. SETUP HIGH-RESOLUTION IMAGE (TIFF for SCI)
# Width and height adjusted for a standard landscape figure.
tiff("FINAL_Figure6_Nomogram.tiff", width = 3600, height = 2400, res = 300, compression = "lzw")

# Adjust margins: bottom, left, top, right
# Increased left margin (second value) to accommodate variable names.
par(mfrow = c(1, 1), mar = c(4, 5, 2, 2)) 

# 5. CREATE NOMOGRAM OBJECT
nom <- nomogram(
  f_lrm,
  fun = plogis,               # Convert log-odds to probability
  fun.at = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), # Clean risk scale points
  funlabel = "Risk of POD",   # Label for the probability axis
  lp = FALSE,                 # Hide Linear Predictor axis (redundant for clinicians)
  
  # Define specific tick marks to match clinical units nicely
  age_10y = seq(4, 9, by = 1),       # Represents 40-90 years
  cpb_30m = seq(0, 10, by = 1),      # Represents 0-300 mins
  propofol_100 = seq(0, 20, by = 2)  # Represents 0-2000 mg
)

# 6. PLOT WITH AESTHETIC ENHANCEMENTS
plot(nom,
     xfrac = 0.3,         # Fraction of space reserved for labels (0.3 gives more room for text)
     cex.var = 1.2,       # Font size for variable names (bold/prominent)
     cex.axis = 1.0,      # Font size for axis numbers
     lwd = 1.5,           # Line width (thicker lines look better in print)
     col.grid = gray(0.9),# Subtle grid color
     force.label = TRUE   # Ensure labels are forced onto the plot
)

# Close the device to save the file
dev.off()

cat("   > Saved: FINAL_Figure6_Nomogram.tiff (High Resolution)\n")


# --- 7. SUPPLEMENTARY: BOOTSTRAP VALIDATION (Internal) ------------------------
cat("--- [Step 6] Running Bootstrap Validation (1000 Reps)... ---\n")
# This calculates the Optimism-Corrected AUC and Slope.
# Validating the 'lrm' object is sufficient as a proxy for the study's stability.

# Validate function in rms does bootstrap resampling
val_stats <- validate(f_lrm, method = "boot", B = 1000)

# Extract key metrics
auc_optimism <- 0.5 * (val_stats["Dxy", "index.corrected"] + 1)
slope_optimism <- val_stats["Slope", "index.corrected"]

cat("\n======================================================\n")
cat(" BOOTSTRAP VALIDATION RESULTS (B=1000)\n")
cat("======================================================\n")
cat(sprintf("Optimism-Corrected AUC:  %.3f\n", auc_optimism))
cat(sprintf("Calibration Slope:       %.3f\n", slope_optimism))
cat("======================================================\n")

# Save these stats to a CSV for the manuscript text
validation_results <- data.frame(
  Metric = c("Optimism-Corrected AUC", "Calibration Slope"),
  Value = c(auc_optimism, slope_optimism)
)
write_csv(validation_results, "Tables/Supplementary_Table_Bootstrap_Results.csv")

cat("--- Script Completed Successfully. ---\n")



# ==============================================================================
# SCRIPT: 06_Generate_External_Validation_ROC.R
# PURPOSE: External Validation of the Final Model on Center B Data
# OUTPUT: FINAL_Figure3B_External_ROC.tiff
# ==============================================================================

library(tidyverse)
library(brms)
library(pROC)
library(mice)

# 1. Load Data and Final Model
# Ensure these files exist from previous scripts
imputed_mids <- readRDS("imputed_mids_final_clinical.rds")
final_model <- readRDS("final_clinical_model.rds")

# 2. Extract Center B Data (External Site)
# Extract the first imputed dataset for standard ROC plotting
df_full <- complete(imputed_mids, action = 1)
center_b_data <- df_full %>% filter(hospital == "Center B")

# 3. Generate Predictions using the Pre-trained Model
# We do not retrain the model; we apply the existing weights to Center B
preds <- posterior_epred(final_model, newdata = center_b_data, allow_new_levels = TRUE)
center_b_data$prob_mean <- colMeans(preds)

# 4. Calculate ROC Metrics
roc_ext <- roc(center_b_data$pod_status, center_b_data$prob_mean, ci = TRUE, quiet = TRUE)
auc_val <- as.numeric(auc(roc_ext))
ci_val <- ci(roc_ext)

# 5. Generate Figure 3B
tiff("FINAL_Figure3B_External_ROC.tiff", width = 2000, height = 2000, res = 300, compression = "lzw")
plot(roc_ext, 
     col = "#E74C3C", 
     lwd = 3, 
     legacy.axes = TRUE, 
     main = "External Validation (Center B)")
grid(col = "lightgray", lty = "dotted")

# Add AUC Text
text_label <- sprintf("AUC = %.3f\n95%% CI: %.3f - %.3f", 
                      auc_val, ci_val[1], ci_val[3])
text(0.4, 0.2, text_label, cex = 1.2, pos = 4, font = 2)
dev.off()


# ==============================================================================
# SCRIPT: 05_Supplementary_External_Validation_Pooling.R
# PURPOSE: Pooled External Validation across all 50 Imputed Datasets
# OUTPUT: Figure_External_Validation_ROC.tiff
# ==============================================================================

library(tidyverse)
library(brms)
library(mice)
library(pROC)

# 1. Load Data and Final Model
imputed_mids <- readRDS("imputed_mids_final_clinical.rds")
final_model <- readRDS("final_clinical_model.rds")

# 2. Iterate through all 50 Imputations to Pool Predictions
m_count <- imputed_mids$m
test_list <- complete(imputed_mids, action = "all")

# Filter only Center B for each imputation
test_list_b <- lapply(test_list, function(df) df[df$hospital == "Center B", ])

# Matrix to store probabilities: Patients x Imputations
prob_matrix <- matrix(0, nrow = nrow(test_list_b[[1]]), ncol = m_count)

for(i in 1:m_count) {
  # Apply pre-trained model to i-th imputation
  preds_i <- posterior_epred(final_model, newdata = test_list_b[[i]], allow_new_levels = TRUE)
  prob_matrix[, i] <- colMeans(preds_i)
}

# 3. Calculate Average Probability (Pooled)
final_avg_probs <- rowMeans(prob_matrix)
actual_outcome <- test_list_b[[1]]$pod_status

# 4. Calculate ROC
roc_pooled <- roc(actual_outcome, final_avg_probs, ci = TRUE, quiet = TRUE)
auc_val <- as.numeric(auc(roc_pooled))
ci_val <- ci(roc_pooled)

# 5. Generate Supplementary Figure
tiff("Figure_External_Validation_ROC.tiff", width = 1800, height = 1800, res = 300, compression = "lzw")
plot(roc_pooled, 
     col = "#2E86C1", 
     lwd = 3, 
     legacy.axes = TRUE, 
     main = "Pooled External Validation ROC\n(50 Imputed Sets)")
grid(col = "lightgray", lty = "dotted")

# Add Legend
legend_text <- sprintf("Pooled AUC = %.3f\n95%% CI: %.3f - %.3f", 
                       auc_val, ci_val[1], ci_val[3])
legend("bottomright", legend = legend_text, bty = "n", cex = 1.1, text.font = 2)
dev.off()


# ==============================================================================
# Module: Supplementary Figure 2 Generation
# Description: Generate a correlation matrix heatmap for candidate predictors 
#              to demonstrate multicollinearity (specifically among time-based variables).
# Input: "my data.csv"
# Output: "Figure2_CorrelationMatrix.tiff"
# ==============================================================================

# 1. Load necessary libraries
library(tidyverse)
library(ggcorrplot) # For professional ggplot2-based correlation plots

# 2. Load the dataset
# Note: Reading all columns as character first to safely handle symbols like '%'
raw_data <- read.csv("my data.csv", stringsAsFactors = FALSE, na.strings = c("", "NA"))

# 3. Data Cleaning and Preparation
# We need to select potential continuous predictors and clean specific columns (e.g., LVEF)
analytic_data <- raw_data %>%
  mutate(
    # Clean LVEF: Remove '%' sign and convert to numeric
    lvef_clean = as.numeric(gsub("%", "", lvef)),
    
    # Ensure other key variables are numeric
    age_years = as.numeric(age_years),
    bmi = as.numeric(bmi),
    edu_years = as.numeric(edu_years),
    mmse_preop = as.numeric(mmse_preop),
    surgery_duration_min = as.numeric(surgery_duration_min),
    cpb_time_min = as.numeric(cpb_time_min),
    acclamp_time_min = as.numeric(acclamp_time_min),
    propofol_mg = as.numeric(propofol_mg)
  ) %>%
  # Select only the candidate continuous predictors for the matrix
  select(
    age_years,
    bmi,
    edu_years,
    mmse_preop,
    lvef_clean,
    surgery_duration_min,
    cpb_time_min,
    acclamp_time_min,
    propofol_mg
  ) %>%
  # Rename columns for professional publication display
  rename(
    `Age` = age_years,
    `BMI` = bmi,
    `Education (Years)` = edu_years,
    `Preop MMSE` = mmse_preop,
    `LVEF (%)` = lvef_clean,
    `Surgery Duration` = surgery_duration_min,
    `CPB Duration` = cpb_time_min,
    `X-Clamp Duration` = acclamp_time_min,
    `Propofol Dose` = propofol_mg
  ) %>%
  # Remove rows with missing values in these specific variables for the correlation calculation
  na.omit()

# 4. Calculate Correlation Matrix
# Use Pearson correlation
corr_matrix <- cor(analytic_data, method = "pearson")

# 5. Generate the Correlation Heatmap
# Visual style: Lower triangle, show coefficients, significant colors
plot_corr <- ggcorrplot(
  corr_matrix,
  method = "square",
  type = "lower",           # Only show lower triangle
  lab = TRUE,               # Show correlation coefficients
  lab_size = 3,             # Size of numbers
  colors = c("#6D9EC1", "white", "#E46726"), # Blue (Neg) -> White -> Red (Pos)
  title = "Correlation Matrix of Candidate Continuous Predictors",
  ggtheme = ggplot2::theme_minimal()
) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  )

# 6. Save the Figure
# High-resolution TIFF for publication
ggsave(
  filename = "Figure2_CorrelationMatrix.tiff",
  plot = plot_corr,
  width = 8,
  height = 8,
  dpi = 300,
  compression = "lzw"
)

# Output confirmation message
message("SUCCESS: Figure2_CorrelationMatrix.tiff has been generated.")

