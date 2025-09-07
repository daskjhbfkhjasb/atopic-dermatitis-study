################### Cox Regression Analysis for Cancer Risk ##################################
# This R script performs a multivariable Cox proportional hazards regression analysis 
# to assess the association between Atopic Dermatitis (AD) and various cancer types 
# using UK Biobank data.
# 
# Analysis includes:
#   1. Overall population analysis
#   2. Gender-stratified analysis
#
# Excluded Cancer Types:
#   - C27, C28, C29, C35, C36, C59, C87, C89
#
# Confounders Adjusted For:
#   1. Age
#   2. Sex
#   3. Smoking status
#   4. Body type distribution
#   5. Townsend deprivation index (TDI) group
#   6. Diastolic blood pressure (DBP)
#
# Output:
#   - Hazard Ratios (HR)
#   - 95% Confidence Intervals (CI)
#   - p-values
#
# Required package: survival (for Cox regression modeling)
#
# Results will be saved as CSV files.
# 
# By: [Zhengdong Wu; Qiaolin Wang]
# Date: [2025-09-07]
##########################################################################################


library(survival)

# Set R to not display scientific notation (only affects display, not actual storage)
options(scipen = 999)

# Define cancer types to be excluded
excluded_cancers <- c("C27", "C28", "C29", "C35", "C36", "C59", "C87", "C89")

# Get all cancer type column names (C00 - C97)
cancer_cols <- grep("^C[0-9]{2}$", colnames(data_cleaned4_sub), value = TRUE)

# Remove excluded cancer types from the list
cancer_cols <- setdiff(cancer_cols, excluded_cancers)

# Create empty list to store results
results <- list()

# Multivariable Cox regression loop
for (cancer in cancer_cols) {
  tryCatch({
    # Build survival object
    surv_obj <- with(data_cleaned4_sub, Surv(FollowUpTime, get(cancer)))
    
    # Build Cox model formula
    cox_formula <- as.formula(paste(
      "surv_obj ~ AD + Age + Sex + Smoking.status + Body.Type.Distribution + TDI_group + DBP"
    ))
    
    # Fit Cox model
    fit <- coxph(cox_formula, data = data_cleaned4_sub)
    
    # Extract results for AD variable
    summary_fit <- summary(fit)
    ad_index <- which(rownames(summary_fit$coefficients) == "AD")
    
    # Extract key metrics
    HR <- exp(coef(fit))["AD"]
    CI <- exp(confint(fit))["AD", ]
    p_value <- summary_fit$coefficients["AD", "Pr(>|z|)"]
    
    # Handle very small P values
    if (!is.na(p_value) && p_value < 0.0001) {
      p_display <- "<0.0001"
    } else {
      p_display <- format(round(p_value, 4), nsmall = 4)
    }
    
    # Store results
    results[[cancer]] <- data.frame(
      CancerType = cancer,
      HR = round(HR, 4),
      CI_lower = round(CI[1], 4),
      CI_upper = round(CI[2], 4),
      p_value = p_display  # Use character format to avoid scientific notation
    )
  }, error = function(e) {
    message(sprintf("Failed to fit model for %s: %s", cancer, e$message))
    # Return row with error marker when fitting fails
    results[[cancer]] <- data.frame(
      CancerType = cancer,
      HR = NA,
      CI_lower = NA,
      CI_upper = NA,
      p_value = "Fitting failed"
    )
  })
}

# Combine all results
final_results <- do.call(rbind, results)
rownames(final_results) <- NULL

# Display results - ensure no scientific notation format is used
print(final_results)

# Save results avoiding scientific notation
write.csv(final_results, "cox_regression_results.csv", row.names = FALSE)

# Gender-stratified Cox regression analysis
library(survival)

# Set R to not display scientific notation
options(scipen = 999)

# Define cancer types to be excluded
excluded_cancers <- c("C27", "C28", "C29", "C35", "C36", "C59", "C87", "C89")

# Get all cancer type column names (C00 - C97)
cancer_cols <- grep("^C[0-9]{2}$", colnames(data_cleaned4_sub), value = TRUE)
cancer_cols <- setdiff(cancer_cols, excluded_cancers)

# Check gender variable format (ensure correct coding)
if (!"Sex" %in% colnames(data_cleaned4_sub)) stop("Sex variable not found in dataset")
if (is.numeric(data_cleaned4_sub$Sex)) {
  message("Note: Sex variable is numeric, ensure 0/1 represents female/male respectively")
} else if (is.factor(data_cleaned4_sub$Sex)) {
  message("Sex variable is factor type")
} else {
  data_cleaned4_sub$Sex <- as.factor(data_cleaned4_sub$Sex)
}

# Create empty lists to store male and female results
male_results <- list()
female_results <- list()

# Helper function: perform stratified analysis
analyze_subgroup <- function(subgroup_name, subgroup_data) {
  results_list <- list()
  
  for (cancer in cancer_cols) {
    tryCatch({
      # Build survival object
      surv_obj <- with(subgroup_data, Surv(FollowUpTime, get(cancer)))
      
      # Build Cox model formula (exclude sex variable)
      cox_formula <- as.formula(paste(
        "surv_obj ~ AD + Age + Smoking.status + Body.Type.Distribution + TDI_group + DBP"
      ))
      
      # Fit Cox model
      fit <- coxph(cox_formula, data = subgroup_data)
      
      # Extract results for AD variable
      HR <- exp(coef(fit))["AD"]
      CI <- exp(confint(fit))["AD", ]
      p_value <- summary(fit)$coefficients["AD", "Pr(>|z|)"]
      
      # Handle very small P values
      if (!is.na(p_value) && p_value < 0.0001) {
        p_display <- "<0.0001"
      } else {
        p_display <- format(round(p_value, 4), nsmall = 4)
      }
      
      # Store results
      results_list[[cancer]] <- data.frame(
        CancerType = cancer,
        HR = round(HR, 4),
        CI_lower = round(CI[1], 4),
        CI_upper = round(CI[2], 4),
        p_value = p_display,
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      message(sprintf("Failed to fit %s in %s group: %s", cancer, subgroup_name, e$message))
      results_list[[cancer]] <- data.frame(
        CancerType = cancer,
        HR = NA,
        CI_lower = NA,
        CI_upper = NA,
        p_value = "Fitting failed",
        stringsAsFactors = FALSE
      )
    })
  }
  
  return(do.call(rbind, results_list))
}

# Gender-stratified analysis
if ("Male" %in% unique(data_cleaned4_sub$Sex) | 1 %in% unique(data_cleaned4_sub$Sex)) {
  male_data <- subset(data_cleaned4_sub, Sex %in% c("Male", 1))
  male_final_results <- analyze_subgroup("Male", male_data)
  message("Male stratified analysis completed!")
} else {
  warning("Male data not found in dataset, please check Sex variable coding")
}

if ("Female" %in% unique(data_cleaned4_sub$Sex) | 0 %in% unique(data_cleaned4_sub$Sex)) {
  female_data <- subset(data_cleaned4_sub, Sex %in% c("Female", 0))
  female_final_results <- analyze_subgroup("Female", female_data)
  message("Female stratified analysis completed!")
} else {
  warning("Female data not found in dataset, please check Sex variable coding")
}

# Add gender identifier to results tables
if (exists("male_final_results")) {
  male_final_results$Gender <- "Male"
}
if (exists("female_final_results")) {
  female_final_results$Gender <- "Female"
}

# Combine results (optional)
if (exists("male_final_results") && exists("female_final_results")) {
  stratified_results <- rbind(male_final_results, female_final_results)
} else {
  warning("Unable to obtain both male and female results simultaneously, please use male or female results tables separately")
}

# Display environment objects
ls()

# Save male and female results separately
write.csv(male_final_results, "cox_male_results.csv", row.names = FALSE)
write.csv(female_final_results, "cox_female_results.csv", row.names = FALSE)