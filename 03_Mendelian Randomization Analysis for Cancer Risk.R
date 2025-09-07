################### Mendelian Randomization Analysis for Cancer Risk ###########################################
# This R script performs comprehensive two-sample Mendelian randomization (MR) analysis to investigate 
# the causal relationship between AD and various cancer types using GWAS summary statistics.
#######################################################################################

# Load required packages
library(TwoSampleMR)
library(MendelianRandomization)
library(MRPRESSO)
library(RadialMR)
library(openxlsx)
library(dplyr)

# Extract exposure data (AD genetic instruments)
exposure_dat <- extract_instruments("ebi-a-GCST90018784",
                                    p1 = 5e-08,
                                    r2 = 0.01,
                                    kb = 10000)
cat(sprintf("Number of exposure instruments: %d\n", nrow(exposure_dat)))

# Define outcome variables (cancer types)
outcomes_info <- data.frame(
  cancer_type = c("Esophageal cancer", "Lung cancer", "Skin cancer", "Breast cancer", "Uterine cancer", 
                  "Prostate cancer", "Kidney cancer", "Brain tumor", "T/NK cell lymphoma", 
                  "Non-Hodgkin lymphoma", "Primary lymphoid hematopoietic neoplasm"),
  outcome_id = c("ebi-a-GCST90018841", "ebi-a-GCST90018655", "finn-b-C3_OTHER_SKIN",
                 "finn-b-C3_BREAST", "finn-b-C3_CORPUS_UTERI", "ukb-b-2160",
                 "ukb-b-1316", "finn-b-C3_BRAIN", "finn-b-CD2_TNK_LYMPHOMA",
                 "finn-b-CD2_NONHODGKIN_NAS", "finn-b-CD2_PRIMARY_LYMPHOID_HEMATOPOIETIC_NAS"),
  code = c("C15", "C34", "C44", "C50", "C54", "C61", "C64", "C71", "C84", "C85", "C96"),
  stringsAsFactors = FALSE
)

# Create lists to store results
all_results <- list()
detailed_results <- list()
removed_snps_summary <- list()
summary_results <- data.frame()
significant_results <- data.frame()

# Custom function: Format p-values (avoid scientific notation, keep two decimal places)
format_pvalue <- function(p) {
  if (is.na(p) || is.null(p)) return("NA")
  if (is.character(p)) return(p)  # If already character, return directly
  if (!is.numeric(p)) return("NA")
  if (p < 0.01) {
    return(sprintf("%.2f", p))
  } else {
    return(sprintf("%.2f", p))
  }
}

# Custom function: Format confidence intervals
format_ci <- function(or, se) {
  if (is.na(or) || is.na(se) || is.null(or) || is.null(se)) return("NA")
  if (!is.numeric(or) || !is.numeric(se)) return("NA")
  lower <- exp(log(or) - 1.96 * se)
  upper <- exp(log(or) + 1.96 * se)
  return(sprintf("%.3f-%.3f", lower, upper))
}

# Start MR analysis loop
cat("\n=== Starting MR Analysis Loop ===\n")

for (i in 1:nrow(outcomes_info)) {
  cancer_type <- outcomes_info$cancer_type[i]
  outcome_id <- outcomes_info$outcome_id[i]
  cancer_code <- outcomes_info$code[i]
  
  cat(sprintf("\n\n--- Analysis %d/%d: %s (%s) ---\n", 
              i, nrow(outcomes_info), cancer_type, cancer_code))
  
  tryCatch({
    # Extract outcome data
    outcome_dat <- extract_outcome_data(snps = exposure_dat$SNP, 
                                        outcomes = outcome_id)
    
    if (nrow(outcome_dat) == 0) {
      cat("No matching outcome data found, skipping this analysis\n")
      next
    }
    
    # Harmonize data
    harmonised_dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
    initial_snp_count <- nrow(harmonised_dat)
    
    cat(sprintf("Initial SNP count: %d\n", initial_snp_count))
    
    # Record removed SNPs
    removed_snps <- data.frame(
      SNP = character(),
      Reason = character(),
      Source = character(),
      stringsAsFactors = FALSE
    )
    
    current_dat <- harmonised_dat
    
    # MR-PRESSO heterogeneity test and outlier removal
    presso_result <- NULL
    presso_global_p_initial <- NA
    presso_global_p_final <- NA
    
    if (nrow(current_dat) >= 3) {
      tryCatch({
        presso_result <- mr_presso(BetaOutcome = "beta.outcome", 
                                   BetaExposure = "beta.exposure", 
                                   SdOutcome = "se.outcome", 
                                   SdExposure = "se.exposure", 
                                   OUTLIERtest = TRUE, 
                                   DISTORTIONtest = TRUE, 
                                   data = current_dat, 
                                   NbDistribution = 1000,
                                   SignifThreshold = 0.05)
        
        presso_global_p_initial <- presso_result$`MR-PRESSO results`$`Global Test`$Pvalue
        
        if (presso_global_p_initial < 0.05 && !is.null(presso_result$`Main MR results`)) {
          outliers <- which(presso_result$`MR-PRESSO results`$`Outlier Test`$Pvalue < 0.05)
          if (length(outliers) > 0) {
            removed_snps_presso <- current_dat[outliers, ]
            removed_snps <- rbind(removed_snps, 
                                  data.frame(SNP = removed_snps_presso$SNP,
                                             Reason = "Outlier (p<0.05)",
                                             Source = "PRESSO",
                                             stringsAsFactors = FALSE))
            current_dat <- current_dat[-outliers, ]
            cat(sprintf("MR-PRESSO outlier removal: %d SNPs\n", length(outliers)))
            
            # Re-run PRESSO
            if (nrow(current_dat) >= 3) {
              presso_result_final <- mr_presso(BetaOutcome = "beta.outcome", 
                                               BetaExposure = "beta.exposure", 
                                               SdOutcome = "se.outcome", 
                                               SdExposure = "se.exposure", 
                                               OUTLIERtest = TRUE, 
                                               DISTORTIONtest = TRUE, 
                                               data = current_dat, 
                                               NbDistribution = 1000,
                                               SignifThreshold = 0.05)
              presso_global_p_final <- presso_result_final$`MR-PRESSO results`$`Global Test`$Pvalue
            }
          }
        }
      }, error = function(e) {
        cat("MR-PRESSO analysis failed:", e$message, "\n")
      })
    }
    
    # RadialMR outlier detection
    if (nrow(current_dat) >= 3) {
      tryCatch({
        radial_dat <- format_radial(current_dat$beta.exposure, 
                                    current_dat$beta.outcome,
                                    current_dat$se.exposure, 
                                    current_dat$se.outcome,
                                    current_dat$SNP)
        
        radial_result <- ivw_radial(radial_dat, alpha = 0.05, weights = 3)
        
        if (nrow(radial_result$outliers) > 0) {
          outlier_snps <- radial_result$outliers$SNP
          outlier_indices <- which(current_dat$SNP %in% outlier_snps)
          
          removed_snps <- rbind(removed_snps,
                                data.frame(SNP = outlier_snps,
                                           Reason = "Radial outlier",
                                           Source = "RadialMR",
                                           stringsAsFactors = FALSE))
          
          current_dat <- current_dat[!current_dat$SNP %in% outlier_snps, ]
          cat(sprintf("RadialMR outlier removal: %d SNPs\n", length(outlier_snps)))
        }
      }, error = function(e) {
        cat("RadialMR analysis failed:", e$message, "\n")
      })
    }
    
    final_snp_count <- nrow(current_dat)
    removed_snp_count <- initial_snp_count - final_snp_count
    
    cat(sprintf("Final SNP count: %d\n", final_snp_count))
    cat(sprintf("Removed SNP count: %d\n", removed_snp_count))
    
    # Display removed SNP details
    if (nrow(removed_snps) > 0) {
      cat("\nRemoved SNP details:\n")
      for (j in 1:nrow(removed_snps)) {
        cat(sprintf("  %s - %s (%s)\n", 
                    removed_snps$SNP[j], 
                    removed_snps$Reason[j], 
                    removed_snps$Source[j]))
      }
    }
    
    # Perform MR analysis
    if (final_snp_count < 3) {
      cat("Insufficient SNP count, skipping MR analysis\n")
      next
    }
    
    mr_results <- mr(current_dat)
    
    # Heterogeneity test
    heterogeneity_results <- mr_heterogeneity(current_dat)
    
    # Egger intercept test (pleiotropy test)
    pleiotropy_test <- mr_pleiotropy_test(current_dat)
    
    # Organize results
    detailed_result <- data.frame(
      Cancer_Type = cancer_type,
      Cancer_Code = cancer_code,
      Initial_SNPs = initial_snp_count,
      Final_SNPs = final_snp_count,
      Removed_SNPs = removed_snp_count,
      stringsAsFactors = FALSE
    )
    
    # Add MR results
    methods <- c("Inverse variance weighted", "MR Egger", "Weighted median", 
                 "Weighted mode", "Simple mode")
    method_short <- c("IVW", "MR_Egger", "Weighted_median", "Weighted_mode", "Simple_mode")
    
    for (k in 1:length(methods)) {
      method_result <- mr_results[mr_results$method == methods[k], ]
      if (nrow(method_result) > 0) {
        or_val <- exp(method_result$b)
        ci_val <- format_ci(or_val, method_result$se)
        p_val <- format_pvalue(method_result$pval)
        
        detailed_result[paste0(method_short[k], "_OR")] <- round(or_val, 3)
        detailed_result[paste0(method_short[k], "_95CI")] <- ci_val
        detailed_result[paste0(method_short[k], "_P")] <- p_val
      } else {
        detailed_result[paste0(method_short[k], "_OR")] <- NA
        detailed_result[paste0(method_short[k], "_95CI")] <- "NA"
        detailed_result[paste0(method_short[k], "_P")] <- "NA"
      }
    }
    
    # Add heterogeneity results
    ivw_het <- heterogeneity_results[heterogeneity_results$method == "Inverse variance weighted", ]
    egger_het <- heterogeneity_results[heterogeneity_results$method == "MR Egger", ]
    
    detailed_result$IVW_Q <- if(nrow(ivw_het) > 0) round(ivw_het$Q, 3) else NA
    detailed_result$IVW_Het_P <- if(nrow(ivw_het) > 0) format_pvalue(ivw_het$Q_pval) else "NA"
    detailed_result$Egger_Q <- if(nrow(egger_het) > 0) round(egger_het$Q, 3) else NA
    detailed_result$Egger_Het_P <- if(nrow(egger_het) > 0) format_pvalue(egger_het$Q_pval) else "NA"
    
    # Add quality control results
    detailed_result$PRESSO_Global_P_Initial <- if(!is.na(presso_global_p_initial)) format_pvalue(presso_global_p_initial) else "NA"
    detailed_result$PRESSO_Global_P_Final <- if(!is.na(presso_global_p_final)) format_pvalue(presso_global_p_final) else "NA"
    detailed_result$Egger_Intercept_P <- if(!is.na(pleiotropy_test$pval)) format_pvalue(pleiotropy_test$pval) else "NA"
    
    # Save detailed results
    detailed_results[[cancer_type]] <- detailed_result
    removed_snps_summary[[cancer_type]] <- removed_snps
    
    # Output current analysis results
    cat("\n=== MR Analysis Results ===\n")
    for (k in 1:length(methods)) {
      method_result <- mr_results[mr_results$method == methods[k], ]
      if (nrow(method_result) > 0) {
        or_val <- exp(method_result$b)
        ci_val <- format_ci(or_val, method_result$se)
        p_val <- format_pvalue(method_result$pval)
        cat(sprintf("%s: OR=%.3f, 95%%CI=%s, P=%s\n", 
                    methods[k], or_val, ci_val, p_val))
      }
    }
    
    cat("\n=== Heterogeneity Test ===\n")
    if(nrow(ivw_het) > 0) {
      cat(sprintf("IVW: Q=%.3f, P=%s\n", ivw_het$Q, format_pvalue(ivw_het$Q_pval)))
    }
    if(nrow(egger_het) > 0) {
      cat(sprintf("MR Egger: Q=%.3f, P=%s\n", egger_het$Q, format_pvalue(egger_het$Q_pval)))
    }
    
    cat("\n=== Quality Control ===\n")
    if (!is.na(presso_global_p_initial)) {
      cat(sprintf("MR-PRESSO global test (initial): P=%s\n", format_pvalue(presso_global_p_initial)))
    }
    if (!is.na(presso_global_p_final)) {
      cat(sprintf("MR-PRESSO global test (final): P=%s\n", format_pvalue(presso_global_p_final)))
    }
    cat(sprintf("Egger intercept test: P=%s\n", format_pvalue(pleiotropy_test$pval)))
    
    # Check for significant results
    ivw_result <- mr_results[mr_results$method == "Inverse variance weighted", ]
    if (nrow(ivw_result) > 0 && !is.na(ivw_result$pval) && ivw_result$pval < 0.05) {
      sig_result <- data.frame(
        Cancer_Type = cancer_type,
        Cancer_Code = cancer_code,
        Method = "IVW",
        OR = round(exp(ivw_result$b), 3),
        CI_95 = format_ci(exp(ivw_result$b), ivw_result$se),
        P_Value = format_pvalue(ivw_result$pval),
        SNP_Count = final_snp_count,
        stringsAsFactors = FALSE
      )
      significant_results <- rbind(significant_results, sig_result)
    }
    
  }, error = function(e) {
    cat(sprintf("Error occurred during analysis of %s: %s\n", cancer_type, e$message))
  })
}

# Combine all detailed results
if (length(detailed_results) > 0) {
  complete_results <- do.call(rbind, detailed_results)
  rownames(complete_results) <- NULL
  
  # Create simplified report table (suitable for publication)
  publication_table <- complete_results[, c("Cancer_Type", "Cancer_Code", "Final_SNPs",
                                            "IVW_OR", "IVW_95CI", "IVW_P",
                                            "MR_Egger_OR", "MR_Egger_95CI", "MR_Egger_P",
                                            "IVW_Het_P", "Egger_Intercept_P")]
  
  # Output summary results
  cat("\n\n=== Complete Results Table ===\n")
  print(complete_results)
  
  cat("\n\n=== Publication Table (Publication Ready) ===\n")
  print(publication_table)
  
  if (nrow(significant_results) > 0) {
    cat("\n\n=== Significant Results Summary ===\n")
    print(significant_results)
    cat(sprintf("Total significant associations found: %d\n", nrow(significant_results)))
  } else {
    cat("\n\n=== Significant Results Summary ===\n")
    cat("No significant causal associations found\n")
  }
  
  # Export Excel file
  wb <- createWorkbook()
  
  # Complete results table
  addWorksheet(wb, "Complete_Results")
  writeData(wb, "Complete_Results", complete_results)
  
  # Publication table
  addWorksheet(wb, "Publication_Table")
  writeData(wb, "Publication_Table", publication_table)
  
  # Significant results
  if (nrow(significant_results) > 0) {
    addWorksheet(wb, "Significant_Results")
    writeData(wb, "Significant_Results", significant_results)
  }
  
  # Removed SNPs summary
  if (length(removed_snps_summary) > 0) {
    addWorksheet(wb, "Removed_SNPs")
    all_removed_snps <- data.frame()
    for (cancer in names(removed_snps_summary)) {
      if (nrow(removed_snps_summary[[cancer]]) > 0) {
        temp_df <- removed_snps_summary[[cancer]]
        temp_df$Cancer_Type <- cancer
        all_removed_snps <- rbind(all_removed_snps, temp_df)
      }
    }
    if (nrow(all_removed_snps) > 0) {
      writeData(wb, "Removed_SNPs", all_removed_snps)
    }
  }
  
  # Save Excel file
  saveWorkbook(wb, "MR_Analysis_Results.xlsx", overwrite = TRUE)
  cat("\nResults exported to 'MR_Analysis_Results.xlsx'\n")
  
} else {
  cat("\nNo valid analysis results obtained\n")
}

cat("\n=== MR Analysis Loop Complete ===\n")
write.xlsx(all_removed_snps, "all_removed_snps.xlsx")