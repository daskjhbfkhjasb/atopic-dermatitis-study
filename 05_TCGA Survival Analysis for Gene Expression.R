################### TCGA Survival Analysis for Gene Expression ###########################################

#######################################################################################

library(TCGAbiolinks)
library(SummarizedExperiment)
library(survival)
library(survminer)
library(dplyr)
library(biomaRt)

# Set TCGA project to lung adenocarcinoma (TCGA-LUAD)
query_expr <- GDCquery(project = "TCGA-LUAD",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       workflow.type = "STAR - Counts")

# Download and prepare expression data
GDCdownload(query_expr)
expr_data <- GDCprepare(query_expr)

# Extract expression matrix
expr_matrix <- assay(expr_data)

# Remove version numbers from gene IDs (e.g., ENSG00000163558.7 -> ENSG00000163558)
rownames(expr_matrix) <- gsub("\\..*", "", rownames(expr_matrix))

# Set up biomaRt connection
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Get Ensembl ID for target gene (HLA-DPB1 in this example)
gene_map <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                  filters = "hgnc_symbol",
                  values = "HLA-DPB1",
                  mart = mart)

print(paste("Gene found:", gene_map$hgnc_symbol))
print(paste("Corresponding Ensembl ID:", gene_map$ensembl_gene_id))

# Select target gene ID (adjust index as needed)
target_gene_id <- gene_map$ensembl_gene_id[8]

# Check if gene exists in expression matrix
if (!target_gene_id %in% rownames(expr_matrix)) {
  stop("Target gene not found in expression matrix")
}

# Extract target gene expression data
target_gene_expr <- expr_matrix[target_gene_id, ]

# Get clinical data
clinical <- colData(expr_data)
clinical_df <- as.data.frame(clinical)

# Extract survival time and status information
surv_data <- clinical_df %>%
  dplyr::select(submitter_id, vital_status, days_to_death, days_to_last_follow_up) %>%
  mutate(
    time = ifelse(is.na(days_to_death), days_to_last_follow_up, days_to_death),
    status = ifelse(vital_status == "Alive", 0, 1)
  ) %>%
  filter(!is.na(time), time > 0) %>%  # Filter out invalid time data
  dplyr::select(submitter_id, time, status)

# Standardize expression matrix sample IDs (remove extra information)
expr_id <- substr(colnames(expr_matrix), 1, 12)
target_gene_df <- data.frame(submitter_id = expr_id, expression = as.numeric(target_gene_expr))

# Merge expression and survival data
merged_data <- merge(target_gene_df, surv_data, by = "submitter_id")

# Remove samples with zero or NA expression values
merged_data <- merged_data %>%
  filter(!is.na(expression), expression > 0)

print(paste("Final sample size for analysis:", nrow(merged_data)))

# Group by median expression level
median_expr <- median(merged_data$expression)
merged_data$group <- ifelse(merged_data$expression >= median_expr, "High", "Low")

# Survival analysis using Kaplan-Meier method
fit <- survfit(Surv(time, status) ~ group, data = merged_data)

# Create publication-ready survival plot
p1 <- ggsurvplot(
  fit, 
  data = merged_data,
  
  # Basic settings
  pval = TRUE,
  pval.size = 5,
  pval.coord = c(100, 0.25),
  conf.int = TRUE,  # Display confidence intervals
  conf.int.alpha = 0.2,  # Confidence interval transparency
  
  # Risk table
  risk.table = TRUE,
  risk.table.col = "strata",
  risk.table.height = 0.3,
  
  # Titles and labels
  title = "HLA-DPB1 Expression and Overall Survival in LUAD Patients",
  subtitle = "Kaplan-Meier Survival Curves with 95% Confidence Intervals",
  caption = "Data source: TCGA-LUAD",
  xlab = "Overall Survival Time (Days)",
  ylab = "Survival Probability",
  
  # Legend
  legend.title = "HLA-DPB1 Expression",
  legend.labs = c("High Expression", "Low Expression"),
  legend = c(0.8, 0.8),
  
  # Colors and theme
  palette = c("#E31A1C", "#1F78B4"),  # Red and blue
  ggtheme = theme_minimal(),
  
  # Axis settings
  break.time.by = 365,  # Display by year
  xlim = c(0, max(merged_data$time, na.rm = TRUE)),
  
  # Other aesthetic options
  fontsize = 4,
  tables.theme = theme_cleantable()
)

# Further plot customization - fix deprecated warnings
p1$plot <- p1$plot + 
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
    plot.caption = element_text(size = 10, color = "gray50"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)  # Use linewidth instead of size
  )

# Display plot
print(p1)

# Save high-quality image
ggsave("HLA-DPB1_survival_analysis_lung.png", plot = p1$plot, 
       width = 12, height = 8, dpi = 300, units = "in")

# Cox proportional hazards regression analysis
cox_fit <- coxph(Surv(time, status) ~ expression, data = merged_data)
cox_summary <- summary(cox_fit)
print("Continuous Cox Regression Results:")
print(cox_summary)

# Grouped Cox regression analysis
cox_group <- coxph(Surv(time, status) ~ group, data = merged_data)
cox_group_summary <- summary(cox_group)
print("Categorical Cox Regression Results:")
print(cox_group_summary)

# Additional statistical information
cat("\n=== Descriptive Statistics ===\n")
cat("Total sample size:", nrow(merged_data), "\n")
cat("High expression group size:", sum(merged_data$group == "High"), "\n")
cat("Low expression group size:", sum(merged_data$group == "Low"), "\n")
cat("Number of death events:", sum(merged_data$status), "\n")
cat("Median expression level:", round(median_expr, 4), "\n")

# Median survival time
surv_summary <- surv_median(fit)
print("Median Survival Time:")
print(surv_summary)

# Log-rank test
logrank_test <- survdiff(Surv(time, status) ~ group, data = merged_data)
cat("\nLog-rank Test Results:\n")
print(logrank_test)
cat("P-value:", 1 - pchisq(logrank_test$chisq, length(logrank_test$n) - 1), "\n")