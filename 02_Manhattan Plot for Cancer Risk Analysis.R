################### Manhattan Plot for Cancer Risk Analysis ###########################################

#######################################################################################

library(ggplot2)
library(ggrepel)
library(dplyr)
library(readr)

# Load Cox regression results data
data1 <- read.csv("cox_regression_results.CSV")
data2 <- read_csv("cox_male_results.csv")
data3 <- read.csv("cox_female_results.CSV")

# Define group labels (15 groups based on anatomical classification)
chr_lab <- c('Lip, oral cavity and pharynx',
             'Digestive organs',
             'Respiratory organs',
             'Bone and articular cartilage',
             'Skin',
             'Mesothelial and soft tissue',
             'Breast',
             'Female genital organs',
             'Male genital organs',
             'Urinary tract',
             'Central nervous system',
             'Endocrine glands',
             'ISUS',
             'Haematopoietic',
             'Primary multiple sites')

# Prepare final_results data
# Ensure data frame has the following columns: CancerType, p_value

# Convert p-values to numeric format (handle "<0.0001" format)
data3$p_value <- as.numeric(gsub("<", "", data3$p_value))
data3 <- data3 %>% filter(!is.na(p_value))

# Manually assign group IDs (15 groups based on ICD-10 classification)
assign_group_id <- function(cancer_type) {
  group_id <- case_when(
    cancer_type %in% c("C00", "C01", "C02", "C03", "C04", "C05", "C06", "C07", "C08", "C09", "C10", "C11", "C12", "C13", "C14") ~ 1,
    cancer_type %in% c("C15", "C16", "C17", "C18", "C19", "C20", "C21", "C22", "C23", "C24", "C25", "C26") ~ 2,
    cancer_type %in% c("C30", "C31", "C32", "C33", "C34", "C35", "C36", "C37", "C38", "C39") ~ 3,
    cancer_type %in% c("C40", "C41") ~ 4,
    cancer_type %in% c("C43", "C44") ~ 5,
    cancer_type %in% c("C45", "C46", "C47", "C48", "C49") ~ 6,
    cancer_type == "C50" ~ 7,
    cancer_type %in% c("C51", "C52", "C53", "C54", "C55", "C56", "C57", "C58") ~ 8,
    cancer_type %in% c("C60", "C61", "C62", "C63") ~ 9,
    cancer_type %in% c("C64", "C65", "C66", "C67", "C68") ~ 10,
    cancer_type %in% c("C69", "C70", "C71", "C72") ~ 11,
    cancer_type %in% c("C73", "C74", "C75") ~ 12,
    cancer_type %in% c("C76", "C77", "C78", "C79", "C80") ~ 13,
    cancer_type %in% c("C81", "C82", "C83", "C84", "C85", "C86", "C87", "C88", "C89", "C90", "C91", "C92", "C93", "C94", "C95", "C96") ~ 14,
    cancer_type == "C97" ~ 15,
    TRUE ~ NA_integer_  # Cancer types not in the above groups
  )
  return(group_id)
}

# Function to prepare Manhattan plot data
prepare_manhattan_data <- function(data) {
  # Assign group IDs
  data$CHR <- sapply(data$CancerType, assign_group_id)
  
  # Remove unrecognized groups
  data <- data %>% filter(!is.na(CHR))
  
  # Add log10 p-values
  data$log10p <- -log10(data$p_value)
  
  # Set upper limit to prevent extreme values
  data$log10p <- pmin(data$log10p, 10)
  
  # Add point positions within each group
  data <- data %>% 
    group_by(CHR) %>% 
    mutate(
      position_in_group = row_number(),
      bp = position_in_group
    ) %>% 
    ungroup()
  
  # Add cumulative positions for x-axis
  data_cum <- data %>% 
    group_by(CHR) %>% 
    summarise(max_bp = max(position_in_group)) %>% 
    mutate(
      bp_add = lag(cumsum(max_bp), default = 0)
    )
  
  data <- data %>% 
    inner_join(data_cum, by = "CHR") %>% 
    mutate(bp_cum = bp + bp_add)
  
  return(data)
}

# Function to create Manhattan plot
manhattan_plot <- function(data) {
  # Prepare data
  tb1 <- prepare_manhattan_data(data)
  
  # Get existing groups in the data
  existing_groups <- unique(tb1$CHR)
  
  # Check if all groups are present
  if (length(existing_groups) < 15) {
    message("Note: Data contains only ", length(existing_groups), "/15 cancer groups")
  }
  
  # Calculate x-axis tick positions
  axis_set <- tb1 %>%
    group_by(CHR) %>%
    summarize(center = mean(bp_cum))
  
  # Calculate y-axis range
  min_pvalue <- min(tb1$p_value, na.rm = TRUE)
  ylim <- max(ceiling(-log10(min_pvalue)) + 1, 5)  # Ensure y-axis shows at least up to 5
  
  # Get labels for existing groups
  group_labels <- chr_lab[existing_groups]
  
  # Ensure groups are ordered properly
  axis_set <- axis_set %>% arrange(CHR)
  
  # Create Manhattan plot
  p1 <- ggplot(tb1, aes(x = bp_cum, y = log10p, color = as.factor(CHR))) +
    geom_hline(yintercept = seq(0, ceiling(ylim), 1), color = "#DCDCDC", 
               linetype = "solid", alpha = 0.5) + 
    geom_hline(yintercept = -log10(0.05), color = "red", 
               linetype = "dashed", alpha = 0.8) +
    geom_point(alpha = 1, shape = 16, size = 4) +
    geom_text_repel(data = tb1 %>% filter(p_value < 0.05), 
                    aes(label = CancerType), 
                    color = "black",
                    size = 3.5,
                    max.overlaps = 20,
                    min.segment.length = 0.3,
                    box.padding = 0.5) +
    scale_x_continuous(
      breaks = axis_set$center,
      labels = group_labels,
      expand = c(0.15, 0.15)  # Control x-axis expansion to avoid crowded labels
    ) +
    scale_y_continuous(
      breaks = seq(0, ceiling(ylim), 1),
      limits = c(0, ceiling(ylim))
    ) + 
    scale_color_manual(values = rep(c("#91CCC0", "#7FABD1", "#EC6E66"), length = 15)) +
    labs(x = "Cancer Groups", y = expression(-log[10](p)), 
         title = "Manhattan Plot of Cancer Risk Associations",
         subtitle = paste("Showing", length(existing_groups), "of 15 cancer groups")) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(family = "mono", face = "bold", angle = 45, hjust = 1, size = 10),  # Fine-tune x-axis labels
      axis.title.x = element_text(size = 13, margin = margin(t = 15)),
      axis.title.y = element_text(size = 13, margin = margin(r = 15)),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray30"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
}

# Generate Manhattan plot
p <- manhattan_plot(data)

# Display the plot
print(p)