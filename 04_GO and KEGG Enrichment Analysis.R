################### GO and KEGG Enrichment Analysis ###########################################

library(readxl)
diff <- read_excel("DATA.xlsx")

# 1. Load required packages
BiocManager::install("AnnotationDbi")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("clusterProfiler")
library(AnnotationDbi)
library(org.Hs.eg.db)  # Human genome annotation package
library(clusterProfiler)  # Enrichment analysis package
library(dplyr)
library(ggplot2)  # Plotting package

# 2. Gene ID conversion (KEGG and GO enrichment require ENTREZID format)
gene.df <- bitr(diff$SYMBOL, fromType="SYMBOL", toType="ENTREZID", OrgDb = org.Hs.eg.db)
gene <- gene.df$ENTREZID

# 3. GO enrichment analysis
## CC = Cellular Component, MF = Molecular Function, BP = Biological Process
## ALL = enrichment of all three categories simultaneously
## Choose based on your needs. Generally, BP, MF, CC are analyzed separately then combined into one dataframe for easier pathway selection and plotting

# Enrichment analysis for all GO categories
ego_ALL <- enrichGO(gene = gene,  # Gene list defined above
                    OrgDb = org.Hs.eg.db,
                    keyType = "ENTREZID",
                    ont = "ALL",  # GO enrichment type
                    pAdjustMethod = "BH",  # Multiple testing correction method (standard: BH)
                    minGSSize = 1,
                    pvalueCutoff = 0.05,  # P-value threshold
                    qvalueCutoff = 0.05,
                    readable = TRUE)

# Cellular Component enrichment
ego_CC <- enrichGO(gene = gene,
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "CC",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)

# Biological Process enrichment
ego_BP <- enrichGO(gene = gene,
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)

# Molecular Function enrichment
ego_MF <- enrichGO(gene = gene,
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "MF",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)

# 4. Save results to current directory
ego_ALL <- as.data.frame(ego_ALL)
ego_result_BP <- as.data.frame(ego_BP)
ego_result_CC <- as.data.frame(ego_CC)
ego_result_MF <- as.data.frame(ego_MF)
ego <- rbind(ego_result_BP, ego_result_CC, ego_result_MF)  # Alternative way to get same result as ego_ALL

# Export enrichment results
write.csv(ego_ALL, file = "ego_ALL.csv", row.names = T)
write.csv(ego_result_BP, file = "ego_result_BP.csv", row.names = T)
write.csv(ego_result_CC, file = "ego_result_CC.csv", row.names = T)
write.csv(ego_result_MF, file = "ego_result_MF.csv", row.names = T)
write.csv(ego, file = "ego.csv", row.names = T)

# 5. Select subset of pathways for visualization (when too many pathways are enriched)
# Skip step 4 and proceed to step 5 for selective plotting
display_number = c(9, 13, 13)  # Number of pathways to display for BP, CC, MF respectively (adjustable)
ego_result_BP <- as.data.frame(ego_BP)[1:display_number[1], ]
ego_result_CC <- as.data.frame(ego_CC)[1:display_number[2], ]
ego_result_MF <- as.data.frame(ego_MF)[1:display_number[3], ]

## Combine selected pathways into a single dataframe
go_enrich_df <- data.frame(
  ID = c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),
  Description = c(ego_result_BP$Description, ego_result_CC$Description, ego_result_MF$Description),
  geneID = c(ego_result_BP$geneID, ego_result_CC$geneID, ego_result_MF$geneID),
  GeneNumber = c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
  type = factor(c(rep("biological process", display_number[1]), 
                  rep("cellular component", display_number[2]),
                  rep("molecular function", display_number[3])), 
                levels = c("biological process", "cellular component", "molecular function")))

## Truncate pathway names (too long for visualization)
## Select first 20 words of pathway description as pathway name
for(i in 1:nrow(go_enrich_df)){
  description_splite = strsplit(go_enrich_df$Description[i], split = " ")
  description_collapse = paste(description_splite[[1]][1:20], collapse = " ")  # 20 words (adjustable)
  go_enrich_df$Description[i] = description_collapse
  go_enrich_df$Description = gsub(pattern = "NA", "", go_enrich_df$Description)
}

write.csv(go_enrich_df, "go_enrich_df_11111.csv")
go_enrich_df <- read.csv("go_enrich_df.csv", row.names = 1)

## Create GO bar chart
### Horizontal bar chart
go_enrich_df$type_order = factor(rev(as.integer(rownames(go_enrich_df))), 
                                 labels = rev(go_enrich_df$Description))  # Essential step for proper bar ordering
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")  # Color scheme

plot1 <- ggplot(data = go_enrich_df, aes(x = type_order, y = GeneNumber, fill = type)) +
  geom_bar(stat = "identity", width = 0.8) +  # Bar width (adjustable)
  scale_fill_manual(values = COLS) +  # Color assignment
  coord_flip() +  # Flip coordinates for horizontal bars
  xlab("GO term") + 
  ylab("Gene_Number") + 
  labs(title = "The Most Enriched GO Terms") +
  geom_text(aes(label = GeneNumber), hjust = -0.2, size = 3) +
  theme_bw()  # Classic black and white theme

# Alternative horizontal bar chart without y-axis values
ggplot(data = go_enrich_df, aes(x = type_order, fill = type)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = COLS) +
  coord_flip() +
  xlab("GO term") + 
  ylab("Gene_Number") + 
  labs(title = "The Most Enriched GO Terms") +
  geom_text(aes(label = GeneNumber), hjust = -0.2, size = 3) +
  theme_bw()

# Save horizontal bar chart
ggsave("GO.pdf", plot = plot1, width = 10, height = 6, dpi = 600)

### Vertical bar chart 
go_enrich_df$type_order = factor(rev(as.integer(rownames(go_enrich_df))), 
                                 labels = rev(go_enrich_df$Description))
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")

ggplot(data = go_enrich_df, aes(x = type_order, y = GeneNumber, fill = type)) + 
  geom_bar(stat = "identity", width = 0.6) + 
  scale_fill_manual(values = COLS) + 
  theme_bw() + 
  xlab("GO term") + 
  ylab("Num of Genes") + 
  labs(title = "The Most Enriched GO Terms") + 
  theme(axis.text.x = element_text(face = "bold", color = "gray50", angle = 70, vjust = 1, hjust = 1))  # angle = rotation angle for x-axis labels (adjustable)

write.csv(go_enrich_df, "go_enrich_df.csv")

# KEGG Pathway Enrichment Analysis
# 1. KEGG enrichment
kk <- enrichKEGG(gene = gene, keyType = "kegg", organism = "human", 
                 qvalueCutoff = 0.05, pvalueCutoff = 0.05)

# 2. Visualization
### Bar chart
hh <- as.data.frame(kk)  # Remember to save results!
rownames(hh) <- 1:nrow(hh)
hh$order = factor(rev(as.integer(rownames(hh))), labels = rev(hh$Description))

ggplot(hh, aes(y = order, x = Count, fill = p.adjust)) +
  geom_bar(stat = "identity", width = 0.7) +  # Bar width
  scale_fill_gradient(low = "red", high = "blue") +  # Color gradient (customizable)
  labs(title = "KEGG Pathways Enrichment",
       x = "Gene numbers", 
       y = "Pathways") +
  theme(axis.title.x = element_text(face = "bold", size = 16),
        axis.title.y = element_text(face = "bold", size = 16),
        legend.title = element_text(face = "bold", size = 16)) +
  theme_bw()

hh <- as.data.frame(kk)
hh <- read.csv("kegg_enrich_df_22222.csv", row.names = 1)

### Bubble plot
rownames(hh) <- 1:nrow(hh)
hh$order = factor(rev(as.integer(rownames(hh))), labels = rev(hh$Description))

ggplot(hh, aes(y = order, x = Count)) +
  geom_point(aes(size = Count, color = -1 * p.adjust)) +  # Point size modification
  scale_color_gradient(low = "green", high = "red") +
  geom_text(aes(label = Count), hjust = -0.7, size = 4, color = "black") +
  labs(color = expression(p.adjust, size = "Count"), 
       x = "Gene Number", y = "Pathways", title = "KEGG Pathway Enrichment") +
  theme_bw()

# Alternative bubble plot
ggplot(hh, aes(y = order, x = Count)) +
  geom_point(aes(size = Count, color = -1 * p.adjust)) +
  scale_color_gradient(low = "green", high = "red") +
  geom_text(aes(label = Count), hjust = -0.5, size = 4) +
  labs(color = expression(p.adjust, size = "Count"), 
       x = "Gene Number", y = "Pathways", title = "KEGG Pathway Enrichment") +
  theme_bw()

# Save KEGG plot (note: plot2 variable needs to be defined)
ggsave("KEGG.pdf", plot = plot2, width = 6, height = 6, dpi = 600)

# Remove specific rows from BP results (optional filtering)
ego_result_BP <- ego_result_BP[-c(2, 4, 6, 13), ]