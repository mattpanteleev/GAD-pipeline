library(ggplot2)
library(ggrepel)
library(dplyr)
library(edgeR)

##### First experiment
# raffinose - black
# galactose - red
count_files <- c(
  "31_red_gene_counts_f.txt",
  "31_black_gene_counts_f.txt"
)
count_files <- c(
  "35_black_gene_counts_f.txt",
  "35_red_gene_counts_f.txt"
)

count_files <- c(
  "2990_black_gene_counts_f.txt",
  "2990_red_gene_counts_f.txt"
)
##### Second experiment
# glucose - black
# raffinose - blue
# galactose - red

# glucse 
count_files <- c(
  "new_90_black_gene_counts_f.txt",
  "new_90_red_gene_counts_f.txt"
)

count_files <- c(
  "new_111_black_gene_counts_f.txt",
  "new_111_red_gene_counts_f.txt"
)
count_files <- c(
  "new_118_black_gene_counts_f.txt",
  "new_118_red_gene_counts_f.txt"
)
count_files <- c(
  "new_70_black_gene_counts_f.txt",
  "new_70_red_gene_counts_f.txt"
)

# raffinocse
count_files <- c(
  "new_90_blue_gene_counts_f.txt",
  "new_90_red_gene_counts_f.txt"
)

count_files <- c(
  "new_111_blue_gene_counts_f.txt",
  "new_111_red_gene_counts_f.txt"
)
count_files <- c(
  "new_118_blue_gene_counts_f.txt",
  "new_118_red_gene_counts_f.txt"
)
count_files <- c(
  "new_70_blue_gene_counts_f.txt",
  "new_70_red_gene_counts_f.txt"
)

# raffinocse vs glucose

# glucose - black
# raffinose - blue

count_files <- c(
  "new_90_blue_gene_counts_f.txt",
  "new_90_black_gene_counts_f.txt"
)

count_files <- c(
  "new_111_blue_gene_counts_f.txt",
  "new_111_black_gene_counts_f.txt"
)
count_files <- c(
  "new_118_blue_gene_counts_f.txt",
  "new_118_black_gene_counts_f.txt"
)
count_files <- c(
  "new_70_blue_gene_counts_f.txt",
  "new_70_black_gene_counts_f.txt"
)


##### Download the data:
count_data_list <- list()
for (file in count_files) {
  data <- read.table(file, quote = "", header = TRUE, sep = "\t", row.names = 1)
  sample_name <- gsub("_gene_counts_f.txt", "", file)
  colnames(data) <- sample_name
  count_data_list[[sample_name]] <- data
}
count_matrix <- do.call(cbind, count_data_list)

# checking the errors in count_matrix:
sum(is.na(count_matrix))


##### Differential Expression:
group <- factor(c("A", "B","A", "B","C","A", "B","A", "B","C","A", "B","A", "B","A", "B","C"))  # 1 sample per group
group <- factor(c("A", "B"))  # 1 sample per group

# Preparation:
y <- DGEList(counts=count_matrix, group=group)

#y_check <- y[["counts"]]
#y_check <- data.frame(y_check)
#y_check$gene <- rownames(y_check)

keep <- filterByExpr(y)
y <- y[keep,]
y <- calcNormFactors(y, method = "TMM")
normalized_cpm <- cpm(y, normalized.lib.sizes = TRUE, log = FALSE)
#Skip estimateDisp — manually assign a reasonable value (e.g. 0.1 or 0.4) - because not enough numbers of repeats
y$common.dispersion <- 0.1

# Main part:
exact_results <- exactTest(y, dispersion = 0.1)  # specify it here too
results <- topTags(exact_results, n = Inf, adjust.method = "BH")$table

results$gene <- rownames(results)
normalized_cpm_df <- data.frame(normalized_cpm)
normalized_cpm_df$gene <- rownames(normalized_cpm_df)

merged_df <- merge(results,normalized_cpm_df, by = "gene")
rownames(merged_df) <- merged_df$gene
merged_df <- merged_df[, -1]
#write.csv(merged_df, file = "3070_2_raf_vs_gluc_exp.csv")


##### alleles to label:

alleles_to_label <- c(
  "SUC2",
  "REV1",
  "REV3rc_delCTD",
  "REV1wt",
  "REV7rc",
  "REV7wt",
  "REV3",
  "REV3rc",
  "REV7",
  "POL31",
  "POL31wt",
  "POL32wt",
  "POL32",
  "GAL3",
  "GAL1",
  "GAL7",
  "GAL10",
  "HXT1",
  "HXT3"
)


# Create a column indicating which alleles to label
results$highlight <- ifelse(
  results$gene %in% alleles_to_label,
  "yes", 
  "no"
)

fdr_threshold <- 0.05
logfc_threshold <- 1  
results <- results %>%
  mutate(
    color_category = case_when(
      # Alleles in your list AND significant
      gene %in% alleles_to_label & FDR < fdr_threshold & abs(logFC) > logfc_threshold ~ "Significant Gene of Interest",
      # Alleles in your list but NOT significant
      gene %in% alleles_to_label ~ "Non-Significant Allele",
      # Genes not in list but significant
      FDR < fdr_threshold & abs(logFC) > logfc_threshold ~ "Significant Gene",
      # All others
      TRUE ~ "Non-Significant"
    )
  )
results <- results %>%
  mutate(
    color_category = case_when(
      gene %in% alleles_to_label & 
        FDR < fdr_threshold & 
        abs(logFC) > logfc_threshold ~ "Significant Gene of Interest",
      gene %in% alleles_to_label ~ "Non-Significant Gene of Interest",
      FDR < fdr_threshold & 
        abs(logFC) > logfc_threshold ~ "Significant Gene",
      TRUE ~ "Non-Significant Gene"
    )
  )

results
# Define colors for each category
color_palette <- c(
  "Significant Gene of Interest" = "red",
  "Non-Significant Gene of Interest" = "orange",
  "Significant Gene" = "blue",
  "Non-Significant Gene" = "gray"
)

unique(results$color_category)

### Visualization of results. Volcano Plot
p <- ggplot(results, aes(x = logFC, y = -log10(PValue))) +
  geom_point(
    data = subset(results, color_category == "Non-Significant Gene"),
    color = "gray", alpha = 0.5, size = 1.5
  ) +
  geom_point(
    data = subset(results, color_category == "Significant Gene"),
    color = "blue", alpha = 0.4, size = 2
  ) +
  geom_point(
    data = subset(results, color_category == "Non-Significant Gene of Interest"),
    color = "orange", alpha = 0.8, size = 3
  ) +
  geom_point(
    data = subset(results, color_category == "Significant Gene of Interest"),
    color = "red", alpha = 1, size = 3
  ) +
  geom_text_repel(
    data = subset(results, gene %in% alleles_to_label),
    aes(label = gene),
    box.padding = 0.5,
    max.overlaps = Inf,
    size = 4
  ) +
  labs(x = "Log2 Fold Change", y = "-Log10(P-value)") +
  theme_minimal()
print(p)




# Check for NA values
sum(is.na(results$logFC))
sum(is.na(results$PValue))



tryCatch(dev.off(), error = function(e) message("No active graphics device"))
while (!is.null(dev.list())) 
  dev.off()
graphics.off()








# Enrichment analysis 

#BiocManager::install(c("clusterProfiler", "org.Sc.sgd.db", "DOSE"))


library(clusterProfiler)
library(org.Sc.sgd.db)
library(DOSE)
library(ggplot2)



significant_genes <- subset(results, FDR < 0.05 & abs(logFC) > 1)
gene_list <- rownames(significant_genes)  # Gene IDs (e.g., "YAL001C")

gene_list <- gsub("REV3rc", "REV3", gene_list)
gene_list <- gsub("REV7rc", "REV7", gene_list)
rename_row_names <- function(df) {
  rownames(df) <- gsub("REV3rc", "REV3", rownames(df))  #
  rownames(df) <- gsub("REV7rc", "REV7", rownames(df))    # Rename REV7 to REV7wt
  return(df)
}
results_new <- rename_row_names(results)


#upregulated
significant_genes <- subset(results_new, FDR < 0.05 & logFC > 1)
gene_list <- rownames(significant_genes)  # Gene IDs (e.g., "YAL001C")

go_enrich <- enrichGO(
  gene          = gene_list,
  universe      = rownames(results_new),  # Background: all tested genes
  OrgDb         = org.Sc.sgd.db,
  keyType       = "GENENAME",              # Use systematic ORF names (e.g., YAL001C)
  ont           = "BP",               # Biological Process (or "MF"/"CC")
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = FALSE  # Convert ORFs to gene names
)

# Visualize top GO terms
dotplot(go_enrich, showCategory = 14) + 
  ggtitle("Top Enriched GO Upregulated Biological Pathways")

#downrgulated
significant_genes <- subset(results, FDR < 0.05 & logFC < -1)
gene_list <- rownames(significant_genes)  # Gene IDs (e.g., "YAL001C")


go_enrich <- enrichGO(
  gene          = gene_list,
  universe      = rownames(results_new),  # Background: all tested genes
  OrgDb         = org.Sc.sgd.db,
  keyType       = "GENENAME",              # Use systematic ORF names (e.g., YAL001C)
  ont           = "BP",               # Biological Process (or "MF"/"CC")
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = FALSE  # Convert ORFs to gene names
)

# Visualize top GO terms
dotplot(go_enrich, showCategory = 7) + 
  ggtitle("Top Enriched GO Downregulated Biological Pathways")

#downrgulated
significant_genes <- subset(results, FDR < 0.05 & logFC < -1)
gene_list <- rownames(significant_genes)  # Gene IDs (e.g., "YAL001C")

head(go_enrich$Description) 
top_go <- head(go_enrich$Description, 10)
top_go
go_enrich$Description

write.csv(go_enrich@result, "GO_enrichment_results.csv")