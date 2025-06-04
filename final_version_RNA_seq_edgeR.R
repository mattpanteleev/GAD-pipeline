library(ggplot2)
library(ggrepel)
library(dplyr)
library(edgeR)

##### Download the data:
count_files <- c(
  "control_strain_gene_counts.txt",
  "case_strain_gene_counts.txt"
) # data after featureCounts
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

########### Differential Expression ############
group <- factor(c("A", "B"))  # 1 sample per group
# Preparation:
y <- DGEList(counts=count_matrix, group=group)
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
write.csv(merged_df, file = "example.csv") # save a table of the gene expression analysis 


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
      # Alleles in your list and significant:
      gene %in% alleles_to_label & FDR < fdr_threshold & abs(logFC) > logfc_threshold ~ "Significant Gene of Interest",
      # Alleles in your list but not significant:
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


########### Enrichment analysis ############

#BiocManager::install(c("clusterProfiler", "org.Sc.sgd.db", "DOSE"))


library(clusterProfiler)
library(org.Sc.sgd.db)
library(DOSE)
library(ggplot2)


#upregulated
significant_genes <- subset(results_new, FDR < 0.05 & logFC > 1)
gene_list <- rownames(significant_genes)  # Gene IDs of upregulated significant genes

go_enrich <- enrichGO(
  gene          = gene_list,
  universe      = rownames(results_new), 
  OrgDb         = org.Sc.sgd.db,
  keyType       = "GENENAME",          
  ont           = "BP",          # Biological Pathways 
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = FALSE 
)

# Visualize top GO terms
dotplot(go_enrich, showCategory = 10) + 
  ggtitle("Top Enriched GO Upregulated Biological Pathways")

#downrgulated
significant_genes <- subset(results, FDR < 0.05 & logFC < -1)
gene_list <- rownames(significant_genes)  # Gene IDs of downregulated significant genes


go_enrich <- enrichGO(
  gene          = gene_list,
  universe      = rownames(results_new),  
  OrgDb         = org.Sc.sgd.db,
  keyType       = "GENENAME",              
  ont           = "BP",               
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvawritelueCutoff  = 0.2,
  readable      = FALSE  
)

# Visualize top GO terms
dotplot(go_enrich, showCategory = 10) + 
  ggtitle("Top Enriched GO Downregulated Biological Pathways")
