library(Seurat) 
library(dplyr)
library(purrr)
library(ggplot2)
library(cowplot)
library(EnhancedVolcano)
library(pheatmap)
library(Matrix)
library(tibble)

# samples with correct paths to my computer
samples <- list(
  Lesional_early = "C:/Users/Brakebusch/Downloads/GSM7488053_2330-AB-1.h5",
  Lesional_late = "C:/Users/Brakebusch/Downloads/GSM7488055_4902-AB-1.h5",
  NonLesional_early = "C:/Users/Brakebusch/Downloads/GSM7488054_2330-AB-3.h5",
  NonLesional_late = "C:/Users/Brakebusch/Downloads/GSM7488056_4902-AB-2.h5",
  Control1 = "C:/Users/Brakebusch/Downloads/GSM7488057_2057-AB-3.h5",
  Control2 = "C:/Users/Brakebusch/Downloads/GSM7488058_2057-AB-5.h5",
  Control3 = "C:/Users/Brakebusch/Downloads/GSM7488059_2057-AB-6.h5"
)

#read H5 file and create Seurat object
read_h5_to_seurat <- function(file_path, sample_name) {
  counts <- Read10X_h5(file_path)
  seurat_obj <- CreateSeuratObject(counts = counts, project = sample_name)
  seurat_obj$sample <- sample_name
  seurat_obj$condition <- ifelse(grepl('Control', sample_name), 'Control', 'Morphea')
  return(seurat_obj)
}

#create object
seurat_objects <- map2(samples, names(samples), ~read_h5_to_seurat(.x, .y))

# Preprocess each Seurat object individually,standard settings used
for (i in 1:length(seurat_objects)) {
  seurat_objects[[i]] <- NormalizeData(seurat_objects[[i]])
  seurat_objects[[i]] <- FindVariableFeatures(seurat_objects[[i]], selection.method = "vst", nfeatures = 2000)
}

# Select features that are repeatedly variable across datasets for integration,standard settings used
features <- SelectIntegrationFeatures(object.list = seurat_objects, nfeatures = 2000)

# Find integration anchors
anchors <- FindIntegrationAnchors(object.list = seurat_objects, anchor.features = features)

# Integrate the data
combined <- IntegrateData(anchorset = anchors)

DefaultAssay(combined) <- "integrated"

# Scaling, PCA, and UMAP
combined <- ScaleData(combined)
combined <- RunPCA(combined)
combined <- RunUMAP(combined, dims = 1:30)
combined <- FindNeighbors(combined, dims = 1:30)
combined <- FindClusters(combined, resolution = 0.5)

# Plot UMAP colored by condition
umap_plot_condition <- DimPlot(combined, reduction = "umap", group.by = "condition", label = TRUE) +
  ggtitle("UMAP plot of cells colored by condition (Integrated Data)")
print(umap_plot_condition)
ggsave("umap_plot_condition_integrated.png", umap_plot_condition, width = 10, height = 8)

# Plot UMAP colored by clusters
umap_plot_clusters <- DimPlot(combined, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("UMAP plot of cells colored by clusters (Integrated Data)")
print(umap_plot_clusters)
ggsave("umap_plot_clusters_integrated.png", umap_plot_clusters, width = 10, height = 8)

# Save the integrated Seurat object
saveRDS(combined, file = "combined_integrated_seurat_object.rds")

# summary and counts here
cell_counts <- table(combined$sample, combined$condition)
print(cell_counts)

total_cells <- sum(cell_counts)
cat("Total number of cells:", total_cells, "\n")


cells_per_condition <- colSums(cell_counts)
cat("Number of cells per condition:\n")
print(cells_per_condition)





#DGE

DefaultAssay(combined) <- "RNA"

combined <- JoinLayers(combined)


Idents(combined) <- "condition"
de_results <- FindMarkers(
  combined,
  ident.1 = "Morphea",
  ident.2 = "Control",
  assay = "RNA",
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# Sort the results by adjusted p-value
de_results <- de_results[order(de_results$p_val_adj), ]

# the top 20 
print(head(de_results, 20))

# Save the full DE results
write.csv(de_results, "de_results_morphea_vs_control.csv", row.names = TRUE)

# Create a volcano plot
volcano_plot <- EnhancedVolcano(
  de_results,
  lab = rownames(de_results),
  x = 'avg_log2FC',
  y = 'p_val_adj',
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 3.0,
  labSize = 6.0,
  title = 'Morphea vs Control',
  subtitle = 'Differential Expression'
)
print(volcano_plot)
ggsave("volcano_plot_morphea_vs_control.png", plot = volcano_plot, width = 12, height = 10)

# Get the top 10 up-regulated and down-regulated genes
top_up <- head(rownames(de_results[de_results$avg_log2FC > 0, ]), 10)
top_down <- head(rownames(de_results[de_results$avg_log2FC < 0, ]), 10)

# Create a heatmap of top DE genes
top_genes <- c(top_up, top_down)
heatmap_data <- GetAssayData(combined, assay = "RNA", layer =  "data")[top_genes, ]
heatmap_data <- as.matrix(heatmap_data)
heatmap_data <- heatmap_data - rowMeans(heatmap_data)

# Create a data frame for annotation
annotation_col <- data.frame(Condition = combined$condition)
rownames(annotation_col) <- colnames(heatmap_data)

heatmap_plot <- pheatmap::pheatmap(
  heatmap_data,
  annotation_col = annotation_col,
  show_colnames = FALSE,
  main = "Top 20 Differentially Expressed Genes",
  fontsize_row = 8,
  cluster_cols = FALSE
)

# Summarize the number of DE genes
de_summary <- data.frame(
  Total_DE_genes = nrow(de_results[de_results$p_val_adj < 0.05, ]),
  Upregulated = nrow(de_results[de_results$p_val_adj < 0.05 & de_results$avg_log2FC > 0, ]),
  Downregulated = nrow(de_results[de_results$p_val_adj < 0.05 & de_results$avg_log2FC < 0, ])
)
print(de_summary)


#MANUAL ANNOTATION!!!!

# Get the list of genes in our dataset
available_genes <- rownames(combined)

# Function to filter marker genes
filter_markers <- function(markers, available_genes) {
  lapply(markers, function(x) intersect(x, available_genes))
}

# Define marker genes for each cell type
cell_markers <- list(
  Fibroblasts = c("COL1A1", "COL1A2", "DCN", "LUM"),
  Keratinocytes = c("KRT5", "KRT14", "KRT1", "KRT10"),
  Th1_cells = c("CD4", "IFNG", "TBX21"),
  Th2_cells = c("CD4", "IL4", "GATA3"),
  Myeloid_cells = c("CD14", "CD68", "LYZ"),
  B_cells = c("CD19", "MS4A1", "CD79A"),
  Endothelial_cells = c("PECAM1", "VWF", "CDH5"),
  Langerhans_cells = c("CD1A", "CD207", "LANGERIN"),
  Red_blood_cells = c("HBA1", "HBB"),
  Melanocytes = c("TYRP1", "DCT", "MLANA"),
  Mast_cells = c("TPSAB1", "CPA3", "KIT"),
  Pericytes = c("PDGFRB", "RGS5", "ACTA2"),
  Smooth_muscle_cells = c("ACTA2", "MYH11", "CNN1"),
  Treg_cells = c("CD4", "FOXP3", "IL2RA"),
  Th17_cells = c("CD4", "RORC", "IL17A", "CCR6"),
  CD8_T_cells = c("CD3D", "CD3E", "CD3G", "CD8A", "CD8B", "GZMA", "GZMB", "PRF1")
)

# Filter marker genes
filtered_markers <- filter_markers(cell_markers, available_genes)

# Print the filtered markers to see which genes are available
print(filtered_markers)

# Function to score cells based on marker genes
score_cells <- function(seurat_obj, markers) {
  scores <- matrix(0, nrow = ncol(seurat_obj), ncol = length(markers))
  colnames(scores) <- names(markers)
  
  for (cell_type in names(markers)) {
    if (length(markers[[cell_type]]) > 0) {
      gene_set <- markers[[cell_type]]
      # Ensure genes are present in the dataset
      gene_set <- intersect(gene_set, rownames(seurat_obj))
      if (length(gene_set) > 0) {
        scores[, cell_type] <- Matrix::colMeans(GetAssayData(seurat_obj, assay = "RNA", slot = "data")[gene_set, ])
      }
    }
  }
  return(scores)
}

# Score cells and add to Seurat object
cell_scores <- score_cells(combined, filtered_markers)
cell_scores_df <- as.data.frame(cell_scores)
combined <- AddMetaData(combined, metadata = cell_scores_df)

# Assign cell types based on highest score
combined$cell_type <- colnames(cell_scores)[max.col(cell_scores, ties.method = "first")]

# Print summary of cell type assignments
print(table(combined$cell_type, combined$condition))

# Visualize UMAP plot with cell type annotations
umap_plot_cell_type <- DimPlot(combined, reduction = "umap", group.by = "cell_type", label = TRUE) +
  ggtitle("UMAP plot of cells colored by annotated cell type")
print(umap_plot_cell_type)
ggsave("umap_plot_annotated_cell_types.png", umap_plot_cell_type, width = 12, height = 10)

# Save the updated Seurat object
saveRDS(combined, file = "combined_annotated_seurat_object.rds")

#DGE ON THE EACH CELL TYPE

# Function 
perform_de_analysis <- function(seurat_obj, cell_type_value) {
  # Subset the Seurat object to include only cells of the specified cell type
  subset_obj <- subset(seurat_obj, subset = cell_type == cell_type_value)
  
  # Check if there are enough cells in each condition
  if (length(unique(subset_obj$condition)) < 2) {
    warning(paste("Not enough conditions for cell type:", cell_type_value))
    return(NULL)
  }
  
  # Set the identity classes to 'condition' for DE analysis
  Idents(subset_obj) <- "condition"
  
  # Perform differential expression analysis between Morphea and Control
  de_results <- FindMarkers(
    subset_obj,
    ident.1 = "Morphea",
    ident.2 = "Control",
    assay = "RNA",
    min.pct = 0.25,
    logfc.threshold = 0.25
  )
  
  return(de_results)
}

# Perform DE analysis for each cell type
cell_types <- unique(combined$cell_type)
de_results_list <- lapply(cell_types, function(ct) perform_de_analysis(combined, ct))
names(de_results_list) <- cell_types

# Remove NULL entries in case some cell types didn't have enough cells
de_results_list <- de_results_list[!sapply(de_results_list, is.null)]

# Function to create a volcano plot for a specific cell type
create_volcano_plot <- function(de_results, cell_type) {
  de_results$gene <- rownames(de_results)
  de_results$significant <- de_results$p_val_adj < 0.05 & abs(de_results$avg_log2FC) > 1
  plot <- ggplot(de_results, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    geom_point(aes(color = significant), alpha = 0.5) +
    scale_color_manual(values = c("grey", "red")) +
    labs(
      title = paste("Volcano Plot -", cell_type),
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-value"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
  return(plot)
}

# Create volcano plots for each cell type
volcano_plots <- lapply(names(de_results_list), function(ct) {
  create_volcano_plot(de_results_list[[ct]], ct)
})

# Combine all volcano plots into a single figure
combined_volcano_plot <- cowplot::plot_grid(plotlist = volcano_plots, ncol = 3)
ggsave("combined_volcano_plot.png", combined_volcano_plot, width = 20, height = 20)

# Function to get top DE genes for a specific cell type
get_top_de_genes <- function(de_results, n = 10) {
  de_results %>%
    tibble::rownames_to_column("gene") %>%
    arrange(p_val_adj) %>%
    slice_head(n = n) %>%
    select(gene, avg_log2FC, p_val_adj)
}

# Get top DE genes for each cell type
top_de_genes_list <- lapply(de_results_list, get_top_de_genes)

# Print top DE genes for each cell type
for (ct in names(top_de_genes_list)) {
  cat("\nTop DE genes for", ct, ":\n")
  print(top_de_genes_list[[ct]])
}

# Save the DE results
saveRDS(de_results_list, file = "de_results_by_cell_type.rds")



#focus on endothelial cells

# Subset the Seurat object to include only endothelial cells
endothelial_cells <- subset(combined, subset = cell_type == "Endothelial_cells")

# Check the number of endothelial cells per condition
table(endothelial_cells$condition)

# UMAP plot of endothelial cells colored by condition
umap_endothelial <- DimPlot(endothelial_cells, reduction = "umap", group.by = "condition", label = TRUE) +
  ggtitle("UMAP plot of Endothelial Cells by Condition")
print(umap_endothelial)
ggsave("umap_endothelial_cells_by_condition.png", umap_endothelial, width = 10, height = 8)


endothelial_cells <- subset(combined, subset = cell_type == "Endothelial_cells")

endothelial_cells <- NormalizeData(endothelial_cells)
endothelial_cells <- FindVariableFeatures(endothelial_cells, selection.method = "vst", nfeatures = 2000)
endothelial_cells <- ScaleData(endothelial_cells)
endothelial_cells <- RunPCA(endothelial_cells)
endothelial_cells <- RunUMAP(endothelial_cells, dims = 1:30)
endothelial_cells <- FindNeighbors(endothelial_cells, dims = 1:30)
endothelial_cells <- FindClusters(endothelial_cells, resolution = 0.5)

umap_plot <- DimPlot(endothelial_cells, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("UMAP of Endothelial Cells")
print(umap_plot)
ggsave("endothelial_clusters_umap.png", umap_plot, width = 10, height = 8)

# Table of cell counts per cluster and condition
cluster_condition_counts <- table(Idents(endothelial_cells), endothelial_cells$condition)

# Print the table
print("Number of cells per cluster in Morphea and Control:")
print(cluster_condition_counts)

# Optionally, save this table to a CSV file for further analysis or review
write.csv(as.data.frame(cluster_condition_counts), file = "cluster_condition_cell_counts.csv", row.names = TRUE)


#focus on 2 clusters that show major differences first (cluster 2 and cluster 7)

# Install clusterProfiler
BiocManager::install("clusterProfiler")

# Install org.Hs.eg.db
BiocManager::install("org.Hs.eg.db")

library(clusterProfiler)
library(org.Hs.eg.db)

# 1. Isolate clusters 2 and 7
clusters_of_interest <- subset(endothelial_cells, idents = c(2, 7))

# 2. Perform differential expression analysis
# Compare each cluster to all other cells
de_results_2 <- FindMarkers(clusters_of_interest, ident.1 = 2, ident.2 = NULL, only.pos = FALSE)
de_results_7 <- FindMarkers(clusters_of_interest, ident.1 = 7, ident.2 = NULL, only.pos = FALSE)

# 3. Conduct pathway enrichment analysis
perform_enrichment <- function(de_results, cluster_id) {
  # Get top differentially expressed genes
  top_genes <- rownames(de_results)[de_results$p_val_adj < 0.05 & abs(de_results$avg_log2FC) > 1]
  
  # Convert gene symbols to Entrez IDs
  entrez_ids <- mapIds(org.Hs.eg.db, keys = top_genes, keytype = "SYMBOL", column = "ENTREZID")
  entrez_ids <- entrez_ids[!is.na(entrez_ids)]
  
  # Perform GO enrichment analysis
  go_enrichment <- enrichGO(gene = entrez_ids,
                            OrgDb = org.Hs.eg.db,
                            ont = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.2)
  
  # Perform KEGG pathway analysis
  kegg_enrichment <- enrichKEGG(gene = entrez_ids,
                                organism = "hsa",
                                pAdjustMethod = "BH",
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.2)
  
  return(list(go = go_enrichment, kegg = kegg_enrichment))
}

enrichment_results_2 <- perform_enrichment(de_results_2, 2)
enrichment_results_7 <- perform_enrichment(de_results_7, 7)


# Volcano plot for cluster 2
volcano_plot_2 <- EnhancedVolcano(de_results_2,
                                  lab = rownames(de_results_2),
                                  x = 'avg_log2FC',
                                  y = 'p_val_adj',
                                  title = 'Cluster 2 vs Others',
                                  pCutoff = 0.05,
                                  FCcutoff = 1,
                                  pointSize = 3.0,
                                  labSize = 6.0)

# Volcano plot for cluster 7
volcano_plot_7 <- EnhancedVolcano(de_results_7,
                                  lab = rownames(de_results_7),
                                  x = 'avg_log2FC',
                                  y = 'p_val_adj',
                                  title = 'Cluster 7 vs Others',
                                  pCutoff = 0.05,
                                  FCcutoff = 1,
                                  pointSize = 3.0,
                                  labSize = 6.0)

# Save plots
ggsave("volcano_plot_cluster2.png", volcano_plot_2, width = 12, height = 10)
ggsave("volcano_plot_cluster7.png", volcano_plot_7, width = 12, height = 10)

# Print top GO terms
print(head(enrichment_results_2$go))
print(head(enrichment_results_7$go))

# Print top KEGG pathways
print(head(enrichment_results_2$kegg))
print(head(enrichment_results_7$kegg))

# Save detailed results
write.csv(de_results_2, "de_results_cluster2.csv")
write.csv(de_results_7, "de_results_cluster7.csv")
write.csv(as.data.frame(enrichment_results_2$go), "go_enrichment_cluster2.csv")
write.csv(as.data.frame(enrichment_results_2$kegg), "kegg_enrichment_cluster2.csv")
write.csv(as.data.frame(enrichment_results_7$go), "go_enrichment_cluster7.csv")
write.csv(as.data.frame(enrichment_results_7$kegg), "kegg_enrichment_cluster7.csv")






#check the up and down genes for pSTAT3 (tocilizumab issue)

# Upstream genes related to STAT3 activation
stat3_upstream_genes <- c(
  "JAK1", "JAK2", "JAK3", "TYK2",    # JAK Kinases
  "EGFR", "PDGFR",                     # Receptor Tyrosine Kinases
  "SRC", "FYN",                        # Src Family Kinases
  "IL6", "IL10", "IFNG", "EGF"         # Cytokines and Growth Factors
)

# Downstream genes regulated by STAT3
stat3_downstream_genes <- c(
  "MYC", "CCND1", "CDK2",              # Cell Cycle Regulators
  "BCL2", "MCL1",                      # Anti-apoptotic Proteins
  "SOCS3", "VEGFA", "IL6",             # Immune Modulators
  "MMP2", "MMP9",                       # Matrix Metalloproteinases
  "HIF1A"                               # Hypoxia-Inducible Factor 1-alpha
)

# Combined list of STAT3-related genes
stat3_genes <- c("STAT3", stat3_upstream_genes, stat3_downstream_genes)


# Function to get average expression of STAT3-related genes for each cluster
get_stat3_expression <- function(seurat_obj) {
  avg_exp <- AverageExpression(seurat_obj, features = stat3_genes, group.by = "seurat_clusters")
  return(as.data.frame(avg_exp$RNA))
}

# Get expression data for all clusters
all_clusters_stat3_exp <- get_stat3_expression(endothelial_cells)
