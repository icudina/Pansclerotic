# Load necessary libraries
library(Seurat)

# file paths
file_paths <- c("C:/Users/Brakebusch/Downloads/GSM7488058_2057-AB-5.h5", "C:/Users/Brakebusch/Downloads/GSM7488057_2057-AB-3.h5", "C:/Users/Brakebusch/Downloads/GSM7488054_2330-AB-3.h5", 
                "C:/Users/Brakebusch/Downloads/GSM7488055_4902-AB-1.h5", "C:/Users/Brakebusch/Downloads/GSM7488059_2057-AB-6.h5", "C:/Users/Brakebusch/Downloads/GSM7488053_2330-AB-1.h5", 
                "C:/Users/Brakebusch/Downloads/GSM7488056_4902-AB-2.h5")

# Load the data from each file into a list of Seurat objects
seurat_objects <- lapply(file_paths, function(file) {
  Read10X_h5(file)
})

# Print the structure of the first Seurat object to understand the data
str(seurat_objects[[1]])


# Create Seurat objects and assign groups
seurat_list <- list()
group_labels <- c("Control", "Control", "Morphea_nonlesional", "Morphea_lesional", 
                  "Control", "Morphea_lesional", "Morphea_nonlesional")

for (i in 1:length(file_paths)) {
  seurat_obj <- CreateSeuratObject(counts = seurat_objects[[i]], project = group_labels[i])
  seurat_obj$group <- group_labels[i]
  seurat_list[[i]] <- seurat_obj
}


# Add filtration step for low-quality cells
# Define mitochondrial genes pattern
mito_genes_pattern <- '^MT-'

# Filter cells based on the given criteria
filtered_seurat_list <- lapply(seurat_list, function(seurat_obj) {
  # Calculate percentage of mitochondrial genes
  seurat_obj[['percent.mt']] <- PercentageFeatureSet(seurat_obj, pattern = mito_genes_pattern)
  
  # Filter cells
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 100 & nCount_RNA > 500 & 
                         nCount_RNA < 1e5 & percent.mt < 15)
  return(seurat_obj)
})

# Normalize and find variable features for each filtered dataset
filtered_seurat_list <- lapply(X = filtered_seurat_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# Select features for integration
filtered_features <- SelectIntegrationFeatures(object.list = filtered_seurat_list)

# Find anchors and integrate filtered data
filtered_anchors <- FindIntegrationAnchors(object.list = filtered_seurat_list, anchor.features = filtered_features)
filtered_integrated_seurat <- IntegrateData(anchorset = filtered_anchors)

# Switch to integrated assay for downstream analysis
DefaultAssay(filtered_integrated_seurat) <- "integrated"

# Scale data, run PCA and UMAP
filtered_integrated_seurat <- ScaleData(filtered_integrated_seurat, verbose = FALSE)
filtered_integrated_seurat <- RunPCA(filtered_integrated_seurat, npcs = 30, verbose = FALSE)
filtered_integrated_seurat <- RunUMAP(filtered_integrated_seurat, reduction = "pca", dims = 1:30)

# Visualize the filtered integrated data
filtered_p1 <- DimPlot(filtered_integrated_seurat, reduction = "umap", group.by = "group")
filtered_p2 <- DimPlot(filtered_integrated_seurat, reduction = "umap", group.by = "orig.ident", label = TRUE, repel = TRUE)

# Print some summary statistics
print(table(filtered_integrated_seurat$group))
print(table(filtered_integrated_seurat$orig.ident))

# Save the plots
pdf("filtered_integrated_umap_plots.pdf", width = 12, height = 6)
print(filtered_p1)
print(filtered_p2)
dev.off()

# Display the plots
print(filtered_p1)
print(filtered_p2)







# Load required libraries
library(ggplot2)
library(dplyr)

# Perform differential gene expression analysis between the groups


# Define a function to perform DGE and plot results
perform_dge <- function(seurat_obj, group1, group2) {
  # Find differentially expressed genes
  dge_results <- FindMarkers(seurat_obj, ident.1 = group1, ident.2 = group2, min.pct = 0.25)
  
  # Print top 10 differentially expressed genes
  print(paste("Top 10 differentially expressed genes between", group1, "and", group2))
  print(head(dge_results, 10))
  
  # Create a volcano plot
  dge_results$gene <- rownames(dge_results)
  dge_results$significant <- dge_results$p_val_adj < 0.05
  
  # Plot
  p <- ggplot(dge_results, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significant)) +
    geom_point(alpha = 0.5) +
    theme_minimal() +
    labs(title = paste("DGE between", group1, "and", group2),
         x = "Average Log2 Fold Change",
         y = "-Log10 Adjusted P-Value") +
    scale_color_manual(values = c("grey", "red"))
  
  return(list(results = dge_results, plot = p))
}

# Perform DGE for Control vs Morphea_lesional
dge_control_vs_lesional <- perform_dge(filtered_integrated_seurat, "Control", "Morphea_lesional")

# Display the plot
print(dge_control_vs_lesional$plot)

# Save the plot
ggsave("dge_control_vs_lesional_volcano.pdf", dge_control_vs_lesional$plot, width = 10, height = 8)




# Perform DGE for Control vs Morphea_nonlesional
dge_control_vs_nonlesional <- perform_dge(filtered_integrated_seurat, "Control", "Morphea_nonlesional")

# Display the plot
print(dge_control_vs_nonlesional$plot)

# Save the plot
ggsave("dge_control_vs_nonlesional_volcano.pdf", dge_control_vs_nonlesional$plot, width = 10, height = 8)

# Print top 10 differentially expressed genes
print("Top 10 differentially expressed genes between Control and Morphea_nonlesional:")
print(head(dge_control_vs_nonlesional$results, 10))

# Perform DGE for Morphea_lesional vs Morphea_nonlesional
dge_lesional_vs_nonlesional <- perform_dge(filtered_integrated_seurat, "Morphea_lesional", "Morphea_nonlesional")

# Display the plot
print(dge_lesional_vs_nonlesional$plot)

# Save the plot
ggsave("dge_lesional_vs_nonlesional_volcano.pdf", dge_lesional_vs_nonlesional$plot, width = 10, height = 8)

# Print top 10 differentially expressed genes
print("Top 10 differentially expressed genes between Morphea_lesional and Morphea_nonlesional:")
print(head(dge_lesional_vs_nonlesional$results, 10))

# Calculate summary statistics
summary_stats <- list(
  control_vs_lesional = sum(dge_control_vs_lesional$results$p_val_adj < 0.05),
  control_vs_nonlesional = sum(dge_control_vs_nonlesional$results$p_val_adj < 0.05),
  lesional_vs_nonlesional = sum(dge_lesional_vs_nonlesional$results$p_val_adj < 0.05)
)

print("Number of significantly differentially expressed genes (adjusted p-value < 0.05):")
print(summary_stats)

# Annotate cell types using the expanded marker list
cell_markers <- list(
  Fibroblasts = c("COL1A1", "COL1A2", "DCN", "LUM", "FAP", "PDGFRA", "PDGFRB", "ACTA2"),
  Keratinocytes = c("KRT5", "KRT14", "KRT1", "KRT10", "DSG1", "FLG", "IVL"),
  Th1_cells = c("CD4", "IFNG", "TBX21", "CXCR3"),
  Th2_cells = c("CD4", "IL4", "GATA3", "IL5", "IL13"),
  Myeloid_cells = c("CD14", "CD68", "LYZ", "ITGAM", "CD33", "CSF1R"),
  B_cells = c("CD19", "MS4A1", "CD79A", "CD79B", "CD20", "CD22"),
  Endothelial_cells = c("PECAM1", "VWF", "CDH5", "FLT1", "KDR", "ENG"),
  Langerhans_cells = c("CD1A", "CD207", "CD1B", "CD1C"),
  Red_blood_cells = c("HBA1", "HBA2", "HBB"),
  Melanocytes = c("TYRP1", "DCT", "MLANA", "MITF", "PMEL"),
  Mast_cells = c("TPSAB1", "CPA3", "KIT", "FCER1A"),
  Pericytes = c("PDGFRB", "RGS5", "ACTA2", "CSPG4", "MCAM"),
  Smooth_muscle_cells = c("ACTA2", "MYH11", "CNN1", "TAGLN"),
  Treg_cells = c("CD4", "FOXP3", "IL2RA", "CTLA4", "IKZF2"),
  Th17_cells = c("CD4", "RORC", "IL17A", "CCR6", "IL23R"),
  CD8_T_cells = c("CD3D", "CD3E", "CD3G", "CD8A", "CD8B", "GZMA", "GZMB", "PRF1")
)

# Function to safely add module score
safe_add_module_score <- function(seurat_obj, cell_type, markers) {
  tryCatch({
    seurat_obj <- AddModuleScore(
      object = seurat_obj,
      features = list(markers),
      name = cell_type,
      ctrl = min(vapply(cell_markers, length, FUN.VALUE = numeric(1))) # Use the minimum number of markers across all cell types
    )
    return(seurat_obj)
  }, error = function(e) {
    message(paste("Error adding module score for", cell_type, ":", e$message))
    return(seurat_obj)
  })
}

# Add cell type scores to the Seurat object
for (cell_type in names(cell_markers)) {
  filtered_integrated_seurat <- safe_add_module_score(filtered_integrated_seurat, cell_type, cell_markers[[cell_type]])
}

# Check which cell types were successfully added
added_cell_types <- grep("^Fibroblasts|Keratinocytes|Th1_cells|Th2_cells|Myeloid_cells|B_cells|Endothelial_cells|Langerhans_cells|Red_blood_cells|Melanocytes|Mast_cells|Pericytes|Smooth_muscle_cells|Treg_cells|Th17_cells|CD8_T_cells", 
                         colnames(filtered_integrated_seurat@meta.data), value = TRUE)

print("Successfully added cell types:")
print(added_cell_types)

# Assign cell type based on the highest score
cell_type_scores <- as.data.frame(filtered_integrated_seurat@meta.data[, added_cell_types])
filtered_integrated_seurat$cell_type <- apply(cell_type_scores, 1, function(x) {
  cell_type <- sub("1$", "", names(which.max(x)))
  return(cell_type)
})

# Visualize the cell type annotations
cell_type_plot <- DimPlot(filtered_integrated_seurat, reduction = "umap", group.by = "cell_type", label = TRUE, repel = TRUE) +
  ggtitle("Cell Type Annotations")

# Save and display the plot
ggsave("cell_type_annotation_umap.pdf", cell_type_plot, width = 12, height = 10)
print(cell_type_plot)

# Print cell type proportions
cell_type_proportions <- prop.table(table(filtered_integrated_seurat$cell_type))
print("Cell type proportions:")
print(cell_type_proportions)

# Visualize cell type proportions
cell_type_prop_plot <- ggplot(data.frame(cell_type = names(cell_type_proportions), proportion = as.numeric(cell_type_proportions)), 
                              aes(x = reorder(cell_type, -proportion), y = proportion)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Cell Type", y = "Proportion", title = "Cell Type Proportions")

# Save and display the proportion plot
ggsave("cell_type_proportions.pdf", cell_type_prop_plot, width = 12, height = 8)
print(cell_type_prop_plot)





# required libraries
library(Seurat)
library(ggplot2)
library(dplyr)

# Subset to endothelial cells based on the cell type annotation
endothelial_cells <- subset(filtered_integrated_seurat, subset = cell_type == "Endothelial_cells")

# Confirm the number of endothelial cells
print(paste("Number of endothelial cells:", ncol(endothelial_cells)))


# Set the default assay to RNA
DefaultAssay(endothelial_cells) <- "RNA"

# Re-normalize the data
endothelial_cells <- NormalizeData(endothelial_cells)

# Find variable features
endothelial_cells <- FindVariableFeatures(endothelial_cells, selection.method = "vst", nfeatures = 2000)

# Scale the data
endothelial_cells <- ScaleData(endothelial_cells, features = rownames(endothelial_cells))

# Run PCA
endothelial_cells <- RunPCA(endothelial_cells, features = VariableFeatures(object = endothelial_cells))

# Determine the number of PCs to use (you can adjust the number based on the ElbowPlot)
ElbowPlot(endothelial_cells)

# Run UMAP
endothelial_cells <- RunUMAP(endothelial_cells, dims = 1:12)

# Find neighbors
endothelial_cells <- FindNeighbors(endothelial_cells, dims = 1:12)

# Find clusters
endothelial_cells <- FindClusters(endothelial_cells, resolution = 0.5)

# Visualize UMAP with clusters
p1 <- DimPlot(endothelial_cells, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE) +
  ggtitle("Endothelial Cells Clustering")

# Visualize UMAP with group labels (e.g., Control, Morphea_lesional, etc.)
p2 <- DimPlot(endothelial_cells, reduction = "umap", group.by = "group") +
  ggtitle("Endothelial Cells by Group")

# Save the plots
ggsave("endothelial_cells_umap_clusters.pdf", p1, width = 8, height = 6)
ggsave("endothelial_cells_umap_groups.pdf", p2, width = 8, height = 6)

# Display the plots
print(p1)
print(p2)


# Join data layers in the Seurat object
endothelial_cells <- JoinLayers(endothelial_cells)

# Now, identify cluster markers
cluster_markers <- FindAllMarkers(
  endothelial_cells,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)


# View top markers for each cluster
top_markers <- cluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)
print(top_markers)

# Save the markers
write.csv(cluster_markers, "endothelial_cells_cluster_markers.csv")


# Load necessary libraries
library(Seurat)
library(ggplot2)
library(dplyr)  # For the arrange() function

# Set the identity classes to 'group' for group-wise comparison
Idents(endothelial_cells) <- endothelial_cells$group

# Differential expression: Morphea_lesional vs Control
dge_lesional_vs_control <- FindMarkers(
  endothelial_cells,
  ident.1 = "Morphea_lesional",
  ident.2 = "Control",
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# Sort the results by adjusted p-value
de_results <- arrange(dge_lesional_vs_control, p_val_adj)

# Display the top 20 differentially expressed genes
print("Top 20 differentially expressed genes between Morphea_lesional and Control in Endothelial Cells:")
print(head(de_results, 20))

# Save the full results to a CSV file
write.csv(de_results, file = "morphea_lesional_vs_control_DGE.csv", row.names = TRUE)

print("Differential expression analysis completed. Results saved to 'morphea_lesional_vs_control_DGE.csv'")

# Create a volcano plot for the DGE results
de_results$gene <- rownames(de_results)
de_results$significant <- de_results$p_val_adj < 0.05

volcano_plot_lesional <- ggplot(de_results, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significant)) +
  geom_point(alpha = 0.5) +
  theme_minimal() +
  labs(
    title = "DGE: Morphea_lesional vs Control in Endothelial Cells",
    x = "Average Log2 Fold Change",
    y = "-Log10 Adjusted P-Value"
  ) +
  scale_color_manual(values = c("grey", "red")) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank()
  )

# Save and display the volcano plot
ggsave("endothelial_cells_dge_lesional_vs_control_volcano.pdf", volcano_plot_lesional, width = 10, height = 8)
print(volcano_plot_lesional)



# Load necessary libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)

# Read the differential expression results
de_results <- read.csv("morphea_lesional_vs_control_DGE.csv", row.names = 1)

# Prepare the gene list for enrichment analysis
gene_list <- de_results$avg_log2FC
names(gene_list) <- rownames(de_results)
gene_list <- sort(gene_list, decreasing = TRUE)

# Convert gene symbols to Entrez IDs
gene_ids <- bitr(names(gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Remove any duplicate IDs
gene_ids <- gene_ids[!duplicated(gene_ids$SYMBOL),]

# Create a named vector of fold changes
gene_list_entrez <- gene_list[gene_ids$SYMBOL]
names(gene_list_entrez) <- gene_ids$ENTREZID

# Perform GO enrichment analysis
go_results <- enrichGO(gene = names(gene_list_entrez),
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

# Display top GO terms
print("Top GO terms:")
print(head(go_results))

# Save GO results
write.csv(as.data.frame(go_results), file = "GO_enrichment_results.csv")

# Perform GSEA with Hallmark gene sets
hallmark <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- split(hallmark$gene_symbol, hallmark$gs_name)

gsea_results <- GSEA(gene_list, TERM2GENE = hallmark[, c("gs_name", "gene_symbol")])

# Display top GSEA results
print("Top GSEA Hallmark results:")
print(head(gsea_results))

# Save GSEA results
write.csv(as.data.frame(gsea_results), file = "GSEA_Hallmark_results.csv")

# Perform KEGG pathway analysis
kegg_results <- enrichKEGG(gene = names(gene_list_entrez),
                           organism = "hsa",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.2)

# Display top KEGG pathways
print("Top KEGG pathways:")
print(head(kegg_results))

# Save KEGG results
write.csv(as.data.frame(kegg_results), file = "KEGG_pathway_results.csv")

print("Enrichment analyses completed. Results saved to CSV files.")

# Generate plots
pdf("enrichment_plots.pdf", width = 12, height = 8)

# GO dotplot
print(dotplot(go_results, showCategory = 20, title = "GO Biological Process"))

# GSEA plot for top pathway
top_pathway <- gsea_results$ID[1]
print(gseaplot2(gsea_results, geneSetID = top_pathway, title = paste("GSEA -", top_pathway)))

# KEGG dotplot
print(dotplot(kegg_results, showCategory = 20, title = "KEGG Pathways"))

dev.off()

print("Plots saved to 'enrichment_plots.pdf'")


# Load necessary libraries
library(Seurat)
library(dplyr)
library(ggplot2)

# Assuming 'endothelial_cells' is your Seurat object containing endothelial cells

# Get the number of cells in each cluster
cluster_counts <- table(endothelial_cells$seurat_clusters)
print("Number of cells in each cluster:")
print(cluster_counts)

# Get the number of cells in each group
group_counts <- table(endothelial_cells$group)
print("Number of cells in each group:")
print(group_counts)

# Get the number of cells in each cluster for each group
cluster_group_counts <- table(endothelial_cells$seurat_clusters, endothelial_cells$group)
print("Number of cells in each cluster for each group:")
print(cluster_group_counts)


