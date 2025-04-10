library(dplyr)
library(Seurat)
library(SeuratData)
library(future)
library(sctransform)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(fgsea)
library(viridis)
library(scales)

setwd("/gpfs/commons/groups/sanjana_lab/Cas13/Liver_Atlas/cellranger")

# Load the data
# Define the file mapping
samples <- data.frame(
  Run = c("SRR17374977", "SRR17374981", "SRR17374985",
          "SRR17374989", "SRR17374997", "SRR17374998", "SRR17374999", "SRR17375000",
          "SRR17375001", "SRR17375002", "SRR17375003", "SRR17375004", "SRR17375005", "SRR17375006",
          "SRR17375007", "SRR17375008", "SRR17375009", "SRR17375010"),
  Filename = c("CS118", "CS120", "CS122",
               "CS124", "CS182", "CS183", "CS184", "CS185",
               "CS188", "CS189", "CS190", "CS191", "CS192", "CS193",
               "CS194", "CS195", "CS196", "CS197")
)

base_dir <- "mouse"

# List to store individual Seurat objects
seurat_objects <- list()

# Loop through the samples to create Seurat objects
for (i in 1:nrow(samples)) {
  run <- samples$Run[i]
  project_name <- samples$Filename[i]
  h5_file <- file.path(base_dir, paste0("cellranger_", run, "/outs/filtered_feature_bc_matrix.h5"))

  # Check if the file exists
  if (file.exists(h5_file)) {
    print(paste("Reading H5 file for run:", run, "with project name:", project_name))
    new_data <- Read10X_h5(h5_file)

    # Create a Seurat object with a unique project name
    new_seurat <- CreateSeuratObject(new_data, project = project_name)

    # Store the Seurat object in the list
    seurat_objects[[run]] <- new_seurat
  } else {
    print(paste("File not found for run:", run, "- Skipping"))
  }
}

# Metadata table for samples
samples <- data.frame(
  Run = c("SRR17374977", "SRR17374981", "SRR17374985",
          "SRR17374989", "SRR17374997", "SRR17374998", "SRR17374999", "SRR17375000",
          "SRR17375001", "SRR17375002", "SRR17375003", "SRR17375004", "SRR17375005", "SRR17375006",
          "SRR17375007", "SRR17375008", "SRR17375009", "SRR17375010"),
  Assay = c("v3", "v3", "v3", "v3",
            "v3.1", "v3.1", "v3.1", "v3.1",
            "v3.1", "v3.1", "v3.1", "v3.1",
            "v3.1", "v3.1", "v3.1", "v3.1",
            "v3.1", "v3.1"),
  Molecule = c("Total RNA", "Total RNA", "Total RNA", "Total RNA",
               "Nuclear RNA", "Nuclear RNA", "Nuclear RNA", "Nuclear RNA",
               "Nuclear RNA", "Nuclear RNA", "Nuclear RNA", "Nuclear RNA",
               "Nuclear RNA", "Nuclear RNA", "Nuclear RNA", "Nuclear RNA",
               "Nuclear RNA", "Nuclear RNA"),
  Description = c(
    "Standard diet fed (36 weeks) Liver CD45- cells_Mouse_001_Sample_001",
    "Western diet fed (36 weeks) Liver CD45- cells_Mouse_001_Sample_001",
    "Standard diet fed (24 weeks) Liver CD45- cells_Mouse_003_Sample_001",
    "Western diet fed (24 weeks) Liver CD45- cells_Mouse_003_Sample_001",
    "Standard diet fed (24 weeks) Liver nuclei_Mouse_004_Sample_001",
    "Western diet fed (24 weeks) Liver nuclei_Mouse_004_Sample_001",
    "Standard diet fed (24 weeks) Liver nuclei_Mouse_005_Sample_001",
    "Western diet fed (24 weeks) Liver nuclei_Mouse_005_Sample_001",
    "Standard diet fed (24 weeks) Liver nuclei_Mouse_006_Sample_001",
    "Standard diet fed (24 weeks) Liver nuclei_Mouse_007_Sample_001",
    "Standard diet fed (24 weeks) Liver nuclei_Mouse_008_Sample_001",
    "Western diet fed (24 weeks) Liver nuclei_Mouse_006_Sample_001",
    "Western diet fed (24 weeks) Liver nuclei_Mouse_007_Sample_001",
    "Western diet fed (24 weeks) Liver nuclei_Mouse_008_Sample_001",
    "Standard diet fed (36 weeks) Liver nuclei_Mouse_003_Sample_001",
    "Standard diet fed (36 weeks) Liver nuclei_Mouse_004_Sample_001",
    "Western diet fed (36 weeks) Liver nuclei_Mouse_003_Sample_001",
    "Western diet fed (36 weeks) Liver nuclei_Mouse_004_Sample_001"
  )
)

# Add metadata to Seurat objects and consolidate diet/duration metadata
for (run in samples$Run) {
  if (run %in% names(seurat_objects)) {
    # Extract metadata for the current run
    sample_data <- samples[samples$Run == run, ]
    description <- sample_data$Description

    # Extract diet and duration from the description
    diet <- ifelse(grepl("^Standard", description), "SD", "WD")
    duration <- ifelse(grepl("\\(36 weeks\\)", description), "36wk", "24wk")

    # Add metadata to the Seurat object
    seurat_objects[[run]]$assay <- sample_data$Assay
    seurat_objects[[run]]$molecule <- sample_data$Molecule
    seurat_objects[[run]]$description <- description
    seurat_objects[[run]]$diet <- diet
    seurat_objects[[run]]$duration <- duration
  }
}

# Combine all Seurat objects into a single object
combined <- Reduce(function(x, y) merge(x, y), seurat_objects)
combined[["RNA"]] <- JoinLayers(combined[["RNA"]])


## Filter cells
DefaultAssay(combined) <- "RNA"
combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("ASD_QC_VlnPlot.pdf", width = 16, height = 10)


# Visualize QC metrics as a scatter plot
plot1 <- FeatureScatter(combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
# save
ggsave("ASD_QC_scatter.pdf", width = 16, height = 10)

# Low quality cell filtering
combined <- subset(combined, subset = nCount_RNA > 100 & nFeature_RNA > 200)
combined <- subset(combined, subset = percent.mt < 20)

####################################################################################################

#split the layer for different assay/molecule
combined[["RNA"]] <- split(combined[["RNA"]], f = combined$molecule)


setwd("/gpfs/commons/home/jameslee/Cas13/liver_atlas/Seurat")
saveRDS(combined, file = "combined_before_preprocessing.rds")

# run standard anlaysis workflow
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined)
combined <- ScaleData(combined)
combined <- RunPCA(combined)
combined <- FindNeighbors(combined, dims = 1:30, reduction = "pca")
combined <- FindClusters(combined, resolution = 2, cluster.name = "unintegrated_clusters")
combined <- RunUMAP(combined, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
#DimPlot(combined, reduction = "umap.unintegrated", group.by = c("stim", "seurat_clusters"))

saveRDS(combined, file = "combined_after_normalization.rds")

# Use Integration function
temp <- IntegrateLayers(object = combined, method = CCAIntegration, normalization.method = "SCT", verbose = F)
temp <- FindNeighbors(temp, reduction = "integrated.dr", dims = 1:30)
temp <- FindClusters(temp, resolution = 0.6)

# Reassign temp to combined
combined <- temp

# Save the new combined object
saveRDS(combined, file = "combined_after_integration.rds")

# Generate UMAPs grouped by sample name, molecule, duration, and diet
group_by_metadata <- c("sample", "molecule", "duration", "diet")

# Save UMAP plots as PDFs
for (group in group_by_metadata) {
  pdf_filename <- paste0("UMAP_", group, ".pdf")
  pdf(pdf_filename)
  DimPlot(combined, reduction = "umap", group.by = group) +
    ggtitle(paste("UMAP grouped by", group))
  dev.off()
}




####### DEG
# combined[["RNA"]] <- JoinLayers(combined[["RNA"]])

DefaultAssay(combined)

# Plot Ldlr expression on UMAP

p1 <- FeaturePlot(combined, features = c("Ldlr"))
ggsave("Ldlr.pdf", p1, width = 10, height = 6)

# Plot Ldlr expression by diet using Violin Plot
pdf("Ldlr_By_Diet.pdf")
VlnPlot(combined, features = c("Ldlr"), group.by = "diet")
dev.off()

# Perform DEG analysis comparing SD and WD
deg_results <- FindMarkers(
  object = combined,
  ident.1 = "SD",
  ident.2 = "WD",
  group.by = "diet",
  test.use = "wilcox",
  logfc.threshold = 0 # No logFC threshold
)

# Save DEG results as CSV for record
write.csv(deg_results, file = "DEG_SD_vs_WD.csv", row.names = TRUE)

# Prepare DEG data for heatmap
library(pheatmap)
top_degs <- rownames(deg_results[order(deg_results$p_val_adj), ][1:50, ])  # Top 50 DEGs
heatmap_data <- FetchData(combined, vars = top_degs)

# Generate the heatmap
pdf("DEG_SD_vs_WD_Heatmap.pdf")
pheatmap(
  mat = heatmap_data,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  main = "Top 50 DEGs (SD vs WD)"
)
dev.off()



# Add significance columns
deg_results$Significant <- with(deg_results, ifelse(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25, "Significant", "Not Significant"))

# Create Volcano Plot
volcano_plot <- ggplot(deg_results, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Significant)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
  theme_minimal() +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  labs(
    title = "Volcano Plot: DEG SD vs WD",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Significance"
  ) +
  theme(
    text = element_text(size = 14),
    legend.position = "top"
  )

# Save the plot
pdf("Volcano_Plot_SD_vs_WD.pdf", width = 8, height = 6)
print(volcano_plot)
dev.off()



