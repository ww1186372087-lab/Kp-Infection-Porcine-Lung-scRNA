##For example ATII and RAS cell

library(monocle3)
library(Seurat)
library(SeuratObject)
library(tidyselect)
library(dplyr)
library(BiocGenerics)


RASATII = scRNA_NEW[,scRNA_NEW@meta.data$seurat_clusters  %in% c(12,15)]
Idents(RASATII) <- RASATII$Group
RASATII_NC= RASATII[, Idents(RASATII) %in% c( "NC")] 
RASATII_KP= RASATII[, Idents(RASATII) %in% c( "KP")] 

data <- GetAssayData(RASATII, assay = 'RNA', layer = 'counts')

cell_metadata <- RASATII@meta.data

gene_annotation <- data.frame(
  gene_short_name = rownames(data), 
  row.names = rownames(data)        
)


cds <- new_cell_data_set(
  expression_data = data,
  cell_metadata = cell_metadata,
  gene_metadata = gene_annotation
)


cds <- preprocess_cds(cds, num_dim = 10)

p_pca <- plot_pc_variance_explained(cds)

cds <- reduce_dimension(cds, preprocess_method = "PCA") 

p_umap <- plot_cells(cds, reduction_method = "UMAP", color_cells_by = "seurat_clusters") +
  ggtitle('Monocle3 UMAP by Seurat Clusters') +
  theme(plot.title = element_text(hjust = 0.5))

cds <- cluster_cells(cds, resolution = 0.5) 


cds <- learn_graph(cds)

p_trajectory <- plot_cells(
  cds, 
  color_cells_by = "Cellstyle", 
  label_groups_by_cluster = FALSE,
  label_leaves = TRUE, 
  label_branch_points = TRUE,
  graph_label_size = 3
) + ggtitle('Cell Trajectory by Cellstyle') +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("results/Cell_Trajectory.pdf", p_trajectory, width = 10, height = 8)

cds <- order_cells(cds)

p_pseudotime <- plot_cells(
  cds, 
  color_cells_by = "pseudotime", 
  label_cell_groups = FALSE, 
  label_leaves = FALSE,  
  label_branch_points = FALSE,
  cell_size = 0.6
) + ggtitle('Cell Pseudotime Distribution') +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("results/Pseudotime.pdf", p_pseudotime, width = 8, height = 6)