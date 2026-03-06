library(AUCell)
library(clusterProfiler)
library(ggplot2)
library(Seurat)
library(GSEABase)
library(Seurat)
library(dplyr)
library(tidyr)


myeloid_celltypes <- c("Alveolar_macrophage", "Monocyte", "pDC", 
                       "Interstitial_macrophage", "Neutrophil")
scRNA_myeloid <- subset(scRNA, subset = Celltype %in% myeloid_celltypes)
table(scRNA_myeloid$Celltype)
DimPlot(scRNA_myeloid, group.by = "Celltype", label = TRUE) 
print(comparison_results)
go_df <- read.csv("C:/Users/2023/Desktop/GO.csv")
all_genes <- rownames(GetAssayData(scRNA_myeloid, layer = "counts"))
gene_list <- unique(go_df$Gene[go_df$Gene %in% all_genes])

cells_rankings <- AUCell_buildRankings(
  GetAssayData(scRNA_myeloid, layer = "counts"), 
  plotStats = FALSE
)
gene_set <- GeneSet(gene_list, setName = "Myeloid_GO_Signature")
cells_AUC <- AUCell_calcAUC(GeneSetCollection(gene_set), cells_rankings)
scRNA_myeloid$AUC_Score <- as.numeric(getAUC(cells_AUC)["Myeloid_GO_Signature", ])

#FDR
library(dplyr)
library(tidyr)

plot_df <- data.frame(
  scRNA_myeloid@meta.data,
  UMAP_1 = scRNA_myeloid@reductions$umap@cell.embeddings[, 1],
  UMAP_2 = scRNA_myeloid@reductions$umap@cell.embeddings[, 2]
)

stat_results <- plot_df %>%
  group_by(Celltype) %>%
  summarise(
    p_value_raw = suppressWarnings(
      wilcox.test(AUC_Score ~ Group, exact = FALSE)$p.value
    )
  ) %>%
  ungroup() %>%
  mutate(
    p_value_adj = p.adjust(p_value_raw, method = "BH"),
    significance = case_when(
      p_value_adj < 0.001 ~ "***",
      p_value_adj < 0.01 ~ "**",
      p_value_adj < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    y_position = sapply(Celltype, function(ct) {
      max(plot_df$AUC_Score[plot_df$Celltype == ct], na.rm = TRUE) * 1.08
    })
  )

print(stat_results)

library(ggplot2)
library(gghalves)

mytheme <- theme(
  axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
  axis.text.y = element_text(size = 12),
  axis.title = element_text(size = 13),
  legend.text = element_text(size = 12),
  legend.title = element_text(size = 12),
  axis.line = element_line(size = 0.7),
  panel.border = element_blank(),
  panel.grid = element_blank()
)

p <- ggplot() +
  geom_half_violin(
    data = plot_df %>% filter(Group == "NC"),
    aes(x = Celltype, y = AUC_Score),
    fill = "#2372A9", side = "l", color = NA, alpha = 0.7
  ) +
  geom_half_violin(
    data = plot_df %>% filter(Group == "KP"),
    aes(x = Celltype, y = AUC_Score),
    fill = "#CA2A28", side = "r", color = NA, alpha = 0.7
  ) +
  geom_point(
    data = plot_df,
    aes(x = Celltype, y = AUC_Score, fill = Group),
    stat = "summary", fun = mean,
    position = position_dodge(width = 0.8), shape = 21, size = 2
  ) +
  stat_summary(
    data = plot_df,
    aes(x = Celltype, y = AUC_Score, group = Group),
    fun.min = function(z) { quantile(z, 0.25) },
    fun.max = function(z) { quantile(z, 0.75) },
    geom = "errorbar", width = 0.1,
    position = position_dodge(width = 0.8)
  ) +
  geom_text(
    data = stat_results,
    aes(x = Celltype, y = y_position, label = significance),
    size = 6, vjust = 0.5
  ) +
  mytheme +
  theme_bw() +
  labs(y = "AUC Score (Myeloid Gene Signature)", x = "") +
  scale_fill_manual(values = c("NC" = "#2372A9", "KP" = "#CA2A28"))

print(p)










