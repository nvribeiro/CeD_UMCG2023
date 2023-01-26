# Cluster indentification - batch A

# Load libraries
library(Seurat)
library(tidyverse)
library(SeuratDisk)
library(patchwork)

# Loading function to plot markers
source('./functions/my_functions.R')

A_obj <- LoadH5Seurat('./data_preprocessing/outputs/A/SeuratObj_A_merged_genotyped_filtered_clustered.h5seurat')

DimPlot(A_obj, reduction = 'umap', label = TRUE)
VlnPlot(A_obj, features = 'percent.mt')
VlnPlot(A_obj, features = 'nFeature_RNA')
VlnPlot(A_obj, features = 'nCount_RNA')

# Finding cluster markers
A.markers <- FindAllMarkers(A_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Save
write.csv(A.markers, './data_preprocessing/outputs/A/A_all_markers_clusters.csv', row.names = F)

# Getting the top 10 markers per cluster
A.top.markers <- A.markers %>% 
  group_by(cluster) %>% 
  arrange(p_val_adj, by.group = TRUE) %>%
  top_n(n = 10, wt = avg_log2FC)

# Reshaping the table
A.top.markers.clean <- A.top.markers %>%
  select(cluster, gene) %>%
  group_by(cluster) %>%
  summarise(markers = str_c(unlist(cur_data()), collapse=', '))

# Saving
write.csv(A.top.markers.clean, './data_preprocessing/outputs/A/A_clusters_top10_markers_clusters.csv', row.names = FALSE)

## Cluster identification using UCell
epithelial_markers <- '../../resources/cell_type_markers_small_intestine_epithelial.csv'
other_markers <- '../../resources/cell_type_markers_small_intestine_immune_and_others.csv'

PlotUCell(A_obj, epithelial_markers, output = './data_preprocessing/outputs/plots/A_epithelial_clusters_Ucell')
PlotUCell(A_obj, other_markers, output = './data_preprocessing/outputs/plots/A_immune_others_clusters_Ucell')

# Manual inspestions
FeaturePlot(A_obj, features = c('MUC3A', 'MUC5B'))


VlnPlot(A_obj, features = c('KCNJ8', 'ABCC9', 'RGS5', 'NCAM1'))
FeaturePlot(A_obj, features = c('JCHAIN'), label = TRUE)
VlnPlot(A_obj, features = 'percent.mt') + NoLegend() + geom_hline(yintercept = 10) + geom_hline(yintercept = 20) + geom_hline(yintercept = 75)
FeaturePlot(A_obj, features = 'percent.mt', label = TRUE, repel = TRUE)

DimPlot(A_obj, reduction = 'umap', group.by = 'genotype') + NoLegend()

tmp <- subset(A_obj, subset = seurat_clusters == 0)
FeaturePlot(tmp, features = 'percent.mt')
VlnPlot(tmp, features = 'percent.mt')


pdf('./data_preprocessing/outputs/plots/A_UMAP_genotyped_filtered.pdf', width = 11, height = 8, paper = 'a4r')
DimPlot(A_obj, reduction = 'umap', label = TRUE, repel = TRUE)
dev.off()

# Annotating the clusters
clusters.ids <- read.csv('./data_preprocessing/outputs/A/A_clusters_top10_markers_clusters.csv')
clusters.ids <- clusters.ids$sub.type.1
names(clusters.ids) <- levels(A_obj)
A_obj <- RenameIdents(A_obj, clusters.ids)

pdf('./data_preprocessing/outputs/plots/A_UMAP_genotyped_filtered_identified.pdf', width = 11, height = 8, paper = 'a4r')
DimPlot(A_obj, reduction = 'umap', label = TRUE, repel = TRUE) + NoLegend()
dev.off()