# Cluster indentification - batch F

# Load libraries
library(Seurat)
library(tidyverse)
library(SeuratDisk)
library(patchwork)

# Loading function to plot markers
source('./functions/my_functions.R')

F_obj <- LoadH5Seurat('./outputs/F/SeuraObj_f_reclustered_res04.h5seurat')
F_obj <- LoadH5Seurat('./outputs/F/SeuratObj_F_singlets_genotyped_clustered_res04_harmony.h5seurat')

# Finding cluster markers
F.markers <- FindAllMarkers(F_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
F.top.markers <- F.markers %>% 
  group_by(cluster) %>% 
  arrange(p_val_adj, by.group = TRUE) %>%
  top_n(n = 10, wt = avg_log2FC)

# Reshaping the table
F.top.markers.clean <- F.top.markers %>%
  select(cluster, gene) %>%
  group_by(cluster) %>%
  summarise(markers = str_c(unlist(cur_data()), collapse=', '))

write.csv(F.top.markers.clean, './outputs/F/F_clusters_top_markers_reclustered.csv', row.names = FALSE)

# Inspecting the clusters with my expected markers and cell types
plotMarkers(F_obj, markers_csv = '../../resources/cell_markers.csv', output = './outputs/plots/F_')

# QC plots
pdf('./outputs/plots/F_QCplots_reclustered.pdf', width = 11, height = 8, paper = 'a4r')
VlnPlot(F_obj, features = 'percent.mt')
VlnPlot(F_obj, features = 'nCount_RNA')
VlnPlot(F_obj, features = 'nFeature_RNA')
FeaturePlot(F_obj, features = 'percent.mt')
dev.off()

# Inspections
VlnPlot(F_obj, features = 'percent.mt') & geom_hline(yintercept = 65) # for cluster 2

VlnPlot(F_obj, features = 'EPCAM')
VlnPlot(F_obj, features = 'TIMD4')
VlnPlot(F_obj, features = 'BEST4')
VlnPlot(F_obj, features = 'NCR3')
VlnPlot(F_obj, features = 'LPP')

FeaturePlot(F_obj, features = c('ITGAE', 'CD8A', 'CD4', 'GZMB', 'PRF1'))
FeaturePlot(F_obj, features = c('TGM2', 'IFNG'), split.by = 'status')
FeaturePlot(F_obj, features = c('VIM', 'APOA4', 'MKI67'))
FeaturePlot(F_obj, features = c( 'LPP'), split.by = 'status')
FeaturePlot(F_obj, features = c('TXNIP', 'IL32', 'CD52', 'S100A4'), split.by = 'status')

VlnPlot(F_obj, features = 'percent.mt') + geom_hline(yintercept = 10) + geom_hline(yintercept = 15) + geom_hline(yintercept = 65) + geom_hline(yintercept = 75) + geom_hline(yintercept = 80) + geom_hline(yintercept = 30) + NoLegend()

# Assigning cluster IDs
clusters.IDs_df <- read.csv('./outputs/F/F_clusters_ID.csv')
clusters.IDs <- clusters.IDs_df$cell.type

names(clusters.IDs) <- levels(F_obj)
F_obj <- RenameIdents(F_obj, clusters.IDs)

png('./outputs/plots/F_clustered_identified.png', width = 600, height = 600)
DimPlot(F_obj, reduction='umap', label = TRUE, repel = TRUE) + NoLegend() + labs(title = 'Batch F clusters - resolution = 0.4')
dev.off()

cells.per.cluster <- as.data.frame(F_obj@active.ident)
p <- ggplot(cells.per.cluster, aes(F_obj@active.ident, fill=F_obj@active.ident)) +
  geom_bar()
p

# Save this object
SaveH5Seurat(F_obj, './outputs/F/SeuratObj_F_clustered_res04_identified')

## New identification with clusters filtered -----------------------------------
F_obj <- LoadH5Seurat('./data_preprocessing/outputs/F/SeuraObj_f_reclustered_res04.h5seurat')

DimPlot(F_obj, reduction = 'umap', label = TRUE)
VlnPlot(F_obj, features = 'percent.mt')
VlnPlot(F_obj, features = 'nFeature_RNA')
VlnPlot(F_obj, features = 'nCount_RNA')

# Finding cluster markers
F.markers <- FindAllMarkers(F_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Save
write.csv(F.markers, './data_preprocessing/outputs/F/F_all_markers_reclustered.csv', row.names = F)

# Getting the top 10 markers per cluster
F.top.markers <- F.markers %>% 
  group_by(cluster) %>% 
  arrange(p_val_adj, by.group = TRUE) %>%
  top_n(n = 10, wt = avg_log2FC)

# Reshaping the table
F.top.markers.clean <- F.top.markers %>%
  select(cluster, gene) %>%
  group_by(cluster) %>%
  summarise(markers = str_c(unlist(cur_data()), collapse=', '))

# Saving
write.csv(F.top.markers.clean, './data_preprocessing/outputs/F/F_clusters_top10_markers_reclustered.csv', row.names = FALSE)

## Cluster identification using UCell
epithelial_markers <- '../../resources/cell_type_markers_small_intestine_epithelial.csv'
other_markers <- '../../resources/cell_type_markers_small_intestine_immune_and_others.csv'

PlotUCell(F_obj, epithelial_markers, output = './data_preprocessing/outputs/plots/F_epithelial_clusters_Ucell')
PlotUCell(F_obj, other_markers, output = './data_preprocessing/outputs/plots/F_immune_others_clusters_Ucell')

# Clusters inspection
FeaturePlot(F_obj, features = c('GJA1', 'SOX10'))
VlnPlot(F_obj, features = c('GJA1', 'SOX10'))
