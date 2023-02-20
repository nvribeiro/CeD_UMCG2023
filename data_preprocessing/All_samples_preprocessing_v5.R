## Checking the new clustering after removing extra doublets -------------------
# Dataset v5

# Load libraries
library(Seurat)
library(tidyverse)
library(SeuratDisk)
library(patchwork)

source('./functions/my_functions.R')

path <- './data_preprocessing/outputs/all/'

# Clustering after removing doublets
# Running in the cluster
obj_v5 <- readRDS(paste0(path, 'All_samples_v5.rds'))

DefaultAssay(obj_v5) <- 'RNA'
obj_v5 <- obj_v5 %>% 
  SCTransform(vars.to.regress = "percent.mt", verbose = F) %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:50) %>%
  FindClusters(resolution = 0.4) %>%
  RunUMAP(dims = 1:50)

SaveH5Seurat(obj_v5, filename = paste0(path, 'All_samples_clustered_v5'), overwrite = T)

# Finding markers
All_markers <- FindAllMarkers(obj_v5, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(All_markers, paste0(path, 'All_samples_v5_markers.csv'))

# Check and cluster identification
All_obj <- LoadH5Seurat(paste0(path, 'All_samples_clustered_v5.h5seurat'),
                    assay = 'SCT', reductions = 'umap', graphs = FALSE, plots = FALSE)

# Checking clusters
p1 <- DimPlot(All_obj, reduction = 'umap', label = TRUE)
p2 <- VlnPlot(All_obj, features = 'percent.mt')

metadata <- All_obj@meta.data

p0 <- ggplot(metadata, aes(lane, fill = lane)) +
  geom_bar() +
  geom_text(stat = 'count', aes(label = after_stat(count)), vjust = -1) +
  theme_classic() +
  labs(title='N cells per lane', y = 'n',
       subtitle = paste0('Total number of cells: ', nrow(metadata)))

p3 <- DimPlot(All_obj, reduction = 'umap', group.by = 'batch', label = TRUE) + NoLegend()
p4 <- DimPlot(All_obj, reduction = 'umap', group.by = 'batch', split.by = 'batch') + NoLegend()

p5 <- ggplot(metadata, aes(genotype, fill = genotype)) +
  geom_bar() +
  geom_text(stat = 'count', aes(label = after_stat(count)), vjust = -1) +
  theme_classic() +
  labs(title='N cells per genotype', y = 'n',
       subtitle = paste0('Total number of cells: ', nrow(metadata)))

p6 <- ggplot(metadata, aes(status, fill = status)) +
  geom_bar() +
  geom_text(stat = 'count', aes(label = after_stat(count)), vjust = -1) +
  theme_classic() +
  coord_flip() +
  labs(title='N cells per status', y = 'n',
       subtitle = paste0('Total number of cells: ', nrow(metadata)))

p7 <- DimPlot(All_obj, reduction = 'umap', split.by = 'genotype') + NoLegend()
p8 <- DimPlot(All_obj, reduction = 'umap', split.by = 'status') + NoLegend()

pdf(paste0(path,'All_samples_checking_v5.pdf'), width = 11, height = 8, paper = 'a4r')
p0
p5
p6
p1
p2
p3
p4
p7
p8
dev.off()

# Identifying clusters
markers <- '../../resources/cell_type_markers_small_intestine_all.csv'
PlotUCell(All_obj, markers, output = paste0(path, 'All_samples_UCell_v5'))
plotMarkers(data = All_obj, markers = markers, output = paste0(path, 'All_samples_v5_markers'))

All_markers <- read.csv(paste0(path, 'All_samples_v5_markers.csv'))

# Getting the top 20 markers per cluster
top.markers <- All_markers %>% 
  group_by(cluster) %>%
  filter(p_val_adj < 0.05) %>%
  arrange(p_val_adj, by.group = TRUE) %>%
  top_n(n = 20, wt = avg_log2FC)

# Reshaping the table
top.markers.clean <- top.markers %>%
  select(cluster, gene) %>%
  group_by(cluster) %>%
  summarise(markers = str_c(unlist(cur_data()), collapse=', '))

write.csv(top.markers.clean, paste0(path, 'All_samples_v5_top20markers.csv'), row.names = F)

## Manual inspections
# Pericytes
FeaturePlot(All_obj, features = c('KCNJ8', 'ABCC9', 'RGS5', 'CSP4G'), label = T)
VlnPlot(All_obj, features = c('KCNJ8', 'ABCC9', 'RGS5', 'CSP4G'))
# Fibroblasts
FeaturePlot(All_obj, features = c('THY1', 'COL1A2', 'VIM'), label = T)
VlnPlot(All_obj, features = c('THY1', 'COL1A2', 'VIM'))
# Myofibroblast and mesothelium
FeaturePlot(All_obj, features = c('VIM', 'FOXF1', 'TAGLN', 'RSPO2'), label = T)
VlnPlot(All_obj, features = c('VIM', 'FOXF1', 'TAGLN', 'RSPO2'))
# Muscle
FeaturePlot(All_obj, features = c('MYH11', 'ACTG2', 'TAGLN'), label = T)
VlnPlot(All_obj, features = c('MYH11', 'ACTG2', 'TAGLN'))

## Checking the dataset with batch correction using Harmony --------------------
####
# Normalization and PCA
obj_v5 <- obj_v5 %>%
  SCTransform(vars.to.regress = "percent.mt", verbose = F) %>%
  RunPCA()

# Run Harmony
obj_v5 <- RunHarmony(obj_v5, "batch")

# Clustering
obj_v5 <- obj_v5 %>%
  RunUMAP(reduction = "harmony", dims = 1:50) %>%
  FindNeighbors(reduction = "harmony", dims = 1:50) %>%
  FindClusters(resolution = 0.4)
####


All_obj <- LoadH5Seurat(paste0(path, 'All_samples_clustered_v5_harmony.h5seurat'),
                        assay = 'SCT', reductions = 'umap', graphs = FALSE, plots = FALSE)

# Checking clusters
p1 <- DimPlot(All_obj, reduction = 'umap', label = TRUE)
p2 <- VlnPlot(All_obj, features = 'percent.mt')

metadata <- All_obj@meta.data

p0 <- ggplot(metadata, aes(lane, fill = lane)) +
  geom_bar() +
  geom_text(stat = 'count', aes(label = after_stat(count)), vjust = -1) +
  theme_classic() +
  labs(title='N cells per lane', y = 'n',
       subtitle = paste0('Total number of cells: ', nrow(metadata)))

p3 <- DimPlot(All_obj, reduction = 'umap', group.by = 'batch', label = TRUE) + NoLegend()
p4 <- DimPlot(All_obj, reduction = 'umap', group.by = 'batch', split.by = 'batch') + NoLegend()

p5 <- ggplot(metadata, aes(genotype, fill = genotype)) +
  geom_bar() +
  geom_text(stat = 'count', aes(label = after_stat(count)), vjust = -1) +
  theme_classic() +
  labs(title='N cells per genotype', y = 'n',
       subtitle = paste0('Total number of cells: ', nrow(metadata)))

p6 <- ggplot(metadata, aes(status, fill = status)) +
  geom_bar() +
  geom_text(stat = 'count', aes(label = after_stat(count)), vjust = -1) +
  theme_classic() +
  coord_flip() +
  labs(title='N cells per status', y = 'n',
       subtitle = paste0('Total number of cells: ', nrow(metadata)))

p7 <- DimPlot(All_obj, reduction = 'umap', split.by = 'genotype') + NoLegend()
p8 <- DimPlot(All_obj, reduction = 'umap', split.by = 'status') + NoLegend()

pdf(paste0(path,'All_samples_checking_v5_harmony.pdf'), width = 11, height = 8, paper = 'a4r')
p0
p5
p6
p1
p2
p3
p4
p7
p8
dev.off()

# Identifying clusters
markers <- '../../resources/cell_type_markers_small_intestine_all.csv'
PlotUCell(All_obj, markers, output = paste0(path, 'All_samples_UCell_v5_harmony'))
plotMarkers(data = All_obj, markers = markers, output = paste0(path, 'All_samples_v5_markers_harmony'))

All_markers <- read.csv(paste0(path, 'All_samples_v5_harmony_markers.csv'))

# Getting the top 20 markers per cluster
top.markers <- All_markers %>% 
  group_by(cluster) %>%
  filter(p_val_adj < 0.05) %>%
  arrange(p_val_adj, by.group = TRUE) %>%
  top_n(n = 20, wt = avg_log2FC)

# Reshaping the table
top.markers.clean <- top.markers %>%
  select(cluster, gene) %>%
  group_by(cluster) %>%
  summarise(markers = str_c(unlist(cur_data()), collapse=', '))

write.csv(top.markers.clean, paste0(path, 'All_samples_v5_top20markers_harmony.csv'), row.names = F)

# Further inspection in unknown clusters: 24, checking for innate lymphoid cells
FeaturePlot(All_obj, features = c('PTPRC', 'KIT', 'IL7R', 'ALDOC', 'PCDH9'))
VlnPlot(All_obj, features = c('nCount_RNA', 'nFeature_RNA'))

## Adding cluster identification to the object ---------------------------------
# Load the whole object first
All_obj <- LoadH5Seurat(paste0(path, 'All_samples_clustered_v5_harmony.h5seurat'))

clusters_ids <- read.csv(paste0(path, 'All_samples_v5_top20markers_harmony.csv'))
clusters_ids <- clusters_ids %>%
  select(cluster, cell.type.1, cell.type.2)
clusters_ids$cluster <- as.factor(clusters_ids$cluster)
new_clusters_ids <- clusters_ids %>%
  pull(cell.type.2)

names(new_clusters_ids) <- levels(All_obj)
All_obj <- RenameIdents(All_obj, new_clusters_ids)

# Also adding to the metadata
metadata <- All_obj@meta.data
metadata <- metadata %>%
  select(-cell.type.1, -cell.type.2)
metadata$barcode <- row.names(metadata)
metadata <- left_join(metadata, clusters_ids, by = c('seurat_clusters' = 'cluster'))
row.names(metadata) <- metadata$barcode
All_obj@meta.data <- metadata

# Saving both as RDS and H5
saveRDS(All_obj, file = paste0(path, 'All_obj_v5_final_indentified.rds'))
SaveH5Seurat(All_obj, file = paste0(path, 'All_obj_v5_final_indentified'))

# Also saving just the metadata
saveRDS(metadata, file = paste0(path, 'All_obj_v5_final_metadata.rds'))

## Removing Unknown cluster and saving  -------------------------------
new_dataset <- subset(my_dataset, subset = seurat_clusters != 8 & seurat_clusters != 25)
saveRDS(new_dataset, paste0(path, 'All_SeuratObj_v5_final_identified_unknownremoved.rds'))

metadata <- new_dataset@meta.data
saveRDS(metadata, file = paste0(path, 'All_SeuratObj_v5_final_metadata_unknownremoved.rds'))
