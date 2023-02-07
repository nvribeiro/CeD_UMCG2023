# Merging and integrating samples A and F after preprocessing
# Running this in the cluster

# Load libraries
library(Seurat)
library(tidyverse)
library(SeuratDisk)
library(patchwork)
library(harmony)

source('./functions/my_functions.R')

path <- './data_preprocessing/outputs/all/'

A_samples <- LoadH5Seurat(paste0(path,'SeuratObj_A_identified_filtered.h5seurat'))
F_samples <- LoadH5Seurat(paste0(path,'SeuratObj_F_identified_filtered.h5seurat'))

DefaultAssay(F_samples) <- 'RNA'
DefaultAssay(A_samples) <- 'RNA'

# Merge datasets
All_samples <- merge(A_samples, F_samples)

# Clustering
All_samples <- All_samples %>% 
  SCTransform(vars.to.regress = "percent.mt", verbose = F) %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:50) %>%
  FindClusters(resolution = 0.4) %>%
  RunUMAP(dims = 1:50)

SaveH5Seurat(All_samples, filename = paste0(path, 'All_samples_merged_v2'))

# QC inspection
All_obj <- LoadH5Seurat(paste0(path, 'All_samples_merged_v3.h5seurat'),
                          assay = 'SCT', reductions = 'umap', graphs = FALSE, plots = FALSE)

p1 <- DimPlot(All_obj, reduction = 'umap', label = TRUE)
p2 <- VlnPlot(All_obj, features = 'percent.mt')

metadata <- All_obj@meta.data

p0 <- ggplot(metadata, aes(lane, fill = lane)) +
  geom_bar() +
  geom_text(stat = 'count', aes(label = after_stat(count)), vjust = -1) +
  theme_classic() +
  labs(title='N cells per lane', y = 'n',
       subtitle = paste0('Total number of cells: 27,427'))

# Renaming some columns for clarity and tidying the metadata a bit
metadata <- metadata %>%
  dplyr::rename(lane = sample) %>%
  dplyr::rename(cell_ID_A = cell_ID_v1) %>%
  dplyr::select(-BARCODE, -BARCODE.UPDATED)

# Updating batch column
metadata <- metadata %>% 
  mutate(batch = case_when(
    lane == 'A2' | lane == 'A3' | lane == 'A4' | lane == 'A5' | lane == 'A6' | lane == 'A7' ~ 'A',
    lane == 'F1' | lane == 'F2' | lane == 'F3' ~ 'F'))

# Saveing back to the object
All_obj@meta.data <- metadata

# Checking for batch effects
p3 <- DimPlot(All_obj, reduction = 'umap', split.by = 'batch', label = TRUE) + NoLegend()
p4 <- DimPlot(All_obj, reduction = 'umap', group.by = 'batch', split.by = 'batch') + NoLegend()

p5 <- ggplot(metadata, aes(genotype, fill = genotype)) +
  geom_bar() +
  geom_text(stat = 'count', aes(label = after_stat(count)), vjust = -1) +
  theme_classic() +
  labs(title='N cells per genotype', y = 'n',
       subtitle = paste0('Total number of cells: 27,427'))

p6 <- ggplot(metadata, aes(status, fill = status)) +
  geom_bar() +
  geom_text(stat = 'count', aes(label = after_stat(count)), vjust = -1) +
  theme_classic() +
  coord_flip() +
  labs(title='N cells per status', y = 'n',
       subtitle = paste0('Total number of cells: 27,427'))

p7 <- DimPlot(All_obj, reduction = 'umap', split.by = 'genotype') + NoLegend()

pdf(paste0(path,'All_samples_checking_v3.pdf'), width = 11, height = 8, paper = 'a4r')
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

DimPlot(All_obj, reduction = 'umap', split.by = 'lane') + NoLegend()
p8 <- DimPlot(All_obj, reduction = 'umap', split.by = 'status') + NoLegend()

# Save the object with updated metadata
SaveH5Seurat(All_obj, paste0(path, 'All_samples_v3_with_metadata'))

## Finding markers - do also for samples from donor A0038 alone
All_markers <- FindAllMarkers(All_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(All_markers, paste0(path,'A_all_markers_clusters.csv'), row.names = F)

# Subseting A0338
A_19 <- subset(All_obj, subset = genotype == '19_CeDNN-A0338')
DimPlot(A_19, reduction = 'umap')

A_19_markers <- FindAllMarkers(A_19, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Getting the top 10 markers per cluster
top.markers <- A_19_markers %>% 
  group_by(cluster) %>% 
  arrange(p_val_adj, by.group = TRUE) %>%
  top_n(n = 10, wt = avg_log2FC)

# Reshaping the table
top.markers.clean <- top.markers %>%
  select(cluster, gene) %>%
  group_by(cluster) %>%
  summarise(markers = str_c(unlist(cur_data()), collapse=', '))

write.csv(top.markers.clean, paste0(path, 'A0338_clusters_top10markers.csv'), row.names = F)

# Checking the clusters
epithelial_markers <- '../../resources/cell_type_markers_small_intestine_epithelial.csv'
other_markers <- '../../resources/cell_type_markers_small_intestine_immune_and_others.csv'

PlotUCell(A_19, epithelial_markers, output = paste0(path, 'A0338_epithelial_Ucell'))
PlotUCell(A_19, other_markers, output = paste0(path, 'A0338_others_Ucell'))

# What happens if I recluster just them?
All_obj <- LoadH5Seurat('./data_preprocessing/outputs/all/All_samples_merged_v2.h5seurat')
A_19 <- subset(All_obj, subset = genotype == '19_CeDNN-A0338')
DefaultAssay(A_19) <- 'RNA'
A_19

A_19 <- A_19 %>% 
  SCTransform(vars.to.regress = "percent.mt", verbose = F) %>%
  RunPCA()

# ElbowPlot(A_19, ndims = 50)

A_19 <- A_19 %>% 
  FindNeighbors(dims = 1:50) %>%
  FindClusters(resolution = 0.4) %>%
  RunUMAP(dims = 1:50)

SaveH5Seurat(A_19, paste0(path, 'A0338_SeuratObj'))

DimPlot(A_19, reduction = 'umap', label = TRUE) + NoLegend()

FeaturePlot(A_19, features = c('MUC6', 'MUC3A', 'MUC2'))
VlnPlot(A_19, features = 'percent.mt') + NoLegend()

A_19_markers <- FindAllMarkers(A_19, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Getting the top 10 markers per cluster
top.markers <- A_19_markers %>% 
  group_by(cluster) %>% 
  arrange(p_val_adj, by.group = TRUE) %>%
  top_n(n = 10, wt = avg_log2FC)

# Reshaping the table
top.markers.clean <- top.markers %>%
  select(cluster, gene) %>%
  group_by(cluster) %>%
  summarise(markers = str_c(unlist(cur_data()), collapse=', '))

write.csv(top.markers.clean, paste0(path, 'A0338_clusters_top10markers_v2.csv'), row.names = F)
PlotUCell(A_19, epithelial_markers, output = paste0(path, 'A0338_epithelial_Ucell_v2'))
PlotUCell(A_19, other_markers, output = paste0(path, 'A0338_others_Ucell_v2'))

## Exclude sample 19_CeDNN-A0338 and recluster -----------------------------------------

# Load dataset
All_samples <- LoadH5Seurat(path, 'All_samples_merged_v2')

# Excluding genotype '19_CeDNN-A0338'
All_samples <- subset(All_samples, subset = genotype != '19_CeDNN-A0338')

# Clustering
DefaultAssay(All_samples) <- 'RNA'

# Clustering
All_samples <- All_samples %>% 
  SCTransform(vars.to.regress = "percent.mt", verbose = F) %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:50) %>%
  FindClusters(resolution = 0.4) %>%
  RunUMAP(dims = 1:50)

SaveH5Seurat(All_samples, filename = paste0(path, 'All_samples_merged_v3'), overwrite = T)

# Run QC and inspection lines again

## Starting clusters identification --------------------------------------------------
All_markers <- FindAllMarkers(All_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(All_markers, paste0(path,'A_all_markers_clusters_v3.csv'), row.names = F)

markers <- '../../resources/cell_type_markers_small_intestine_all.csv'

PlotUCell(All_obj, markers, output = paste0(path, 'All_samples_UCell_v3'))
plotMarkers(data = All_obj, markers = markers, output = paste0(path, 'All_samples_v3'))

All_markers <- read.csv('./data_preprocessing/outputs/all/All_markers_clusters_v3.csv')

# Getting the top 10 markers per cluster
top.markers <- All_markers %>% 
  group_by(cluster) %>% 
  arrange(p_val_adj, by.group = TRUE) %>%
  top_n(n = 10, wt = avg_log2FC)

# Reshaping the table
top.markers.clean <- top.markers %>%
  select(cluster, gene) %>%
  group_by(cluster) %>%
  summarise(markers = str_c(unlist(cur_data()), collapse=', '))

write.csv(top.markers.clean, paste0(path, 'All_clusters_top10markers.csv'), row.names = F)

FeaturePlot(All_obj, features = c('KCNJ8', 'ABCC9', 'RGS5', 'CSP4G'), label = T)
VlnPlot(All_obj, features = c('KCNJ8', 'ABCC9', 'RGS5', 'CSP4G'))

FeaturePlot(All_obj, features = c('LYVE1', 'PROX1', 'RELN'))

VlnPlot(All_obj, features = c('ELAVL4', 'CHRNA3', 'GAP43'))

#IEL
FeaturePlot(All_obj, features = c('CCR9', 'ITGAE', 'KLRK1'))

## Saving clusters IDs to metadata ----------------------------------------------
clusters.IDs <- read.csv('./data_preprocessing/outputs/all/All_clusters_top10markers.csv')
clusters.IDs <- clusters.IDs %>% select(cluster, cell.type.1, cell.type.2)
clusters.IDs$cluster <- as.factor(clusters.IDs$cluster)

metadata$barcode <- row.names(metadata)
metadata <- metadata %>% select(-cell_ID_A, -cell_ID_F)
metadata <- left_join(metadata, clusters.IDs, by = c('seurat_clusters' = 'cluster'))
row.names(metadata) <- metadata$barcode
metadata <- metadata %>% select(-barcode)

# Load the whole obj to save the metadata
All_obj <- LoadH5Seurat(paste0(path, 'All_samples_merged_v3.h5Seurat'))
All_obj@meta.data <- metadata
SaveH5Seurat(All_obj, paste0(path, 'All_samples_v3_identified'))

# Saving just the metadata
write.csv(metadata, paste0(path, 'All_samples_v3_metadata.csv'))

## Checking for where A7 cells are ----------------------------------------------
obj <- LoadH5Seurat('./data_preprocessing/outputs/all/All_samples_v3_identified.h5seurat',
                    assay = 'SCT', reductions = 'umap', graphs = FALSE, plots = FALSE)


A7_cells <- WhichCells(obj, expression = lane == 'A7')
DimPlot(obj, split.by = 'batch', label = T) + NoLegend()
FeaturePlot(obj, features = 'percent.mt')
# Not a problem

## Checking the new clustering after removing extra doublets -------------------
