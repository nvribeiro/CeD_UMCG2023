# Reclustering dataset v5 after removing unknown cluster
set.seed(1234)
library(Seurat)
library(tidyverse)
library(harmony)
library(patchwork)
library(UCell)

source('functions/my_functions.R')

out_path <- '/groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/NR03_scRNAseq/ongoing/outputs/data_preprocessing/'
#obj_path <- '/groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/NR03_scRNAseq/ongoing/All_samples_preprocessing/outputs/All_SeuratObj_v5_final_identified_unknownremoved.rds'

#obj <- readRDS(obj_path)

DefaultAssay(obj) <- 'RNA'

# Normalization and PCA
obj <- obj %>%
  SCTransform(vars.to.regress = "percent.mt", verbose = T) %>%
  RunPCA()

# Run Harmony
obj <- RunHarmony(obj, "batch")

# Clustering
obj <- obj %>%
  RunUMAP(reduction = "harmony", dims = 1:50) %>%
  FindNeighbors(reduction = "harmony", dims = 1:50) %>%
  FindClusters(resolution = 0.4)

# Saving
saveRDS(obj, paste0(out_path, 'SeuratObj_AllCells_SCT_RES04_clustered_final.rds'))
obj <- readRDS(paste0(out_path, 'SeuratObj_AllCells_SCT_RES04_clustered_final.rds'))
  
# Find markers
All_markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(All_markers, paste0(out_path, 'AllCells_SCT_RES04_clustered_final_markers.csv'))

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

write.csv(top.markers.clean, paste0(out_path, 'AllCells_SCT_RES04_clustered_final_top20markers.csv'), row.names = F)

# Plotting markers and identifying clusters
markers <- '../resources/cell_type_markers_small_intestine_all.csv'
PlotUCell(obj, markers, output = paste0(out_path, 'AllCells_SCT_RES04_clustered_final_Ucell'))
plotMarkers(obj, markers, output = paste0(out_path, 'AllCells_SCT_RES04_clustered_final_markers'))

# Adding clusters id to the metadata
clusters_ids <- read.csv(paste0(out_path, 'AllCells_SCT_RES04_clustered_final_top20markers_identity.csv'))
clusters_ids <- clusters_ids %>%
  select(cluster, cell.type.1, cell.type.2)
clusters_ids$cluster <- as.factor(clusters_ids$cluster)

metadata <- obj@meta.data
metadata <- metadata %>%
  select(-cell.type.1, -cell.type.2)
metadata <- left_join(metadata, clusters_ids, by = c('seurat_clusters' = 'cluster'))
row.names(metadata) <- metadata$barcode
obj@meta.data <- metadata

# Saving the annotated obj
saveRDS(obj, paste0(out_path, 'SeuratObj_AllCells_SCT_RES04_clustered_final_annotated.rds'))

# Saving just the metadata
saveRDS(metadata, paste0(out_path, 'AllCells_SCT_RES04_clustered_final_onlymetadata.rds'))

# UMAP plots
pdf(paste0(out_path, 'AllCells_SCT_RES04_clustered_final_UMAPs.pdf'), width = 11, height = 8, paper = 'a4r')
DimPlot(obj, reduction = 'umap', label = T) + NoLegend()
DimPlot(obj, reduction = 'umap', split.by = 'batch') + NoLegend()
DimPlot(obj, reduction = 'umap', split.by = 'status') + NoLegend()
(VlnPlot(obj, features = 'percent.mt') + NoLegend()) / (VlnPlot(obj, features = 'nFeature_RNA') + NoLegend()) / (VlnPlot(obj, features = 'nCount_RNA') + NoLegend())
DimPlot(obj, reduction = 'umap', group.by = 'cell.type.1')
DimPlot(obj, reduction = 'umap', group.by = 'cell.type.2', label = T, repel = T) + NoLegend()
dev.off()

## Solving the problem with A7 barcodes (see build_full_metadata.R) -------------
obj <- readRDS(paste0(out_path, 'AllCells_SCT_RES04_clustered_final_annotated.rds'))
tmp_A7 <- filter(obj@meta.data, lane == 'A7')

# Finding the A7 intra-sample doublets to remove
obj@meta.data <- obj@meta.data %>%
  mutate(A7_remove = if_else(barcode %in% A7_intra_doublets == TRUE, TRUE, FALSE))

# Now removing the A7 intra-doublets
obj <- subset(obj, subset = A7_remove == FALSE)

# Remove this column
obj@meta.data <- select(obj@meta.data, -A7_remove)

## Re-run the normalization and clustering from beginning ----------------------
