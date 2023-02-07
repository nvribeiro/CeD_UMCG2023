# Removing intra-sample doublets (doublets with same genotype, undetected by demuxlet)
# Using recoverDoublets to find the intra doublets based on the known inter doublets from demuxlet

library(Seurat)
library(tidyverse)
library(SeuratDisk)
library(scDblFinder)
library(scater)

# In-/Output path
path <- './data_preprocessing/outputs/all/'

# Defining function to find the doublets in each lane of each batch
findDoublets <- function(lane, obj) {
  start.time <- Sys.time()
  
  lane_obj <- subset(obj, subset = sample == lane)
  lane_obj <- as.SingleCellExperiment(lane_obj)
  doublets <- recoverDoublets(lane_obj, doublets = lane_obj$doublet, samples = table(lane_obj$genotype))
  
  doublets <- as.data.frame(doublets)
  doublets$barcode <- lane_obj$BARCODE.UPDATED
  doublets$lane <- lane
  
  return(doublets)
  
  end.time() <- Sys.time()
  time.taken <- end.time - start.time
  
  print(paste0('Lane ', lane, 'done. Time taken: ', time.taken))
}

# Batch F
F_obj <- readRDS(paste0(path, 'F_for_doublets.rds'))
lanes <- unique(F_obj$sample)

all_doublets <- map(lanes, ~findDoublets(obj = F_obj))
saveRDS(all_doublets, paste0(path, 'F_doublets.rds'))

# Batch A
A_obj <- readRDS(paste0(path, 'A_for_doublets.rds'))
lanes <- unique(A_obj$sample)

all_doublets <- map(lanes, ~findDoublets(obj = A_obj))
saveRDS(all_doublets, paste0(path, 'A_doublets.rds'))


# Checking results
# Getting barcodes for doublets to check them in the final dataset
# Loading dataset (v3)
obj <- LoadH5Seurat(paste0(path, 'All_samples_v3_identified.h5seurat'),
                    assay = 'SCT', reductions = 'umap', graphs = FALSE, plots = FALSE)

A_doublets <- readRDS(paste0(path, 'A_doublets.rds'))
A_doublets <- bind_rows(A_doublets)
F_doublets <- readRDS(paste0(path, 'F_doublets.rds'))
F_doublets <- bind_rows(F_doublets)

F_doublets <- F_doublets %>%
  filter(known == TRUE | predicted == TRUE) %>%
  mutate(barcode = paste0(barcode, '_2')) %>%
  pull(barcode)

A_doublets <- A_doublets %>%
  filter(known == TRUE | predicted == TRUE) %>%
  mutate(barcode = paste0(barcode, '_1')) %>%
  pull(barcode)

DimPlot(obj, reduction = 'umap', cells.highlight = c(A_doublets, F_doublets), label = T, repel = T)

barcodes <- as.vector(row.names(obj@meta.data))
table(c(A_doublets, F_doublets) %in% barcodes)

# Doublets per cluster
metadata <- obj@meta.data
metadata$barcode <- row.names(metadata)
cluster_count <- metadata %>%
  group_by(seurat_clusters) %>%
  summarise(total_cells = n())

doublets_per_cluster <- metadata %>%
  filter(barcode %in% c(A_doublets, F_doublets)) %>%
  group_by(seurat_clusters) %>%
  summarise(count_doublets = n()) %>%
  left_join(cluster_count, by = 'seurat_clusters') %>%
  mutate(pct = count_doublets / total_cells)

ggplot(doublets_per_cluster, aes(seurat_clusters, pct, label = scales::percent(pct, accuracy = 0.1))) +
  geom_col() +
  geom_text(nudge_y = 0.01) +
  labs(title = '% of intra-doublets per cluster') +
  theme_classic()

write.csv(doublets_per_cluster, paste0(path, './doublets_per_cluster_v3.csv'), row.names = F)

## Removing these doublets from the final dataset and reclustering -------------
barcodes_to_keep <- barcodes[barcodes %in% c(A_doublets, F_doublets) == FALSE]

obj <- LoadH5Seurat(paste0(path, 'All_samples_v3_identified.h5seurat'))
obj_v5 <- subset(obj, cells = barcodes_to_keep)

saveRDS(obj_v5, paste0(path, 'All_samples_v5.rds'))

# Running in the cluster
obj_v5 <- readRDS(paste0(path, 'All_samples_v5.rds'))

DefaultAssay(obj_v5) <- 'RNA'
obj_v5 <- All_samples %>% 
  SCTransform(vars.to.regress = "percent.mt", verbose = F) %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:50) %>%
  FindClusters(resolution = 0.4) %>%
  RunUMAP(dims = 1:50)

SaveH5Seurat(obj_v5, filename = paste0(path, 'All_samples_clustered_v5'), overwrite = T)

# Finding markers
All_markers <- FindAllMarkers(obj_v5, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(All_markers, paste0(path, 'All_samples_v5_markers.csv'))