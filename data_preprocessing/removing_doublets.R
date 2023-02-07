## Removing doublets (same genotype) from filtered clustered objects

library(Seurat)
library(tidyverse)
library(patchwork)
library(SeuratDisk)
library(scDblFinder)
library(scater)

path <- './data_preprocessing/outputs/all/'

## Investigating doublets -------------------------------------------------------

# Using the doblets from A7 as a 'benchmark'
# Loading the dataset
A7_obj <- Read10X('./data_preprocessing/data/A7/') %>%
  CreateSeuratObject(project = 'NR03')

A7_obj$percent.mt <- PercentageFeatureSet(A7_obj, pattern = '^MT-')

# Starting demultiplex
A7_demuxlet.best <- read.delim('./data_preprocessing/demuxlet/A7_demuxlet.best')

A7_demuxlet.best <- A7_demuxlet.best %>% select(BARCODE, BEST, SNG.1ST)

tmp <- A7_obj@meta.data
tmp$BARCODE <- row.names(tmp)
tmp <- left_join(tmp, A7_demuxlet.best, by="BARCODE")

tmp <- tmp %>% mutate(
  droplet = case_when(
    grepl("SNG", BEST) == TRUE ~ 'singlet',
    grepl("DBL", BEST) == TRUE ~ 'doblet',
    grepl("AMB", BEST) == TRUE ~ 'ambiguous')) %>% 
  dplyr::rename(genotype = SNG.1ST)

samples_info <- read.csv('../samples_info/samples_info.csv')
tmp <- left_join(tmp, samples_info, by = c('genotype' = 'sample_ID'))

row.names(tmp) <- tmp$BARCODE

# Replacing 'wrong' genotypes (doblets and ambiguous) with NA
tmp <- tmp %>% 
  mutate(genotype = replace(genotype, droplet != 'singlet', NA))

# Excluding BEST column
tmp <- tmp %>% select(!BEST)
A7_data <- tmp

# Checking the distribution of RNA count and features in doblets and singlets
p1 <- ggplot(filter(A7_data, droplet != 'ambiguous'), aes(nCount_RNA, color = droplet, fill = droplet)) +
  geom_density(alpha = 0.3) +
  xlim(c(0, 10000))

p2 <- ggplot(filter(A7_data, droplet != 'ambiguous'), aes(nFeature_RNA, color = droplet, fill = droplet)) +
  geom_density(alpha = 0.3) +
  xlim(c(0, 5000))

p1 | p2

## Now checking the pattern of a 'good' lane -----------------------------------
F1_obj <- Read10X('./data_preprocessing/data/A2/') %>%
  CreateSeuratObject(project = 'NR03')

F1_obj$percent.mt <- PercentageFeatureSet(F1_obj, pattern = '^MT-')

# Starting demultiplex
F1_demuxlet.best <- read.delim('./data_preprocessing/demuxlet/A2_demuxlet.best')

F1_demuxlet.best <- F1_demuxlet.best %>% select(BARCODE, BEST, SNG.1ST)

tmp <- F1_obj@meta.data
tmp$BARCODE <- row.names(tmp)
tmp <- left_join(tmp, F1_demuxlet.best, by="BARCODE")

tmp <- tmp %>% mutate(
  droplet = case_when(
    grepl("SNG", BEST) == TRUE ~ 'singlet',
    grepl("DBL", BEST) == TRUE ~ 'doblet',
    grepl("AMB", BEST) == TRUE ~ 'ambiguous')) %>% 
  dplyr::rename(genotype = SNG.1ST)

samples_info <- read.csv('../samples_info/samples_info.csv')
tmp <- left_join(tmp, samples_info, by = c('genotype' = 'sample_ID'))

row.names(tmp) <- tmp$BARCODE

# Replacing 'wrong' genotypes (doblets and ambiguous) with NA
tmp <- tmp %>% 
  mutate(genotype = replace(genotype, droplet != 'singlet', NA))

# Excluding BEST column
tmp <- tmp %>% select(!BEST)
F1_data <- tmp

# Checking the distribution of RNA count and features in doblets and singlets
p3 <- ggplot(filter(F1_data, droplet != 'ambiguous'), aes(nCount_RNA, color = droplet, fill = droplet)) +
  geom_density(alpha = 0.3) +
  xlim(c(0, 10000))

p4 <- ggplot(filter(F1_data, droplet != 'ambiguous'), aes(nFeature_RNA, color = droplet, fill = droplet)) +
  geom_density(alpha = 0.3) +
  xlim(c(0, 5000))

(p1 | p2) / (p3 | p4)

## Trying scDblFinder ----------------------------------------------------------
# Trying just with A7 singles
A7_obj@meta.data <- A7_data
A7_sing <- subset(A7_obj, subset = droplet == 'singlet' & nCount_RNA > 0)
A7_sing <- as.SingleCellExperiment(A7_sing)
A7_sing <- scDblFinder(A7_sing)
table(A7_sing$scDblFinder.class)

# Some visual inspection
A7_sing <-  as.Seurat(A7_sing)
tmp <- A7_sing@meta.data

p1 <- ggplot(tmp, aes(nCount_RNA, color = scDblFinder.class, fill = scDblFinder.class)) +
  geom_density(alpha = 0.3) +
  xlim(c(0, 10000))

p2 <- ggplot(tmp, aes(nFeature_RNA, color = scDblFinder.class, fill = scDblFinder.class)) +
  geom_density(alpha = 0.3) +
  xlim(c(0, 5000))

p1 / p2


# Trying using the known doublets
A7_data <- A7_data %>%
  mutate(doublet = if_else(droplet == 'doblet', TRUE, FALSE))
A7_obj@meta.data <- A7_data

# Filtering 'empty' cells
A7_obj <- subset(A7_obj, subset = nCount_RNA > 1000 & nFeature_RNA > 300)
A7_obj <- as.SingleCellExperiment(A7_obj)

A7_obj <- scDblFinder(A7_obj, knownDoublets = A7_obj$doublet, knownUse = 'positive')
table(A7_obj$scDblFinder.class)
A7_obj <- as.Seurat(A7_obj)

# Checking the cells that are left
A7_singlets <- subset(A7_obj, subset = (scDblFinder.class == 'singlet' & droplet == 'singlet'))
VlnPlot(A7_singlets, features = c('percent.mt', 'nFeature_RNA', 'nCount_RNA'))

A7_singlets <- A7_singlets %>% 
  SCTransform(vars.to.regress = "percent.mt", verbose = T) %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:50) %>%
  FindClusters(resolution = 0.4) %>%
  RunUMAP(dims = 1:50)

DimPlot(A7_singlets, reduction = 'umap')

## Finding doublets in batch A -------------------------------------------------
# Using as input my final object already filtered, clustered and identified
A_obj <- LoadH5Seurat('./data_preprocessing/outputs/A/SeuraObj_A_identified_filtered_reclustered.h5seurat')
A_obj <- as.SingleCellExperiment(A_obj)

# Testing with no clusters
A_obj <- scDblFinder(A_obj)
table(A_obj$scDblFinder.class)

# With clustering
A_obj <- scDblFinder(A_obj, clusters = 'seurat_clusters', samples = 'sample')
table(A_obj$scDblFinder.class)

plotUMAP(A_obj, colour_by = 'scDblFinder.score')


## Finding doublets in batch F -------------------------------------------------
F_obj <- LoadH5Seurat('./data_preprocessing/outputs/F/SeuratObj_F_identified_filtered_reclustered.h5Seurat')
F_obj <- as.SingleCellExperiment(F_obj)

# Testing with no clusters
F_obj <- scDblFinder(F_obj)
table(F_obj$scDblFinder.class)

# With clustering
DefaultAssay(F_obj) <- 'RNA'

F_obj <- F_obj %>% 
  SCTransform(vars.to.regress = "percent.mt", verbose = F) %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:50) %>%
  FindClusters(resolution = 0.4) %>%
  RunUMAP(dims = 1:50)

DimPlot(F_obj, reduction = 'umap')

F_obj <- as.SingleCellExperiment(F_obj)
F_obj <- scDblFinder(F_obj, clusters = 'seurat_clusters', samples = 'sample')
table(F_obj$scDblFinder.class)

plotUMAP(F_obj, colour_by = 'scDblFinder.score')
hist(F_obj$scDblFinder.score)

F_obj <- as.Seurat(F_obj)
SaveH5Seurat(F_obj, './data_preprocessing/outputs/F/SeuratObj_F_identified_filtered_reclustered')

############################

## Trying with the whole dataset

obj <- LoadH5Seurat('./data_preprocessing/outputs/all/All_samples_v3_identified.h5seurat')

obj <- as.SingleCellExperiment(obj)
set.seed(1234)
obj <- scDblFinder(obj, clusters = 'seurat_clusters', samples = 'lane')
table(obj$scDblFinder.class)
obj <- as.Seurat(obj)

obj_doublets <- WhichCells(obj, expression = scDblFinder.class == 'doublet')
DimPlot(obj, reduction = 'UMAP', cells.highlight = obj_doublets)

doublets_info <- obj@meta.data %>%
  filter(scDblFinder.class == 'doublet')

## Trying findDoubletClusters() -------------------------------------------------
obj <- as.SingleCellExperiment(obj)
set.seed(1234)

dbl.out <- findDoubletClusters(obj, clusters = obj$seurat_clusters)

## Trying recoverDoublets() - run in the cluster --------------------------------
# Testing with F1 after initial filter

F1 <- subset(F_filtered, subset = sample == 'F1')
saveRDS(F1, './data_preprocessing/outputs/F/F1_for_doublet.rds')

F1_known_doublets <- F_filtered@meta.data %>%
  filter(sample == 'F1' & doublet == TRUE)

F1 <- as.SingleCellExperiment(F1)
F1_doublets <- recoverDoublets(F1, doublets = F1$doublet, samples = c(0.2243512, 0.2586657, 0.2093712, 0.1781920, 0.1294200))
F1 <- as.Seurat(F1)

# Comparing with scDblFinder results
set.seed(1234)
F1_scDbl <- scDblFinder(F1, knownDoublets = F1$doublet, knownUse = 'discard')

F1_scDbl <- as.Seurat(F1_scDbl)

# Get the barcodes for the doublets found
F1_doublets <- as.data.frame(readRDS('./data_preprocessing/outputs/F/F1_revoverDoublets.rds'))

tmp <- F1@meta.data
barcodes <- tmp$BARCODE.UPDATED
F1_doublets$barcode <- barcodes
tmp <- left_join(tmp, F1_doublets, by = c('BARCODE.UPDATED' = 'barcode'))
row.names(tmp) <- barcodes
F1@meta.data <- tmp

# Get the intra-sample doublets
intra_barcodes <- F1_doublets %>%
  filter(predicted == TRUE) %>%
  pull(BARCODE.UPDATED)

intra_barcodes <- paste0(intra_barcodes, '_2')

# Finding them in the final dataset
obj <- LoadH5Seurat('./data_preprocessing/outputs/all/All_samples_v3_identified.h5seurat',
                    assay = 'SCT', reductions = 'umap', graphs = FALSE, plots = FALSE)

DimPlot(obj, reduction = 'umap', cells.highlight = intra_barcodes)

barcodes <- as.vector(row.names(obj@meta.data))
table(intra_barcodes %in% barcodes)

# Clustering F1 to visualize if the results make sense
F1 <- F1 %>% 
  SCTransform(vars.to.regress = "percent.mt", verbose = F) %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:50) %>%
  FindClusters(resolution = 0.4) %>%
  RunUMAP(dims = 1:50)

F1@meta.data <- F1@meta.data %>% 
  mutate(genotype = replace_na(genotype, 'unknown'))

p1 <- FeaturePlot(F1, features = 'proportion', pt.size = 1, cols = c('darkblue', 'yellow'))
p2 <- DimPlot(F1, reduction = 'umap', group.by = 'known', pt.size = 1, cols = c(rgb(0,0,0,0.05), 'red'))
p3 <- DimPlot(F1, reduction = 'umap', group.by = 'predicted', pt.size = 1, cols = c(rgb(0,0,0,0.05), 'red'))

p1 / (p2 | p3)

####
# Testing with A7
A7 <- readRDS('./data_preprocessing/outputs/A/A7_for_doublets.rds')
A7 <- as.SingleCellExperiment(A7)
A7_doublets <- recoverDoublets(A7, doublets = A7$doublet, samples = table(A7$genotype))

saveRDS(A7_doublets, paste0(path, 'A7_recoverDoublets.rds'))


## Using recoverDoublets() in all lanes ----------------------------------------
# Run in the cluster
library(Seurat)
library(tidyverse)
library(scDblFinder)
library(scater)

# batch F
F_obj <- readRDS('./outputs/F_for_doublets.rds')

lanes <- unique(F_obj$sample)
all_doublets <- data.frame()

for (i in lanes) {
  lane_obj <- subset(F_obj, subset = sample == i)
  lane_obj <- as.SingleCellExperiment(lane_obj)
  doublets <- recoverDoublets(lane_obj, doublets = lane_obj$doublet, samples = table(lane_obj$genotype))
  
  doublets <- as.data.frame(doublets)
  
  doublets$barcode <- lane_obj$BARCODE.UPDATED
  doublets$lane <- i
  
  all_doublets <- bind_rows(all_doublets, doublets)
}

saveRDS(all_doublets, './outputs/F_doublets.rds')

# Same for batch A - but with a function
library(Seurat)
library(tidyverse)
library(scDblFinder)
library(scater)

# Defining my function
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

A_obj <- readRDS('./outputs/A_for_doublets.rds')
lanes <- unique(A_obj$sample)

all_doublets <- map(lanes, ~findDoublets(obj = A_obj))

saveRDS(all_doublets, './outputs/A_doublets.rds')

# Checking results
A_doublets <- readRDS('./data_preprocessing/outputs/A/A_doublets.rds')
A_doublets <- bind_rows(A_doublets)

##########
## Getting barcodes for doublets to check them in the final dataset
# Loading dataset (v3)
obj <- LoadH5Seurat('./data_preprocessing/outputs/all/All_samples_v3_identified.h5seurat',
                    assay = 'SCT', reductions = 'umap', graphs = FALSE, plots = FALSE)

F_doublets <- readRDS('./data_preprocessing/outputs/F/F_doublets.rds')
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

write.csv(doublets_per_cluster, './data_preprocessing/outputs/all/doublets_per_cluster_v3.csv', row.names = F)

## Removing these doublets from the final dataset and reclustering ----
barcodes_to_keep <- barcodes[barcodes %in% c(A_doublets, F_doublets) == FALSE]

obj <- LoadH5Seurat('./data_preprocessing/outputs/all/All_samples_v3_identified.h5seurat')
obj_v5 <- subset(obj, cells = barcodes_to_keep)

saveRDS(obj_v5, './data_preprocessing/outputs/all/All_samples_v5.rds')

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