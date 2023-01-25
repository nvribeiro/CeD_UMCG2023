# Preprocessing of A samples without A7 after demultiplexing
# Color-blind friendly palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")

# Load libraries
library(Seurat)
library(tidyverse)
library(SeuratDisk)
library(patchwork)

## Preprocessing A withouht A7 after demultiplexing ---------------------------------------------------------
# Load the dataset
A_obj <- LoadH5Seurat('./data_preprocessing/outputs/A/SeuraObj_A_filtered_clustered_woA7_genotyped.h5seurat')

# Excluding ambiguous and doblets
A_obj <- subset(A_obj, subset = droplet == 'singlet')

# Change the active assay back ro RNA
DefaultAssay(A_obj) <- 'RNA'

# Normalizing and clustering
A_obj <- SCTransform(A_obj, vars.to.regress = "percent.mt", verbose = T) %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:50) %>%
  FindClusters(dims = 1:50, resolution = 0.4) %>%
  RunUMAP(dims = 1:50)

# Save
SaveH5Seurat(A_obj, './data_preprocessing/outputs/A/SeuraObj_A_filtered_reclustered_woA7_genotyped.h5seurat')

# Inspecting clusters
A_obj <- LoadH5Seurat('./data_preprocessing/outputs/A/SeuraObj_A_filtered_reclustered_woA7_genotyped.h5seurat',
                      assays = 'SCT', reductions = 'umap', graphs = FALSE, images = FALSE)
View(A_obj@meta.data)

pdf('./data_preprocessing/outputs/plots/QC_plots_A_merged_filtered_clustered_genotyped.pdf', width = 11, height = 8, paper = 'a4r')
  DimPlot(A_obj, reduction = 'umap', label = T)
  DimPlot(A_obj, reduction = 'umap', label = T, split.by = 'genotype')
  DimPlot(A_obj, reduction = 'umap', label = T, split.by = 'sample')
  DimPlot(A_obj, reduction = 'umap', label = T, split.by = 'status')
  VlnPlot(A_obj, features = 'percent.mt')
  VlnPlot(A_obj, features = 'nFeature_RNA')
  VlnPlot(A_obj, features = 'nCount_RNA')
dev.off()

## Preprocessing A7 after demultiplexing ----------------------------------------
# Load the dataset
A7_obj <- LoadH5Seurat('./data_preprocessing/outputs/A/SeuraObj_A7_genotyped.h5seurat')

A7_obj <- subset(A7_obj, subset = nCount_RNA > 1000 &
                   nCount_RNA < 30000 &
                   nFeature_RNA > 300 &
                   nFeature_RNA < 6000 &
                   percent.mt < 90 &
                   droplet == 'singlet')

A7_obj <- A7_obj %>%
  SCTransform(vars.to.regress = "percent.mt", verbose = T) %>%
  RunPCA()

A7_obj <- A7_obj %>%
  FindNeighbors(dims = 1:50) %>%
  FindClusters(dims = 1:50, resolution = 0.4) %>%
  RunUMAP(dims = 1:50)

pdf('./data_preprocessing/outputs/plots/QC_plots_A7_filtered_clustered_genotyped.pdf', width = 11, height = 8, paper = 'a4r')
  DimPlot(A7_obj, reduction = 'umap', label = T)
  DimPlot(A7_obj, reduction = 'umap', label = T, split.by = 'genotype')
  DimPlot(A7_obj, reduction = 'umap', label = T, split.by = 'status')
  VlnPlot(A7_obj, features = 'percent.mt')
  VlnPlot(A7_obj, features = 'nFeature_RNA')
  VlnPlot(A7_obj, features = 'nCount_RNA')
dev.off()

SaveH5Seurat(A7_obj, './data_preprocessing/outputs/A/SeuraObj_A7_genotyped_filtered_clustered.h5seurat')

## Merging A and A7 after demultiplex and QC ------------------------------------
# Path outputs 
path <- './outputs/' # in the cluster
#path <- './data_preprocessing/outputs/A/' # locally

# Load datasets
A_obj <- LoadH5Seurat(paste0(path, 'SeuraObj_A_filtered_reclustered_woA7_genotyped.h5seurat'))
A7_obj <- LoadH5Seurat(paste0(path, 'SeuraObj_A7_genotyped_filtered_clustered.h5seurat'))

# Changing the active assay
DefaultAssay(A_obj) <- 'RNA'
DefaultAssay(A7_obj) <- 'RNA'

# Adding lane information to A7 before merging
A7_obj$sample <- 'A7'

# Merging
A_obj <- merge(A_obj, A7_obj)

# Clustering
A_obj <- A_obj %>%
  SCTransform(vars.to.regress = "percent.mt", verbose = T) %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:50) %>%
  FindClusters(dims = 1:50, resolution = 0.4) %>%
  RunUMAP(dims = 1:50)

# Saving
SaveH5Seurat(A_obj, paste0(path, 'SeuratObj_A_merged_genotyped_filtered_clustered.h5seurat'))

# Checking the results
A_obj <- LoadH5Seurat('./data_preprocessing/outputs/A/SeuratObj_A_merged_genotyped_filtered_clustered.h5seurat',
                      assays = 'SCT', reductions = 'umap', graphs = FALSE, images = FALSE)
A_obj

View(A_obj@meta.data)

# Removing extra BARCODE (A7) column from meta.data
A_obj@meta.data <- select(A_obj@meta.data, !BARCODE)
tmp <- A_obj@meta.data

# QC plots 
p <- tmp %>%
  ggplot(aes(x=sample, fill=sample)) +
  geom_bar() +
  geom_text(stat='count', aes(label=after_stat(count)), vjust=-1) +
  theme_classic() +
  scale_fill_manual(values=cbbPalette) +
  labs(title='N cells per lane', y = 'n')

p1 <- tmp %>% 
  ggplot(aes(x=nCount_RNA, color=sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  scale_fill_manual(values=cbbPalette) +
  scale_colour_manual(values=cbbPalette) +
  labs(title='RNA counts per cell', y = 'Cell density')

p2 <- tmp %>% 
  ggplot(aes(x=nFeature_RNA, color=sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  scale_fill_manual(values=cbbPalette) +
  scale_colour_manual(values=cbbPalette) +
  labs(title='Features per cell', y = 'Cell density')

p3 <- tmp %>% 
  ggplot(aes(x=percent.mt, color=sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  scale_fill_manual(values=cbbPalette) +
  scale_colour_manual(values=cbbPalette) +
  labs(title='% MT RNA per cell', y = 'Cell density')

# Checking clusters
p4 <- VlnPlot(A_obj, features = 'percent.mt')
p5 <- VlnPlot(A_obj, features = 'nCount_RNA')
p6 <- VlnPlot(A_obj, features = 'nFeature_RNA')


p7 <- DimPlot(A_obj, reduction = 'umap', label = TRUE, repel = TRUE)
p8 <- DimPlot(A_obj, reduction = 'umap', label = TRUE, repel = TRUE, split.by = 'sample')
p9 <- DimPlot(A_obj, reduction = 'umap', label = TRUE, repel = TRUE, split.by = 'genotype')
p10 <- DimPlot(A_obj, reduction = 'umap', label = TRUE, repel = TRUE, split.by = 'status')
p11 <- DimPlot(A_obj, reduction = 'umap', split.by = 'genotype', group.by = 'sample') + scale_colour_manual(values=cbbPalette)

pdf('./data_preprocessing/outputs/plots/A_merged_filtered_clusteres_genotyped.pdf', width = 11, height = 8, paper = 'a4r')
p / (p1 | p2 | p3)
p7
p8
p9
p10
p11
p4 / p5 / p6 & theme(legend.position = 'none')
dev.off()
