# Merging and integrating samples A and F after preprocessing
# Running this in the cluster

# Load libraries
library(Seurat)
library(tidyverse)
library(SeuratDisk)
library(patchwork)
library(harmony)

path <- './outputs/'

A_samples <- LoadH5Seurat(paste0(path,''))
F_samples <- LoadH5Seurat(paste0(path,''))

DefaultAssay(F_samples) <- 'RNA'
DefaultAssay(A_samples) <- 'RNA'

# Merge datasets
All_samples <- merge(A_samples, F_samples)

# Clustering
All_samples <- All_samples %>% 
  SCTransform(All_samples, vars.to.regress = "percent.mt", verbose = F) %>%
  RunPCA(All_samples) %>%
  FindNeighbors(All_samples, dims = 1:50) %>%
  FindClusters(All_samples, resolution = 0.4) %>%
  RunUMAP(All_samples, dims = 1:50)

SaveH5Seurat(All_samples, filename = paste0(path, 'All_samples_merged_v2'))

# QC inspection
All_clustered <- Connect('./outputs/All_samples_clustered.h5Seurat')
All_clustered$index()
All_clustered$close_all()

UMAP_data <- LoadH5Seurat('./outputs/All_samples_clustered.h5Seurat', assay = 'SCT', reductions = 'umap')
UMAP_data

metadata <- UMAP_data@meta.data

# Cells per sample
metadata %>% count(sample)

# Grouping samples from the same batch
UMAP_data@meta.data <- UMAP_data@meta.data %>% mutate(batch = case_when(
  sample == 'A2' | sample == 'A3' | sample == 'A4' | sample == 'A5' | sample == 'A6' | sample == 'A7' ~ 'A',
  sample == 'F1' | sample == 'F2' | sample == 'F3' ~ 'F'
))

# Checking the clusters
png('./outputs/plots/UMAP_allsamples_batch.png', width = 700, height = 700)
DimPlot(UMAP_data, reduction = 'umap', group.by = 'batch')
dev.off()

png('./outputs/plots/UMAP_allsamples_batch_split.png', width = 700, height = 400)
DimPlot(UMAP_data, reduction = 'umap', group.by = 'batch', split.by = 'batch')
dev.off()

png('./outputs/plots/UMAP_allsamples.png', width = 700, height = 700)
DimPlot(UMAP_data, reduction = 'umap')
dev.off()

png('./outputs/plots/percentmt_Allsamples.png', width = 800, height = 500)
VlnPlot(UMAP_data,  features = 'percent.mt')
dev.off()

png('./outputs/plots/nFeature_Allsamples.png', width = 800, height = 500)
VlnPlot(UMAP_data,  features = 'nFeature_RNA')
dev.off()

png('./outputs/plots/nCount_Allsamples.png', width = 800, height = 500)
VlnPlot(UMAP_data,  features = 'nCount_RNA')
dev.off()

####### Trying integration with Harmony - run in the cluster
# Dataset with all cells
All_samples <- LoadH5Seurat('./outputs/All_samples_merged.h5Seurat')

# Add the batch column to this dataset
All_samples@meta.data <- All_samples@meta.data %>% mutate(batch = case_when(
  sample == 'A2' | sample == 'A3' | sample == 'A4' | sample == 'A5' | sample == 'A6' | sample == 'A7' ~ 'A',
  sample == 'F1' | sample == 'F2' | sample == 'F3' ~ 'F'
))

# Normalization and PCA
# Clustering
All_samples <- SCTransform(All_samples, vars.to.regress = "percent.mt", verbose = F)
All_samples <- RunPCA(All_samples)

# Run Harmony
All_samples <- RunHarmony(All_samples, "batch")

# Clustering
All_samples <- All_samples %>% 
  RunUMAP(reduction = "harmony", dims = 1:50) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:50) %>% 
  FindClusters(resolution = 0.4)

# Save the object
SaveH5Seurat(All_samples, filename = './outputs/All_samples_integrated_harmony.h5Seurat')

#####

# Inspecting Harmony results
harmony <- Connect('./outputs/All_samples_integrated_harmony.h5Seurat')
harmony$index()
harmony$close_all()

harmony_obj <- LoadH5Seurat('./outputs/All_samples_integrated_harmony.h5Seurat',
                            assays = list(SCT = c("data", "scale.data")),
                            reductions = 'umap',
                            graphs = FALSE,)

View(harmony_obj@meta.data)

png('./outputs/plots/UMAP_all_harmony_bybatch.png', width = 700, height = 400)
DimPlot(harmony_obj, reduction = 'umap', group.by = 'batch', split.by = 'batch')
dev.off()

png('./outputs/plots/UMAP_all_harmony_bysample.png', width = 700, height = 400)
DimPlot(harmony_obj, reduction = 'umap', group.by = 'sample', split.by = 'batch')
dev.off()

png('./outputs/plots/UMAP_all_harmony_clusters.png', width = 700, height = 600)
DimPlot(harmony_obj, reduction = 'umap')
dev.off()

VlnPlot(harmony_obj, features = 'percent.mt')
VlnPlot(harmony_obj, features = 'nFeature_RNA')
VlnPlot(harmony_obj, features = 'nCount_RNA')
