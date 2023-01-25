# Preprocessing of F samples after demultiplexing
setwd('./Desktop/IMIM/Groningen/RPII/Projects/NR03_scRNAseq_main_analysis/data_preprocessing/')
source('C:/Users/nvrib/Desktop/IMIM/Groningen/RPII/Projects/scripts/function_plot_markers.R')

# Color blind palette
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Load libraries
library(Seurat)
library(tidyverse)
library(gridExtra)
library(SeuratDisk)

# Load the dataset
F_obj <- readRDS('./outputs/F/SeuratObj_F_filtered_genotyped.rds')

# Excluding ambiguous and doblets
F_obj <- subset(F_obj, subset = droplet == 'singlet')

# Save this one
SaveH5Seurat(F_obj, './outputs/F/SeuratObj_F_onlysinglets_genotyped', overwrite = T)

# Normalizing and clustering
F_obj <- SCTransform(F_obj, vars.to.regress = "percent.mt", verbose = T)
F_obj <- RunPCA(F_obj)
F_obj <- FindNeighbors(F_obj, dims = 1:50)
F_obj <- FindClusters(F_obj, resolution = 0.4)
F_obj <- RunUMAP(F_obj, dims = 1:50)

SaveH5Seurat(F_obj, './outputs/F/SeuratObj_F_singlets_genotyped_clustered_res04', overwrite = T)
#F_obj <- LoadH5Seurat('./outputs/F/SeuratObj_F_singlets_genotyped_clustered_res04.h5seurat')

# Plots
png('./outputs/plots/UMAP_F_genotypes.png', width = 900, height = 400)
DimPlot(F_obj, reduction = 'umap', split.by = 'genotype', label = TRUE)
dev.off()

png('./outputs/plots/UMAP_F_genotyped_clusters.png', width = 600, height = 500)
DimPlot(F_obj, reduction = 'umap', label = TRUE)
dev.off()

VlnPlot(F_obj, features = 'percent.mt')
VlnPlot(F_obj, features = 'nCount_RNA')
VlnPlot(F_obj, features = 'nFeature_RNA')

### Trying to figure out what is happening with samples from 19_CeDNN-A0338
plotMarkers(F_obj, markers_csv = '../../resources/cell_markers.csv', output = './outputs/plots/F_')

FeaturePlot(F_obj, features = 'EPCAM', split.by = 'genotype')

# Integrating based on genotype using Harmony
library(harmony)
F_obj <- RunHarmony(F_obj, "genotype")
F_obj <- F_obj %>% 
  RunUMAP(reduction = "harmony", dims = 1:50) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:50) %>% 
  FindClusters(resolution = 0.4)

SaveH5Seurat(F_obj, './outputs/F/SeuratObj_F_singlets_genotyped_clustered_res04_harmony', overwrite = T)

# Ploting again
png('./outputs/plots/UMAP_F_genotypes_harmony.png', width = 900, height = 400)
DimPlot(F_obj, reduction = 'umap', split.by = 'genotype', label = TRUE)
dev.off()

png('./outputs/plots/UMAP_F_clustered_harmony.png', width = 600, height = 600)
DimPlot(F_obj, reduction = 'umap', label = TRUE)
dev.off()

FeaturePlot(F_obj, features = 'EPCAM', split.by = 'genotype')
