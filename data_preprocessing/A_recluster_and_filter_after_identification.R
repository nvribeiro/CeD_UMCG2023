# Re-filtering and clustering A batch applying the new %MT filters per cluster
source('C:/Users/nvrib/Desktop/IMIM/Groningen/RPII/Projects/scripts/function_plot_markers.R')

# Load libraries
library(Seurat)
library(tidyverse)
library(SeuratDisk)

# Load the data
A_obj <- LoadH5Seurat('./data_preprocessing/outputs/A/SeuratObj_A_merged_genotyped_filtered_clustered.h5seurat')

# Grouping all immune cells
immune <- subset(A_obj, subset = (((seurat_clusters == 3 | seurat_clusters == 6 | seurat_clusters == 7 | seurat_clusters == 20 |
                                    seurat_clusters == 26 | seurat_clusters == 27) & percent.mt < 10) |
                                  (seurat_clusters == 8 | seurat_clusters == 9 | seurat_clusters == 10 | seurat_clusters == 12 | 
                                     seurat_clusters == 13 | seurat_clusters == 21 | seurat_clusters == 30) & percent.mt < 25))

VlnPlot(immune, features = 'percent.mt') + NoLegend()
immune

# Grouping all epithelial cells
epithelial <- subset(A_obj, subset =  (seurat_clusters == 1 | seurat_clusters == 2 |
                                      seurat_clusters == 4 | seurat_clusters == 14 | seurat_clusters == 16 | seurat_clusters == 17 |
                                      seurat_clusters == 23 |  seurat_clusters == 24 |  seurat_clusters == 25 | 
                                      seurat_clusters == 29) & percent.mt < 75)
epithelial
 
VlnPlot(epithelial, features='percent.mt') + NoLegend()

# Other cells 
others <- subset(A_obj, subset = (((seurat_clusters == 5 | seurat_clusters == 15 | seurat_clusters == 18 | seurat_clusters == 19 | seurat_clusters == 22) & percent.mt < 15) |
                                  ((seurat_clusters == 11 | seurat_clusters == 28) & percent.mt < 75)))
others

VlnPlot(others, features = 'percent.mt') + NoLegend()

# Merging everything
A_new <- merge(x = immune, y = list(epithelial, others))
VlnPlot(A_new, features = 'percent.mt') + NoLegend()

A_new

# Adding cluster names to the metadata
clusters.ids <- read.csv('./data_preprocessing/outputs/A/cell_IDs_A.csv')
tmp <- A_new@meta.data
tmp$seurat_clusters <- as.integer(tmp$seurat_clusters)
tmp <- left_join(tmp, clusters.ids, by = "seurat_clusters")

A_new@meta.data <- tmp

# Saving
SaveH5Seurat(A_new, './data_preprocessing/outputs/A/SeuratObj_A_identified_filtered')
