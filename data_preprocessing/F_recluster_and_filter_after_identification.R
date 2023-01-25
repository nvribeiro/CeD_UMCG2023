# Re-filtering and clustering F batch applying the new %MT filters per cluster
setwd('C:/Users/nvrib/Desktop/IMIM/Groningen/RPII/Projects/NR03_scRNAseq_main_analysis/data_preprocessing/')
source('C:/Users/nvrib/Desktop/IMIM/Groningen/RPII/Projects/scripts/function_plot_markers.R')

# Load libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(SeuratDisk)

# Load the data - using the dataset without annotation just to make it easier
F_obj <- LoadH5Seurat('./outputs/F/SeuratObj_F_singlets_genotyped_clustered_res04_harmony.h5seurat')

# Subseting and applying the filters for each cluster
filter_10 <- subset(F_obj, subset = (seurat_clusters == 11 | seurat_clusters == 12 | seurat_clusters == 18)
                                    & percent.mt < 10)

filter_15 <- subset(F_obj, subset = (seurat_clusters == 0 | seurat_clusters == 1 | seurat_clusters == 3 |
                                    seurat_clusters == 5 | seurat_clusters == 7 | seurat_clusters == 14 |
                                    seurat_clusters == 17 | seurat_clusters == 19) & percent.mt < 15)

filter_30 <- subset(F_obj, subset = seurat_clusters == 16 & percent.mt < 30)

filter_65 <- subset(F_obj, subset = (seurat_clusters == 2 | seurat_clusters == 4)
                                    & percent.mt < 65)

filter_75 <- subset(F_obj, subset = (seurat_clusters == 8 | seurat_clusters == 9 | seurat_clusters == 13 |
                                    seurat_clusters == 15)
                                    & percent.mt < 75)

filter_80 <- subset(F_obj, subset = (seurat_clusters == 6 | seurat_clusters == 10)
                                    & percent.mt < 80)

# Merging them back together
F_obj_new_filters <- merge(x=filter_10, y=list(filter_15, filter_30, filter_65, filter_65, filter_75, filter_80))
SaveH5Seurat(F_obj_new_filters, './outputs/F/SeuraObj_f_reclustered_res04')

VlnPlot(F_obj_new_filters, features = 'percent.mt')

# Re-clustering
library(harmony)
F_obj_new_filters <- SCTransform(F_obj_new_filters, vars.to.regress = "percent.mt", verbose = T)
F_obj_new_filters <- RunPCA(F_obj_new_filters)
F_obj_new_filters <- RunHarmony(F_obj_new_filters, "genotype")

F_obj_new_filters <- F_obj_new_filters %>% 
  RunUMAP(reduction = "harmony", dims = 1:50) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:50) %>% 
  FindClusters(resolution = 0.4)

DimPlot(F_obj_new_filters, reduction = 'umap', label = TRUE) + NoLegend()
VlnPlot(F_obj_new_filters, features = 'percent.mt')

plotMarkers(F_obj_new_filters, markers_csv = '../../resources/cell_markers.csv', output = './outputs/plots/F_reclustered')

## Trying a different approach --------------------------------------------------
hfile <- Connect('./outputs/F/SeuratObj_F_clustered_res04_identified.h5seurat')
hfile$index()
hfile[["assays/SCT"]]
hfile$close_all()

#F_obj <- LoadH5Seurat('./outputs/F/SeuratObj_F_clustered_res04_identified.h5seurat', assays = list(SCT = c("data", "scale.data")))

F_obj <- LoadH5Seurat('./outputs/F/SeuratObj_F_clustered_res04_identified.h5seurat')

# Grouping all immune cells
immune <- subset(F_obj, subset = (seurat_clusters == 0 | seurat_clusters == 1 | seurat_clusters == 5 |
                                    seurat_clusters == 11 | seurat_clusters == 12 | seurat_clusters == 17 | seurat_clusters == 18))

# MT filter
VlnPlot(immune, features = 'percent.mt') + geom_hline(yintercept = 10)
immune <- subset(immune, subset = percent.mt < 10)
immune

VlnPlot(immune, features = 'percent.mt')

# Grouping all epithelial cells
epithelial <- subset(F_obj, subset = (seurat_clusters == 4 | seurat_clusters == 6 | seurat_clusters == 16 |
                                      seurat_clusters == 8 | seurat_clusters == 9 | seurat_clusters == 10 | seurat_clusters == 13) |
                                      (seurat_clusters == 7 & percent.mt < 20))
epithelial
 
VlnPlot(epithelial, features='percent.mt') + geom_hline(yintercept = 20)

# Other cells 
others <- subset(F_obj, subset = ((seurat_clusters == 3 | seurat_clusters == 14 | seurat_clusters == 19) & percent.mt < 15) |
                                  seurat_clusters == 15 | seurat_clusters == 2)
others

VlnPlot(others, features = 'percent.mt') + NoLegend() + geom_hline(yintercept = 15)

# Merging everything
F_new <- merge(x = immune, y = list(epithelial, others))
VlnPlot(F_new, features = 'percent.mt') + NoLegend()

# Re-cluster again
DefaultAssay(F_new) <- 'RNA'
F_new

F_new <- SCTransform(F_new, vars.to.regress = "percent.mt", verbose = T)
F_new <- RunPCA(F_new)

library(harmony)
F_new <- RunHarmony(F_new, "genotype")
F_new <- F_new %>% 
  RunUMAP(reduction = "harmony", dims = 1:50) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:50) %>% 
  FindClusters(resolution = 0.4)

DimPlot(F_new, reduction = 'umap', label = TRUE)
VlnPlot(F_new, features = 'percent.mt')

SaveH5Seurat(F_new, './outputs/F/SeuraObj_f_reclustered_res04.h5seurat', overwrite = T)

plotMarkers(F_new, markers_csv = '../../resources/cell_markers.csv', output = './outputs/plots/NEW2001_F_')
