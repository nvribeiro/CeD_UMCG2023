# Preprocessing of A samples - doing some QC, checking for technical/batch effects

# Color-blind friendly palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")

# Load libraries
library(Seurat)
library(tidyverse)
library(SeuratDisk)

## Preprocessing A batch without lane A7 -------------------------------------

# Get data location
dirs <- str_subset(list.dirs(path = 'data', recursive = F, full.names = T), 'A')

# Create Seurat objects
for (x in dirs) {
  name <- gsub('data/','',x)
  
  data <- Read10X(paste0(x,'/'))
  assign(name, CreateSeuratObject(counts = data, project = 'NR03'))
}

obj_list <- list(A2 = A2, A3 = A3, A4 = A4, A5 = A5, A6 = A6)

# Calculating percent.mt
for (i in c(1:5)) {
  obj_list[[i]]$percent.mt <- PercentageFeatureSet(obj_list[[i]], pattern = '^MT-')
  obj_list[[i]]$sample <- paste0('A',i+1)
}
list2env(obj_list, .GlobalEnv)

# Merge the objects
A_merged <- merge(A2, y = c(A3,A4,A5,A6))
A_merged
# n cells = 77531
SaveH5Seurat(A_merged, file = './outputs/A/SeuratObj_A_merged_without_A7', overwrite = T)

# QC check
metadata <- A_merged@meta.data

pdf(file = './outputs/plots/QCPlots_A_wo_A7.pdf', width = 11, height = 8, paper = 'a4r')
VlnPlot(A_merged, features = 'nFeature_RNA', group.by = 'sample') + geom_hline(yintercept = 6000) + geom_hline(yintercept = 300)
VlnPlot(A_merged, features = 'nCount_RNA', group.by = 'sample') + geom_hline(yintercept = 50000) + geom_hline(yintercept = 1000)
VlnPlot(A_merged, features = 'percent.mt', group.by = 'sample') + geom_hline(yintercept = 90)

metadata %>%
  ggplot(aes(x=sample, fill=sample)) +
  geom_bar() +
  geom_text(stat='count', aes(label=after_stat(count)), vjust=-1) +
  theme_classic() +
  scale_fill_manual(values=cbbPalette) +
  labs(title='N cells per lane', y = 'n')

p1 <- metadata %>% 
  ggplot(aes(x=nCount_RNA, color=sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  scale_fill_manual(values=cbbPalette) +
  scale_colour_manual(values=cbbPalette) +
  geom_vline(xintercept = 1000) +
  geom_vline(xintercept = 50000) +
  labs(title='RNA counts per cell', y = 'Cell density')

p2 <- metadata %>% 
  ggplot(aes(x=nFeature_RNA, color=sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  scale_fill_manual(values=cbbPalette) +
  scale_colour_manual(values=cbbPalette) +
  geom_vline(xintercept = 200) +
  geom_vline(xintercept = 6000) +
  labs(title='Features per cell', y = 'Cell density')

p3 <- metadata %>% 
  ggplot(aes(x=percent.mt, color=sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  scale_fill_manual(values=cbbPalette) +
  scale_colour_manual(values=cbbPalette) +
  geom_vline(xintercept = 90) +
  labs(title='% MT RNA per cell', y = 'Cell density')

grid.arrange(p1, p2, p3, ncol = 3)

metadata %>%
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) +
  geom_point(size = 0.5) +
  scale_colour_gradient(low = "gray90", high = "blue") +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 1000) +
  geom_vline(xintercept = 50000) +
  geom_hline(yintercept = 200) +
  geom_hline(yintercept = 6000) +
  facet_wrap(~sample)

dev.off()

# Filtering
# percent.mt < 90%
# nFeature_RNA > 200 and < 6000
# nCount_RNA > 1000 and < 60000

A_merged <- LoadH5Seurat('./outputs/A/SeuratObj_A_merged_without_A7.h5seurat')

A_filtered <- subset(A_merged, subset = nCount_RNA > 1000 &
                       nCount_RNA < 60000 &
                       nFeature_RNA > 200 &
                       nFeature_RNA < 6000 &
                       percent.mt < 90)

# Normalizing and clustering
A_filtered <- SCTransform(A_filtered, vars.to.regress = "percent.mt")
A_filtered <- RunPCA(A_filtered)
A_filtered <- FindNeighbors(A_filtered, dims = 1:50)
A_filtered <- FindClusters(A_filtered, resolution = 0.4)
A_filtered <- RunUMAP(A_filtered, dims = 1:50)

# Saving
SaveH5Seurat(A_filtered, filename = './outputs/SeuratObj_A_filtered_clustered_woA7.h5Seurat', overwrite = T)

# Checking clusters
A_obj <- LoadH5Seurat('./outputs/A/SeuratObj_A_filtered_clustered_woA7.h5Seurat',
                      assays = 'SCT', reductions = 'umap', graphs = FALSE)

DimPlot(A_obj, reduction = 'umap')
VlnPlot(A_obj, features = 'percent.mt')
VlnPlot(A_obj, features = 'nCount_RNA')
VlnPlot(A_obj, features = 'nFeature_RNA')

##
## Preprocessing A7 alone --------------------------------------------------------
##

A7_obj <- Read10X('./data_preprocessing/data/A7/') %>%
  CreateSeuratObject(project = 'NR03')

A7_obj$percent.mt <- PercentageFeatureSet(A7_obj, pattern = '^MT-')

# QC
metadata <- A7_obj@meta.data

pdf(file = './data_preprocessing/outputs/plots/A7_QC_plots.pdf', width = 11, height = 8, paper = 'a4r')

VlnPlot(A7_obj, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'))

p1 <- metadata %>% 
  ggplot(aes(x=nCount_RNA)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  geom_vline(xintercept = 1000) +
  geom_vline(xintercept = 30000) +
  labs(title='RNA counts per cell', y = 'Cell density')

p2 <- metadata %>% 
  ggplot(aes(x=nFeature_RNA)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  geom_vline(xintercept = 300) +
  geom_vline(xintercept = 6000) +
  labs(title='Features per cell', y = 'Cell density')

p3 <- metadata %>% 
  ggplot(aes(x=percent.mt)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  geom_vline(xintercept = 90) +
  labs(title='% MT RNA per cell', y = 'Cell density')

grid.arrange(p1, p2, p3, ncol = 3)

metadata %>%
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) +
  geom_point(size = 0.5) +
  scale_colour_gradient(low = "gray90", high = "blue") +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 1000) +
  geom_vline(xintercept = 30000) +
  geom_hline(yintercept = 300) +
  geom_hline(yintercept = 6000)

dev.off()

## Archive - tryouts all A together - do not run -------------------------------
# Get data location
dirs <- str_subset(list.dirs(path = 'data/', recursive = F, full.names = T), 'A')

# Create Seurat objects
for (x in dirs) {
  name <- gsub('data//','',x)
  
  data <- Read10X(paste0(x,'/'))
  assign(name, CreateSeuratObject(counts = data, project = 'NR03'))
}

# QC check individually
obj_list <- list(A2 = A2, A3 = A3, A4 = A4, A5 = A5, A6 = A6, A7 = A7)

# Calculating percent.mt
for (i in c(1:6)) {
  obj_list[[i]]$percent.mt <- PercentageFeatureSet(obj_list[[i]], pattern = '^MT-')
  obj_list[[i]]$sample <- paste0('A',i+1)
}
list2env(obj_list, .GlobalEnv)

# Merge the objects
A_merged <- merge(A2, y = c(A3,A4,A5,A6,A7))
A_merged
# n cells = 77531
saveRDS(A_merged, file = './outputs/SeuratObj_A_merged.rds')


# QC check
metadata <- A_merged@meta.data

png('./outputs/plots/QC_plots_A_merged.png', width = 800, height = 600)
VlnPlot(A_merged, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3) & geom_jitter(alpha=0.1, size = 0.2)
dev.off()

VlnPlot(A_merged, features = 'nFeature_RNA', group.by = 'sample')
VlnPlot(A_merged, features = 'nCount_RNA', group.by = 'sample')
VlnPlot(A_merged, features = 'percent.mt', group.by = 'sample')

png('./outputs/plots/Ncells_A_samples.png', width = 800, height = 700)
metadata %>%
  ggplot(aes(x=sample, fill=sample)) +
  geom_bar() +
  geom_text(stat='count', aes(label=..count..), vjust=-1) +
  theme_classic() +
  scale_fill_manual(values=cbbPalette) +
  labs(title='N cells per lane', y = 'n')
dev.off()

p1 <- metadata %>% 
  ggplot(aes(x=nCount_RNA, color=sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  scale_fill_manual(values=cbbPalette) +
  scale_colour_manual(values=cbbPalette) +
  #geom_vline(xintercept = 1000) +
  #geom_vline(xintercept = 80000) +
  labs(title='RNA counts per cell', y = 'Cell density')

p2 <- metadata %>% 
  ggplot(aes(x=nFeature_RNA, color=sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  scale_fill_manual(values=cbbPalette) +
  scale_colour_manual(values=cbbPalette) +
  #geom_vline(xintercept = 300) +
  #geom_vline(xintercept = 10000) +
  labs(title='Features per cell', y = 'Cell density')

p3 <- metadata %>% 
  ggplot(aes(x=percent.mt, color=sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  scale_fill_manual(values=cbbPalette) +
  scale_colour_manual(values=cbbPalette) +
  #geom_vline(xintercept = 90) +
  labs(title='% MT RNA per cell', y = 'Cell density')

png('./outputs/plots/QC_distribution_A.png', width = 800, height = 400)
grid.arrange(p1, p2, p3, ncol = 3)
dev.off()

png('./outputs/plots/Feature_vs_Count_vs_MT_A.png', width = 800, height = 600)
metadata %>%
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) +
  geom_point(size = 0.5) +
  scale_colour_gradient(low = "gray90", high = "blue") +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  #geom_vline(xintercept = 1000) +
  #geom_vline(xintercept = 80000) +
  #geom_hline(yintercept = 300) +
  #geom_hline(yintercept = 10000) +
  facet_wrap(~sample)
dev.off()

## Standard workflow with no filtering - run in the cluster
A_merged <- SCTransform(A_merged, vars.to.regress = "percent.mt", verbose = T)
A_merged <- RunPCA(A_merged)

png('./outputs/plots/ElbowPlot_A_merged.png', width = 500, height = 500)
ElbowPlot(F_merged, ndims = 50)
dev.off()

A_merged <- FindNeighbors(A_merged, dims = 1:50)
A_merged <- FindClusters(A_merged, resolution = 0.4)
A_merged <- RunUMAP(A_merged, dims = 1:50)

saveRDS(A_merged, file = './outputs/SeuratObj_A_merged_clustered_res04.rds')
##

# Inspecting the clusters
A_clustered <- readRDS('./outputs/SeuratObj_A_merged_clustered_res04.rds')

png('./outputs/plots/UMAP_A_merged_nofilters.png', width = 800, height = 700)
DimPlot(A_clustered, reduction = 'umap', group.by = 'sample', cols = cbbPalette)
dev.off()

# Checking QC and deciding filters
png('./outputs/plots/percent.mt_A_merged_nofilters.png', width = 800, height = 500)
VlnPlot(A_clustered, features = 'percent.mt') & geom_hline(yintercept = 90)
dev.off()

png('./outputs/plots/nFeatures_A_merged_nofilters.png', width = 800, height = 500)
VlnPlot(A_clustered, features = 'nFeature_RNA') & geom_hline(yintercept = 300) & geom_hline(yintercept = 7500)
dev.off()

png('./outputs/plots/nCount_A_merged_nofilters.png', width = 800, height = 500)
VlnPlot(A_clustered, features = 'nCount_RNA') & geom_hline(yintercept = 2500) & geom_hline(yintercept = 60000)
dev.off()

metadata <- A_clustered@meta.data

p1 <- metadata %>% 
  ggplot(aes(x=nCount_RNA, color=sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  geom_vline(xintercept = 2500) +
  geom_vline(xintercept = 60000) +
  labs(title='RNA counts per cell', y = 'Cell density')

p2 <- metadata %>% 
  ggplot(aes(x=nFeature_RNA, color=sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  geom_vline(xintercept = 300) +
  geom_vline(xintercept = 7500) +
  labs(title='Features per cell', y = 'Cell density')

p3 <- metadata %>% 
  ggplot(aes(x=percent.mt, color=sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  geom_vline(xintercept = 90) +
  labs(title='% MT RNA per cell', y = 'Cell density')

png('./outputs/plots/QC_distribution_A_clustered_nofilter.png', width = 800, height = 400)
grid.arrange(p1, p2, p3, ncol = 3)
dev.off()

png('./outputs/plots/Feature_vs_Count_vs_MT_A_clustered_nofilter.png', width = 800, height = 600)
metadata %>%
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) +
  geom_point(size = 0.5) +
  scale_colour_gradient(low = "gray90", high = "blue") +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 2500) +
  geom_vline(xintercept = 60000) +
  geom_hline(yintercept = 300) +
  geom_hline(yintercept = 7500) +
  facet_wrap(~sample)
dev.off()


##### Now working in the cluster 
# Load the merged dataset
A_merged <- readRDS('./outputs/SeuratObj_A_merged.rds')

# Applying filters: percent.mt < 90%
# nFeature_RNA > 300 and < 75000
# nCount_RNA > 2500 and < 60000

A_filtered <- subset(A_merged, subset = nCount_RNA > 2500 &
                       nCount_RNA < 60000 &
                       nFeature_RNA > 300 &
                       nFeature_RNA < 7500 &
                       percent.mt < 90)

saveRDS(A_filtered, file = './outputs/SeuratObj_A_filtered.rds')

# Normalizing and clustering again
A_filtered <- SCTransform(A_filtered, vars.to.regress = "percent.mt", verbose = F)
A_filtered <- RunPCA(A_filtered)
A_filtered <- FindNeighbors(A_filtered, dims = 1:50)
A_filtered <- FindClusters(A_filtered, resolution = 0.4)
A_filtered <- RunUMAP(A_filtered, dims = 1:50)

# Saving as RDS and H5file
saveRDS(A_filtered, file = './outputs/SeuratObj_A_filtered_clustered_dim04.rds')
SaveH5Seurat(A_filtered, filename = './outputs/SeuratObj_A_filtered_clustered_dim04.h5Seurat', overwrite = T)

##
# for h5 files:
SaveH5Seurat(A_merged, filename = './outputs/SeuratObj_A_merged.h5Seurat', overwrite = T) # saves Seurar as H5
A_merged_h <- Connect('./outputs/SeuratObj_A_merged.h5Seurat') # connect and reads the file
A_merged_h$close_all() # closes connection to the file
##

# Checking the filtered clusters
A_filtered_h <- Connect('./outputs/SeuratObj_A_filtered_clustered_dim04.h5Seurat')
A_filtered_h
A_filtered_h$index()

# Getting just the data I need
A_data <- LoadH5Seurat(A_filtered_h, assays = 'SCT')

png('./outputs/plots/UMAP_A_merged_filtered_res04.png', width = 700, height = 700)
DimPlot(A_data, reduction = 'umap')
dev.off()

png('./outputs/plots/percentmt_A_merged_filtered.png', width = 800, height = 500)
VlnPlot(A_data,  features = 'percent.mt')
dev.off()

png('./outputs/plots/nFeature_A_merged_filtered.png', width = 800, height = 500)
VlnPlot(A_data,  features = 'nFeature_RNA')
dev.off()

png('./outputs/plots/nCount_A_merged_filtered.png', width = 800, height = 500)
VlnPlot(A_data,  features = 'nCount_RNA')
dev.off()

# How many cells I have left?
A_filtered <- readRDS('./outputs/SeuratObj_A_filtered.rds')
A_merged <- readRDS('./outputs/SeuratObj_A_merged.rds')