# Preprocessing of F samples - doing some QC, checking for technical/batch effects
setwd('./Desktop/IMIM/Groningen/RPII/Projects/NR03_scRNAseq_main_analysis/data_preprocessing/')

# Load libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)

# Get data location
dirs <- str_subset(list.dirs(path = 'data/', recursive = F, full.names = T), 'F')

# Create Seurat objects
for (x in dirs) {
  name <- gsub('data//','',x)
  
  data <- Read10X(paste0(x,'/'))
  assign(name, CreateSeuratObject(counts = data, project = 'NR03'))
}

F1
F2
F3

# QC check individually
F1$percent.mt <- PercentageFeatureSet(F1, pattern = '^MT-')
F2$percent.mt <- PercentageFeatureSet(F2, pattern = '^MT-')
F3$percent.mt <- PercentageFeatureSet(F3, pattern = '^MT-')

VlnPlot(F1, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3) & geom_jitter(alpha=0.1, size = 0.5)
VlnPlot(F2, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3) & geom_jitter(alpha=0.1, size = 0.5)
VlnPlot(F3, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3) & geom_jitter(alpha=0.1, size = 0.5)

# Merging the datasets
F1$sample <- 'F1'
F2$sample <- 'F2'
F3$sample <- 'F3'

F_merged <- merge(F1, y = c(F2, F3))
View(F_merged@meta.data)
saveRDS(F_merged, file = './outputs/SeuratObj_F_merged.rds')

# QC merged
png('./outputs/plots/QC_plots_F_merged.png', width = 800, height = 600)
VlnPlot(F_merged, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3) & geom_jitter(alpha=0.1, size = 0.2)
dev.off()

# Inspecting nFeatures to decide filter
VlnPlot(F_merged, features = 'nFeature_RNA') & geom_jitter(alpha=0.1, size = 0.1) & scale_y_continuous(limits = c(0,5000)) & geom_hline(yintercept = 300)

# Standard workflow with no filtering
F_merged <- SCTransform(F_merged, vars.to.regress = "percent.mt", verbose = T)
F_merged <- RunPCA(F_merged)

png('./outputs/plots/ElbowPlot_F_merged.png', width = 500, height = 500)
ElbowPlot(F_merged, ndims = 50)
dev.off()

F_merged <- FindNeighbors(F_merged, dims = 1:50)
F_merged <- FindClusters(F_merged, resolution = 0.4)
F_merged <- RunUMAP(F_merged, dims = 1:50)

saveRDS(F_merged, file = './outputs/SeuratObj_F_merged_clustered_res04.rds')

F_clustered <- readRDS('./outputs/SeuratObj_F_merged_clustered_res04.rds')

# Ploting UMAP
png('./outputs/plots/UMAP_F_merged_nofilters.png', width = 600, height = 600)
DimPlot(F_clustered, reduction = 'umap', group.by = 'sample', cols = c('red', 'yellow', 'blue'))
dev.off()

# Checking QC and deciding filters
png('./outputs/plots/percent.mt_F_merged_nofilters.png', width = 800, height = 500)
VlnPlot(F_clustered, features = 'percent.mt') & geom_hline(yintercept = 90)
dev.off()

png('./outputs/plots/nFeatures_F_merged_nofilters.png', width = 800, height = 500)
VlnPlot(F_clustered, features = 'nFeature_RNA') & geom_hline(yintercept = 300) & geom_hline(yintercept = 10000)
dev.off()

png('./outputs/plots/nCount_F_merged_nofilters.png', width = 800, height = 500)
VlnPlot(F_clustered, features = 'nCount_RNA') & geom_hline(yintercept = 2500) & geom_hline(yintercept = 80000)
dev.off()

FeaturePlot(F_clustered, features = c('percent.mt','nFeature_RNA'))

# Loading the dataset before clustering again
F_merged <- readRDS('./outputs/SeuratObj_F_merged.rds')
metadata <- F_merged@meta.data

# Some other QC check
png('./outputs/plots/Ncells_F_samples.png', width = 800, height = 700)
metadata %>%
  ggplot(aes(x=sample, fill=sample)) +
  geom_bar() +
  geom_text(stat='count', aes(label=..count..), vjust=-1) +
  theme_classic() +
  labs(title='N cells per lane', y = 'n')
dev.off()

p1 <- metadata %>% 
  ggplot(aes(x=nCount_RNA, color=sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  geom_vline(xintercept = 1000) +
  geom_vline(xintercept = 80000) +
  labs(title='RNA counts per cell', y = 'Cell density')

p2 <- metadata %>% 
  ggplot(aes(x=nFeature_RNA, color=sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  geom_vline(xintercept = 300) +
  geom_vline(xintercept = 10000) +
  labs(title='Features per cell', y = 'Cell density')

p3 <- metadata %>% 
  ggplot(aes(x=percent.mt, color=sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  geom_vline(xintercept = 90) +
  labs(title='% MT RNA per cell', y = 'Cell density')

png('./outputs/plots/QC_distribution_F.png', width = 800, height = 400)
grid.arrange(p1, p2, p3, ncol = 3)
dev.off()

png('./outputs/plots/Feature_vs_Count_vs_MT_F.png', width = 800, height = 600)
metadata %>%
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) +
  geom_point(size = 0.5) +
  scale_colour_gradient(low = "gray90", high = "blue") +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 1000) +
  geom_vline(xintercept = 80000) +
  geom_hline(yintercept = 300) +
  geom_hline(yintercept = 10000) +
  facet_wrap(~sample)
dev.off()

# Applying filters percent.mt < 90%
# nFeature_RNA > 300 and < 10000
# nCount_RNA > 2500 and < 80000
F_filtered <- subset(F_merged, subset = nCount_RNA > 1000 &
                       nCount_RNA < 80000 &
                       nFeature_RNA > 300 &
                       nFeature_RNA < 10000 &
                       percent.mt < 90)

saveRDS(F_filtered, file = './outputs/SeuratObj_F_filtered.rds')

## Running the code below in the cluster
F_filtered <- readRDS('./outputs/SeuratObj_F_filtered.rds')

# Normalizing and clustering again
F_filtered <- SCTransform(F_filtered, vars.to.regress = "percent.mt", verbose = F)
F_filtered <- RunPCA(F_filtered)
F_filtered <- FindNeighbors(F_filtered, dims = 1:50)
F_filtered <- FindClusters(F_filtered, resolution = 0.4)
F_filtered <- RunUMAP(F_filtered, dims = 1:50)

saveRDS(F_filtered, file = './outputs/F_filtered_clustered_dim04.rds')
##

# Checking the new clusters
F_filtered_clustered <- readRDS('./outputs/SeuratObj_F_filtered_clustered_res04.rds')

png('./outputs/plots/UMAP_F_merged_filtered_res04.png', width = 700, height = 700)
DimPlot(F_filtered_clustered, reduction = 'umap')
dev.off()

png('./outputs/plots/percentmt_F_merged_filtered.png', width = 800, height = 500)
VlnPlot(F_filtered_clustered,  features = 'percent.mt')
dev.off()

png('./outputs/plots/nFeature_F_merged_filtered.png', width = 800, height = 500)
VlnPlot(F_filtered_clustered,  features = 'nFeature_RNA')
dev.off()

png('./outputs/plots/nCount_F_merged_filtered.png', width = 800, height = 500)
VlnPlot(F_filtered_clustered,  features = 'nCount_RNA')
dev.off()

# How many cells I have left?
F_merged <- readRDS('./outputs/SeuratObj_F_merged.rds')
F_filtered <- readRDS('./outputs/SeuratObj_F_filtered.rds')
