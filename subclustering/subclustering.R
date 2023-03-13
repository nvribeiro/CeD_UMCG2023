set.seed(1234)
library(Seurat)
library(tidyverse)
library(harmony)
source('functions/my_functions.R')

out_path <- '/groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/NR03_scRNAseq/ongoing/outputs/subclustering/'
in_path <- '/groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/NR03_scRNAseq/processed/data_preprocessing/'

## Subclustering epithelial cells ----------------------------------------------
# Subseting
obj <- readRDS(paste0(in_path, 'SeuratObj_AllCells_SCT_RES04_clustered_final_annotated.rds'))
epithelia <- subset(obj, subset = cell.type.1 == 'Epithelial')

DefaultAssay(epithelia) <- 'RNA'

# Normalizing
epithelia <- SCTransform(epithelia, vars.to.regress = "percent.mt", verbose = T)

# Reclustering
epithelia <- epithelia %>%
  RunPCA() %>%
  RunHarmony("batch") %>%
  RunUMAP(reduction = "harmony", dims = 1:50) %>%
  FindNeighbors(reduction = "harmony", dims = 1:50)

epithelia <- FindClusters(epithelia, resolution = 0.7)

saveRDS(epithelia, paste0(out_path, 'epithelial_SeuratObj_res07.rds'))
#epithelia <- readRDS(paste0(out_path, 'epithelial_SeuratObj_res07.rds'))

pdf(paste0(out_path, 'only_epithelial_res07_umap.pdf'), width = 11, height = 8, paper = 'a4r')
DimPlot(epithelia, reduction = 'umap', label = T) + NoLegend()
DimPlot(epithelia, reduction = 'umap', split.by = 'status') + NoLegend()
dev.off()

markers <- '../resources/epithelial_markers_detailed.csv'
PlotUCell(epithelia, markers, output = paste0(out_path, 'only_epithelial_reclustered_UCell'))
plotMarkers(epithelia, markers, output = paste0(out_path, 'only_epithelial_reclustered'))

# Finding markers
All_markers <- FindAllMarkers(epithelia, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(All_markers, paste0(out_path, 'only_epithelial_allmarkers.csv'), row.names = F)

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

write.csv(top.markers.clean, paste0(out_path, 'only_epithelial_top20markers.csv'), row.names = F)

# Inspections
FeaturePlot(epithelia, features = c('EPCAM', 'PTPRC'))

# Cluster 13 apperars to be doublets of epithelial and immune cells - removed but saved for future inspection
epi_immune_doublets <- subset(epithelia, subset = seurat_clusters == 13)
saveRDS(epi_immune_doublets, paste0(out_path, 'epithelial_immune_doublets.rds'))

epithelia <- subset(epithelia, subset = seurat_clusters != 13)

## Running normalization, clustering and findmarkers again...
DefaultAssay(epithelia) <- 'RNA'

epithelia <- SCTransform(epithelia, vars.to.regress = "percent.mt", verbose = T)

epithelia <- epithelia %>%
  RunPCA() %>%
  RunHarmony("batch") %>%
  RunUMAP(reduction = "harmony", dims = 1:50) %>%
  FindNeighbors(reduction = "harmony", dims = 1:50)

epithelia <- FindClusters(epithelia, resolution = 0.7)

saveRDS(epithelia, paste0(out_path, 'epithelial_SeuratObj_res07_v2.rds'))

pdf(paste0(out_path, 'only_epithelial_res07_umap_v2.pdf'), width = 11, height = 8, paper = 'a4r')
DimPlot(epithelia, reduction = 'umap', label = T) + NoLegend()
DimPlot(epithelia, reduction = 'umap', split.by = 'status') + NoLegend()
dev.off()

markers <- '../resources/epithelial_markers_detailed.csv'
PlotUCell(epithelia, markers, output = paste0(out_path, 'only_epithelial_reclustered_UCell_v2'))
plotMarkers(epithelia, markers, output = paste0(out_path, 'only_epithelial_reclustered_v2'))

# Finding markers
All_markers <- FindAllMarkers(epithelia, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(All_markers, paste0(out_path, 'only_epithelial_allmarkers_v2.csv'), row.names = F)

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

write.csv(top.markers.clean, paste0(out_path, 'only_epithelial_top20markers_v2.csv'), row.names = F)

# Inspections
# EEC Precursors Neurog3 Neurod1 Sox4
FeaturePlot(epithelia, features = c('NEUROG3', 'NEUROD1', 'SOX4'), label = T)
FeaturePlot(epithelia, features = c('NEUROG3', 'NEUROD1', 'SOX4'), label = T)
VlnPlot(epithelia, features = c('NEUROG3', 'NEUROD1', 'SOX4'))
FeaturePlot(epithelia, features = c('NEUROG3', 'CLPS', 'MDK'), label = T)

## EEC L
FeaturePlot(epithelia, features = c('PYY', 'GCG', 'UCN3'), label = T)

DimPlot(epithelia, reduction = 'umap', label = T) + NoLegend()
DimPlot(epithelia, reduction = 'umap', split.by = 'genotype', label = T) + NoLegend()

# Checking different enterocytes clusters
enterocytes <- subset(epithelia, subset = seurat_clusters == 0 | seurat_clusters == 2 | seurat_clusters == 3 | seurat_clusters == 5 | seurat_clusters == 6)

# Getting the marker genes from Moor et al 2018 - don't run
#library(readxl)
#enterocytes_markers <- read_excel('../../resources/gene_markers_enterocytes_Moor2018.xlsx')
#mouse_genes <- c(enterocytes_markers$`Cluster 1`, enterocytes_markers$`Cluster 2`, enterocytes_markers$`Cluster 3`, enterocytes_markers$`Cluster 4`, enterocytes_markers$`Cluster 5`)
#mouse_genes <- unique(mouse_genes)
#write_clip(cat(paste(mouse_genes, "\n")))

# Convert these genes to the human equivalent
#human_genes <- read.csv('../../resources/human_orthologs_moor2018.txt') %>% pull(Human.gene.name)
human_genes_bottom <- read.csv('../resources/human_orthologs_moor2018_bottom_landmark.txt') %>% pull(Human.gene.name)
human_genes_top <- read.csv('../resources/human_orthologs_moor2018_top_landmark.txt') %>% pull(Human.gene.name)
burclaff_markers <- c('CA2', 'AKR1B10', 'SLC4A4', 'PCK1', 'PFKFB2', 'APOB', 'APOA4', 'SI', 'FABP2', 'ANPEP', 'PRAP1', 'FANP6', 'TMIGD1', 'PBL1', 'CUBN', 'SLC10A2', 'CLCA4', 'SLC13A1', 'FAM151A', 'SLC13A1', 'SLC17A8')
enterocytes_avg <- AverageExpression(enterocytes, assay = 'SCT', slot = 'scale.data', return.seurat =  T)

pdf(paste0(out_path, 'enterocytes_markers_topbottom.pdf'), height = 11, width = 8, paper = 'a4')
DoHeatmap(enterocytes, features = human_genes_bottom, raster = FALSE) + labs(title = 'Bottom landmark genes')
DoHeatmap(enterocytes, features = human_genes_top, raster = FALSE) + labs(title = 'Top landmark genes')
DoHeatmap(enterocytes_avg, features = human_genes_bottom, raster = FALSE, draw.lines = F) + labs(title = 'Bottom landmark genes - avg cluster expression')
DoHeatmap(enterocytes_avg, features = human_genes_top, raster = FALSE, draw.lines = F) + labs(title = 'Top landmark genes - avg cluster expression')
DoHeatmap(enterocytes, features = burclaff_markers, raster = FALSE) + labs(title = 'Burclaff marker genes')
DoHeatmap(enterocytes_avg, features = burclaff_markers, raster = FALSE, draw.lines = F) + labs(title = 'Burclaff marker genes - avg cluster expression')
dev.off()

# Plotting my markers from FindAllMarkers
#all_markers <- read_csv('data_preprocessing/outputs/subclusters/only_epithelial_markers.csv')
#enterocytes_markers <- all_markers %>%
#  filter(cluster == 0 | cluster == 1 | cluster == 2 | cluster == 3 | cluster == 4 | cluster == 6) %>%
#  distinct(gene) %>%
#  pull(gene)

# Adding annotation to metadata
clusters.ids <- read.csv(paste0(out_path, 'only_epithelial_top20markers_v2_identity.csv'))
clusters.ids <- clusters.ids %>% select(-markers)
clusters.ids$cluster <- as.factor(clusters.ids$cluster)

metadata <- epithelia@meta.data
metadata <- left_join(metadata, clusters.ids, by = c('seurat_clusters' = 'cluster'))
row.names(metadata) <- metadata$barcode

epithelia@meta.data <- metadata

saveRDS(epithelia, paste0(out_path, 'only_epithelial_res07_indentified_v2.rds'))
#epithelia <- readRDS(paste0(out_path, 'only_epithelial_res07_indentified.rds'))

pdf(paste0(out_path, 'only_epithelial_res07_umap_identified_v2.pdf'), width = 11, height = 8, paper = 'a4r')
DimPlot(epithelia, reduction = 'umap', label = T, repel = T, group.by = 'cell.type.3') + NoLegend()
DimPlot(epithelia, reduction = 'umap', group.by = 'cell.type.3', split.by = 'status') + NoLegend()
dev.off()

# Saving table for scDODA
# Add sample info
epithelia@meta.data <- epithelia@meta.data %>%
  mutate(sample = case_when(
    genotype == '30_219-0607' ~ 'CeD1',
    genotype == '51_219-1056' ~ 'CeD2',
    genotype == '13_CeDNN-A0325' ~ 'Ctrl1',
    genotype == '15_CeDNN-A0326' ~ 'Ctrl2',
    genotype == '20_CeDNN-A0339' ~ 'CeD3',
    genotype == '23_CeDNN-A0357' ~ 'Ctrl4',
    genotype == '24_CeDNN-A0381' ~ 'Ctrl5',
    genotype == '21_CeDNN-A0340' ~  'CeD4'
  ), .after = genotype)

cell_counts <- epithelia@meta.data %>%
  group_by(sample, age, sex, batch, status, cell.type.3) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  pivot_wider(names_from = cell.type.3, values_from = count)

cell_counts <- replace(cell_counts, is.na(cell_counts), 0)
write.csv(cell_counts, paste0(out_path, 'epithelial_cells_counts.csv'), row.names = F)

## Reordering enterocytes clusters based on trajectory analysis
epithelia <- readRDS(paste0(out_path, 'only_epithelial_res07_indentified_v2.rds'))

clusters.ids <- read.csv(paste0(out_path, 'only_epithelial_top20markers_v2_identity.csv'))
clusters.ids <- clusters.ids %>% select(-markers)
clusters.ids$cluster <- as.factor(clusters.ids$cluster)

metadata <- epithelia@meta.data
# Deleting old cell.type.3
metadata <- select(metadata, -cell.type.3)
metadata <- left_join(metadata, clusters.ids, by = c('seurat_clusters' = 'cluster'))
row.names(metadata) <- metadata$barcode

epithelia@meta.data <- metadata
DimPlot(epithelia, reduction = 'umap', group.by = 'cell.type.3', label= T) + NoLegend()

# Saving
saveRDS(epithelia, paste0(out_path, 'only_epithelial_res07_indentified_v3.rds'))

# Saving counts for scCODA
epithelia@meta.data <- epithelia@meta.data %>%
  mutate(sample = case_when(
    genotype == '30_219-0607' ~ 'CeD1',
    genotype == '51_219-1056' ~ 'CeD2',
    genotype == '13_CeDNN-A0325' ~ 'Ctrl1',
    genotype == '15_CeDNN-A0326' ~ 'Ctrl2',
    genotype == '20_CeDNN-A0339' ~ 'CeD3',
    genotype == '23_CeDNN-A0357' ~ 'Ctrl4',
    genotype == '24_CeDNN-A0381' ~ 'Ctrl5',
    genotype == '21_CeDNN-A0340' ~  'CeD4'
  ), .after = genotype)

cell_counts <- epithelia@meta.data %>%
  group_by(sample, age, sex, batch, status, cell.type.3) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  pivot_wider(names_from = cell.type.3, values_from = count)

cell_counts <- replace(cell_counts, is.na(cell_counts), 0)
write.csv(cell_counts, paste0(out_path, 'epithelial_cells_counts_v3.csv'), row.names = F)


# Saving SeuratObj as AnnData for scanpy - already with PCA
library(SeuratDisk)
epithelia <- readRDS(paste0(out_path, 'only_epithelial_res07_indentified_v3.rds'))

## Preparing dataset for PAGA - merging EEC clusters
DefaultAssay(epithelia) <- 'RNA'
epithelia$cell.type.4 <- epithelia$cell.type.3
epithelia$cell.type.4[epithelia$seurat_clusters %in% c(12, 14, 15, 16, 17, 18, 19)] <- 'EEC'


epithelia <- NormalizeData(epithelia, normalization.method = "LogNormalize")
epithelia <- FindVariableFeatures(epithelia, nfeatures = 2000)
epithelia <- ScaleData(epithelia)
epithelia <- RunPCA(epithelia, features = VariableFeatures(epithelia))

SaveH5Seurat(epithelia, filename = paste0(out_path, 'epithelia_v3.h5Seurat'), overwrite = T)
Convert(paste0(out_path, 'epithelia_v3.h5Seurat'), dest = 'h5ad', assay = 'RNA', overwrite = T)

# Just raw data
epithelia <- readRDS(paste0(out_path, 'only_epithelial_res07_indentified_v3.rds'))
DefaultAssay(epithelia) <- 'RNA'

SaveH5Seurat(epithelia, filename = paste0(out_path, 'epithelia_v3_RNAraw.h5Seurat'), overwrite = T)
Convert(paste0(out_path, 'epithelia_v3_RNAraw.h5Seurat'), dest = 'h5ad', assay = 'RNA', overwrite = T)

##

FeaturePlot(epithelia, features = c('ORC6', 'TOP2A'), label = T)
FeaturePlot(epithelia, features = c('CDX1', 'CDX2'), label = T)

## Increasing resolution of epithelial cells clustering --------------------------
epithelia <- readRDS(paste0(out_path, 'only_epithelial_res07_indentified_v3.rds'))

epithelia <- FindClusters(epithelia, resolution = 1.1)
DimPlot(epithelia, reduction = 'umap', label = T, repel = T) + NoLegend()

FeaturePlot(epithelia, features = c('ATOH1', 'SOX9', 'KLF4', 'HES1'))

markers <- '../resources/epithelial_markers_detailed.csv'
PlotUCell(epithelia, markers, output = paste0(out_path, 'only_epithelial_reclustered_UCell_V4'))
plotMarkers(epithelia, markers, output = paste0(out_path, 'only_epithelial_reclustered_V4'))

All_markers <- FindAllMarkers(epithelia, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(All_markers, paste0(out_path, 'only_epithelial_allmarkers_v4.csv'), row.names = F)

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
  summarise(markers = str_c(unlist(pick(gene)), collapse=', '))

write.csv(top.markers.clean, paste0(out_path, 'only_epithelial_top20markers_v4.csv'), row.names = F)

# Inspections
DimPlot(epithelia, reduction = 'umap', label = T, split.by = 'status') + NoLegend()
FeaturePlot(epithelia, features = c('LGR5', 'ASCL2', 'MYC'), label = T)
FeaturePlot(epithelia, features = c('DLL1', 'DLL4', 'LGR5'), label = T)
FeaturePlot(epithelia, features = c('FOXA2', 'INSM1', 'HES6'), label = T)
FeaturePlot(epithelia, features = c('NEUROG3'), label = T)
FeaturePlot(epithelia, features = c('BMI1', 'HOPX', 'TERT', 'LRIG1', 'LGR5'), label = T)
FeaturePlot(epithelia, features = c('HES1', 'HNF4A', 'HNF4G'), label = T)
FeaturePlot(epithelia, features = c('FABP2', 'APOE', 'FABP1'), label = T)
FeaturePlot(epithelia, features = c('ATOH1'), label = T)
FeaturePlot(epithelia, features = c('SPIB', 'GP2'), label = T)

# Cluster 27 is immune and 19 is endothelial - removing
FeaturePlot(epithelia, features = c('EPCAM', 'PTPRC', 'PECAM1'))
epithelia <- subset(epithelia, subset = seurat_clusters != 27 & seurat_clusters != 19)
DimPlot(epithelia, reduction = 'umap', label = T) + NoLegend()

# Standard workflow again
DefaultAssay(epithelia) <- 'RNA'
epithelia <- SCTransform(epithelia, vars.to.regress = "percent.mt", verbose = T)

epithelia <- epithelia %>%
  RunPCA() %>%
  RunHarmony("batch") %>%
  RunUMAP(reduction = "harmony", dims = 1:40) %>%
  FindNeighbors(reduction = "harmony", dims = 1:40)

epithelia <- FindClusters(epithelia, resolution = 0.9)
DimPlot(epithelia, reduction = 'umap', label = T, repel = T) + NoLegend()

markers <- '../resources/epithelial_markers_detailed.csv'
PlotUCell(epithelia, markers, output = paste0(out_path, 'only_epithelial_reclustered_UCell_V4'))
plotMarkers(epithelia, markers, output = paste0(out_path, 'only_epithelial_reclustered_V4'))

All_markers <- FindAllMarkers(epithelia, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(All_markers, paste0(out_path, 'only_epithelial_allmarkers_v4.csv'), row.names = F)

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
  summarise(markers = str_c(unlist(pick(gene)), collapse=', '))

write.csv(top.markers.clean, paste0(out_path, 'only_epithelial_top20markers_v4.csv'), row.names = F)

DimPlot(epithelia, reduction = 'umap', split.by = 'status', label = T) + NoLegend()

FeaturePlot(epithelia, features = c('ATOH1', 'SOX9', 'KLF4', 'HES1', 'LGR5'))
FeaturePlot(epithelia, features = c('MUC2', 'NEUROG3'), label = T)
FeaturePlot(epithelia, features = c('NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4'))
FeaturePlot(epithelia, features = c('FABP1', 'SI', 'APOA1'), label = T)
FeaturePlot(epithelia, features = c('EPHB2', 'SOX9'), label = T)
FeaturePlot(epithelia, features = c('SPIB', 'TRAF6', 'CLU'), label = T)

epithelia.avg <- AverageExpression(epithelia, assay = 'SCT', return.seurat = T)
DoHeatmap(epithelia.avg, features = c('EPHB2', 'SOX9'))
DoHeatmap(epithelia.avg, features = c('PIGR', 'SI'))
DoHeatmap(epithelia.avg, features = c('FABP1', 'SI', 'APOA1', 'REG1A', 'RPS23', 'RPL3'))
DoHeatmap(epithelia.avg, features = c('SPIB', 'TRAF6', 'CLU'))

saveRDS(epithelia, paste0(out_path, 'SeuratObj_epithelia_v4_res09.rds'))

# Adding id
clusters.ids <- read.csv(paste0(out_path, 'only_epithelial_top20markers_v4_identity.csv'))
clusters.ids <- clusters.ids %>% select(cluster, cell.type.3)
clusters.ids$cluster <- as.factor(clusters.ids$cluster)

metadata <- epithelia@meta.data
# Deleting old cell.type.3
metadata <- select(metadata, -cell.type.3)
metadata <- left_join(metadata, clusters.ids, by = c('seurat_clusters' = 'cluster'))
row.names(metadata) <- metadata$barcode

epithelia@meta.data <- metadata

# Cleaning the metadata a little
epithelia@meta.data <- dplyr::select(epithelia@meta.data, -grep('SCT_', colnames(metadata), value = T))

saveRDS(epithelia, paste0(out_path, 'SeuratObj_epithelia_v4_res09_identified.rds'))
epithelia <- readRDS(paste0(out_path, 'SeuratObj_epithelia_v4_res09_identified.rds'))

# Reorder conditions
epithelia$status <- as.factor(epithelia$status)
epithelia$status <- relevel(epithelia$status, 'Ctrl')
saveRDS(epithelia, paste0(out_path, 'SeuratObj_epithelia_v4_res09_identified.rds'))

DimPlot(epithelia, reduction = 'umap', group.by = 'cell.type.3', split.by = 'status', label = T, repel = T) + NoLegend()

## Preparing dataset for PAGA
library(SeuratDisk)

SaveH5Seurat(epithelia, filename = paste0(out_path, 'epithelia_v4_TAsplit.h5Seurat'), overwrite = T)
Convert(paste0(out_path, 'epithelia_v4_TAsplit.h5Seurat'), dest = 'h5ad', assay = 'SCT', overwrite = T)

# Cell counts for scCODA
cell_counts <- epithelia@meta.data %>%
  group_by(sample, age, sex, batch, status, cell.type.3) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  pivot_wider(names_from = cell.type.3, values_from = count)

cell_counts <- replace(cell_counts, is.na(cell_counts), 0)
write.csv(cell_counts, paste0(out_path, 'epithelial_cells_counts_v4_TAsplit.csv'), row.names = F)





## Subclustering immune cells --------------------------------------------------
# Get barcodes of immune group
obj <- readRDS(paste0(in_path, 'SeuratObj_AllCells_SCT_RES04_clustered_final_annotated.rds'))
immune <- subset(obj, subset = cell.type.1 == 'Immune')

DefaultAssay(immune) <- 'RNA'

# Normalizing
immune <- SCTransform(immune, vars.to.regress = "percent.mt", verbose = T)

# Reclustering
immune <- immune %>%
  RunPCA() %>%
  RunHarmony("batch") %>%
  RunUMAP(reduction = "harmony", dims = 1:50) %>%
  FindNeighbors(reduction = "harmony", dims = 1:50)

immune <- FindClusters(immune, resolution = 0.8)

saveRDS(immune, paste0(out_path, 'immune_SeuratObj_res08.rds'))
immune <- readRDS(paste0(out_path, 'immune_SeuratObj_res08.rds'))

# Checking
DimPlot(immune, reduction = 'umap', label = T) + NoLegend()
DimPlot(immune, reduction = 'umap', split.by = 'status', label = T) + NoLegend()

FeaturePlot(immune, features = c('PTPRC', 'EPCAM'), label = T)
VlnPlot(immune, features = 'percent.mt')

# Finding markers
All_markers <- FindAllMarkers(immune, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(All_markers, paste0(out_path, 'only_immune_markers.csv'), row.names = F)

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
  summarise(markers = str_c(unlist(pick(gene)), collapse=', '))

write.csv(top.markers.clean, paste0(out_path, 'only_immune_markers_top20.csv'), row.names = F)

markers <- '../resources/immune_markers.csv'
PlotUCell(immune, markers, output = paste0(out_path, 'only_immune_reclustered_UCell'))
plotMarkers(immune, markers, output = paste0(out_path, 'only_immune_reclustered'))

## Removing dead and reclustering
immune <- subset(immune, subset = seurat_clusters != 16 & seurat_clusters != 17)

DefaultAssay(immune) <- 'RNA'

# Normalizing
immune <- SCTransform(immune, vars.to.regress = "percent.mt", verbose = T)

# Reclustering
immune <- immune %>%
  RunPCA() %>%
  RunHarmony("batch") %>%
  RunUMAP(reduction = "harmony", dims = 1:50) %>%
  FindNeighbors(reduction = "harmony", dims = 1:50)

immune <- FindClusters(immune, resolution = 0.8)

saveRDS(immune, paste0(out_path, 'immune_SeuratObj_res08_v2.rds'))
#immune <- readRDS(paste0(out_path, 'immune_SeuratObj_res08_v2.rds'))

# Checking
DimPlot(immune, reduction = 'umap', label = T) + NoLegend()
DimPlot(immune, reduction = 'umap', split.by = 'status', label = T) + NoLegend()

FeaturePlot(immune, features = c('PTPRC', 'EPCAM'), label = T)
VlnPlot(immune, features = 'percent.mt')

# Finding markers
All_markers <- FindAllMarkers(immune, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(All_markers, paste0(out_path, 'only_immune_markers_v2.csv'), row.names = F)

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
  summarise(markers = str_c(unlist(pick(gene)), collapse=', '))

write.csv(top.markers.clean, paste0(out_path, 'only_immune_markers_top20_v2.csv'), row.names = F)

# Adding identity
clusters.ids <- read_csv(paste0(out_path, 'only_immune_markers_top20_v2_identity.csv')) %>% dplyr::select(cluster, cell.type.3)
clusters.ids$cluster <- as.factor(clusters.ids$cluster)

metadata <- immune@meta.data
metadata <- left_join(metadata, clusters.ids, by = join_by('seurat_clusters' == 'cluster'))
row.names(metadata) <- row.names(immune@meta.data)
immune@meta.data <- metadata

DimPlot(immune, reduction = 'umap', group.by = 'cell.type.3', label = T, repel = T) + NoLegend()
saveRDS(immune, paste0(out_path, 'immune_SeuratObj_res08_identified_v2.rds'))

immune <- readRDS(paste0(out_path, 'immune_SeuratObj_res08_identified_v2.rds'))
immune$status <- as.factor(immune$status)
immune$status <- relevel(immune$status, 'Ctrl')
DimPlot(immune, reduction = 'umap', group.by = 'cell.type.3', split.by = 'status', label = T, repel = T) + NoLegend()



# Saving counts for scCODA
immune@meta.data <- immune@meta.data %>%
  mutate(sample = case_when(
    genotype == '30_219-0607' ~ 'CeD1',
    genotype == '51_219-1056' ~ 'CeD2',
    genotype == '13_CeDNN-A0325' ~ 'Ctrl1',
    genotype == '15_CeDNN-A0326' ~ 'Ctrl2',
    genotype == '20_CeDNN-A0339' ~ 'CeD3',
    genotype == '23_CeDNN-A0357' ~ 'Ctrl4',
    genotype == '24_CeDNN-A0381' ~ 'Ctrl5',
    genotype == '21_CeDNN-A0340' ~  'CeD4'
  ), .after = genotype)

cell_counts <- immune@meta.data %>%
  group_by(sample, age, sex, batch, status, cell.type.3) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  pivot_wider(names_from = cell.type.3, values_from = count)

cell_counts <- replace(cell_counts, is.na(cell_counts), 0)
write.csv(cell_counts, paste0(out_path, 'immune_cells_counts.csv'), row.names = F)

## T cell compartment only -----
tcell <- subset(immune, idents = c(0, 1, 2, 3, 7, 8, 12, 14, 18, 19))

# Re-normalize and cluster
DefaultAssay(tcell) <- 'RNA'

tcell <- SCTransform(tcell, vars.to.regress = "percent.mt", verbose = T)

tcell <- tcell %>%
  RunPCA() %>%
  RunHarmony("batch") %>%
  RunUMAP(reduction = "harmony", dims = 1:50) %>%
  FindNeighbors(reduction = "harmony", dims = 1:50)

tcell <- FindClusters(tcell, resolution = 0.7)

saveRDS(tcell, paste0(out_path, 'tcells_SeuratObj_res07.rds'))

DimPlot(tcell, reduction = 'umap', label = T) + NoLegend()
DimPlot(tcell, reduction = 'umap', split.by = 'status', label = T) + NoLegend()
DimPlot(tcell, reduction = 'umap', group.by = 'status')

# Finding markers
All_markers <- FindAllMarkers(tcell, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(All_markers, paste0(out_path, 'only_tcell_markers.csv'), row.names = F)

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
  summarise(markers = str_c(unlist(pick(gene)), collapse=', '))

write.csv(top.markers.clean, paste0(out_path, 'only_tcell_markers_top20.csv'), row.names = F)

## Inspections
markers <- '../resources/tcell_markers.csv'
PlotUCell(tcell, markers, output = paste0(out_path, 'only_tcells_UCell'))
plotMarkers(tcell, markers, output = paste0(out_path, 'only_tcells'))

FeaturePlot(tcell, features = c('CD4', 'CD8A'), label = T)
FeaturePlot(tcell, features = c('KLRB1', 'DPP4'), label = T)

FeaturePlot(tcell, features = c('E2F2'), label = T)
FeaturePlot(tcell, features = c('CD8A', 'TRDC', 'TRGC2', 'IL2RB'), label = T)

FeaturePlot(tcell, features = c('CD69', 'ITGAE', 'ITGA1', 'CXCR3', 'KLRG1', 'CCR7', 'SELL', 'IL7R', 'CD8A'), label = T)

FeaturePlot(tcell, features = c('GZMK', 'AHNAK', 'KLRG1', 'PLEK', 'A2M-AS1', 'CD44'), label = T)
FeaturePlot(tcell, features = c('EOMES', 'CCL4', 'LYST', 'SAMD3', 'TC2N'), label = T)

# C0 Atlasy
FeaturePlot(tcell, features = c('CD8A', 'ITGAE', 'ENTPD1', 'ID2', 'RBPJ', 'GBP5', 'PRRG3'), label = T)

# C1 Atlasy
FeaturePlot(tcell, features = c('CD8A', 'B2M', 'BCL11B', 'KLRB1', 'TPT1'), label = T)