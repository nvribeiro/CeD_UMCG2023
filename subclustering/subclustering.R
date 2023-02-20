set.seed(1234)
library(Seurat)
library(tidyverse)
library(harmony)
source('./functions/my_functions.R')

path <- './subclustering/outputs/'

## Subclustering epithelial cells ----------------------------------------------
# Subseting
obj <- readRDS('./data_preprocessing/outputs/all/All_SeuratObj_v5_final_identified_unknownremoved.rds')
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

saveRDS(epithelia, paste0(path, 'subclusters/epithelial_SeuratObj_res07.rds'))
epithelia <- readRDS(paste0(path, 'subclusters/epithelial_SeuratObj_res07.rds'))

pdf(paste0(path, 'subclusters/only_epithelial_res07_umap.pdf'), width = 11, height = 8, paper = 'a4r')
DimPlot(epithelia, reduction = 'umap', label = T) + NoLegend()
dev.off()

markers <- '../../resources/epithelial_markers_detailed.csv'
PlotUCell(epithelia, markers, output = paste0(path, 'subclusters/only_epithelial_reclustered'))
plotMarkers(epithelia, markers, output = paste0(path, 'subclusters/only_epithelial_reclustered'))

# Finding markers
All_markers <- FindAllMarkers(epithelia, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(All_markers, paste0(path, 'subclusters/only_epithelial_markers.csv'), row.names = F)

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

write.csv(top.markers.clean, paste0(path, 'subclusters/only_epithelial_markers_top20.csv'), row.names = F)

# Inspections
FeaturePlot(epithelia, features = c('EPCAM', 'PTPRC'))
# Enterocytes compartment
FeaturePlot(epithelia, features = c('REG1A', 'SLC2A2', 'ADA'), label = T)
FeaturePlot(epithelia, features = c('GCG', 'PYY'), label = T)
VlnPlot(epithelia, features = c('GCG', 'PYY'))

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
enterocytes <- subset(epithelia, subset = seurat_clusters == 0 | seurat_clusters == 1 | seurat_clusters == 3 | seurat_clusters == 4 | seurat_clusters == 6)

# Getting the marker genes from Moor et al 2018
library(readxl)
enterocytes_markers <- read_excel('../../resources/gene_markers_enterocytes_Moor2018.xlsx')
mouse_genes <- c(enterocytes_markers$`Cluster 1`, enterocytes_markers$`Cluster 2`, enterocytes_markers$`Cluster 3`, enterocytes_markers$`Cluster 4`, enterocytes_markers$`Cluster 5`)
mouse_genes <- unique(mouse_genes)

write_clip(cat(paste(mouse_genes, "\n")))

# Convert these genes to the human equivalent
human_genes <- read.csv('../../resources/human_orthologs_moor2018.txt') %>% pull(Human.gene.name)
human_genes_bottom <- read.csv('../../resources/human_orthologs_moor2018_bottom_landmark.txt') %>% pull(Human.gene.name)
human_genes_top <- read.csv('../../resources/human_orthologs_moor2018_top_landmark.txt') %>% pull(Human.gene.name)
burclaff_markers <- c('CA2', 'AKR1B10', 'SLC4A4', 'PCK1', 'PFKFB2', 'APOB', 'APOA4', 'SI', 'FABP2', 'ANPEP', 'PRAP1', 'FANP6', 'TMIGD1', 'PBL1', 'CUBN', 'SLC10A2', 'CLCA4', 'SLC13A1', 'FAM151A', 'SLC13A1', 'SLC17A8')
enterocytes_avg <- AverageExpression(enterocytes, assay = 'SCT', slot = 'scale.data', return.seurat =  T)

pdf(paste0(path, 'subclusters/enterocytes_markers_topbottom.pdf'), height = 11, width = 8, paper = 'a4')
DoHeatmap(enterocytes, features = human_genes_bottom, raster = FALSE) + labs(title = 'Bottom landmark genes')
DoHeatmap(enterocytes, features = human_genes_top, raster = FALSE) + labs(title = 'Top landmark genes')
DoHeatmap(enterocytes_avg, features = human_genes_bottom, raster = FALSE, draw.lines = F) + labs(title = 'Bottom landmark genes - avg cluster expression')
DoHeatmap(enterocytes_avg, features = human_genes_top, raster = FALSE, draw.lines = F) + labs(title = 'Top landmark genes - avg cluster expression')
DoHeatmap(enterocytes, features = burclaff_markers, raster = FALSE) + labs(title = 'Burclaff marker genes')
DoHeatmap(enterocytes_avg, features = burclaff_markers, raster = FALSE, draw.lines = F) + labs(title = 'Burclaff marker genes - avg cluster expression')
dev.off()

# Plotting my markers from FindAllMarkers
all_markers <- read_csv('data_preprocessing/outputs/subclusters/only_epithelial_markers.csv')
enterocytes_markers <- all_markers %>%
  filter(cluster == 0 | cluster == 1 | cluster == 2 | cluster == 3 | cluster == 4 | cluster == 6) %>%
  distinct(gene) %>%
  pull(gene)

# Adding anotation to metadata
clusters.ids <- read.csv('data_preprocessing/outputs/subclusters/only_epithelial_markers_top20.csv')
clusters.ids <- clusters.ids %>% select(-markers)
clusters.ids$cluster <- as.factor(clusters.ids$cluster)

metadata <- epithelia@meta.data
metadata <- left_join(metadata, clusters.ids, by = c('seurat_clusters' = 'cluster'))
row.names(metadata) <- metadata$barcode

epithelia@meta.data <- metadata

saveRDS(epithelia, paste0(path, 'subclusters/only_epithelial_res07_indentified.rds'))
epithelia <- readRDS(paste0(path, 'subclusters/only_epithelial_res07_indentified.rds'))

pdf(paste0(path, 'subclusters/only_epithelial_res07_umap_identified.pdf'), width = 11, height = 8, paper = 'a4r')
DimPlot(epithelia, reduction = 'umap', label = T, repel = T, group.by = 'cell.type.3') + NoLegend()
dev.off()

# Check unknown cluster - can't figure out
FeaturePlot(epithelia, features = 'EPCAM')
FeaturePlot(epithelia, features = 'PTPRC')
FeaturePlot(epithelia, features = c('GJA4', 'IGFBP3', 'UNC5B')) # arterial markers
FeaturePlot(epithelia, features = c('THY1', 'COL1A2', 'VIM')) # fibroblasts
FeaturePlot(epithelia, features = c('PECAM1', 'CDH5', 'CLDN5')) # endothelial
FeaturePlot(epithelia, features = c('KCNJ8', 'ABCC9', 'RGS5', 'CSP4G'))
FeaturePlot(epithelia, features = c('MYH11', 'ACTG2', 'TAGLN'))

# Removing the immune cluster
epithelia <- subset(epithelia, subset = seurat_clusters != 13)
saveRDS(epithelia, 'data_preprocessing/outputs/subclusters/only_epithelial_res07_identified_immuneremoved.rds')

p1 <- DimPlot(epithelia, reduction = 'umap', group.by = 'cell.type.3', label = T, repel = T) + NoLegend()
p2 <- DimPlot(epithelia, reduction = 'umap', group.by = 'cell.type.3', split.by = 'status') + NoLegend()
p3 <- DimPlot(epithelia, reduction = 'umap', group.by = 'cell.type.3', split.by = 'sample') + NoLegend()

pdf(paste0(path, 'subclusters/only_epithelial_res07_umaps_immuneremoved.pdf'), width = 11, height = 8, paper = 'a4r')
p1
p2
p3
dev.off()

# Saving table for scDODA
cell_counts <- epithelia@meta.data %>%
  group_by(sample, age, sex, batch, status, cell.type.3) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  pivot_wider(names_from = cell.type.3, values_from = count)

cell_counts <- replace(cell_counts, is.na(cell_counts), 0)
write.csv(cell_counts, paste0(path, 'subclusters/epithelial_cells_counts.csv'), row.names = F)

# Saving the immune cluster to append in the immune group
extra_immune <- subset(epithelia, subset = seurat_clusters == 13)
saveRDS(extra_immune, paste0(path, 'subclusters/extra_immune_cluster.rds'))

# Some more enterocytes inspection
epithelia <- readRDS(paste0(path, 'subclusters/only_epithelial_res07_identified_immuneremoved.rds'))
enterocytes <- subset(epithelia, subset = seurat_clusters == 0 | seurat_clusters == 1 | seurat_clusters == 3 | seurat_clusters == 4 | seurat_clusters == 6)
all_markers <- read_csv('data_preprocessing/outputs/subclusters/only_epithelial_markers.csv')
enterocytes_markers <- all_markers %>%
  filter(cluster == 0 | cluster == 1 | cluster == 2 | cluster == 3 | cluster == 4 | cluster == 6) %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)

p1 <- DoHeatmap(enterocytes, features = enterocytes_markers$gene, raster = FALSE)

enterocytes_avg <- AverageExpression(enterocytes, assay = 'SCT', slot = 'scale.data', return.seurat =  T)
p2 <- DoHeatmap(enterocytes_avg, features = enterocytes_markers$gene, raster = FALSE, draw.lines = F) + labs(title = 'Avg expression enterocytes clusters')

FeaturePlot(enterocytes, features = burclaff_markers)

pdf('data_preprocessing/outputs/subclusters/enterocytes_cluster_markers.pdf', width = 11, height = 8, paper = 'a4r')
p1
p2
dev.off()


## Subclustering epithelial cells v2 - adding the cells that were in the immune object -----------------------
epithelia <- readRDS(paste0(path, 'only_epithelial_res07_identified_immuneremoved.rds'))
extra_epithelia <- subset(obj, cells = extra_epithelial)

new_epithelia <- merge(epithelia, extra_epithelia)

# Reclustering and repeating everything
DefaultAssay(new_epithelia) <- 'RNA'

# Normalizing
new_epithelia <- SCTransform(new_epithelia, vars.to.regress = "percent.mt", verbose = T)

# Reclustering
new_epithelia <- new_epithelia %>%
  RunPCA() %>%
  RunHarmony("batch") %>%
  RunUMAP(reduction = "harmony", dims = 1:50) %>%
  FindNeighbors(reduction = "harmony", dims = 1:50)

new_epithelia <- FindClusters(new_epithelia, resolution = 0.7)

DimPlot(new_epithelia, reduction = 'umap', label = T) + NoLegend()

markers <- '../../resources/epithelial_markers_detailed.csv'
PlotUCell(new_epithelia, markers, output = paste0(path, 'only_epithelial_reclustered_v2'))

# Finding markers
All_markers <- FindAllMarkers(new_epithelia, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(All_markers, paste0(path, 'only_epithelial_markers_v2.csv'), row.names = F)

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

write.csv(top.markers.clean, paste0(path, 'only_epithelial_markers_top20_v2.csv'), row.names = F)

# Now I have the immune cluster back...
FeaturePlot(new_epithelia, features = c('EPCAM', 'PTPRC'))
extra_immune <- WhichCells(new_epithelia, idents = 12)

# Checking if 'extra immune cells' in epithelia dataset and 'extra epithelia cells' in immune dataset are the same
# Could be doublets messing things up
table(extra_epithelial %in% extra_immune)
# 194 cells are the same (TRUE) and they express EPCAM and PTPRC, probably doublets
# Visualizing the FALSE and TRUE cells
barcodes_true <- extra_epithelial[extra_epithelial %in% extra_immune == TRUE]
barcodes_false <- extra_epithelial[extra_epithelial %in% extra_immune == FALSE]

DimPlot(new_epithelia, reduction = 'umap', cells.highlight = barcodes_true)
DimPlot(new_immune, reduction = 'umap', cells.highlight = barcodes_false)

# Where these cells were in my original dataset?
DimPlot(obj, reduction = 'umap', cells.highlight = extra_epithelial)
DimPlot(obj, reduction = 'umap', cells.highlight = extra_immune)
# Comparing to where I had doublets, these cells are in clusters where I found most of the doublets
# Removing this cluster but saving these barcodes because they can be useful for other analysis

barcodes <- unique(c(extra_epithelial, extra_immune))
saveRDS(barcodes, 'barcodes_doubletes_immune_epithelial.rds')

## Subclustering immune cells --------------------------------------------------
# Get barcodes of immune group
extra_immune_barcodes <- WhichCells(epithelia, idents = 13)

immune <- subset(obj, subset = cell.type.1 == 'Immune' | barcode %in% extra_immune_barcodes)

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

saveRDS(immune, paste0(path, 'only_immune_res08.rds'))
immune <- readRDS(paste0(path, 'only_immune_res08.rds'))

DimPlot(immune, reduction = 'umap', label = T) + NoLegend()
FeaturePlot(immune, features = c('PTPRC', 'EPCAM'))

markers <- '../../resources/immune_markers_detailed.csv'
PlotUCell(immune, markers, output = paste0(path, 'immune_res08_Ucell'))

# Finding markers
All_markers <- FindAllMarkers(immune, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(All_markers, paste0(path, 'only_immune_markers.csv'), row.names = F)

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

write.csv(top.markers.clean, paste0(path, 'only_immune_markers_top20.csv'), row.names = F)

# Cluster 13 is high in markers for enterocytes - removing
extra_epithelial <- WhichCells(immune, idents = 13)
new_immune <- subset(immune, subset = seurat_clusters != 13)
