## Comparing my cell types and annotation with the data from Elmentaite et al 2021 -------------
library(Seurat)
library(tidyverse)
library(patchwork)
library(pals)
library(ggrepel)
library(SeuratDisk)

path <- './data_preprocessing/outputs/all/'
path2 <- 'D:/CeDNN_scRNAseq_2023/elmentaite_2021_data/'
output <- './clusters_refcheck/output/'

elmentaite_adult <- readRDS(paste0(path2, 'elmentaite_2021_adult_smallInt_sct_clus.rds'))
#elmentaite_ped <- readRDS(paste0(path2 , 'elmentaite_2021_pediatric_smallInt_sct_clus.rds'))

my_dataset <- readRDS(paste0(path, 'All_obj_v5_final_indentified.rds'))

# Testing with adult
anchors <- FindTransferAnchors(reference = elmentaite_adult,
                               query = my_dataset,
                               dims = 1:50,
                               normalization.method = 'SCT',
                               reference.reduction = 'pca')

predictions <- TransferData(anchorset = anchors,
                            refdata = elmentaite_adult$cell_type,
                            dims = 1:50)

my_dataset <- AddMetaData(my_dataset, metadata = predictions)

# MapQuery method
my_dataset_2 <- MapQuery(anchorset = anchors,
                       reference = elmentaite_adult,
                       query = my_dataset,
                       refdata = list(celltype = 'cell_type'),
                       reference.reduction = 'pca',
                       reduction.model = 'umap')

p1 <- DimPlot(elmentaite_adult, reduction = 'umap', group.by = 'cell_type', label = T, repel = T) + NoLegend()
p2 <- DimPlot(my_dataset_2, reduction = 'ref.umap', group.by = 'predicted.celltype', label = T, repel = T) + NoLegend()
p1 | p2

# Checking the quality of the querying
metadata <- my_dataset_2@meta.data

# Get the number of cells in each cluster and the median and mean for the score
tmp <- metadata %>%
  select(cell.type.2, predicted.celltype.score, predicted.celltype) %>%
  group_by(cell.type.2, predicted.celltype) %>%
  summarise(count = n(),
            mean.score = mean(predicted.celltype.score),
            median.score = median(predicted.celltype.score),
            sd.score = sd(predicted.celltype.score)) %>%
  ungroup()

# Trying with a corrected score, considering the number of cells in each group and the score
# Corrected score = (N cells in the group / N total cells in the group) * mean score of the group
# High corrected score means a high proportion of cells with a high score
# Low corrected score means high proportion of cells with low score or many groups with similar number of cells and score
tmp2 <- tmp %>%
  group_by(cell.type.2) %>%
  mutate(corrected.score = (count/sum(count)) * mean.score)

# visualizing the results 
p <- ggplot(tmp, aes(x = cell.type.2, y = predicted.celltype, size = count, color = mean.score)) +
  geom_point() +
  geom_text_repel(aes(label = count), size = 3, color = 'black', min.segment.length = 0) +
  labs(title = 'Comparison between my annotation of clusters and annotation in Elmentaite 2021 adult samples',
       x = 'Annotated cell type',
       y = 'Predicted cell type by Elmentaite') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_colour_gradient(low = 'lightgrey', high = 'darkred')

p_corrected <- ggplot(tmp2, aes(x = cell.type.2, y = predicted.celltype, color = corrected.score)) +
  geom_point(size = 5) +
  geom_text_repel(aes(label = count), size = 3, color = 'black', min.segment.length = 0) +
  labs(title = 'Comparison between my annotation of clusters and annotation in Elmentaite 2021 adult samples',
       x = 'Annotated cell type',
       y = 'Predicted cell type by Elmentaite') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_colour_gradient(low = 'lightgrey', high = 'darkred')

pdf(paste0(output, 'comparison_annotation_elmentaite2021.pdf'), width = 11, height = 8, paper = 'a4r')
p
p_corrected
dev.off()
  
write.csv(tmp, paste0(output, 'summary_celltype_prediction_Elmentaite.csv'), row.names = F)

## Comparing only epithelial cells with Elmentaite ---------------------------------------
my_dataset <- epithelia

anchors <- FindTransferAnchors(reference = elmentaite_adult,
                               query = my_dataset,
                               dims = 1:50,
                               normalization.method = 'SCT',
                               reference.reduction = 'pca')

my_dataset <- MapQuery(anchorset = anchors,
                         reference = elmentaite_adult,
                         query = my_dataset,
                         refdata = list(celltype = 'cell_type'),
                         reference.reduction = 'pca',
                         reduction.model = 'umap')

metadata <- my_dataset@meta.data

# Get the number of cells in each cluster and the median and mean for the score
tmp <- metadata %>%
  select(cell.type.3, predicted.celltype.score, predicted.celltype) %>%
  group_by(cell.type.3, predicted.celltype) %>%
  summarise(count = n(),
            mean.score = mean(predicted.celltype.score),
            median.score = median(predicted.celltype.score),
            sd.score = sd(predicted.celltype.score)) %>%
  ungroup() %>%
  group_by(cell.type.3) %>%
  mutate(corrected.score = (count/sum(count)) * mean.score)

p_corrected <- ggplot(tmp, aes(x = cell.type.3, y = predicted.celltype, color = corrected.score)) +
  geom_point(size = 5) +
  geom_text_repel(aes(label = count), size = 3, color = 'black', min.segment.length = 0) +
  labs(title = 'Comparison between my annotation of clusters and annotation in Elmentaite 2021 adult samples',
       x = 'Annotated cell type',
       y = 'Predicted cell type by Elmentaite') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_colour_gradient(low = 'lightgrey', high = 'darkred')

pdf(paste0(output, 'comparison_annotation_elmentaite2021_onlyepithelial.pdf'), width = 11, height = 8, paper = 'a4r')
p_corrected
dev.off()

## Comparing with Burclaff et al 2022 (epithelial cells only) -----------------------------
library(SeuratDisk)
Convert(paste0(path2, 'GSE185224_clustered_annotated_adata_k10_lr0.92_v1.7.h5ad'), dest = "h5seurat", overwrite = T)
burclaff_data <- LoadH5Seurat(paste0(path2, 'GSE185224_clustered_annotated_adata_k10_lr0.92_v1.7.h5Seurat'))
# Not working


