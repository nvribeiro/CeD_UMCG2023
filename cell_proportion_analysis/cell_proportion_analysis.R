# Checking differences in cell proportions
library(tidyverse)
library(Seurat)
library(patchwork)
library(speckle)
library(SingleCellExperiment)
library(CellBench)
library(limma)
library(scater)
library(edgeR)
library(statmod)

in_path <- '/groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/NR03_scRNAseq/ongoing/outputs/subclustering/'
out_path <- '/groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/NR03_scRNAseq/ongoing/outputs/cell_proportion_analysis/'

## Propeller for epithelial cells -------------------------------------------------------------
epithelia <- readRDS(paste0(in_path, 'SeuratObj_epithelia_v4_res09_identified.rds'))

# Adding a more user-friendly name for the samples
epithelia@meta.data <- epithelia@meta.data %>%
  mutate(sample = case_when(
    genotype == '30_219-0607' ~ 'CeD1',
    genotype == '51_219-1056' ~ 'CeD2',
    genotype == '13_CeDNN-A0325' ~ 'Ctrl1',
    genotype == '15_CeDNN-A0326' ~ 'Ctrl2',
    genotype == '20_CeDNN-A0339' ~ 'CeD3',
    genotype == '23_CeDNN-A0357' ~ 'Ctrl4',
    genotype == '24_CeDNN-A0381' ~ 'Ctrl5',
    genotype == '21_CeDNN-A0340' ~  'CeD4',
  ))

# Adding sub-classification of TAs
epithelia$cell.type.3[epithelia$seurat_clusters == 1] <- 'TA G1'
epithelia$cell.type.3[epithelia$seurat_clusters == 3] <- 'TA-Ent'
epithelia$cell.type.3[epithelia$seurat_clusters == 8] <- 'TA S'
epithelia$cell.type.3[epithelia$seurat_clusters == 12] <- 'TA G2M'

saveRDS(epithelia, paste0(out_path, 'SeuratObj_epithelia_v4_res09_identified_TAsplit.rds'))
  
epithelia@meta.data <- epithelia@meta.data %>%
  mutate(sample_full = paste(sample, batch, sex, age, sep = '_'))

props <- getTransformedProps(clusters = epithelia$cell.type.3, sample = epithelia$sample_full, transform = 'logit')

plotCellTypeMeanVar(props$Counts)
plotCellTypePropsMeanVar(props$Counts)

names(props)
props$TransformedProps

group_info <- read_csv(paste0(in_path, 'epithelial_cells_counts_v4.csv')) %>% dplyr::select(1:5)
condition <- group_info$status
batch <- group_info$batch
sex <- group_info$sex
age <- group_info$age

design <- model.matrix(~ 0 + condition + batch + sex)

mycontr <- makeContrasts(conditionCtrl-conditionCeD, levels = design)
res <- propeller.ttest(props, design, contrasts = mycontr, robust = T, trend = F, sort = T)
significant <- filter(res, P.Value < 0.05)

# Saving
write.csv(res, paste0(out_path, 'propeller_epithelia_results_v4.csv'))

## Ploting proportions ---------------------------------------------------------
summary <- epithelia@meta.data %>% 
  group_by(status, sample, cell.type.3) %>% 
  summarise(count = n()) %>%
  mutate(total.cells = sum(count)) %>%
  mutate(pct = count/total.cells)

summary$status <- as.factor(summary$status)
summary$status <- relevel(summary$status, 'Ctrl')

library(RColorBrewer)
p <- ggplot(summary, aes(x = cell.type.3, y = pct, fill = status, color = status)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 1.5, position = position_jitterdodge(jitter.width = 0)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_color_manual(values = c('#1f78b4', '#33a02c')) +
  scale_fill_manual(values = c('#a6cee3', '#b2df8a')) +
  labs(title = 'Difference in cell proportions in epithelium',
       y = 'Proportion',
       x = '')

pdf(paste0(out_path, 'epithelia_cell_proportions_v4.pdf'), width = 11, height = 8, paper = 'a4r')
p
dev.off()


## Checking the cell-cycle


## Propeller for immune cells ---------------------------------------------------
immune <- readRDS(paste0(in_path, 'immune_SeuratObj_res08_identified_v2.rds'))

# Adding a more user-friendly name for the samples
immune@meta.data <- immune@meta.data %>%
  mutate(sample = case_when(
    genotype == '30_219-0607' ~ 'CeD1',
    genotype == '51_219-1056' ~ 'CeD2',
    genotype == '13_CeDNN-A0325' ~ 'Ctrl1',
    genotype == '15_CeDNN-A0326' ~ 'Ctrl2',
    genotype == '20_CeDNN-A0339' ~ 'CeD3',
    genotype == '23_CeDNN-A0357' ~ 'Ctrl4',
    genotype == '24_CeDNN-A0381' ~ 'Ctrl5',
    genotype == '21_CeDNN-A0340' ~  'CeD4',
  ))

## Ploting
summary <- immune@meta.data %>% 
  group_by(status, sample, cell.type.3) %>% 
  summarise(count = n()) %>%
  mutate(total.cells = sum(count)) %>%
  mutate(pct = count/total.cells)

summary$status <- as.factor(summary$status)
summary$status <- relevel(summary$status, 'Ctrl')

library(RColorBrewer)
p <- ggplot(summary, aes(x = cell.type.3, y = pct, fill = status, color = status)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 1.5, position = position_jitterdodge(jitter.width = 0)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
  scale_color_manual(values = c('#1f78b4', '#33a02c')) +
  scale_fill_manual(values = c('#a6cee3', '#b2df8a')) +
  labs(title = 'Difference in cell proportions in immune compartment',
       y = 'Proportion',
       x = '')

pdf(paste0(out_path, 'immune_cell_proportions.pdf'), width = 11, height = 8, paper = 'a4r')
p
dev.off()

## Propeller
immune@meta.data <- immune@meta.data %>%
  mutate(sample_full = paste(sample, batch, sex, age, sep = '_'))

props <- getTransformedProps(clusters = immune$cell.type.3, sample = immune$sample_full, transform = 'logit')

plotCellTypeMeanVar(props$Counts)
plotCellTypePropsMeanVar(props$Counts)

names(props)
props$TransformedProps

group_info <- read_csv(paste0(in_path, 'immune_cells_counts.csv')) %>% dplyr::select(1:5)
condition <- group_info$status
batch <- group_info$batch
sex <- group_info$sex
age <- group_info$age

design <- model.matrix(~ 0 + condition + batch + sex + age)

mycontr <- makeContrasts(conditionCtrl-conditionCeD, levels = design)
res <- propeller.ttest(props, design, contrasts = mycontr, robust = T, trend = F, sort = T)

significant <- filter(res, P.Value < 0.05)

# Saving
write.csv(res, paste0(out_path, 'propeller_immune_results.csv'))