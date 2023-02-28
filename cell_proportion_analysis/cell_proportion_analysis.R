# Checking differences in cell proportions

library(tidyverse)
library(SeuratDisk)
library(patchwork)
library(pals)

in_path <- '/groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/NR03_scRNAseq/ongoing/outputs/subclustering/'
out_path <- '/groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/NR03_scRNAseq/ongoing/outputs/cell_proportion_analysis/'

## Propeller for epithelial cells -------------------------------------------------------------
library(speckle)
library(SingleCellExperiment)
library(CellBench)
library(limma)
library(ggplot2)
library(scater)
library(patchwork)
library(edgeR)
library(statmod)

epithelia <- readRDS(paste0(in_path, 'only_epithelial_res07_indentified_v3.rds'))

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

epithelia@meta.data <- epithelia@meta.data %>%
  mutate(sample_full = paste(sample, batch, sex, age, sep = '_'))

props <- getTransformedProps(clusters = epithelia$cell.type.3, sample = epithelia$sample_full, transform = 'logit')

plotCellTypeMeanVar(props$Counts)
plotCellTypePropsMeanVar(props$Counts)

names(props)
props$TransformedProps

group_info <- read_csv(paste0(in_path, 'epithelial_cells_counts_v3.csv')) %>% dplyr::select(1:5)
condition <- group_info$status
batch <- group_info$batch
sex <- group_info$sex
age <- group_info$age

design <- model.matrix(~ 0 + condition + batch + sex + age)

mycontr <- makeContrasts(conditionCtrl-conditionCeD, levels = design)
res <- propeller.ttest(props, design, contrasts = mycontr, robust = T, trend = F, sort = T)

# Saving
write.csv(res, paste0(out_path, 'propeller_epithelia_results.csv'))

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
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
  scale_color_manual(values = c('#1f78b4', '#33a02c')) +
  scale_fill_manual(values = c('#a6cee3', '#b2df8a')) +
  labs(title = 'Difference in cell proportions in epithelium',
       y = 'Proportion',
       x = '')

pdf(paste0(out_path, 'epithelia_cell_proportions.pdf'), width = 11, height = 8, paper = 'a4r')
p
dev.off()