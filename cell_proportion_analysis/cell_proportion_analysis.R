# Checking differences in cell proportions

library(tidyverse)
library(SeuratDisk)
library(patchwork)
library(pals)

path <- './cell_proportion_analysis/output/'
path2 <- './data_preprocessing/outputs/'

## Some visualization and inspection
# Loading data
metadata <- readRDS('./data_preprocessing/outputs/all/All_obj_v5_final_metadata.rds')

# Calculating percentage of cells
summary <- metadata %>%
  group_by(status, cell.type.1, cell.type.2) %>%
  summarise(count.2 = n()) %>%
  ungroup() %>%
  group_by(status, cell.type.1) %>%
  mutate(count.1 = sum(count.2)) %>%
  mutate(pct = count.2/count.1) %>%
  ungroup()

summary1 <- metadata %>% 
  group_by(status, cell.type.1) %>% 
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(status) %>%
  mutate(total.cells = sum(count)) %>%
  mutate(pct = count/total.cells)

summary2 <- metadata %>% 
  group_by(status, cell.type.2) %>% 
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(status) %>%
  mutate(total.cells = sum(count)) %>%
  mutate(pct = count/total.cells)


palette <- as.vector(polychrome(24))
palette2 <- as.vector(tol())

p1 <- ggplot(filter(summary, cell.type.1 != 'Unknown'), aes(x = cell.type.1, y= pct, fill = cell.type.2, label = scales::percent(pct, accuracy = 0.1))) +
  geom_col(position = 'stack') +
  geom_text(position = position_stack(vjust = 0.5),
            size = 4) +
  scale_fill_manual(values = palette, name = 'Cell type') +
  theme_classic() +
  labs(title = 'Distribution of cell types in CeD vs. Controls (%)') +
  xlab('Cell group') +
  ylab('%') +
  facet_wrap(~status)

p2 <- ggplot(summary1, aes(x = status, y = pct, fill = cell.type.1, label = scales::percent(pct, accuracy = 0.1))) +
  geom_col(position = 'stack') +
  geom_text(position = position_stack(vjust = 0.5), size = 4) +
  theme_classic() +
  scale_fill_discrete(name = 'Cell group') +
  labs(title = 'Distribution of cell groups in CeD vs. Controls (%)') +
  xlab('Group') +
  ylab('%')

pdf(paste0(path, 'cell_composition.pdf'), width = 11, height = 8, paper = 'a4r')
p2
p1
dev.off()

## Checking each sample and batches
# Adding a more user-friendly name for the samples
metadata <- metadata %>%
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

metadata$batch <- as.factor(metadata$batch)
metadata$sample <- as.factor(metadata$sample)

p <- ggplot(filter(metadata, cell.type.1 == 'Epithelial'), aes(sample, fill = cell.type.2)) +
  geom_bar(position = 'fill', color = 'white') +
  scale_fill_manual(values = palette, name = 'Cell type') +
  facet_wrap(~status) +
  theme_classic()

summary_batch <- metadata %>%
  group_by(sample, batch, cell.type.1, cell.type.2) %>%
  summarise(count.2 = n()) %>%
  ungroup() %>%
  group_by(sample, batch, cell.type.1) %>%
  mutate(count.1 = sum(count.2)) %>%
  mutate(pct = count.2/count.1) %>%
  ungroup()

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

epithelia <- readRDS(paste0(path, 'subclusters/only_epithelial_res07_identified_immuneremoved.rds'))

prop <- propeller(epithelia, clusters = epithelia$cell.type.3, sample = epithelia$sample, group = epithelia$status, transform = 'asin')

### other way
epithelia@meta.data <- epithelia@meta.data %>%
  mutate(sample_full = paste(sample, batch, sex, age, sep = '_'))

props <- getTransformedProps(clusters = epithelia$cell.type.3, sample = epithelia$sample_full, transform = 'logit')
#props <- getTransformedProps(clusters = epithelia$cell.type.3, sample = epithelia$sample, transform = 'logit')

plotCellTypeMeanVar(props$Counts)
plotCellTypePropsMeanVar(props$Counts)

names(props)
props$TransformedProps

group_info <- read_csv(paste0(path2, 'subclusters/epithelial_cells_counts.csv')) %>% select(1:5)
condition <- group_info$status
batch <- group_info$batch
sex <- group_info$sex
age <- group_info$age

design <- model.matrix(~ 0 + condition + batch + sex + age)

mycontr <- makeContrasts(conditionCtrl-conditionCeD, levels = design)
res <- propeller.ttest(props, design, contrasts = mycontr, robust = T, trend = F, sort = T)

# Saving
write.csv(res, paste0(path, 'propeller_results_complete_model.csv'))

# Ploting
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

pdf(paste0(path, 'epithelium_cell_proportions.pdf'), width = 11, height = 8, paper = 'a4r')
p
dev.off()

  