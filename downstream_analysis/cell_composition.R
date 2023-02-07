# Starting some downstream processing with all samples

library(tidyverse)
library(SeuratDisk)
library(patchwork)
library(pals)

path <- './downstream_analysis/output/'

# Loading data
metadata <- read.csv('./data_preprocessing/outputs/all/All_samples_v3_metadata.csv')

# Changing a few things to better visualization
tmp <- read.csv('./data_preprocessing/outputs/all/All_clusters_top10markers_clustersrenamed.csv')
tmp <- tmp %>% select(cluster, cell.type.1, cell.type.2)
tmp$cluster <- as.factor(tmp$cluster)
metadata$seurat_clusters <- as.factor(metadata$seurat_clusters)
metadata <- left_join(metadata, tmp, by = c('seurat_clusters' = 'cluster'))

metadata <- metadata %>%
  dplyr::rename(cell.type.1 = cell.type.1.y) %>%
  dplyr::rename(cell.type.2 = cell.type.2.y) %>% 
  select(seurat_clusters:cell.type.2, -cell.type.1.x, -cell.type.2.x)


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


palette <- as.vector(polychrome(19))
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
    genotype == '30_219-0607' ~ 'CeD_1',
    genotype == '51_219-1056' ~ 'CeD_2',
    genotype == '13_CeDNN-A0325' ~ 'Ctrl_1',
    genotype == '15_CeDNN-A0326' ~ 'Ctrl_2',
    genotype == '20_CeDNN-A0339' ~ 'CeD_3',
    genotype == '23_CeDNN-A0357' ~ 'Ctrl_3',
    genotype == '24_CeDNN-A0381' ~ 'Ctrl_4',
    genotype == '21_CeDNN-A0340' ~  'CeD_4',
  ))

metadata$batch <- as.factor(metadata$batch)
metadata$sample <- as.factor(metadata$sample)

p <- ggplot(filter(metadata, cell.type.1 != 'Unknown'), aes(sample, fill = cell.type.2)) +
  geom_bar(position = 'fill') +
  scale_fill_manual(values = palette, name = 'Cell type') +
  facet_wrap(~cell.type.1) +
  theme_classic()

summary_batch <- metadata %>%
  group_by(sample, batch, cell.type.1, cell.type.2) %>%
  summarise(count.2 = n()) %>%
  ungroup() %>%
  group_by(sample, batch, cell.type.1) %>%
  mutate(count.1 = sum(count.2)) %>%
  mutate(pct = count.2/count.1) %>%
  ungroup()
