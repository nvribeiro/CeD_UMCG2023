# Demultiplexing F samples
setwd('./Desktop/IMIM/Groningen/RPII/Projects/NR03_scRNAseq_main_analysis/data_preprocessing/')

# Color blind palette
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

library(Seurat)
library(tidyverse)
library(patchwork)

# Load demuxlet data
files <- str_subset(list.files(path = './data_preprocessing/demuxlet', recursive = F, full.names = T), 'F._demuxlet.best')

for (x in files) {
  name <- gsub('./data_preprocessing/demuxlet/','',x)
  
  data <- read.delim(x)
  assign(name, as.data.frame(data))
}

demux_list <- list(F1 = F1_demuxlet.best, F2 = F2_demuxlet.best, F3 = F3_demuxlet.best)

# Updating the barcodes to match the merged dataset. Lane F1 = suffix -1, F2 = -2, F3 = -3
demux_list <- imap(demux_list, ~mutate(.x, BARCODE.UPDATED = case_when(
  .y == 'F1' ~ paste0(BARCODE,'_1'),
  .y == 'F2' ~ paste0(BARCODE,'_2'),
  .y == 'F3' ~ paste0(BARCODE,'_3')), .after = BARCODE))

# Finding singlets (SNG), doublets (DBL), ambigouos (AMB)
SNG.list <- map(demux_list, ~filter(.x, grepl("SNG", BEST)))
DBL.list <- map(demux_list, ~filter(.x, grepl("DBL", BEST)))
AMB.list <- map(demux_list, ~filter(.x, grepl("AMB", BEST)))

# Creating a final merged demuxlet df just with the barcodes and assigned genotype
full_demuxlet <- bind_rows(demux_list) %>% select(BARCODE.UPDATED, BEST, SNG.1ST)

# Load F dataset
F_filtered <- readRDS('./data_preprocessing/outputs/F/SeuratObj_F_filtered.rds')

# Adding the genotype to the object metadata
tmp <- F_filtered@meta.data
tmp$BARCODE.UPDATED <- row.names(tmp)
tmp <- left_join(tmp, full_demuxlet, by="BARCODE.UPDATED")

# Adding some extra information in the table
tmp <- tmp %>% mutate(
  droplet = case_when(
    grepl("SNG", BEST) == TRUE ~ 'singlet',
    grepl("DBL", BEST) == TRUE ~ 'doublet',
    grepl("AMB", BEST) == TRUE ~ 'ambiguous')) %>%
  dplyr::rename(genotype = SNG.1ST)

samples_info <- read.csv('../samples_info/samples_info.csv')
tmp <- left_join(tmp, samples_info, by = c('genotype' = 'sample_ID'))

row.names(tmp) <- tmp$BARCODE.UPDATED

# Replacing 'wrong' genotypes (doblets and ambiguous) with NA
tmp <- tmp %>% 
  mutate(genotype = replace(genotype, droplet != 'singlet', NA))

# Excluding BEST column
tmp <- tmp %>% select(!BEST)

# Saving this info back to the Seurat object
F_filtered@meta.data <- tmp

# Saving the updated SeuratObject
saveRDS(F_filtered, file = './outputs/F/SeuratObj_F_filtered_genotyped.rds')

# Getting the number of singles, doblets and ambiguous
summary <- tmp %>% group_by(droplet) %>% summarise(n()) %>% rename(count = 'n()')

p1 <- ggplot(tmp, aes(sample, fill=droplet)) +
  geom_bar(position = 'fill') +
  labs(x = 'lane', y = 'fraction of cells', title = 'Proportion of singlets, doblets and ambiguous droplets in batch F') +
  theme_classic()

p2 <- ggplot(tmp, aes(sample, fill=genotype)) +
  geom_bar(position = 'fill') +
  labs(x = 'lane', y = 'fraction of cells', title = 'Genotype distribution in each lane') +
  scale_fill_manual(values=cbPalette) +
  theme_classic()

p3 <- ggplot(filter(tmp, !is.na(genotype)), aes(sample, fill=status)) +
  geom_bar(position = 'fill') +
  labs(x = 'lane', y = 'fraction of cells', title = 'Status distribution in each lane') +
  scale_fill_manual(values=cbPalette) +
  theme_classic()

pdf('./data_preprocessing/outputs/plots/F_demultiplexing.pdf', width = 11, height = 8, paper = 'a4r')
p1
p2 | p3
dev.off()

tmp %>% group_by(status) %>% summarise(n())

## Finding more doublets using scDblFinder -------------------------------------
library(scDblFinder)
tmp <- tmp %>%
  mutate(doublet = if_else(droplet == 'doublet', TRUE, FALSE))
F_filtered@meta.data <- tmp

saveRDS(F_filtered, './data_preprocessing/outputs/F/F_for_doublets.rds')