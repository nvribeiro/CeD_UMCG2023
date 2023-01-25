# Demultiplexing A samples

# Color blind palette
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

library(Seurat)
library(tidyverse)
library(gridExtra)
library(SeuratDisk)

## Demultiplexing without A7 --------------------
# Load demuxlet data
files <- str_subset(list.files(path = './data_preprocessing/demuxlet/', recursive = F, full.names = T), 'A._demuxlet.best')

for (x in files) {
  name <- gsub('./data_preprocessing/demuxlet/','',x)
  
  data <- read.delim(x)
  assign(name, as.data.frame(data))
}

demux_list <- list(A2 = A2_demuxlet.best,
                   A3 = A3_demuxlet.best,
                   A4 = A4_demuxlet.best,
                   A5 = A5_demuxlet.best,
                   A6 = A6_demuxlet.best)

# Updating the barcodes to match the merged dataset
demux_list <- imap(demux_list, ~mutate(.x, BARCODE.UPDATED = case_when(
  .y == 'A2' ~ paste0(BARCODE,'_1'),
  .y == 'A3' ~ paste0(BARCODE,'_2'),
  .y == 'A4' ~ paste0(BARCODE,'_3'),
  .y == 'A5' ~ paste0(BARCODE,'_4'),
  .y == 'A6' ~ paste0(BARCODE,'_5')),
  .after = BARCODE))

# Creating a final merged demuxlet df just with the barcodes and assigned genotype
full_demuxlet <- bind_rows(demux_list) %>% select(BARCODE.UPDATED, BEST, SNG.1ST)

# Load A dataset
A_obj <- LoadH5Seurat('./data_preprocessing/outputs/A/SeuratObj_A_filtered_clustered_woA7.h5Seurat')

# Adding the genotype to the object metadata
tmp <- A_obj@meta.data
tmp$BARCODE.UPDATED <- row.names(tmp)
tmp <- left_join(tmp, full_demuxlet, by="BARCODE.UPDATED")

# Adding some extra information in the table
tmp <- tmp %>% mutate(
  droplet = case_when(
    grepl("SNG", BEST) == TRUE ~ 'singlet',
    grepl("DBL", BEST) == TRUE ~ 'doblet',
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

# Getting the number of singles, doblets and ambiguous
summary <- tmp %>% group_by(droplet) %>% summarise(n()) %>% dplyr::rename(count = 'n()')

png('./outputs/plots/A_droplets.png', width = 600, height = 500)
p <- ggplot(tmp, aes(sample, fill=droplet)) +
  geom_bar(position = 'fill') +
  labs(x = 'lane', y = 'fraction of cells', title = 'Proportion of singlets, doblets and ambiguous droplets in batch A (excl. A7)') +
  theme_classic()
plot(p)
dev.off()

# Checking the genotype and status distribution across the lanes
png('./outputs/plots/A_genotypes_samples.png', width = 700, height = 400)
p1 <- ggplot(tmp, aes(sample, fill=genotype)) +
  geom_bar(position = 'fill') +
  labs(x = 'lane', y = 'fraction of cells', title = 'Genotype distribution in each lane') +
  scale_fill_manual(values=cbPalette) +
  theme_classic()
p2 <- ggplot(filter(tmp, !is.na(genotype)), aes(sample, fill=status)) +
  geom_bar(position = 'fill') +
  labs(x = 'lane', y = 'fraction of cells', title = 'Status distribution in each lane') +
  scale_fill_manual(values=cbPalette) +
  theme_classic()
grid.arrange(p1,p2, ncol=2)
dev.off()

tmp %>% group_by(status) %>% summarise(n())

# Saving this info back to the Seurat object
A_obj@meta.data <- tmp

# Saving the updated SeuratObject
SaveH5Seurat(A_obj, file = './outputs/A/SeuraObj_A_filtered_clustered_woA7_genotyped')

##
## Demultiplexing A7 alone -------------------------------------------------------------
##

# Loading the dataset
A7_obj <- Read10X('./data_preprocessing/data/A7/') %>%
  CreateSeuratObject(project = 'NR03')

A7_obj$percent.mt <- PercentageFeatureSet(A7_obj, pattern = '^MT-')

# Starting demultiplex
A7_demuxlet.best <- A7_demuxlet.best %>% select(BARCODE, BEST, SNG.1ST)

tmp <- A7_obj@meta.data
tmp$BARCODE <- row.names(tmp)
tmp <- left_join(tmp, A7_demuxlet.best, by="BARCODE")

tmp <- tmp %>% mutate(
  droplet = case_when(
    grepl("SNG", BEST) == TRUE ~ 'singlet',
    grepl("DBL", BEST) == TRUE ~ 'doblet',
    grepl("AMB", BEST) == TRUE ~ 'ambiguous')) %>% 
  dplyr::rename(genotype = SNG.1ST)

samples_info <- read.csv('../samples_info/samples_info.csv')
tmp <- left_join(tmp, samples_info, by = c('genotype' = 'sample_ID'))

row.names(tmp) <- tmp$BARCODE

# Replacing 'wrong' genotypes (doblets and ambiguous) with NA
tmp <- tmp %>% 
  mutate(genotype = replace(genotype, droplet != 'singlet', NA))

# Excluding BEST column
tmp <- tmp %>% select(!BEST)

# Getting the number of singles, doblets and ambiguous
summary <- tmp %>% group_by(droplet) %>% summarise(n()) %>% dplyr::rename(count = 'n()')

png('./outputs/plots/A_droplets.png', width = 600, height = 500)
p <- ggplot(tmp, aes(orig.ident, fill=droplet)) +
  geom_bar(position = 'fill') +
  labs(x = 'lane', y = 'fraction of cells', title = 'Proportion of singlets, doblets and ambiguous droplets in A7') +
  coord_flip() +
  theme_classic()
plot(p)
dev.off()

# Checking the genotype and status distribution across the lanes
png('./outputs/plots/A_genotypes_samples.png', width = 700, height = 400)
p1 <- ggplot(tmp, aes(orig.ident, fill=genotype)) +
  geom_bar(position = 'fill') +
  labs(x = 'lane', y = 'fraction of cells', title = 'Genotype distribution in A7') +
  scale_fill_manual(values=cbPalette) +
  theme_classic()
p2 <- ggplot(filter(tmp, !is.na(genotype)), aes(orig.ident, fill=status)) +
  geom_bar(position = 'fill') +
  labs(x = 'lane', y = 'fraction of cells', title = 'Status distribution A7') +
  scale_fill_manual(values=cbPalette) +
  theme_classic()
grid.arrange(p1,p2, ncol=2)
dev.off()

tmp %>% group_by(status) %>% summarise(n())

# Saving this info back to the Seurat object
A7_obj@meta.data <- tmp

# Saving the updated SeuratObject
SaveH5Seurat(A7_obj, file = './data_preprocessing/outputs/A/SeuraObj_A7_genotyped', overwrite = T)

