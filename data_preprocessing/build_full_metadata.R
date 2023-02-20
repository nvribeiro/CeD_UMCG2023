# Building full metadata table

library(Seurat)
library(SeuratDisk)
library(tidyverse)

path_hdd <- 'D:/CeDNN_scRNAseq_2023/data_preprocessing_archive/'

A_metadata <- readRDS('data_preprocessing/outputs/A_merged_metadata.rds')
F_metadata <- readRDS('data_preprocessing/outputs/F_merged_metadata.rds')
final_obj_metadata <- readRDS('data_preprocessing/outputs/all/All_SeuratObj_v5_final_metadata_unknownremoved.rds')

# Add column with barcode
A_metadata$barcode <- row.names(A_metadata)
F_metadata$barcode <- row.names(F_metadata)

# Set the barcodes to the same format as they are in the final dataset
# End _1 -> A, _2 _> F
A_metadata <- A_metadata %>%
  mutate(barcode = paste0(barcode, '_1'))

F_metadata <- F_metadata %>%
  mutate(barcode = paste0(barcode, '_2'))

metadata <- bind_rows(A_metadata, F_metadata)

# Add demuxlet information
demuxlet <- read_csv(paste0(path_hdd, 'all/full_result_demuxlet_allcells.csv'))
metadata <- left_join(metadata, demuxlet, by = join_by(barcode == BARCODE.UPDATED))
metadata <- metadata %>%
  select(barcode, sample, BEST, SNG.1ST) %>%
  dplyr::rename(lane = sample) %>%
  mutate(
    droplet = case_when(
      grepl("SNG", BEST) == TRUE ~ 'singlet',
      grepl("DBL", BEST) == TRUE ~ 'doublet',
      grepl("AMB", BEST) == TRUE ~ 'ambiguous')) %>% 
  dplyr::rename(genotype = SNG.1ST) %>%
  mutate(genotype = replace(genotype, droplet != 'singlet', NA)) %>%
  select(!BEST)

metadata <- dplyr::rename(metadata, droplet_demuxlet = droplet)

# Add intra-sample doublet information
A_doublets <- readRDS(paste0(path_hdd, 'A/A_doublets.rds'))
A_doublets <- bind_rows(A_doublets)
F_doublets <- readRDS(paste0(path_hdd, 'F/F_doublets.rds'))
F_doublets <- bind_rows(F_doublets)

F_doublets <- F_doublets %>%
  filter(predicted == TRUE) %>%
  mutate(barcode = paste0(barcode, '_2')) %>%
  pull(barcode)

A_doublets <- A_doublets %>%
  filter(predicted == TRUE) %>%
  mutate(barcode = paste0(barcode, '_1')) %>%
  pull(barcode)

intra_doublets <- c(A_doublets, F_doublets)

metadata <- metadata %>%
  mutate(inter_sample_doublet = if_else(droplet_demuxlet == 'doublet', TRUE, FALSE),
         intra_sample_doublet = if_else(barcode %in% intra_doublets == TRUE, TRUE, FALSE))

# Merge with final obj metadata
final_obj_metadata <- final_obj_metadata %>% 
  select(barcode, SCT_snn_res.0.4, cell.type.1, cell.type.2)

metadata <- left_join(metadata, final_obj_metadata, by = 'barcode')

# Add filtered out column
tmp <- metadata %>%
  mutate(filtered_out = if_else(is.na(SCT_snn_res.0.4) == TRUE, TRUE, FALSE))
