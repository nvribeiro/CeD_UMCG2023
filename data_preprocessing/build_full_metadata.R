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

# Add filter information

# Add doublet information

# Merge everything