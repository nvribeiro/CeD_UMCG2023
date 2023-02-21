# Building full metadata table

library(Seurat)
library(tidyverse)

out_path <- '/groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/NR03_scRNAseq/ongoing/outputs/data_preprocessing/'

A_metadata <- readRDS(paste0(out_path, 'A_merged_metadata.rds'))
F_metadata <- readRDS(paste0(out_path, 'F_merged_metadata.rds'))
final_obj_metadata <- readRDS(paste0(out_path, 'AllCells_SCT_RES04_clustered_final_onlymetadata.rds'))

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
demuxlet <- read_csv(paste0(out_path, 'full_result_demuxlet_allcells.csv'))
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
A_doublets <- readRDS(paste0(out_path, 'A_doublets.rds'))
A_doublets <- bind_rows(A_doublets)
F_doublets <- readRDS(paste0(out_path, 'F_doublets.rds'))
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

# Correcting A7 barcodes before merging
metadata <- metadata %>%
  mutate(barcode = if_else(lane == 'A7', paste0(str_split_fixed(barcode, "_", n = 3)[,1], '_1'), barcode))

# Merge with final obj metadata
final_obj_metadata <- final_obj_metadata %>% 
  select(barcode, SCT_snn_res.0.4, cell.type.1, cell.type.2)

metadata <- left_join(metadata, final_obj_metadata, by = 'barcode')

# Checking if all cells are there
table(final_obj_metadata$barcode %in% metadata$barcode) # all good

# Add filtered out column
metadata <- metadata %>%
  mutate(filtered_out = if_else(is.na(SCT_snn_res.0.4) == TRUE, TRUE, FALSE))

# Saving
write.csv(metadata, paste0(out_path, 'METADATA_NR03_AllCells.csv'))

## Correcting A7 barcodes ------------------------------------------------------
# Missing cells from A7 - or are they duplicated?
tmp <- filter(metadata, lane == 'A7' & inter_sample_doublet == FALSE & intra_sample_doublet == FALSE & droplet_demuxlet == 'singlet' )

# A7 cells in my final dataset have different format of barcode because the way I processed with them
# Trying to convert the barcodes so they match
tmp <- filter(metadata, lane == 'A7')
A7_barcodes_metadata <- tmp$barcode

#final_obj_metadata <- readRDS(paste0(out_path, 'AllCells_SCT_RES04_clustered_final_onlymetadata.rds'))
tmp_A7 <- filter(final_obj_metadata, lane == 'A7')
A7_barcodes_finalObj <- tmp_A7$barcode

# Removing the suffixes from metadata barcodes and getting only the barcode seq
A7_barcodes_metadata_core <- str_split_fixed(A7_barcodes_metadata, "_", n = 3)[,1]
A7_barcodes_finalObj_core <- str_split_fixed(A7_barcodes_finalObj, "_", n = 2)[,1]

# Comparing them - which barcodes in metadata are in my final dataset?
table(A7_barcodes_metadata_core %in% A7_barcodes_finalObj_core) # Number ok, 1645

# Building a dataframe to correct the barcodes in my final obj
correct_barcodes <- tibble(A7_barcodes_finalObj_core)

# Adding the right suffix: _6_1
correct_barcodes <- correct_barcodes %>%
  mutate(corrected_A7_barcode = paste0(A7_barcodes_finalObj_core, '_6_1'))

# Checking these cells, I noticed also that the intra-sample doublets weren't filtered,
# since I used the barcode_suffix to filter them... filtering them out now
A7_check <- filter(metadata, barcode %in% correct_barcodes$corrected_A7_barcode)

# Getting the barcodes of intra-sample doublets
A7_intra_doublets <- filter(A7_check, intra_sample_doublet == TRUE) %>% pull(barcode)

# Converting these barcodes to match in my final obj: suffix _6_1 -> _1
A7_intra_doublets <- paste0(str_split_fixed(A7_intra_doublets, "_", n = 3)[,1], '_1')

# Now I can go back to my obj and filter them out

# Since changing the barcodes in the SeuratObj is more tricky, it's easier to just update them here in the metadata
# Replace _6_1 for _1 in A7 cells
# Running the script from start and changing the barcode


######
