# Building DEG table
library(tidyverse)

input <- 'D:/CeDNN_scRNAseq_2023/backup_cluster_ongoing/outputs/DEG_and_GSE_analysis/DEG_epithelial_clusters.csv'
out_path <- 'DEG_and_GSE_analysis/'

# Loading the data
deg <- read_csv(input)

# Adding cluster info
cluster.ids <- read.csv('D:/CeDNN_scRNAseq_2023/backup_cluster_ongoing/outputs/subclustering/only_epithelial_top20markers_v2_identity.csv')
cluster.ids <- cluster.ids %>% dplyr::select(-markers)

deg <- left_join(deg, cluster.ids, by = 'cluster')

# Cleaning the table a little bit
deg <- deg %>%
  filter(p_val_adj < 0.05 & (avg_log2FC > 0.25 | avg_log2FC < -0.25)) %>%
  dplyr::select(gene, avg_log2FC, p_val_adj, DE, cluster, cell.type.3) %>%
  dplyr::rename(direction = DE)

# Getting info about the genes
library(biomaRt)
ensembl <- useEnsembl(biomart = 'ensembl', datase = 'hsapiens_gene_ensembl')
gene_list <- unique(deg$gene)
attributes <- c('hgnc_symbol', 'ensembl_gene_id', 'description', 'chromosome_name', 'start_position',
                'end_position', 'go_id')

gene_info <- getBM(attributes = attributes, filters = 'hgnc_symbol',
                   values = gene_list, mart = ensembl)

# Filter duplicates for chromosomes
gene_info <- filter(gene_info, grepl('CHR_', gene_info$chromosome_name) == FALSE)

# Getting GO information
library(GO.db)
keys <- unique(gene_info$go_id)
GO_info <- select(GO.db, keys = keys, keytype = 'GOID', columns = c('TERM', 'ONTOLOGY'))

# Filtering just biological process (BP)
GO_info <- filter(GO_info, ONTOLOGY == 'BP')

# Now filtering the gene_info to keep only BP and adding the TERM column there
gene_info <- filter(gene_info, go_id %in% GO_info$GOID)
gene_info <- left_join(gene_info, GO_info, by = join_by('go_id' == 'GOID'))
gene_info <- gene_info %>%
  group_by(hgnc_symbol) %>%
  mutate(GO.BP = str_c(unlist(pick(TERM)), collapse=' | ')) %>%
  dplyr::select(-go_id, -TERM, -ONTOLOGY) %>%
  distinct()

# Now adding this info to my DEG dataframe
deg <- left_join(deg, gene_info, by = join_by('gene' == 'hgnc_symbol'), multiple = 'all')

# Saving
write_csv(deg, paste0(out_path, 'DEG_epithelia_with_genedetails.csv'))
