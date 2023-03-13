# DEG analysis in epithelial cells CeD vs Ctrl
library(Seurat)
library(MAST)
library(clusterProfiler)
library(ReactomePA)
library(tidyverse)
library(ggrepel)
library(org.Hs.eg.db)
library(patchwork)
library(enrichplot)

## Defining functions ------------------------------------------------------------
findDEG <- function(obj, cluster) {
  
  # Subseting SeuratObj
  obj <- subset(obj, subset = cell.type.3 == cluster)
  Idents(obj) <- "status"
  
  # Finding DEGs
  obj_DEG <- FindMarkers(obj, ident.1 = 'CeD', ident.2 = 'Ctrl', 
                         test.use = 'MAST', logfc.threshold = 0,
                         min.pct = 0.25)
  obj_DEG$gene <- rownames(obj_DEG)
  obj_DEG$cluster = cluster
  
  # Add info about the DE
  obj_DEG$DE <- NA
  obj_DEG <- obj_DEG %>%
    mutate(DE = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0.25 ~ 'up',
      p_val_adj < 0.05 & avg_log2FC < -0.25 ~ 'down'))
  
  # Add column with DEG names
  obj_DEG <- obj_DEG %>%
    mutate(DE_gene = if_else(is.na(DE), NA, gene))
  
  return(obj_DEG)
}

plotVolcano <- function(data, cluster_name) {
  if (all(is.na(data$DE_gene)) == TRUE) {
    return('No DEG found.')
  } else {
    p <- ggplot(data, aes(x = avg_log2FC, y = -log10(p_val_adj), color = DE, label = DE_gene)) +
      geom_jitter() +
      geom_text_repel() +
      theme_classic() +
      labs(title = paste0('Cluster: ', cluster_name))
    
    return(p)
  }
}

## Initial setup ----------------------------------------------------------------
in_path <- '/groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/NR03_scRNAseq/ongoing/outputs/subclustering/'
out_path <- '/groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/NR03_scRNAseq/ongoing/outputs/DEG_and_GSE_analysis/'

epithelia <- readRDS(paste0(in_path, 'SeuratObj_epithelia_v4_res09_identified.rds'))

# Setting the object back to RNA and normalizing the counts
DefaultAssay(epithelia) <- 'RNA'
epithelia <- NormalizeData(epithelia, normalization.method = 'LogNormalize')

## Finding DEGs for all clusters all at once -----------------------------------
clusters <- unique(epithelia$cell.type.3)

complete_DEG <- map_dfr(clusters, ~findDEG(obj = epithelia, cluster = .x))
complete_DEG$cluster <- as.factor(complete_DEG$cluster)

# Saving
write_csv(complete_DEG, paste0(out_path, 'DEG_epithelial_clusters_v4.csv'))

## Plotting all volcanos
plots <- imap(clusters, ~plotVolcano(data = filter(complete_DEG, cluster == .x), cluster_name = .x))

# Saving
pdf(paste0(out_path, 'VolcanoPlots_epithelial_clusters_v4.pdf'), width = 11, height = 8, paper = 'a4r')
print(plots)
dev.off()

## Gene set enrichment analysis ------------------------------------------------
complete_DEG <- read.csv(paste0(out_path, 'DEG_epithelial_clusters_v4.csv'))

# Adding ENTREZ ID
ENTREZID <- select(org.Hs.eg.db, keys = complete_DEG$gene, 
                   columns = 'ENTREZID', keytype = 'SYMBOL')

complete_DEG <- left_join(complete_DEG, ENTREZID, by = join_by('gene' == 'SYMBOL'))

# Using compareCluster
complete_DEG <- filter(complete_DEG, p_val_adj < 0.05 & is.na(DE) == FALSE)

write.csv(complete_DEG, paste0(out_path, 'DEG_epithelial_clusters_with_ENTREZID_filtered_p005_v4.csv'))

# Reactome
GSE.Pathway <- compareCluster(ENTREZID~cluster+DE, fun = 'enrichPathway',
                              data = complete_DEG, pvalueCutoff = 0.05)

saveRDS(GSE.Pathway, paste0(out_path, 'GSE_reactome_results.rds'))

# GO: Biological process
GSE.BP <- compareCluster(ENTREZID~cluster+DE, fun = 'enrichGO', OrgDb = org.Hs.eg.db,
                              ont = 'BP', data = complete_DEG, pvalueCutoff = 0.05)

saveRDS(GSE.BP, paste0(out_path, 'GSE_GOBP_results.rds'))

# CC: Cellular compartment
GSE.CC <- compareCluster(ENTREZID~cluster+DE, fun = 'enrichGO', OrgDb = org.Hs.eg.db,
                         ont = 'CC', data = complete_DEG, pvalueCutoff = 0.05)

saveRDS(GSE.CC, paste0(out_path, 'GSE_GOCC_results.rds'))

# MF: Molecular function
GSE.MF <- compareCluster(ENTREZID~cluster+DE, fun = 'enrichGO', OrgDb = org.Hs.eg.db,
                         ont = 'MF', data = complete_DEG, pvalueCutoff = 0.05)

saveRDS(GSE.MF, paste0(out_path, 'GSE_GOMF_results.rds'))

# KEGG
GSE.Kegg <- compareCluster(ENTREZID~cluster+DE, fun = 'enrichKEGG',
                           data = complete_DEG, pvalueCutoff = 0.05)

saveRDS(GSE.Kegg, paste0(out_path, 'GSE_KEGG_results.rds'))

## Ploting ----------------------------------------------------
# Getting cluster names to add in the plot
clusters_ids <- epithelia@meta.data %>%
  dplyr::select(seurat_clusters, cell.type.3) %>%
  distinct() %>%
  arrange(seurat_clusters)

# Getting my files
GSE_files <- list.files(path = out_path, pattern = 'GSE_', full.names = T)
GSE_names <- gsub(pattern = "_results.rds", replacement = "", x = basename(GSE_files))
GSE_list <- lapply(GSE_files, readRDS)
names(GSE_list) <- GSE_names

# Defining my plot function
plotGSE <- function(data, output) {
  
  # Calculate gene ratio as numeric
  tmp <- sapply(data$GeneRatio, function(x) eval(str2expression(x)))
  data$GeneRatio_n <- tmp
  
  # Define plot function
  plots <- function(data, cluster_n) {
    tmp2 <- filter(data, cluster == cluster_n)
    
    plot_title <- clusters_ids$cell.type.3[clusters_ids$seurat_clusters == cluster_n]
    
    ggplot(tmp2, aes(x = cluster, y = Description, color = -log10(p.adjust),
                     size = GeneRatio_n)) +
    geom_point() +
    theme_classic() +
    theme(axis.text = element_text(size = 5)) +
    labs(title = plot_title,
         y = '') +
    facet_wrap(~DE, scales = 'free')
  }
  
  # Plot everything
  clusters <- unique(data$cluster)
  plots_list <- map(clusters, ~plots(data, cluster_n = .x))
  
  # Save as pdf
  pdf(output, width = 11, height = 8, paper = 'a4r')
  print(plots_list)
  dev.off()
  
}

# Plot all GSEs
imap(GSE_list, ~plotGSE(data = as.data.frame(.x, col.names = NULL), output = paste0(out_path, 'plots_', .y, '.pdf')))

## GSE with different filters and reducing redudancy ----------------------------
complete_DEG <- filter(complete_DEG, (avg_log2FC > 0.5 | avg_log2FC < 0.5) & p_val_adj < 0.01 & is.na(DE) == FALSE)

# Reactome
GSE.Pathway <- compareCluster(ENTREZID~cluster+DE, fun = 'enrichPathway',
                              data = complete_DEG, pvalueCutoff = 0.05, readable = T)

saveRDS(GSE.Pathway, paste0(out_path, 'GSE_reactome_results_v2.rds'))

# GO: Biological process
GSE.BP <- compareCluster(ENTREZID~cluster+DE, fun = 'enrichGO', OrgDb = org.Hs.eg.db,
                         ont = 'BP', data = complete_DEG, pvalueCutoff = 0.05)
GSE.BP <- pairwise_termsim(GSE.BP)
GSE.BP <- clusterProfiler::simplify(GSE.BP, cutoff=0.7, by="p.adjust", select_fun=min)
saveRDS(GSE.BP, paste0(out_path, 'GSE_GOBP_results_v2.rds'))

# CC: Cellular compartment
GSE.CC <- compareCluster(ENTREZID~cluster+DE, fun = 'enrichGO', OrgDb = org.Hs.eg.db,
                         ont = 'CC', data = complete_DEG, pvalueCutoff = 0.05)
GSE.CC <- pairwise_termsim(GSE.CC)
GSE.CC <- clusterProfiler::simplify(GSE.CC, cutoff=0.7, by="p.adjust", select_fun=min)
saveRDS(GSE.CC, paste0(out_path, 'GSE_GOCC_results_v2.rds'))

# MF: Molecular function
GSE.MF <- compareCluster(ENTREZID~cluster+DE, fun = 'enrichGO', OrgDb = org.Hs.eg.db,
                         ont = 'MF', data = complete_DEG, pvalueCutoff = 0.05)
GSE.MF <- pairwise_termsim(GSE.MF)
GSE.MF <- clusterProfiler::simplify(GSE.MF, cutoff=0.7, by="p.adjust", select_fun=min)
saveRDS(GSE.MF, paste0(out_path, 'GSE_GOMF_results_v2.rds'))

# KEGG
GSE.Kegg <- compareCluster(ENTREZID~cluster+DE, fun = 'enrichKEGG',
                           data = complete_DEG, pvalueCutoff = 0.05)
GSE.Kegg <- pairwise_termsim(GSE.Kegg)
saveRDS(GSE.Kegg, paste0(out_path, 'GSE_KEGG_results_v2.rds'))

GSE_files <- list.files(path = out_path, pattern = 'results_v2', full.names = T)
GSE_names <- gsub(pattern = "_results_v2.rds", replacement = "", x = basename(GSE_files))
GSE_list <- lapply(GSE_files, readRDS)
names(GSE_list) <- GSE_names

# Plot all GSEs
imap(GSE_list, ~plotGSE(data = as.data.frame(.x, col.names = NULL), output = paste0(out_path, 'plots_', .y, '_v2.pdf')))


## Checking stem and TA cells --------------------------------------------------
epithelia <- readRDS(paste0(in_path, 'only_epithelial_res07_indentified_v3.rds'))

# Setting the object back to RNA and normalizing the counts
DefaultAssay(epithelia) <- 'RNA'
epithelia <- NormalizeData(epithelia)

# Comparing Stem 1 and stem 2 in control
control <- subset(epithelia, subset = (seurat_clusters == 1 | seurat_clusters == 4) & status == 'Ctrl')
control_DEG <- FindMarkers(control, ident.1 = 4, ident.2 = 1, test.use = 'MAST', logfc.threshold = 0)

control_DEG$gene <- rownames(control_DEG)

# Add info about the DE
control_DEG$DE <- NA
control_DEG <- control_DEG %>%
  mutate(DE = case_when(
    p_val_adj < 0.05 & avg_log2FC > 0.25 ~ 'up',
    p_val_adj < 0.05 & avg_log2FC < -0.25 ~ 'down'))

# Add column with DEG names
control_DEG <- control_DEG %>%
  mutate(DE_gene = if_else(is.na(DE), NA, gene))

p1 <- plotVolcano(control_DEG, cluster_name = 'Stem cells 2 vs 1 in Control')

# Checking GSE
control_DEG_f <- filter(control_DEG, p_val_adj < 0.05 & is.na(DE) == FALSE)

ENTREZID <- select(org.Hs.eg.db, keys = control_DEG_f$gene, 
                   columns = 'ENTREZID', keytype = 'SYMBOL')
control_DEG_f <- left_join(control_DEG_f, ENTREZID, by = join_by('gene' == 'SYMBOL'))

control.Pathway <- compareCluster(ENTREZID~DE, fun = 'enrichPathway',
                              data = control_DEG_f, pvalueCutoff = 0.05)
p2 <- dotplot(control.Pathway)

control.GOBP <- compareCluster(ENTREZID~DE, fun = 'enrichGO', OrgDb = org.Hs.eg.db,
                               ont = 'BP', data = control_DEG_f, pvalueCutoff = 0.05)
p3 <- dotplot(control.GOBP)

# Comparing Stem in CeD
ced <- subset(epithelia, subset = (seurat_clusters == 1 | seurat_clusters == 4) & status == 'CeD')
ced_DEG <- FindMarkers(ced, ident.1 = 4, ident.2 = 1, test.use = 'MAST', logfc.threshold = 0)

ced_DEG$gene <- rownames(ced_DEG)

# Add info about the DE
ced_DEG$DE <- NA
ced_DEG <- ced_DEG %>%
  mutate(DE = case_when(
    p_val_adj < 0.05 & avg_log2FC > 0.25 ~ 'up',
    p_val_adj < 0.05 & avg_log2FC < -0.25 ~ 'down'))

# Add column with DEG names
ced_DEG <- ced_DEG %>%
  mutate(DE_gene = if_else(is.na(DE), NA, gene))

p4 <- plotVolcano(ced_DEG, cluster_name = 'Stem cells 2 vs 1 in CeD')

# Checking GSE
ced_DEG_f <- filter(ced_DEG, p_val_adj < 0.05 & is.na(DE) == FALSE)

ENTREZID <- select(org.Hs.eg.db, keys = ced_DEG_f$gene, 
                   columns = 'ENTREZID', keytype = 'SYMBOL')
ced_DEG_f <- left_join(ced_DEG_f, ENTREZID, by = join_by('gene' == 'SYMBOL'))

ced.Pathway <- compareCluster(ENTREZID~DE, fun = 'enrichPathway',
                                  data = ced_DEG_f, pvalueCutoff = 0.05)
p5 <- dotplot(ced.Pathway)

ced.GOBP <- compareCluster(ENTREZID~DE, fun = 'enrichGO', OrgDb = org.Hs.eg.db,
                               ont = 'BP', data = ced_DEG_f, pvalueCutoff = 0.05)
p6 <- dotplot(ced.GOBP)

# Saving
pdf(paste0(out_path, 'comparion_stem_cells.pdf'), width = 11, height = 8, paper = 'a4r')
p1
p2
p3
p4
p5
p6
dev.off()

write_csv(control_DEG, file = paste0(out_path, 'DEG_stemcells_control.csv'))
write_csv(ced_DEG, file = paste0(out_path, 'DEG_stemcells_CeD.csv'))

## Doing the same for TA
# Comparing TA 1 and stem 2 in control
control <- subset(epithelia, subset = (seurat_clusters == 7 | seurat_clusters == 9) & status == 'Ctrl')
control_DEG <- FindMarkers(control, ident.1 = 9, ident.2 = 7, test.use = 'MAST', logfc.threshold = 0)

control_DEG$gene <- rownames(control_DEG)

# Add info about the DE
control_DEG$DE <- NA
control_DEG <- control_DEG %>%
  mutate(DE = case_when(
    p_val_adj < 0.05 & avg_log2FC > 0.25 ~ 'up',
    p_val_adj < 0.05 & avg_log2FC < -0.25 ~ 'down'))

# Add column with DEG names
control_DEG <- control_DEG %>%
  mutate(DE_gene = if_else(is.na(DE), NA, gene))

p1 <- plotVolcano(control_DEG, cluster_name = 'TA cells 2 vs 1 in Control')

# Checking GSE
control_DEG_f <- filter(control_DEG, p_val_adj < 0.05 & is.na(DE) == FALSE)

ENTREZID <- select(org.Hs.eg.db, keys = control_DEG_f$gene, 
                   columns = 'ENTREZID', keytype = 'SYMBOL')
control_DEG_f <- left_join(control_DEG_f, ENTREZID, by = join_by('gene' == 'SYMBOL'))

control.Pathway <- compareCluster(ENTREZID~DE, fun = 'enrichPathway',
                                  data = control_DEG_f, pvalueCutoff = 0.05)
p2 <- dotplot(control.Pathway)

control.GOBP <- compareCluster(ENTREZID~DE, fun = 'enrichGO', OrgDb = org.Hs.eg.db,
                               ont = 'BP', data = control_DEG_f, pvalueCutoff = 0.05)
p3 <- dotplot(control.GOBP)

# Comparing TA in CeD
ced <- subset(epithelia, subset = (seurat_clusters == 7 | seurat_clusters == 9) & status == 'CeD')
ced_DEG <- FindMarkers(ced, ident.1 = 9, ident.2 = 7, test.use = 'MAST', logfc.threshold = 0)

ced_DEG$gene <- rownames(ced_DEG)

# Add info about the DE
ced_DEG$DE <- NA
ced_DEG <- ced_DEG %>%
  mutate(DE = case_when(
    p_val_adj < 0.05 & avg_log2FC > 0.25 ~ 'up',
    p_val_adj < 0.05 & avg_log2FC < -0.25 ~ 'down'))

# Add column with DEG names
ced_DEG <- ced_DEG %>%
  mutate(DE_gene = if_else(is.na(DE), NA, gene))

p4 <- plotVolcano(ced_DEG, cluster_name = 'TA cells 2 vs 1 in CeD')

# Checking GSE
ced_DEG_f <- filter(ced_DEG, p_val_adj < 0.05 & is.na(DE) == FALSE)

ENTREZID <- select(org.Hs.eg.db, keys = ced_DEG_f$gene, 
                   columns = 'ENTREZID', keytype = 'SYMBOL')
ced_DEG_f <- left_join(ced_DEG_f, ENTREZID, by = join_by('gene' == 'SYMBOL'))

ced.Pathway <- compareCluster(ENTREZID~DE, fun = 'enrichPathway',
                              data = ced_DEG_f, pvalueCutoff = 0.05)
p5 <- dotplot(ced.Pathway)

ced.GOBP <- compareCluster(ENTREZID~DE, fun = 'enrichGO', OrgDb = org.Hs.eg.db,
                           ont = 'BP', data = ced_DEG_f, pvalueCutoff = 0.05)
p6 <- dotplot(ced.GOBP)

# Saving
pdf(paste0(out_path, 'comparison_TA_cells.pdf'), width = 11, height = 8, paper = 'a4r')
p1
p2
p3
p4
p5
p6
dev.off()

write_csv(control_DEG, file = paste0(out_path, 'DEG_TAcells_control.csv'))
write_csv(ced_DEG, file = paste0(out_path, 'DEG_TAcells_CeD.csv'))

## GSE with msigdb hallmarks -------------------------------------------
library(msigdbr)

# Getting human hallmark dataset
h_gene_set <- msigdbr(species = "Homo sapiens", category = "H")
head(h_gene_set)

msigdbr_t2g <- h_gene_set %>% 
  dplyr::distinct(gs_name, entrez_gene) %>% 
  as.data.frame()

# Testing with cluster 1 up
gse_fun <- function(gene_id = ENTREZID, term = msigdbr_t2g) {
  enricher(gene = gene_id, TERM2GENE = term, pvalueCutoff = 0.05)
}

gse <- compareCluster(ENTREZID~cluster+DE, fun = gse_fun, data = complete_DEG)

pdf(paste0(out_path, 'GSEA_hallmarks_all_clusters_v4.pdf'), width = 11, height = 8, paper = 'a4r')
dotplot(gse, font.size = 8) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()


## Archive -----------
## Stem cells 1 - cluster 1 -----------------------------------------------------
stem1_DEG <- findDEG(epithelia, cluster = 1)
plotVolcano(stem1_DEG)

## Stem cells 2 - cluster 4 -----------------------------------------------------
stem2_DEG <- findDEG(epithelia, cluster = 4)
plotVolcano(stem2_DEG)

stem2_GSE <- runGSE(stem2_DEG, mode = 'GO_BP', FC_cutoff = 0.1)
dotplot(stem2_GSE)

## Enterocytes zone 1 - cluster 3 -----------------------------------------------
ent1_DEG <- findDEG(epithelia, cluster = 3)
plotVolcano(ent1_DEG)

## Enterocytes zone 2 - cluster 0 ---------------
ent2_DEG <- findDEG(epithelia, cluster = 0)
plotVolcano(ent2_DEG)

## Enterocytes zone 3 - cluster 6 ----------------------------------------------
ent3_DEG <- findDEG(epithelia, cluster = 6)
plotVolcano(ent3_DEG)

## Enterocytes zone 4 - cluster 5 ----------------------------------------------
ent4_DEG <- findDEG(epithelia, cluster = 5)
plotVolcano(ent4_DEG)

## Enterocytes zone 5 - cluster 2 ----------------------------------------------
ent5_DEG <- findDEG(epithelia, cluster = 2)
plotVolcano(ent5_DEG)

## TA cells 1 - cluster 7 -----------------------------------------------------
TA1_DEG <- findDEG(epithelia, cluster = 7)
plotVolcano(TA1_DEG)

## TA cells 2 - cluster 9 ----------------------------------------------------
TA2_DEG <- findDEG(epithelia, cluster = 9)
plotVolcano(TA2_DEG)

## Tuft cells - cluster 8 ----------
tuft_DEG <- findDEG(epithelia, cluster = 8)
plotVolcano(tuft_DEG)

## Goblet cells - cluster 10 -------
goblet_DEG <- findDEG(epithelia, cluster = 10)
plotVolcano(goblet_DEG)

## Paneth cells - cluster 11 ----------------
paneth_DEG <- findDEG(epithelia, cluster = 11)
plotVolcano(paneth_DEG)

## Enterochromaffin cells - cluster 12 --------------
ecf_DEG <- findDEG(epithelia, cluster = 12)
plotVolcano(ecf_DEG)

## BEST4 cells - cluster 13 ----------------
best4_DEG <- findDEG(epithelia, cluster = 13)
plotVolcano(best4_DEG)

## EEC A/M - cluster 14 ---------------
eecam_DEG <- findDEG(epithelia, cluster = 14)
plotVolcano(eecam_DEG)

## EEC precursor - cluster 15 -----------------
eecpre_DEG <- findDEG(epithelia, cluster = 15)
plotVolcano(eecpre_DEG)

## EEC D - cluster 16 -----------------------
eecd_DEG <- findDEG(epithelia, cluster = 16)
plotVolcano(eecd_DEG)

## EEC I/L - cluster 17 ----------------------
eecil_DEG <- findDEG(epithelia, cluster = 17)
plotVolcano(eecil_DEG)

## EEC K - cluster 18 ------------------------
eeck_DEG <- findDEG(epithelia, cluster = 18)
plotVolcano(eeck_DEG)

## EEC - cluster 19 ----------------------------
eecx_DEG <- findDEG(epithelia, cluster = 19)
plotVolcano(eecx_DEG)
