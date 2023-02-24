# DEG analysis in epithelial cells CeD vs Ctrl
s

## Defining functions ------------------------------------------------------------
findDEG <- function(obj, cluster) {
  
  # Subseting SeuratObj
  obj <- subset(obj, idents = cluster)
  Idents(obj) <- "status"
  
  # Finding DEGs
  obj_DEG <- FindMarkers(obj, ident.1 = 'CeD', ident.2 = 'Ctrl', test.use = 'MAST', logfc.threshold = 0)
  obj_DEG$gene <- rownames(obj_DEG)
  
  # Add info about the DE
  obj_DEG$DE <- NA
  obj_DEG <- obj_DEG %>%
    mutate(DE = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0.5 ~ 'up',
      p_val_adj < 0.05 & avg_log2FC < -0.5 ~ 'down'))
  
  # Add column with DEG names
  obj_DEG <- obj_DEG %>%
    mutate(DE_gene = if_else(is.na(DE), NA, gene))
  
  return(obj_DEG)
}

plotVolcano <- function(data) {
  p <- ggplot(data, aes(x = avg_log2FC, y = -log10(p_val_adj), color = DE, label = DE_gene)) +
    geom_jitter() +
    geom_text_repel() +
    theme_classic()
  
  return(p)
}

runGSE <- function(data, mode, FC_cutoff) {
  
  data <- filter(data, p_val_adj < 0.05 & avg_log2FC > FC_cutoff)
  
  geneList <- data$avg_log2FC
  names(geneList) <- as.character(row.names(data))
  geneList <- sort(geneList, decreasing = TRUE)
  
  if (mode == 'GO_BP') {
    GSE <- gseGO(geneList, OrgDb = org.Hs.eg.db, ont = 'BP', keyType = 'SYMBOL')
    return(GSE)
  }
  
  if (mode == 'reactome') {
    entrez_id <- bitr(names(geneList), fromType = "SYMBOL", toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')
    names(geneList) <- entrez_id$ENTREZID
    GSE <- gsePathway(geneList, pvalueCutoff = 0.05)
    return(GSE)
  }
  
  if (mode != 'GO_BP' | mode != 'reactome') {
    print('Choose a mode: GO_BP or reactome.')
  }
  
}

## Initial setup ----------------------------------------------------------------
in_path <- '/groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/NR03_scRNAseq/ongoing/outputs/subclustering/'
out_path <- '/groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/NR03_scRNAseq/ongoing/outputs/geneset_enrichment'

epithelia <- readRDS(paste0(in_path, 'only_epithelial_res07_indentified_v2.rds'))

# Setting the object back to RNA and normalizing the counts
DefaultAssay(epithelia) <- 'RNA'
epithelia <- NormalizeData(epithelia)

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
