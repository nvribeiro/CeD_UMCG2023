# Gene Set Enrichment Analysis with ClusterProfiler
set.seed(1234)
library(tidyverse)
library(Seurat)
library(MAST)
library(NMF)
library(scater)
library(GGally)
library(limma)
library(GSEABase)
library(reshape2)
library(data.table)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(rsvd)
library(RColorBrewer)
library(clusterProfiler)
library(ReactomePA)

in_path <- '/groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/NR03_scRNAseq/ongoing/outputs/subclustering/'
out_path <- '/groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/NR03_scRNAseq/ongoing/outputs/geneset_enrichment'

## Epithelial cells - checking differences in enterocytes, TA, and stem clusters --------------

epithelia <- readRDS(paste0(in_path, 'only_epithelial_res07_indentified_v2.rds'))

## Trying with FindMarkers ------------------------------------------------------------------
# Checking genes in TA
TA_markers <- FindMarkers(epithelia, ident.1 = 7, ident.2 = 9)

DoHeatmap(TA, features = row.names(TA_markers), slot = 'data')

geneList <- TA_markers$avg_log2FC
names(geneList) <- as.character(row.names(TA_markers))
geneList <- sort(geneList, decreasing = TRUE)

GSE <- gseGO(geneList, OrgDb = org.Hs.eg.db, ont = 'BP', keyType = 'SYMBOL')
dotplot(GSE)

entrez_id <- bitr(names(geneList), fromType = "SYMBOL", toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')
names(geneList) <- entrez_id$ENTREZID
GSE2 <- gsePathway(geneList, pvalueCutoff = 0.05)
dotplot(GSE2)

## Trying with MAST --------------------------------------------------------------------------
# TA cells (clusters 7 and 9)
# Getting the cells and covert to SingleCellExperiment
TA <- subset(epithelia, subset = seurat_clusters == 7 | seurat_clusters == 9)
DefaultAssay(TA) <- 'RNA'

TA <- NormalizeData(TA)
x <- FindMarkers(TA, ident.1 = 7, ident.2 = 9, test.use = 'MAST')

TA <- as.SingleCellExperiment(TA)
TA <- logNormCounts(TA)

sca <- SceToSingleCellAssay(TA)

# Tutorial: https://nbisweden.github.io/excelerate-scRNAseq/session-de/session-de-methods.html
# Tutorial 2: https://bioconductor.org/packages/release/bioc/vignettes/MAST/inst/doc/MAITAnalysis.html

# count number of detected genes and scale
cdr2 <- colSums(assay(sca) > 0)
colData(sca)$cngeneson <- scale(cdr2)

# Adding threshold
thres <- thresholdSCRNACountMatrix(assay(sca), nbins = 20, min_per_bin = 30)
assays(sca, withDimnames = FALSE) <- list(thresh=thres$counts_threshold, tpm=assay(sca))

# Fit a hurdle model, modeling the condition and (centered) ngeneson factor
cond <- factor(colData(sca)$cell.type.3)
cond <- relevel(cond, 'TA cells 1')
colData(sca)$condition <- cond

zlmCond <- zlm(~condition + cngeneson, sca)

# Run likelihood ratio test for the condition coefficient.
summaryCond <- summary(zlmCond, doLRT='conditionTA cells 2') 

# Get the datatable with all outputs
summaryDt <- summaryCond$datatable

# Make a datatable with all the relevant output - combined Hurdle model
fcHurdle <- merge(summaryDt[contrast=='conditionTA cells 2' & component=='H',.(primerid, `Pr(>Chisq)`)], 
                  summaryDt[contrast=='conditionTA cells 2' & component=='logFC', 
                            .(primerid, coef, ci.hi, ci.lo)], by='primerid') 

# change name of Pr(>Chisq) to fdr and sort by fdr
fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
fcHurdle <- fcHurdle[order(fdr),]
head(fcHurdle)

# Save output
write_csv(fcHurdle, paste0(out_path, "DEG_TA1vsTA2.csv", row.names = FALSE))

id_to_plot <- fcHurdle[1:50, primerid]
flat_dat <- as(sca[id_to_plot, ], 'data.table')
ggbase <- ggplot(flat_dat, aes(x=condition, y=thresh, color=condition)) + geom_jitter()+facet_wrap(~primerid, scale='free_y')+ggtitle("Top 50 DE Genes in TA cells")
ggbase+geom_violin()

## Following MAIT tutorial -----------------------------------------------------
# Exploratory Data Analysis
aheatmap(as.matrix(assay(sca[1:1000,])), labRow='', annCol = colData(sca)[, 'cell.type.3'])

plotPCA <- function(sca_obj){
  projection <- rpca(t(assay(sca_obj)), retx=TRUE, k=4)$x
  colnames(projection)=c("PC1","PC2","PC3","PC4")
  pca <- data.table(projection,  as.data.frame(colData(sca_obj)))
  print(ggpairs(pca, columns=c('PC1', 'PC2', 'PC3', 'status', 'batch', 'sex'),
                mapping=aes(color=cell.type.3), upper=list(continuous='blank')))
  invisible(pca)
}

plotPCA(sca)

# Adaptive thresholding
scaSample <- sca[sample(which(freq(sca)>.1), 20),]
flat <- as(scaSample, 'data.table')
ggplot(flat, aes(x=logcounts))+geom_density() +facet_wrap(~primerid, scale='free_y')

thres <- thresholdSCRNACountMatrix(assay(sca), nbins = 20, min_per_bin = 30)
par(mfrow=c(3,4))
plot(thres)

dev.off()

assays(sca, withDimnames = FALSE) <- list(thresh=thres$counts_threshold, tpm=assay(sca))

# Differential expression (between TA1 and TA2)
cond <- factor(colData(sca)$cell.type.3)
cond <- relevel(cond, 'TA cells 1')
colData(sca)$cell.type.3 <- cond

# Model considering the influence of the covariates in the cell type
zlmCond <- zlm(~cell.type.3 + status + batch + sex + age, sca)
