# Slingshot: trajectory analysis between conditions
library(Seurat)
library(slingshot)
library(SingleCellExperiment)
library(RColorBrewer)
library(scales)
library(viridis)
library(UpSetR)
library(pheatmap)
library(msigdbr)
library(fgsea)
library(gridExtra)
library(tradeSeq)
library(grDevices)
library(scater)

in_path <- '/groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/NR03_scRNAseq/ongoing/outputs/subclustering/'
out_path <- '/groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/NR03_scRNAseq/ongoing/outputs/trajectory_analysis'

epithelia <- readRDS(paste0(in_path, 'only_epithelial_res07_indentified_v2.rds'))

# Trying just with stem cells + TA + enterocytes
enterocytes <- subset(epithelia, idents = c(0, 1, 2, 3, 4, 5, 6, 7, 9))

# Convert to sce
enterocytes <- as.SingleCellExperiment(enterocytes, assay = 'RNA')

# Trajectory inference, starting in stem cells 1 and ending in enterocytes 5
enterocytes <- slingshot(enterocytes, clusterLabels = colData(enterocytes)$seurat_clusters, start.clus = '1',
                         end.clus = '2', reducedDim = 'UMAP')
summary(enterocytes$slingPseudotime_1)

colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(enterocytes$slingPseudotime_1, breaks=100)]

plot(reducedDims(enterocytes)$UMAP, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(enterocytes), lwd=2, type = 'lineages', col='black')

plotReducedDim(enterocytes, dimred = 'UMAP', colour_by = 'cell.type.3')
plotColData(enterocytes, x = 'cell.type.3', y = 'slingPseudotime_1')

# Identifying temporally dynamic genes using tradeSeq
enterocytes <- fitGAM(enterocytes)
