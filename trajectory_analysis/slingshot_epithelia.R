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
library(tidyverse)

in_path <- '/groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/NR03_scRNAseq/ongoing/outputs/subclustering/'
out_path <- '/groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/NR03_scRNAseq/ongoing/outputs/trajectory_analysis/'

epithelia <- readRDS(paste0(in_path, 'SeuratObj_epithelia_v4_res09_identified.rds'))
#epithelia$cell.type.4 <- epithelia$cell.type.3
#epithelia$cell.type.4[epithelia$seurat_clusters %in% c(12, 14, 15, 16, 17, 18, 19)] <- 'EEC'

# Convert to sce
sce <- as.SingleCellExperiment(epithelia)
plotUMAP(control, colour_by = 'cell.type.3')

# Checkin imbalance score
scores <- bioc2020trajectories::imbalance_score(
  rd = reducedDims(sce)$UMAP, 
  cl = colData(sce)$status,
  k = 20, smooth = 40)

grad <- viridis::plasma(10, begin = 0, end = 1)
names(grad) <- levels(cut(scores$scaled_scores, breaks = 10))

pdf(paste0(out_path, 'imbalance_score_epithelia.pdf'), width = 11, height = 8, paper = 'a4r')
plot(reducedDims(sce)$UMAP, col = grad[cut(scores$scaled_scores, breaks = 10)],
     asp = 1, pch = 16, xlab = "UMAP-1", ylab = "UMAP-2", cex = .8)
legend("topleft", legend = names(grad), col = grad, pch = 16, bty = "n", cex = 2 / 3)
dev.off()

# Spliting my dataset into control and CeD
control <- subset(sce, , status == 'Ctrl')
ced <- subset(sce, , status == 'CeD')

# Trajectory inference - Ctrl
control <- slingshot(control, clusterLabels = colData(control)$cell.type.3,
                     reducedDim = 'UMAP', start.clus = 'Stem cells')

sling <- SlingshotDataSet(control)

# Check lineages
lineages_control <- slingLineages(control)
pseudo.paths <- slingPseudotime(control)
head(pseudo.paths)

# Check trajectory IDs
curve.assignments <- slingBranchID(control)
table(curve.assignments)

# Get share pseudotime
shared.pseudo <- rowMeans(pseudo.paths, na.rm=TRUE)

# Plotting the lineages in the UMAP
# Need to loop over the paths and add each one separately.
gg_control <- plotUMAP(control, colour_by=I(shared.pseudo))
embedded <- embedCurves(control, "UMAP")
embedded <- slingCurves(embedded)

# Setting a different color for each trajectory
colors <- brewer.pal(6, 'Set1')

i = 1
for (path in embedded) {
  embedded <- data.frame(path$s[path$ord,])
  gg_control <- gg_control + geom_path(data=embedded, aes(x=UMAP_1, y=UMAP_2), size=1.2, color = colors[i])
  i = i+1
}
gg_control <- gg_control + labs(title = 'Pseudotime and trajectories - Control')
gg_control

# Ploting pseudotimes without the curves
pseudo.plot_control <- plotUMAP(control, colour_by=I(shared.pseudo)) + labs(title = 'Pseudotime - Control')
pseudo.plot_control

# Trajectory inference - CeD
ced <- slingshot(ced, clusterLabels = colData(ced)$cell.type.3,
                     reducedDim = 'UMAP', start.clus = 'Stem cells')

# Check lineages
lineages_ced <- slingLineages(ced)
pseudo.paths <- slingPseudotime(ced)
head(pseudo.paths)

# Check trajectory IDs
curve.assignments <- slingBranchID(ced)
table(curve.assignments)

# Get share pseudotime
shared.pseudo <- rowMeans(pseudo.paths, na.rm=TRUE)

# Plotting the lineages in the UMAP
# Need to loop over the paths and add each one separately.
gg_ced <- plotUMAP(ced, colour_by=I(shared.pseudo))
embedded <- embedCurves(ced, "UMAP")
embedded <- slingCurves(embedded)
colors <- brewer.pal(5, 'Set1')

i = 1
for (path in embedded) {
  embedded <- data.frame(path$s[path$ord,])
  gg_ced <- gg_ced + geom_path(data=embedded, aes(x=UMAP_1, y=UMAP_2), size=1.2, color = colors[i])
  i = i+1
}
gg_ced <- gg_ced + labs(title = 'Pseudotime and trajectories - CeD')
gg_ced

pseudo.plot_ced <- plotUMAP(ced, colour_by=I(shared.pseudo)) + labs(title = 'Pseudotime - CeD')
pseudo.plot_ced

library(patchwork)
pseudo.plot_control + pseudo.plot_ced

# Checking how the pseudotime changes depending on the condition
shared.pseudo_control <- tibble(
  pseudotime = rowMeans(slingPseudotime(control), na.rm=TRUE),
  condition = 'Ctrl')
shared.pseudo_ced <- tibble(
  pseudotime = rowMeans(slingPseudotime(ced), na.rm=TRUE),
  condition = 'CeD')
  
shared.pseudo_df <- rbind(shared.pseudo_control, shared.pseudo_ced)

p_pt <- ggplot(shared.pseudo_df, aes(x = pseudotime, color = condition, fill = condition)) +
  geom_density(alpha = 0.5) +
  theme_classic() +
  labs(title = 'Distribution of average pseudotime per condition')
p_pt

# Differential progression
# Kolmogorov-Smirnov Test to assess whether the two groups of pseudotime values 
# are derived from the same distribution
ks.test(shared.pseudo_control$pseudotime, shared.pseudo_ced$pseudotime)
# p-value < 2.2e-16

# Saving things
pdf(paste0(out_path, 'slingshot_trajectory_plots_epithelia.pdf'), width = 11, height = 8, paper = 'a4r')
gg_control
gg_ced
p_pt
dev.off()

saveRDS(control, paste0(out_path, 'control_epithelia_sce_slingshot.rds'))
saveRDS(ced, paste0(out_path, 'ced_epithelia_sce_slingshot.rds'))

## Differential expression along the trajectory - tradeSeq --------------------
# Control - trying just with lineage 2 (stem -> enterocyte mature)
pseudotime <- slingPseudotime(control)[,2]
cell.weigths <- slingCurveWeights(control)[,2]
# get cells with 0 weigth
valid.weigths <- cell.weigths[cell.weigths > 0]
# and exclude them from pseudotime and cell weights
pseudotime <- pseudotime[names(pseudotime) %in% names(valid.weigths)]
pseudotime <- as.matrix(pseudotime)
valid.weigths <- as.matrix(valid.weigths)
# and subset only these cells from the sce object
lin2 <- control[, colData(control)$barcode %in% names(valid.weigths)]
counts <- as.matrix(assays(lin2)$counts)

# Deciding k nodes
BPPARAM <- BiocParallel::bpparam()
set.seed(5)
icMat <- evaluateK(counts = counts,
                   pseudotime = pseudotime,
                   cellWeights = valid.weigths,
                   nGenes = 200,
                   k = 3:10,
                   verbose = T,
                   BPPARAM = BPPARAM)

# Fitting the GAM model
set.seed(6)

# Matrix with covariates
batch <- as.factor(lin2$batch)
sex <- as.factor(lin2$sex)
age <- as.factor(lin2$age)

U <- model.matrix(~batch + sex + age)
valid.weigths[,1] <- 1

lin2 <- fitGAM(counts = counts,
               pseudotime = pseudotime,
               cellWeights = valid.weigths,
               nknots = 6,
               verbose = TRUE)

saveRDS(lin2, paste0(out_path, 'lin2_sce_fitGAM.rds'))

table(rowData(lin2)$tradeSeq$converged)

assoRes <- associationTest(lin2)
head(assoRes)

startRes <- startVsEndTest(lin2)

gene <- 'LGR5'
plotSmoothers(lin2, counts, gene = gene) + labs(title = gene)

## Archive ----

# Differential expression
# Checking k nodes to use
set.seed(3)
icMat <- evaluateK(counts = as.matrix(assays(control)$counts),
                   pseudotime = slingPseudotime(control),
                   cellWeights = slingCurveWeights(control),
                   nGenes = 300,
                   k = 3:10,
                   verbose = T)

# Checking how the pseudotime changes depending on the condition
Pseudotime <- sce$slingPseudotime_1
Status <- sce$status
xx <- data.frame(Status, Pseudotime)
p_pt <- ggplot(xx, aes(x = Pseudotime, color = Status, fill = Status)) +
  geom_density(alpha = 0.5) +
  theme_classic()

pdf(paste0(out_path, 'slingshot_epithelia_plots.pdf'), width = 11, height = 8, paper= 'a4r')
plotUMAP(sce, colour_by = 'slingPseudotime_1')
plot(SlingshotDataSet(sce), type = 'lineages', show.constraints = T)
plot(SlingshotDataSet(sce), type = 'curves')
plotGeneCount(sce, clusters = colData(sce)$cell.type.4)
pairs(SlingshotDataSet(sce), type = 'lineages', show.constraints = T)
p_pt
dev.off()

saveRDS(sce, paste0(out_path, 'sce_epithelia_with_slingshot.rds'))
sce <- readRDS(paste0(out_path, 'sce_epithelia_with_slingshot.rds'))

# Differential expression
# Checking k nodes to use
set.seed(3)
icMat <- evaluateK(counts = as.matrix(assays(sce)$counts),
                   pseudotime = slingPseudotime(sce),
                   cellWeights = slingCurveWeights(sce),
                   conditions = factor(colData(sce)$status),
                   nGenes = 300,
                   k = 3:10,
                   verbose = T)

# Checking trajectories in control
control <- readRDS(paste0(out_path, 'control_epithelia_sce_slingshot.rds'))

lineages_ctrl <- slingLineages(control)
pseudo.paths <- slingPseudotime(control)
head(pseudo.paths)

print(lineages_ctrl)

library(ggraph)
library(igraph)

## Trying to build the lineage in a DF to build a tree plot -------

Lineages.dataframe <- function(lineage, lineage_index) {
  lin_df <- tibble()
  
  for (i in c(1:length(lineage))) {
    # Stop when reaches the last element
    if (i == length(lineage)) {} else {
      lin_df[i,1] <- lineage[i]
      lin_df[i,2] <- lineage[i+1]
      
      i <- i+1
    }
  }
  
  colnames(lin_df) <- c('from', 'to')
  lin_df$lineage <- as.factor(lineage_index)
  
  return(lin_df)
}

lin_df <- imap_dfr(lineages_ctrl, ~Lineages.dataframe(.x, .y))

test <- select(lin_df, -lineage)
graph <- graph_from_data_frame(lin_df)

ggraph(graph) +
  geom_edge_diagonal()
