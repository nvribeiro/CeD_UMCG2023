# Cell-cycle inspection
library(Seurat)
library(tidyverse)
library(patchwork)
library(ggrepel)
library(RColorBrewer)

out_path <- '/groups/umcg-wijmenga/tmp01/users/umcg-aramirezsanchez/umcg-nribeiro/NR03_scRNAseq/ongoing/outputs/subclustering/'

# Loading dataset
epithelia <- readRDS(paste0(out_path, 'SeuratObj_epithelia_v4_res09_identified_TAsplit.rds'))

# List of cell-cycle genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Assign cell-cycle score
epithelia <- CellCycleScoring(epithelia, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

# Inspecting
DimPlot(epithelia, reduction = 'umap', group.by = 'Phase', split.by = 'status')

p1 <- ggplot(epithelia@meta.data, aes(x = cell.type.3, fill = Phase)) +
  geom_bar(position = 'fill') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = '', y = 'proportion', title = 'Proportion of cells in each phase of cell-cycle')

p2 <- ggplot(epithelia@meta.data, aes(x = cell.type.3, fill = Phase)) +
  geom_bar(position = 'fill') +
  facet_wrap(~status) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = '', y = 'proportion', title = 'Proportion of cells in each phase of cell-cycle (per condition)')

# Inspecting stem alone
p3 <- ggplot(filter(epithelia@meta.data, seurat_clusters == 1 | seurat_clusters == 4), aes(x = cell.type.3, fill = Phase)) +
  geom_bar(position = 'fill') +
  facet_wrap(~status) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = '', y = 'proportion', title = 'Proportion of stem cells in each phase of cell-cycle (per condition)')

stem <- subset(epithelia, subset = seurat_clusters == 1 | seurat_clusters == 4)
DefaultAssay(stem) <- 'RNA'
stem <- stem %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()
stem <- CellCycleScoring(stem, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

stem <- RunPCA(stem, features = c(s.genes, g2m.genes))
DimPlot(stem, dims = c(1,2), group.by = 'Phase', shape.by = 'cell.type.3', split.by = 'status', pt.size = 2)
DimPlot(stem, dims = c(2,3), group.by = 'Phase', shape.by = 'cell.type.3', split.by = 'status', pt.size = 2)

# Correcting for cell cycle
stem <- ScaleData(stem, vars.to.regress = c('S.Score', 'G2M.Score'), features = row.names(stem))
stem <- RunPCA(stem, features = VariableFeatures(stem))
DimPlot(stem, dims = c(1,2), group.by = 'cell.type.3') + labs(title = 'PCA after cell-cycle correction')

# Now let's check the DEG using the corrected data
# Checking only the control
stem_ctrl <- subset(stem, subset = status == 'Ctrl')
stem_deg <- FindMarkers(stem_ctrl, slot = 'scale.data', ident.1 = 4, ident.2 = 1, test.use = 'MAST', min.pct = 0.25, logfc.threshold = 0)

stem_deg$gene <- rownames(stem_deg)

# Add info about the DE - being very stricted so I can get only the most significant changes
stem_deg$DE <- NA
stem_deg <- stem_deg %>%
  mutate(DE = case_when(
    p_val_adj < 0.01 & avg_diff> 0.5 ~ 'up',
    p_val_adj < 0.01 & avg_diff < -0.5 ~ 'down'))

# Add column with DEG names
stem_deg <- stem_deg %>%
  mutate(DE_gene = if_else(is.na(DE), NA, gene))

# Plot
ggplot(stem_deg, aes(x = avg_diff, y = -log10(p_val_adj), color = DE, label = DE_gene)) +
  geom_jitter() +
  geom_text_repel() +
  theme_classic() +
  labs(title = 'DEGs between stem cell 2 and 1 in control')

# Getting the genes that are not RP
`%!in%` = Negate(`%in%`)
stem_deg_f <- filter(stem_deg, DE_gene %!in% grep('RP.', stem_deg$DE_gene, value = T))

ggplot(stem_deg_f, aes(x = avg_diff, y = -log10(p_val_adj), color = DE, label = DE_gene)) +
  geom_jitter() +
  geom_text_repel() +
  theme_classic() +
  labs(title = 'DEGs between stem cell 2 and 1 in control')

# GSE on these genes
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ReactomePA)

ENTREZID <- select(org.Hs.eg.db, keys = stem_deg_f$gene, 
                   columns = 'ENTREZID', keytype = 'SYMBOL')

stem_deg_f <- left_join(stem_deg_f, ENTREZID, by = join_by('gene' == 'SYMBOL'))

GSE.BP <- compareCluster(ENTREZID~DE, fun = 'enrichGO', OrgDb = org.Hs.eg.db,
                         ont = 'BP', data = filter(stem_deg_f, is.na(DE) == FALSE), pvalueCutoff = 0.05)

dotplot(GSE.BP)

GSE.Pathway <- compareCluster(ENTREZID~DE, fun = 'enrichPathway',
                              data = filter(stem_deg_f, is.na(DE) == FALSE), pvalueCutoff = 0.05)

dotplot(GSE.Pathway)

# Expression inspections
FeaturePlot(epithelia, features = c('RBP2', 'FBP1', 'S100A6'))

FeaturePlot(epithelia, features = c('DACH1', 'EPHB2', 'SLC38A11'))

## Checking distribution of PCS per cluster and cycle -----------------
# PCA on cell cycle
stem_check <- RunPCA(stem, features = c(s.genes, g2m.genes))
DimPlot(stem_check, group.by = 'Phase')

# Building my PCA data frame
# Gettin cells info
metadata <- dplyr::select(stem_check@meta.data, barcode, cell.type.3, Phase)
pca <- as.data.frame(stem_check@reductions[["pca"]]@cell.embeddings)
pca$barcode <- row.names(pca)
pca <- left_join(pca, metadata, by = 'barcode')

pc1 <- ggplot(pca, aes(x = PC_1, color = Phase, fill = Phase)) +
  geom_density(alpha = 0.5) +
  theme_classic() +
  facet_wrap(~cell.type.3)

pc2 <- ggplot(pca, aes(x = PC_2, color = Phase, fill = Phase)) +
  geom_density(alpha = 0.5) +
  theme_classic() +
  facet_wrap(~cell.type.3)

pc1 | pc2


## Checking TAs, what genes pop up when I correct for the cell cycle? -------------
ta <- subset(epithelia, subset = seurat_clusters == 7 | seurat_clusters == 9)
DefaultAssay(ta) <- 'RNA'

ta <- ta %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()

ta <- CellCycleScoring(ta, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
ta <- RunPCA(ta, features = c(s.genes, g2m.genes))
DimPlot(ta, dims = c(1,2), group.by = 'Phase', shape.by = 'cell.type.3', pt.size = 2)

# Regressing out cell cycle
ta <- ScaleData(ta, vars.to.regress = c('S.Score', 'G2M.Score'), features = rownames(ta))
ta <- RunPCA(ta, features = VariableFeatures(ta))

DimPlot(ta, dims = c(1,2), group.by = 'cell.type.3')
DimPlot(ta, dims = c(2,3), group.by = 'cell.type.3')

ta_ctrl <- subset(ta, subset = status == 'Ctrl')
ta_deg <- FindMarkers(ta_ctrl, slot = 'scale.data', ident.1 = 9, ident.2 = 7, test.use = 'MAST', min.pct = 0.25, logfc.threshold = 0)

ta_deg$gene <- rownames(ta_deg)

# Add info about the DE - being very stricted so I can get only the most significant changes
ta_deg$DE <- NA
ta_deg <- ta_deg %>%
  mutate(DE = case_when(
    p_val_adj < 0.01 & avg_diff> 0.5 ~ 'up',
    p_val_adj < 0.01 & avg_diff < -0.5 ~ 'down'))

# Add column with DEG names
ta_deg <- ta_deg %>%
  mutate(DE_gene = if_else(is.na(DE), NA, gene))

# Plot
ggplot(ta_deg, aes(x = avg_diff, y = -log10(p_val_adj), color = DE, label = DE_gene)) +
  geom_jitter() +
  geom_text_repel() +
  theme_classic() +
  labs(title = 'DEGs between TA cell 2 and 1 in control')

# GSE on these genes
ENTREZID <- select(org.Hs.eg.db, keys = ta_deg$gene, 
                   columns = 'ENTREZID', keytype = 'SYMBOL')

ta_deg <- left_join(ta_deg, ENTREZID, by = join_by('gene' == 'SYMBOL'))

## Putting all stem and ta cells together -------------------
prolif <- subset(epithelia, subset = seurat_clusters == 1 | seurat_clusters == 4 | seurat_clusters == 7 | seurat_clusters == 9)

# Check the PCA
DefaultAssay(prolif) <- 'RNA'

prolif <- prolif %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()

# Check PCs based on cell cycle genes
prolif <- CellCycleScoring(prolif, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
prolif <- RunPCA(prolif, features = c(s.genes, g2m.genes))
DimPlot(prolif, dims = c(1,2), group.by = 'Phase', shape.by = 'cell.type.3', pt.size = 2)

metadata <- dplyr::select(prolif@meta.data, barcode, cell.type.3, Phase)
pca <- as.data.frame(prolif@reductions[["pca"]]@cell.embeddings)
pca$barcode <- row.names(pca)
pca <- left_join(pca, metadata, by = 'barcode')

pc1 <- ggplot(pca, aes(x = PC_1, color = cell.type.3, fill = cell.type.3)) +
  geom_density(alpha = 0.5) +
  theme_classic() +
  facet_wrap(~Phase)

pc2 <- ggplot(pca, aes(x = PC_2, color = cell.type.3, fill = cell.type.3)) +
  geom_density(alpha = 0.5) +
  theme_classic() +
  facet_wrap(~Phase)

pc1 | pc2

# Regressing out cell cycle - but using the difference in this case because I want to mantain the proliferation signal
prolif$CC.Difference <- prolif$S.Score - prolif$G2M.Score
prolif_corrected <- ScaleData(prolif, vars.to.regress = 'CC.Difference', features = rownames(prolif))
prolif_corrected <- RunPCA(prolif_corrected, features = VariableFeatures(prolif))

saveRDS(prolif_corrected, paste0(out_path, 'prolif_corrected.rds'))
prolif_corrected <- readRDS(paste0(out_path, 'prolif_corrected.rds'))

# Inspecting the PCA of the corrected object
# Pull info for 20 features, balanced positive/negative, for all PCs in given object
feature_list <- list()

for (i in 1:length(x = prolif_corrected@reductions$pca)) {
  feature_list[[i]] <- TopFeatures(object = prolif_corrected[["pca"]], dim = i, nfeatures = 40, balanced = T)
}

ElbowPlot(prolif_corrected)
colors <- brewer.pal(4, 'Paired')

p1 <- DimPlot(prolif_corrected, dims = c(1,2), group.by = 'cell.type.3', pt.size = 1, cols = colors)
p2 <- DimPlot(prolif_corrected, dims = c(2,3), group.by = 'cell.type.3', pt.size = 1, cols = colors)
p3 <- DimPlot(prolif_corrected, dims = c(3,4), group.by = 'cell.type.3', pt.size = 1, cols = colors)
p1 + p2 + p3 + plot_layout(guides = 'collect')

# Run DEG stem 2 vs TA 1 in control
prolif_ctrl <- subset(prolif_corrected, subset = status == 'Ctrl')
prolif_deg <- FindMarkers(prolif_corrected, slot = 'scale.data', ident.1 = 4, ident.2 = 7, test.use = 'MAST', min.pct = 0.25, logfc.threshold = 0)

prolif_deg$gene <- rownames(prolif_deg)

# Add info about the DE - being very stricted so I can get only the most significant changes
prolif_deg$DE <- NA
prolif_deg <- prolif_deg %>%
  mutate(DE = case_when(
    p_val_adj < 0.01 & avg_diff> 0.5 ~ 'up',
    p_val_adj < 0.01 & avg_diff < -0.5 ~ 'down'))

# Add column with DEG names
prolif_deg <- prolif_deg %>%
  mutate(DE_gene = if_else(is.na(DE), NA, gene))

# Plot
ggplot(prolif_deg, aes(x = avg_diff, y = -log10(p_val_adj), color = DE, label = DE_gene)) +
  geom_jitter() +
  geom_text_repel() +
  theme_classic() +
  labs(title = 'DEGs between Stem cell 2 and TA 1 in control')

# Getting the genes that are not RP
`%!in%` = Negate(`%in%`)
prolif_deg_f <- filter(prolif_deg, DE_gene %!in% grep('RP.', prolif_deg$DE_gene, value = T))

ggplot(prolif_deg_f, aes(x = avg_diff, y = -log10(p_val_adj), color = DE, label = DE_gene)) +
  geom_jitter() +
  geom_text_repel() +
  theme_classic() +
  labs(title = 'DEGs between Stem cell 2 and TA 1 in control')