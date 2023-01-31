# Loading required packages
library(Seurat)
library(UCell)
library(tidyverse)
library(patchwork)

## Function plotMarkers() ------------------------------------------------------
# This function plots a FeaturePlot and a VlnPlot (optional) for each marker in a provided dataset.
#

# Usage: plotMarkers(data, markers, output)
# data: a clustered Seurat object
# markers: a .csv file containing two columns: 1) cell type name (Cell.Type) 2) markers for that cell type separated by ', ' (Markers)
# vln: TRUE plots also the VlnPlots for the markers, defaults is FALSE (will plot only FeaturePlots)
# output: file name to save the plots (saved as .pdf by default)

plotMarkers <- function(data, markers, vln = FALSE, output){
  
  # Loading and preparing the markers list
  signatures <- read.csv(markers, na.strings = '')
  
  signatures <- signatures %>% separate_rows(2, sep = ', ')
  signatures$Cell.Type <- str_replace_all(signatures$Cell.Type, '\\W', '.')
  cell_types <- unique(signatures$Cell.Type)
  
  signatures.list <- list()
  
  for (i in cell_types) {
    tmp <- signatures %>%
      filter(Cell.Type == i) %>%
      pull(Markers)
    
    signatures.list[[i]] <- tmp
  }
  
  # Ploting
  if (vln == TRUE) {
    Vln <- imap(signatures.list, ~VlnPlot(data, features = .x) + NoLegend())
    
    pdf(paste0(output,'_VlnPlots.pdf'), width = 11, height = 8, paper = 'a4r')
    print(Feat)
    dev.off()
    
    return(paste0('Violin plots saved to: ', output, '_VlnPlots.pdf'))
  }
  
  Feat <- imap(signatures.list, ~FeaturePlot(data, features = .x))
  
  # Saving
  pdf(paste0(output,'_FeaturePlots.pdf'), width = 11, height = 8, paper = 'a4r')
  print(Feat)
  dev.off()
  
  return(paste0('Feature plots saved to: ' ,output, '_FeaturePlots.pdf'))
  
}

## Function plotUcell() -------------------------------------------------------
# This function is used to identify clusters using the Ucell (https://github.com/carmonalab/UCell) signature scoring package.
# Given a set of markers for a cell type, Ucell will calculate the score of cluster X being the cell type you gave the gene signature for.
# The output are Feature Plots and Violin Plots for each cell signature.
#
# Usage: plotUcell(data, markers, output, reduction = 'umap')
# data: a clustered Seurat object
# markers: a .csv file containing two columns: 1) cell type name (Cell.Type) 2) markers for that cell type separated by ', ' (Markers)
# output: file name to save the plots (saved as .pdf by default)
# reduction: reduction to use in the UCell score calculations. Default = 'umap'

PlotUCell <- function(data, markers, output = getwd(), reduction = 'umap') {
  
  signatures <- read.csv(markers, na.strings = '')
  
  # Preparing the signatures list based on the input table
  signatures <- signatures %>% separate_rows(2, sep = ', ')
  signatures$Cell.Type <- str_replace_all(signatures$Cell.Type, '\\W', '.')
  cell_types <- unique(signatures$Cell.Type)
  
  signatures.list <- list()
  
  for (i in cell_types) {
    tmp <- signatures %>%
      filter(Cell.Type == i) %>%
      pull(Markers)
    
    signatures.list[[i]] <- tmp
  }
  
  # Calculate UCell score on a clustered SeuratObject
  obj <- AddModuleScore_UCell(data, features = signatures.list)
  
  # Use SmoothKNN to reduce sparsity of the UCell scores
  signature_names <- imap_chr(signatures.list, ~paste0(.y, '_UCell'))
  obj <- SmoothKNN(obj, reduction = reduction, signature.names = signature_names)
  
  # Ploting the results
  feat_names <- imap_chr(signatures.list, ~paste0(.y, '_UCell_kNN'))
  Feat <- map(feat_names, ~FeaturePlot(obj, features = .x, label = TRUE) + labs(title = ''))
  Vln <- map(feat_names, ~VlnPlot(obj, features = .x) + labs(title = '') + theme(legend.position = 'none'))
  
  pdf(paste0(output,'.pdf'), width = 11, height = 8, paper = 'a4r')
  #print(FeaturePlot(obj, features = feat_names, label = TRUE))
  #print(VlnPlot(obj, features = feat_names))
  map(cell_types, ~print(Feat[[.x]] / Vln[[.x]] + plot_annotation(title = paste0(names(feat_names[.x]), '_UCell_kNN_score'), subtitle = paste0('Markers used: ', str_c(signatures.list[[.x]], collapse = ', ')))))
  dev.off()
  
  return((paste0('Plots saved to ', output, '.pdf')))
  
}