# Loading required packages
library(Seurat)
library(UCell)
library(tidyverse)
library(patchwork)

## Function plotMarkers() ------------------------------------------------------
# Function to plot VlnPlots and FeaturePlots for all markers in a given dataset
#

## This function requires the following arguments:
# obj: A clustered Seurat object.
# markers: A .csv file with two colums: 1 - cell.type, 2 - marker. One marker per row, multiple markers per cell are allowed.
# output: Path where the files should be saved. Plots are saved as VlnPlots.pdf and FeaturePlots.pdf. Default is current working directory
##

## Output
# Two pdf files containing all the plots for all markers and cell types provided
##

plotMarkers <- function(obj, markers_csv, output = getwd()){
  
  # Loading the markers list
  markers_list <- read.csv(markers_csv)
  
  # Reformating the input table into a list with cell types and their respective markers
  cell.type <- unique(markers_list$cell.type)
  names(cell.type) <- cell.type
  cell.markers.list <- map(cell.type, ~markers_list$marker[markers_list$cell.type == .])
  
  # Ploting the VlnPlots
  Vln <- imap(cell.markers.list, ~VlnPlot(obj, features = .x) & labs(subtitle=paste0('marker for: ',.y)))
  pdf(file=paste0(output,'VlnPlots.pdf'), width = 11, height = 8, paper = 'a4r')
  print(Vln)
  dev.off()
  
  # Ploting the FeaturePlots
  Feat <- imap(cell.markers.list, ~FeaturePlot(F_obj, features = .x) & labs(subtitle=paste0('marker for: ',.y)))
  pdf(file=paste0(output,'FeaturePlots.pdf'), width = 11, height = 8, paper = 'a4r')
  print(Feat)
  dev.off()
  
  return(paste0('Plots saved to: ',output,'VlnPlots.pdf and ',output,'FeaturePlots.pdf'))
  
}

## Function plot_Ucell() -------------------------------------------------------
# This function is used to identify clusters using the Ucell (https://github.com/carmonalab/UCell) signature scoring package.
# Given a set of markers for a cell type, Ucell will calculate the score of cluster X being the cell type you gave the gene signature for.
# The output are Feature Plots and Violin Plots for each cell signature.
#
# Usage: plot_Ucell(data, markers, output, reduction = 'umap')
# data: a clustered Seurat object
# markers: a .csv file containing two columns: 1) cell type name (Cell.Type) 2) markers for that cell type separated by ', ' (Markers)
# output: file name to save the plots (saved as .pdf by default)
# reduction: reduction to use in the UCell score calculations. Default = 'umap'

PlotUCell <- function(data, markers, output, reduction = 'umap') {
  
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