##
# Function to plot VlnPlots and FeaturePlots for all markers in a given dataset
##

## This function requires the following arguments:
# obj: A clustered Seurat object.
# markers: A .csv file with two colums: 1 - cell.type, 2 - marker. One marker per row, multiple markers per cell are allowed.
# output: Path where the files should be saved. Plots are saved as VlnPlots.pdf and FeaturePlots.pdf. Default is current working directory
##

## Output
# Two pdf files containing all the plots for all markers and cell types provided
##

library(Seurat)
library(ggplot2)
library(tidyverse)

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