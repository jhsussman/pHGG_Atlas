library(Seurat)
library(tidyverse)
library(patchwork)
library(readr)
library(dplyr)
library(ggrepel)
library(ggplot2)
library(pheatmap)
library(future)
library(dittoSeq)
library(gridExtra)
library(readr)
library(imcRtools)
library(SingleCellExperiment)
library(scales)
library(viridis)
library(ggpubr)
library(ggrastr)
library(sf)


#Function to filter Seurat 
filter_xenium_seurat <- function(xenium_seurat, coords_csv){
  # Read in Selection csv file 
  csv_data <- read.csv(coords_csv, row.names = NULL, header = T)
  polygon_list <- list()
  
  i_max = as.numeric(strsplit(csv_data$row.names[nrow(csv_data)], " ")[[1]][-1])
  for(i in 1:i_max){
    selection = paste0("Selection ", as.character(i))
    csv_data_i = csv_data[csv_data$row.names==selection, ]
    polygon_data <- data.frame(x = as.numeric(csv_data_i[,2]), y = as.numeric(csv_data_i[,3]))
    polygon_list[[length(polygon_list) + 1]] <- polygon_data
  }
  
  coord_table = xenium_seurat@images$fov@boundaries$centroids@coords
  xenium_seurat$x.coord <- coord_table[,1]
  xenium_seurat$y.coord <- coord_table[,2]
  gc_csd_raw = xenium_seurat@meta.data
  
  #Apply to image
  combined_filtered_data <- data.frame()
  # Loop through each polygon and filter the data frame
  for (i in 1:length(polygon_list)) {
    print(i)
    polygon_coords <- polygon_list[[i]]
    # Check if the polygon is closed (first and last point are the same)
    polygon_1_list <- list(polygon_coords)
    polygon_1_list[[1]] = data.matrix(polygon_1_list[[1]])
    x=polygon_1_list[[1]][,1]
    y=polygon_1_list[[1]][,2]
    plot(c(min(x), max(x)), c(min(-y), max(-y)), type = "n", xlab = "X", ylab = "Y")
    polygon(x, -y, col = "blue")
    polygon <- st_polygon(polygon_1_list)
    gc_csd_raw_sf <- st_as_sf(gc_csd_raw, coords = c("x.coord", "y.coord"))
    filtered_data <- st_intersection(gc_csd_raw_sf, polygon)
    # Append the filtered data to the combined_filtered_data data frame
    combined_filtered_data <- rbind(combined_filtered_data, as.data.frame(filtered_data))
  }
  
  filtered_coordinates <- st_coordinates(combined_filtered_data$geometry)
  combined_filtered_data$x.coord <- filtered_coordinates[, "X"]
  combined_filtered_data$y.coord <- filtered_coordinates[, "Y"]
  combined_filtered_data$geometry <- NULL
  head(combined_filtered_data)
  
  filtered_seurat = subset(xenium_seurat, cells = rownames(combined_filtered_data))
  dim(filtered_seurat)
  
  return(filtered_seurat)
}

#Function to plot image 
plot_points <- function(xenium_seurat, alpha = 0.01, dpi = 300){
  coord_table = xenium_seurat@images$fov@boundaries$centroids@coords
  xenium_seurat$x.coord <- coord_table[,1]
  xenium_seurat$y.coord <- coord_table[,2]
  plot <- ggplot(xenium_seurat@meta.data, aes(x=`x.coord`, y = -`y.coord`)) +
    geom_point(alpha = alpha) + theme_bw() + coord_fixed()
  plot_raster = rasterize(plot, layers = 'Point', dpi=dpi)
  print(plot_raster)
  return(plot_raster)
}

#Function to plot metadata
plot_metadata <- function(xenium_seurat, col = "orig.ident", alpha = 0.01, dpi = 300){
  coord_table = xenium_seurat@images$fov@boundaries$centroids@coords
  xenium_seurat$x.coord <- coord_table[,1]
  xenium_seurat$y.coord <- coord_table[,2]
  plot <- ggplot(xenium_seurat@meta.data, aes(x=`x.coord`, y = -`y.coord`)) +
    geom_point(alpha = alpha, aes(color = !!rlang::sym(col))) + theme_bw() + coord_fixed()
  plot_raster = rasterize(plot, layers = 'Point', dpi=dpi)
  print(plot_raster)
  return(plot_raster)
}

#Function to export mask
export_mask <- function(xenium_seurat, category = "orig.ident", filename){
  barcode_ids <- colnames(xenium_seurat)
  cluster_assignments <- xenium_seurat[[category]]
  output_df <- data.frame(cell_id = barcode_ids, group = cluster_assignments)
  write.csv(output_df,
            file = filename,
            row.names = FALSE
  )
}