library(Seurat)
library(tidyverse)
library(patchwork)
library(readr)
library(dplyr) 
library(ggrepel)
library(ggplot2)
library(RColorBrewer)
library(sf)
library(jsonlite)
library(sfheaders)
library(ggrastr)

sample_names <- c("7622", "6477", "5928", "4740", "4337", "3058", "942", "371", "339", "161", "5335")
##################
i=10
##################
name = sample_names[i]
pix = "Pix4"
print(name)

options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = "v5")
setwd("/mnt/isilon/tan_lab/sussmanj/CODEX/Seurat_Analysis/CPTCA_Full/")
source("/mnt/isilon/tan_lab/sussmanj/CODEX/Seurat_Analysis/CODEX_Functions.R")

#Read in channel names
channel_names <- read_csv(file = "MarkerList.txt", col_names = FALSE)

#Specify the directory of your data
file_path <- paste0("/mnt/isilon/tan_lab/parvaresha/Pediatric_Glioma/pHGG_Oldridge_Samples/",name,"/mesmer/",pix,"/combined_markers.csv")
gc_csd_raw <- read_csv(file.path(file_path))

colnames(gc_csd_raw)[5:61] <- channel_names$X1
colnames(gc_csd_raw)[2] <- "Size"
colnames(gc_csd_raw)[1] <- "CellID" 
colnames(gc_csd_raw)[3:4] <- c("x.coord", "y.coord")
gc_csd_raw[,3:61] <- gc_csd_raw[,3:61] / gc_csd_raw$Size 
hist(gc_csd_raw$Size[2:length(gc_csd_raw$Size)], breaks = 200)

#Check QC for DAPI 
gc_csd_raw %>% ggplot(aes(x=`DAPI`)) + geom_histogram(bins = 200) + xlim(c(0, 300))
gc_csd_raw <- gc_csd_raw %>% dplyr::filter(`DAPI`>10, `DAPI`<250)
gc_csd_raw <- gc_csd_raw %>% dplyr::filter(`CellID`!=0)

#View Image
plot <- ggplot(gc_csd_raw, aes(x=`x.coord`, y = -`y.coord`)) + geom_point(alpha = 0.01) 
rasterize(plot, layers = 'Point', dpi=300)

#Read in QuPath annotations 
polygon_list <- c()
geojson_data <- fromJSON(paste0("Circles/ROI_7316-", name, ".geojson"))
features <- geojson_data$features
num_features <- dim(features)[1]
for (i in 1:num_features) {
  coordinates <- features$geometry$coordinates[[i]]
  x_coords <- coordinates[,,1]
  y_coords <- coordinates[,,2]
  polygon_data <- data.frame(
    x = x_coords,
    y = y_coords
  )
  polygon_list[[i]] <- polygon_data
}

#Function to reorder the polygon, when needed will go in circle 
reorder_polygon <- function(coords_matrix) {
  centroid <- colMeans(coords_matrix)
  #angles <- atan2(coords_matrix[, 2] - centroid[2], coords_matrix[, 1] - centroid[1])
  #sorted_coords <- coords_matrix[order(angles), ]
  sorted_coords <- coords_matrix
  if (!identical(sorted_coords[1, ], sorted_coords[nrow(sorted_coords), ])) {
    sorted_coords <- rbind(sorted_coords, sorted_coords[1, ])
  }
  return(sorted_coords)
}

#Apply to image
combined_filtered_data <- data.frame()
# Loop through each polygon and filter the data frame
for (i in 1:length(polygon_list)) {
  print(i)
  polygon_coords <- polygon_list[[i]]
  # Check if the polygon is closed (first and last point are the same)
  polygon_1_list <- list(polygon_coords)
  polygon_1_list[[1]] = data.matrix(reorder_polygon(polygon_1_list[[1]]))
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
gc_csd_raw <- combined_filtered_data

#Inspect output
plot <- ggplot(combined_filtered_data, aes(x=`x.coord`, y = -`y.coord`)) + geom_point(alpha = 0.01) 
rasterize(plot, layers = 'Point', dpi=300)

#Save output to variable and continue 
gc_csd_raw <- combined_filtered_data

# Remove the size and x/y coords and other erroneous values 
gc_csd_raw_clean <- gc_csd_raw %>% select(-contains("isilon")) 
rownames(gc_csd_raw_clean) <- rownames(gc_csd_raw_clean)

gc_csd2 <- gc_csd_raw %>% select(-starts_with('x.coord'), -contains("Empty"), -starts_with('y.coord'),-contains('size'),-contains("Std.Dev"), -contains("indexes"), -contains('RowSum'), -contains('Blank'), -contains("isilon"), -starts_with("Ch")) 
gc_csd2 <- gc_csd2 %>% select(-contains("Index"),-contains("CellID"), -matches("Python_Index"))

#Rename objects
gc_csd2 = gc_csd2 %>%
  rename_with(~str_remove_all(.x, 'Nucleus Intensity')) %>%
  rename_with(~str_remove_all(.x, '\\(.*\\) ')) %>%
  rename_with(~str_remove_all(.x, ' Mean'))
gc_csd2<-na.omit(gc_csd2)

# Create a Seurat object
t_gc_csd = t(as.matrix(gc_csd2))
colnames(t_gc_csd) <- rownames(gc_csd2)
rownames(t_gc_csd) <- colnames(gc_csd2)
Seurat <- CreateSeuratObject(t_gc_csd, project = paste0("7316-",name,"_",pix), assay="CODEX", meta.data = gc_csd_raw_clean)
saveRDS(Seurat, paste0("Individual_Seurat_Objects/CPTCA_7316-",name,"_",pix,"_revised.RDS"))       


