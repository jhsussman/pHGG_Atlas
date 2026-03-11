library(spatstat)
library(png)
library(stringr)
library(foreach)
library(doParallel)
library(parallel)

setwd("/mnt/isilon/tan_lab/sussmanj/CODEX/Seurat_Analysis/CPTCA_Full")
integrated_seurat <- readRDS("/mnt/isilon/tan_lab/sussmanj/CODEX/Seurat_Analysis/CPTCA_Full/Annotated_Seurat_Objects/Merged_Annotated_Filtered_2_10162023.RDS")

table(integrated_seurat$annotations_lv2)
DefaultAssay(integrated_seurat) <- "CODEX"
sample_names <- c("7622", "6477", "5928", "4740", "4337", "3058", "942", "371", "339", "161", "5335")
pix = "Pix4"


target_cell_name = "Mature Neurons"
target_cell_prefix = "Neurons"

sample_list <- list()
celltype_list <- list()
distance_median_list <- list()
distance_sd_list <- list()
pop_list <- list()
pvalue_list <- list()
iqr_list <- list()

n_permutations <- 100

for (i in 1:length(sample_names)) {
  sample_name <- sample_names[i]
  print(sample_name)
  sample_id <- paste0("7316-",sample_name,"_",pix)
  sample_data <- subset(integrated_seurat, subset=orig.ident==sample_id)
  unique_cell_types <- unique(sample_data$annotations_lv2)
  cell_target <- subset(sample_data, subset=annotations_lv2==target_cell_name)
  cts <- cell_target@meta.data

  xmin <- min(cts$x.coord)
  xmax <- max(cts$x.coord)
  ymin <- min(cts$y.coord)
  ymax <- max(cts$y.coord)
  sample_data_crop <- subset(sample_data, subset = x.coord > xmin & x.coord < xmax & y.coord > ymin & y.coord < ymax)
  sample_data_crop_xy <- sample_data_crop@meta.data[, c('x.coord', 'y.coord')]
  for (j in 1:length(unique_cell_types)) {
    cell_type <- unique_cell_types[j]
    print(cell_type)
    coords_crop <- subset(sample_data_crop, subset = annotations_lv2 == cell_type)
    num_cells <- ncol(coords_crop)
    print(num_cells)
    if (num_cells > 50000) {
      coords_crop <- subset(coords_crop, downsample = max(0.25 * num_cells, 50000))
    } else {
      message("Number of cells for this type is too small to downsample.")
    }
    print(length(colnames(coords_crop)))
    x <- coords_crop$x.coord
    y <- coords_crop$y.coord
    population <- length(x)
    if (length(x) < 10) {
      next
    }
    
    cells <- ppp(x, y, c(xmin-10, xmax+10), c(ymin-10, ymax+10))
    cell_target_cts <- psp(cts$x.coord, cts$y.coord, cts$x.coord, cts$y.coord, window = cells$window)
    
    cm <-list(cell=cells,target=cell_target_cts)
    cm$dtarget <- distfun(cm$target)
        
    distance <- cm$dtarget(x,y)
    # remove NA values
    distance <- distance[!is.na(distance)]
    # remove infinite values
    distance <- distance[!is.infinite(distance)]
    distance_median <- median(distance)
    # distance_sd <- sd(distance)
    # compute the std of the distance from 25th to 75th percentile
    distance_sd <- sd(distance[which(distance > quantile(distance, 0.25) & distance < quantile(distance, 0.75))])
    
    # compute the interquartile range
    distance_iqr <- quantile(distance, 0.75) - quantile(distance, 0.25)
    
    # perform permutation test
    p_distance_median <- list() 
    for (k in 1:n_permutations) {
      print(k)
      # permute the cell coordinates
      permuted_xy <- sample_data_crop_xy[sample(nrow(sample_data_crop_xy), population), ]
      permuted_x <- permuted_xy[, 1]
      permuted_y <- permuted_xy[, 2]
      p_distance <- cm$dtarget(permuted_x, permuted_y)
      p_distance <- p_distance[!is.na(p_distance)]
      p_distance <- p_distance[!is.infinite(p_distance)]
      p_distance_median <- append(p_distance_median, median(p_distance))
    }
    
    sample_list <- append(sample_list, sample_name)
    celltype_list <- append(celltype_list, as.character(cell_type))
    pop_list <- append(pop_list, population)
    distance_median_list <- append(distance_median_list, distance_median)
    distance_sd_list <- append(distance_sd_list, distance_sd)
    pvalue_list <- append(pvalue_list, sum(p_distance_median < distance_median) / n_permutations)
    iqr_list <- append(iqr_list, distance_iqr)
  }
}

sample_vector <- as.character(sample_list)
writeLines(sample_vector,paste0("Distance_Analysis/sample_list_", target_cell_prefix, '.txt'))
celltype_vector <- as.character(celltype_list)
writeLines(celltype_vector, paste0("Distance_Analysis/celltype_list_", target_cell_prefix, '.txt'))
pop_vector <- as.character(pop_list)
writeLines(pop_vector, paste0("Distance_Analysis/pop_list_", target_cell_prefix,  '.txt'))
distance_median_vector <- as.character(distance_median_list)
writeLines(distance_median_vector, paste0("Distance_Analysis/distance_median_list_", target_cell_prefix, '.txt'))
distance_sd_vector <- as.character(distance_sd_list)
writeLines(distance_sd_vector, paste0("Distance_Analysis/distance_sd_list_", target_cell_prefix, '.txt'))
pvalue_vector <- as.character(pvalue_list)
writeLines(pvalue_vector, paste0("Distance_Analysis/pvalue_list_", target_cell_prefix, '.txt'))
iqr_vector <- as.character(iqr_list)
writeLines(iqr_vector, paste0("Distance_Analysis/iqr_list_", target_cell_prefix, '.txt'))