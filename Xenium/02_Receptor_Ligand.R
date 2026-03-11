library(liana)
library(Seurat)
library(stringr)
library(tidyverse)
library(magrittr)
library(circlize)
library(gridExtra) 
library(ComplexHeatmap) 
library(dplyr)
library(spatstat)
library(Matrix)
library(imcRtools)
library(SingleCellExperiment)
library(data.table)

setwd("/mnt/isilon/tan_lab/sussmanj/Xenium/pHGG")
options(future.globals.maxSize = 80000 * 1024^2)
options(Seurat.object.assay.version = "v5")

#Get overlapping receptor ligand gene pairs from LIANA
lr_list = select_resource(c('Consensus'))$Consensus
receptor_ligand_df <- lr_list[, c("source_genesymbol", "target_genesymbol")]

split_df <- data.frame(source_genesymbol = character(), target_genesymbol = character(), stringsAsFactors = FALSE)
for (i in 1:nrow(receptor_ligand_df)) {
  source_split <- unlist(strsplit(receptor_ligand_df$source_genesymbol[i], "_"))
  target_split <- unlist(strsplit(receptor_ligand_df$target_genesymbol[i], "_"))
  for (source_gene in source_split) {
    for (target_gene in target_split) {
      split_df <- rbind(split_df, data.frame(source_genesymbol = source_gene, target_genesymbol = target_gene, stringsAsFactors = FALSE))
    }
  }
}
split_df <- unique(split_df)
colnames(split_df) = c("Ligand", "Receptor")
head(split_df)

xenium.obj <- readRDS("/mnt/isilon/tan_lab/sussmanj/Xenium/pHGG/Seurat_Objects2/Merged_Seurat_Process_RPCA_C10_F5_Log.rds")
xenium_genes = rownames(xenium.obj)

lr_df <- split_df[split_df$Ligand %in% xenium_genes & split_df$Receptor %in% xenium_genes, ]
head(lr_df)
unique(c(lr_df$Ligand, lr_df$Receptor))

#Differential expression 
#markers = FindAllMarkers(xenium.obj, only.pos = T)
#markers$diff.pct = markers$pct.1 - markers$pct.2

#Define Clusters
cluster_counts <- table(xenium.obj$Xenium_snn_res.0.3)
xenium.obj@meta.data <- xenium.obj@meta.data %>%
  mutate(Xenium_snn_res.0.3 = ifelse(Xenium_snn_res.0.3 %in% names(cluster_counts[cluster_counts < 10]), 
                                     "Other", Xenium_snn_res.0.3))

#DimPlot(xenium.obj, group.by = "Xenium_snn_res.0.3", reduction = "umap.rpca", 
#        alpha = 0.3, raster = F, label = T, label.box = T) + coord_fixed()
#VlnPlot(xenium.obj, assay = "Xenium", features = "SOX2", pt.size = 0, group.by = "Xenium_snn_res.0.3")
table(xenium.obj$Xenium_snn_res.0.1)

Idents(xenium.obj) = "Xenium_snn_res.0.3"
xenium.obj = RenameIdents(xenium.obj, 
                          "1" = "Tumor", 
                          "13" = "Tumor",
                          "2" = "Tumor", 
                          "24" = "Other",
                          "28" = "Other",
                          "29" = "Other",
                          "3" = "Other",
                          "30" = "Other",
                          "31" = "Myeloid",
                          "32" = "Other",
                          "33" = "Other",
                          "4" = "Other",
                          "5" = "Other",
                          "6" = "Other",
                          "7" = "Other")
xenium.obj$cell_type_coarse = Idents(xenium.obj)
p1 = DimPlot(xenium.obj, group.by = "cell_type_coarse", reduction = "umap.rpca", 
        alpha = 0.9, raster = T, raster.dpi = c(2400, 2400), label = T, label.box = T, pt.size = 2.8) + coord_fixed()
ggsave(plot = p1, filename = "Figures/UMAP_Global.pdf", width = 9, height = 9)

p1 = FeaturePlot(xenium.obj, features = "EGFR", reduction = "umap.rpca", 
             raster = T, raster.dpi = c(2400, 2400), pt.size = 2.8, max.cutoff = 'q99') + coord_fixed()
ggsave(plot = p1, filename = "Figures/UMAP_EGFR.pdf", width = 9, height = 9)

#Split object 
obj.list = SplitObject(xenium.obj, split.by = "orig.ident")

#Functions for empirical shuffling test 
get_connectivity_chart <- function(seurat_object, cell_type_A, cell_type_B){
  sce.object  <- SingleCellExperiment(seurat_object@assays$Xenium@layers$data, colData=seurat_object@meta.data)
  sce.object <- buildSpatialGraph(sce.object, img_id = "orig.ident", type = "expansion", 
                                  coords = c("x.coord", "y.coord"), threshold = distance_threshold)
  connectivity_chart = as.data.frame(colPair(sce.object, "expansion_interaction_graph"))
  all_cells = colnames(sce.object)
  
  print("Generating connectivity chart.")
  colnames(connectivity_chart) = c("Cell_Number_A", "Cell_Number_B")
  connectivity_chart$Cell_A = as.character(all_cells[connectivity_chart$Cell_Number_A])
  connectivity_chart$Cell_B = as.character(all_cells[connectivity_chart$Cell_Number_B])
  connectivity_chart$Cell_Type_A_Identity = as.character(seurat_object@meta.data$cell_type_coarse[connectivity_chart$Cell_Number_A])
  connectivity_chart$Cell_Type_B_Identity = as.character(seurat_object@meta.data$cell_type_coarse[connectivity_chart$Cell_Number_B])
  connectivity_chart = connectivity_chart[connectivity_chart$Cell_Type_A_Identity == cell_type_A & connectivity_chart$Cell_Type_B_Identity == cell_type_B, ]
  
  return(connectivity_chart)
}

calculate_ligand_receptor_score <- function(seurat_object, ligand_gene, 
                                            receptor_gene, 
                                            distance_threshold, 
                                            cell_type_A, cell_type_B, connectivity_chart) {
  
  normalized_data = GetAssayData(seurat_object, assay = "Xenium", layer = "data")
  ligand_expr_A <- normalized_data[ligand_gene, connectivity_chart$Cell_A]
  receptor_expr_B <- normalized_data[receptor_gene, connectivity_chart$Cell_B]
  
  lr_score <- sum(ligand_expr_A * receptor_expr_B)
  return(lr_score)
}

get_non_connecting_chart <-  function(seurat_object,
                                      distance_threshold, 
                                      cell_type_A, cell_type_B, connectivity_chart) {
  all_cells_A <- WhichCells(seurat_object, ident = cell_type_A)
  all_cells_B <- WhichCells(seurat_object, ident = cell_type_B)
  
  print("Generating all possible non-connecting  pairs.")
  connectivity_chart_dt <- as.data.table(connectivity_chart)
  connecting_pairs <- unique(connectivity_chart_dt[, .(Cell_A, Cell_B)])
  all_possible_pairs <- CJ(Cell_A = all_cells_A, Cell_B = all_cells_B)
  non_connecting_pairs <- all_possible_pairs[!connecting_pairs, on = .(Cell_A, Cell_B)]
  return(non_connecting_pairs)
}

empirical_shuffling <- function(seurat_object, ligand_gene, receptor_gene, 
                                distance_threshold, 
                                cell_type_A, cell_type_B, 
                                n_shuffles = 1000, 
                                connectivity_chart, non_connecting_pairs) {
  print("Getting LR score for contacting pairs.")
  actual_lr_score <- calculate_ligand_receptor_score(seurat_object, 
                                                     ligand_gene, receptor_gene, 
                                                     distance_threshold, cell_type_A, cell_type_B, connectivity_chart)
  shuffled_scores <- numeric(n_shuffles)
  normalized_data = GetAssayData(seurat_object, assay = "Xenium", layer = "data")
  
  print("Calculating shuffled distribution")
  for (i in seq_len(n_shuffles)) {
    if (i %% 100 == 0) {
      print(paste0("Shuffles: ", as.character(i)))
    }
    
    num_pairs = dim(connectivity_chart)[1]
    sampled_pairs <- non_connecting_pairs[sample(nrow(non_connecting_pairs), num_pairs)]
    
    ligand_expr_A_shuffled <- normalized_data[ligand_gene, sampled_pairs$Cell_A]
    receptor_expr_B_shuffled <- normalized_data[receptor_gene, sampled_pairs$Cell_B]
    
    shuffled_scores[i] <- sum(as.numeric(ligand_expr_A_shuffled) * as.numeric(receptor_expr_B_shuffled))
  }
  
  print("Computing statistical significance.")
  mean_shuffled_score <- mean(shuffled_scores)
  sd_shuffled_score <- sd(shuffled_scores)
  
  z_score <- (actual_lr_score - mean_shuffled_score) / sd_shuffled_score
  p_value <- pnorm(z_score, lower.tail = FALSE)  # One-tailed Z-test for enrichment
  fold_change <- actual_lr_score / mean_shuffled_score
  return(list(actual_score = actual_lr_score, 
              mean_shuffled = mean_shuffled_score, 
              sd_shuffled = sd_shuffled_score, 
              z_score = z_score, 
              p_value = p_value, 
              fold_change = fold_change))
}

#Now run it 
distance_threshold <- 30
cell_type_A <- "Myeloid"  
cell_type_B <- "Tumor"  
number_tests = 1000
results_file <- "Empirical_shuffling_results_Myeloid_Tumor_50_Log3_partial3.csv"

results_df <- data.frame(Sample = character(),
                         Ligand = character(),
                         Receptor = character(),
                         Actual_Score = numeric(),
                         Mean_Shuffled = numeric(),
                         SD_Shuffled = numeric(),
                         Z_Score = numeric(),
                         P_Value = numeric(),
                         Fold_Change = numeric(),
                         stringsAsFactors = FALSE)
for(j in c(3)){
  print(j)
  Idents(obj.list[[j]])= "cell_type_coarse"
  sample_name = unique(obj.list[[j]]$orig.ident)
  fov_name <- names(obj.list[[j]]@images)[1]
  coord_table = obj.list[[j]]@images[[fov_name]]@boundaries$centroids@coords
  obj.list[[j]]$x.coord <- coord_table[,1]
  obj.list[[j]]$y.coord <- coord_table[,2]
  
  connectivity_chart_input = get_connectivity_chart(obj.list[[j]], cell_type_A, cell_type_B)
  non_connecting_pairs_input = get_non_connecting_chart(obj.list[[j]], distance_threshold, cell_type_A, cell_type_B, connectivity_chart_input)
  
  for(k in 60:70){
  #for(k in 1:dim(lr_df)[1]){
    row = lr_df[k, ]
    print(row)
    ligand_gene <- row$Ligand
    receptor_gene <- row$Receptor
    
    result <- empirical_shuffling(obj.list[[j]], 
                                  ligand_gene, receptor_gene, 
                                  distance_threshold, 
                                  cell_type_A, cell_type_B, number_tests, 
                                  connectivity_chart_input, non_connecting_pairs_input)
    current_result <- data.frame(Sample = sample_name,
                                 Ligand = ligand_gene,
                                 Receptor = receptor_gene,
                                 Actual_Score = result$actual_score,
                                 Mean_Shuffled = result$mean_shuffled,
                                 SD_Shuffled = result$sd_shuffled,
                                 Z_Score = result$z_score,
                                 P_Value = result$p_value,
                                 Fold_Change = result$fold_change,
                                 stringsAsFactors = FALSE)
    results_df <- rbind(results_df, current_result)
    write.table(results_df, file = results_file, sep = ",", quote = F, col.names = !file.exists(results_file), 
                row.names = FALSE, append = TRUE)
  }
}

