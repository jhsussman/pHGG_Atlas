#Library Prep-----

library(ggplot2)
library(RColorBrewer)
library(infotheo)
library(aricode)
library(mclust)
library(SpatialExperiment)
library(pheatmap)
library(viridis)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(patchwork)
library(readr)
library(ggrepel)
library(stringr)
library(Polychrome)
library(purrr)
library(lsa)


#Analysis ------
setwd("/mnt/isilon/tan_lab/yanga11/CytoIntermediates/HyperGlioma")
processed<-read_csv("NolanInstancesProcessed.csv")
data<-read_csv("Nolan/NolanInstances/NolanInstances.csv")
data<-data%>%mutate(full=paste0(Sample_Name, "_", Instance))
data <- data %>%left_join(select(processed, full, cluster), by = "full")
data <- data %>%left_join(select(processed, full, treatment), by = "full")
data <- data %>% mutate(cluster = ifelse(cluster == 1, 2, cluster))
data <- data %>% mutate(cluster = ifelse(cluster == 0, 1, cluster))
data <- data %>%mutate(cluster = ifelse(is.na(cluster), 3, cluster))

for_plot <- table(as.character(data$treatment), data$cluster)
for_plot<-na.omit(for_plot)
neighborhood_order <- 1:3
neighborhood_counts <- as.numeric(table(factor(data$cluster)))
neighborhoods_names <- unique(data$cluster)
annotations <- unique(data$treatment)

p_value_matrix <- matrix(NA, nrow = length(neighborhoods_names), ncol = length(annotations),
                         dimnames = list(neighborhoods_names, annotations))
n <- sum(for_plot)
for (i in 1:length(neighborhoods_names)) {
  for (j in 1:length(annotations)) {
    obs_count <- sum(data$cluster == neighborhoods_names[i] & 
                       data$treatment == annotations[j])
    p_value <- phyper(obs_count - 1, 
                      sum(data$treatment == annotations[j]),
                      n-sum(data$treatment == annotations[j]),
                      sum(data$cluster == neighborhoods_names[i]),
                      lower.tail = FALSE)
    p_value_matrix[i, j] <- p_value
  }
}

p_value_matrix<-as.data.frame(p_value_matrix)
p_value_matrix<- p_value_matrix[order(as.numeric(row.names(p_value_matrix))), ]
rows<-rownames(p_value_matrix)
cols<-colnames(p_value_matrix)
for (i in 1:nrow(p_value_matrix)) {
  for (j in 1:ncol(p_value_matrix)) {
    #if (p_value_matrix[i,j]>0)
    if (p_value_matrix[i, j] < 0.05) {
      print(paste("Row:", rows[i], ", Column:", cols[j], ", Value:", p_value_matrix[i, j]))
    }
  }
}

#Reading in sub clusters of clusters----

data<-data%>%mutate(cluster=5)
filter<-read_csv("/home/yanga11/Desktop/Jupyter/cluster1_0.csv")
filter<-filter$Sample
filter<-sub( "\\.","-", filter)
data<-data%>%mutate(cluster= case_when(full %in% filter ~ 1, TRUE ~ cluster))
filter1<-read_csv("/home/yanga11/Desktop/Jupyter/cluster1_2.csv")
filter1<-filter1$Sample
filter1<-sub( "\\.","-", filter1)
data<-data%>%mutate(cluster= case_when(full %in% filter1 ~ 2, TRUE ~ cluster))
filter2<-read_csv("/home/yanga11/Desktop/Jupyter/cluster0_0.csv")
filter2<-filter2$Sample
filter2<-sub( "\\.","-", filter2)
data<-data%>%mutate(cluster= case_when(full %in% filter2 ~ 3, TRUE ~ cluster))
filter3<-read_csv("/home/yanga11/Desktop/Jupyter/cluster0_1.csv")
filter3<-filter3$Sample
filter3<-sub( "\\.","-", filter3)
data<-data%>%mutate(cluster= case_when(full %in% filter3 ~ 4, TRUE ~ cluster))


for_plot <- table(as.character(data$Cell_Type), data$cluster)
for_plot<-na.omit(for_plot)
neighborhood_order <- 1:5
neighborhood_counts <- as.numeric(table(factor(data$cluster)))
neighborhoods_names <- unique(data$cluster)
annotations <- unique(data$Cell_Type)

p_value_matrix <- matrix(NA, nrow = length(neighborhoods_names), ncol = length(annotations),
                         dimnames = list(neighborhoods_names, annotations))
n <- sum(for_plot)
for (i in 1:length(neighborhoods_names)) {
  for (j in 1:length(annotations)) {
    obs_count <- sum(data$cluster == neighborhoods_names[i] & 
                       data$Cell_Type == annotations[j])
    p_value <- phyper(obs_count - 1, 
                      sum(data$Cell_Type == annotations[j]),
                      n-sum(data$Cell_Type == annotations[j]),
                      sum(data$cluster == neighborhoods_names[i]),
                      lower.tail = FALSE)
    p_value_matrix[i, j] <- p_value
  }
}

p_value_matrix<-as.data.frame(p_value_matrix)
p_value_matrix<- p_value_matrix[order(as.numeric(row.names(p_value_matrix))), ]
rows<-rownames(p_value_matrix)
cols<-colnames(p_value_matrix)
for (i in 1:nrow(p_value_matrix)) {
  for (j in 1:ncol(p_value_matrix)) {
    #if (p_value_matrix[i,j]>0)
    if (p_value_matrix[i, j] < 0.05) {
      print(paste("Row:", rows[i], ", Column:", cols[j], ", Value:", p_value_matrix[i, j]))
    }
  }
}

# Initialize a matrix to hold the odds ratios
odds_ratio_matrix <- matrix(NA, nrow = length(neighborhoods_names), ncol = length(annotations),
                            dimnames = list(neighborhoods_names, annotations))
for_plot<-for_plot/1000
n <- sum(for_plot)
for (i in 1:length(neighborhoods_names)) {
  for (j in 1:length(annotations)) {
    # Observed count of the cell type in the neighborhood
    a <- for_plot[annotations[j], neighborhoods_names[i]]  # Cell type in cluster
    b <- sum(for_plot[, neighborhoods_names[i]], na.rm = TRUE) - a  # Other cell types in cluster
    c <- sum(for_plot[annotations[j], ], na.rm = TRUE) - a  # Cell type in other clusters
    d <- n - (a + b + c)  # Other cell types in other clusters
    
    # Check for zero or NA values that would lead to division errors
    if (is.na(b) || is.na(c) || b == 0 || c == 0) {
      odds_ratio <- NA  # Set odds ratio to NA if b or c is zero or NA
      message("Skipping (b == 0 or c == 0 or NA) for: ", neighborhoods_names[i], ", ", annotations[j])
    } else {
      # Calculate the odds ratio
      odds_ratio <- (a * d) / (b * c)
    }
    
    # Store the odds ratio in the matrix
    odds_ratio_matrix[i, j] <- odds_ratio
  }
}

odds_ratio_matrix[is.infinite(odds_ratio_matrix)] <- NA  # Handle infinite values
odds_ratio_matrix[is.na(odds_ratio_matrix)] <- 0# Set NA values to 0 if needed
odds_ratio_matrix_log<-log(odds_ratio_matrix)
pheatmap(odds_ratio_matrix_log, 
         color = colorRampPalette(c("dark blue", "white", "dark red"))(100), 
         scale = "column",
         cluster_rows = FALSE,        
         cluster_cols = FALSE,              
         cellwidth = 20,                       
         cellheight = 20,                      
         main = "Odds Ratio Heatmap (Cell Types in Clusters)")


processed<-processed%>%mutate(cluster=0)
processed<-processed%>%mutate(cluster= case_when(full %in% filter ~ 1, TRUE ~ cluster))
processed<-processed%>%mutate(cluster= case_when(full %in% filter1 ~ 2, TRUE ~ cluster))
processed<-processed%>%mutate(cluster= case_when(full %in% filter2 ~ 3, TRUE ~ cluster))
processed<-processed%>%mutate(cluster= case_when(full %in% filter3 ~ 4, TRUE ~ cluster))
processed%>%group_by(sample)%>%summarize(count=n())%>%arrange(count)
processed$treatment <- factor(processed$treatment, levels = c("Pre-Treatment", "Post-Treatment"))
processed$cluster <- factor(processed$cluster, levels = c("0", "1", "2", "3", "4"))
ggplot(processed, aes(x = `UMAP 1`, y = `UMAP 2`)) +
  geom_point(aes(color = factor(cluster)), size = 0.3) +  # Color points by 'cluster'
  labs(x = "UMAP 1", y = "UMAP 2", title = "UMAP 1 and UMAP 2 colored by Cluster") +  
  scale_color_manual(values = c("0" = "gray", "1" = "purple", "2" = "green", "3"="red","4"="cyan"),  # Custom colors for clusters
                     name = "Cluster",  # Title for legend
                     labels = c("Cluster 0", "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4")) +  # Custom labels for the legend
  theme_minimal() +  # Use a clean theme
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)  # Add black box around the plot
  ) + 
  coord_fixed()  # Set equal aspect ratio


ggplot(processed, aes(x = `UMAP 1`, y = `UMAP 2`)) +
  geom_point(aes(color = factor(treatment)), size=0.3) +  # Color points by 'cluster'
  labs(x = "UMAP 1", y = "UMAP 2", title = "UMAP 1 and UMAP 2 colored by Cluster") +  
  scale_color_manual(values = c("Pre-Treatment" = "#3399FF", "Post-Treatment" = "#FF6666"),  # Set custom colors
                     name = "Timepoint",  # Title for legend
                     labels = c("Pre-Treatment", "Post-Treatment")) +  # Custom labels for the legend
  theme_minimal() +  # Use a clean theme
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    panel.border = element_rect(color="black", fill=NA, linewidth=1)  # Add black box around the plot
  )+coord_fixed()


t_test_PC1 <- t.test(processed$`PCA 1` ~ processed$treatment)
print(t_test_PC1)
wilcox_PC1 <- wilcox.test(processed$`PCA 1` ~ treatment, data = processed)
print(wilcox_PC1)


wilcox_results <- list()
for (cluster in 1:4) {
  # Subset data for the current cluster and for cluster 5
  cluster_data <- subset(data, cluster == cluster)
  cluster_5_data <- subset(data, cluster == 5)
  
  # Perform the Wilcoxon test comparing 'PCA 1' in current cluster vs. cluster 5
  wilcox_test <- wilcox.test(cluster_data$`contribution`, cluster_5_data$`contribution`)
  
  # Store the result with the cluster label
  wilcox_results[[paste("Cluster", cluster, "vs Cluster 5")]] <- wilcox_test
}

# Print all results
print(wilcox_results)



# Ensure the cluster column is a factor
processed$cluster <- factor(processed$cluster)

# Initialize an empty list to store results
wilcox_results <- list()

# Loop through clusters 1 to 4 and perform the Wilcoxon test against cluster 0
for (i in 1:4) {
  # Subset the data for cluster 0 and the current cluster
  cluster_0_data <- processed$`PCA 1`[processed$cluster == "0"]
  current_cluster_data <- processed$`PCA 1`[processed$cluster == as.character(i)]
  
  # Perform Wilcoxon test
  test_result <- wilcox.test(current_cluster_data, cluster_0_data)
  
  # Store result in the list
  wilcox_results[[paste("Cluster", i, "vs Cluster 0")]] <- test_result
}

# Display results
wilcox_results



library(pROC)

# Assuming df is your dataframe with PC1 and treatment columns
# Convert treatment column to a binary outcome for ROC calculation
processed$treatment_binary <- ifelse(processed$treatment == "Post-Treatment", 1, 0)

# Calculate ROC curve
roc_result <- roc(processed$treatment_binary, processed$`PCA 1`, direction=">")

# Flip x-axis and ensure y-axis is restricted to relevant range
plot.roc(
  
  roc_result, add=FALSE,
  main = "ROC Curve for PC1 predicting Pre/Post Treatment", 
  xlim = c(1, 0),        # Flip x-axis from 1 to 0 to follow ROC convention
  ylim = c(0, 1),        # Set y-axis from 0 to 1
  print.auc = TRUE,      # Print the AUC on the plot
  asp = NA               # Maintain aspect ratio of 1
)







# Calculate AUC
auc_value <- auc(roc_result)
print(paste("AUC:", auc_value))


ggplot(processed, aes(x = treatment, y = `PCA 1`, fill = treatment)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, position = position_dodge(0.9), outlier.shape = NA) +  # Adding a boxplot for more detail
  labs(title = "Comparison of PC1 values by Treatment",
       x = "Treatment",
       y = "PC1 Value") +
  theme_bw() +
  scale_fill_manual(values = c("Pre-Treatment" = "#3399FF", "Post-Treatment" = "#FF6666"))

ggplot(processed, aes(x = cluster, y = `PCA 1`, fill = cluster)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, position = position_dodge(0.9), outlier.shape = NA) +  # Adding a boxplot for more detail
  labs(title = "Comparison of PC1 values by Cluster",
       x = "Treatment",
       y = "PC1 Value") +
  theme_bw() +
  scale_fill_manual(
    values = c("0" = "gray", "1" = "#800080", "2" = "#008000", "3" = "#FF0000", "4" = "#00FFFF"),  
    name = "Cluster",  
    labels = c("Cluster 0", "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4")
  )
