library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(reshape2)
library(gridExtra)

#Barplots 
two_timepoint_barplot <- function(data, x, y, width, height, colors, filename) {
  data[[y]] <- factor(data[[y]], levels = names(colors))
  plot  <- ggplot(data, aes(x = !!rlang::sym(x), fill = !!rlang::sym(y))) +
    geom_bar(position = "fill") +
    scale_y_continuous(expand = c(0, 0),
                       breaks = seq(0, 1, 0.2),  
                       labels = function(x) paste0(x * 100)) +  
    labs(y = "Cell type percentage (%)") + 
    scale_x_discrete(labels = c("Initial Resection", "Post Therapy")) +  
    scale_fill_manual(values = colors) +
    theme_bw() +
    coord_fixed(ratio = 8) + 
    RotatedAxis() +
    theme(
      axis.text.x = element_text(
        size = 6,           
        family = "sans",   
        color = "black", 
        margin = margin (r=1)
      ),
      axis.text.y = element_text(
        size = 6,           
        family = "sans",   
        color = "black", 
        margin = margin (r=1)
      ),
      axis.title.x = element_blank(),  
      axis.title.y = element_text(size = 7, family = "sans", color = "black"), 
      axis.ticks = element_line(size = 0.3, lineend = "butt"), 
      panel.grid.major = element_blank(),  
      panel.grid.minor = element_blank(), 
      panel.border = element_blank(), 
      axis.line = element_line()) #+
  #guides(fill = "none")
  print(plot)
  ggsave(plot, filename=filename, device="pdf", 
         width = width, height = height, units = "cm")
  
}

general_barplot <- function(data, x, y, width, height, labels, colors, filename) {
  data$timepoint = factor(data$timepoint)
  data[[y]] <- factor(data[[y]], levels = names(colors))
  plot  <- ggplot(data, aes(x = !!rlang::sym(x), fill = !!rlang::sym(y))) +
    geom_bar(position = "fill") +
    scale_y_continuous(expand = c(0, 0),
                       breaks = seq(0, 1, 0.2),  
                       labels = function(x) paste0(x * 100)) +  
    labs(y = "Cell type percentage (%)") + 
    scale_x_discrete(labels = labels) +  
    scale_fill_manual(values = colors) +
    theme_bw() +
    coord_fixed(ratio = 8) + 
    RotatedAxis() +
    theme(
      axis.text.x = element_text(
        size = 6,           
        family = "sans",   
        color = "black", 
        margin = margin (r=1)
      ),
      axis.text.y = element_text(
        size = 6,           
        family = "sans",   
        color = "black", 
        margin = margin (r=1)
      ),
      axis.title.x = element_blank(),  
      axis.title.y = element_text(size = 7, family = "sans", color = "black"), 
      axis.ticks = element_line(size = 0.3, lineend = "butt"), 
      panel.grid.major = element_blank(),  
      panel.grid.minor = element_blank(), 
      panel.border = element_blank(), 
      axis.line = element_line()) #+
    #guides(fill = "none")
  print(plot)
  ggsave(plot, filename=filename, device="pdf", 
         width = width, height = height, units = "cm")
  
}

samples_barplot <- function(data, x, y, width, height, colors, filename) {
  data[[y]] <- factor(data[[y]], levels = names(colors))
  plot  <- ggplot(data, aes(x = !!rlang::sym(x), fill = !!rlang::sym(y))) +
    geom_bar(position = "fill", show.legend = TRUE) +
    scale_y_continuous(expand = c(0, 0),
                       breaks = seq(0, 1, 0.2),  
                       labels = function(x) paste0(x * 100)) +  
    labs(y = "Cell type percentage (%)") + 
    scale_fill_manual(values = colors) +
    theme_bw() +
    coord_fixed(ratio = 8) + 
    RotatedAxis() +
    theme(
      axis.text.x = element_text(
        size = 6,           
        family = "sans",   
        color = "black", 
        margin = margin (r=1)
      ),
      axis.text.y = element_text(
        size = 6,           
        family = "sans",   
        color = "black", 
        margin = margin (r=1)
      ),
      axis.title.x = element_blank(),  
      axis.title.y = element_text(size = 7, family = "sans", color = "black"), 
      axis.ticks = element_line(size = 0.3, lineend = "butt"), 
      panel.grid.major = element_blank(),  
      panel.grid.minor = element_blank(), 
      panel.border = element_blank(), 
      axis.line = element_line()) 
  print(plot)
  ggsave(plot, filename=filename, device="pdf", 
         width = width, height = height, units = "cm")
  
}

#Boxplots
get_fraction_table <- function(data, label_category, label, time, patient){
  patients_with_two_timepoints <- data %>%
    group_by(!!rlang::sym(patient)) %>%
    summarise(num_timepoints = n_distinct(!!rlang::sym(time))) %>%
    filter(num_timepoints == 2) %>%
    pull(!!rlang::sym(patient))
  
  filtered_data <- data %>%
    filter(!!rlang::sym(patient) %in% patients_with_two_timepoints)
  fraction_table <- prop.table(table(filtered_data[[time]], 
                                     filtered_data[[label_category]], 
                                     filtered_data[[patient]]), margin = c(3,1))
  fraction_data <- as.data.frame(fraction_table)
  colnames(fraction_data) <- c("timemerge", "label", "patient", "fraction")
  fraction_data <- fraction_data[fraction_data$label == label, ]
  return(fraction_data)
}
custom_labels <- function(x) {
  x * 100
}
create_paired_box_plot <- function(fractions, label, time, patient) {
  library(car)
  #fractions$fraction = logit(fractions$fraction)
  plot <- ggpaired(
    data = fractions,
    x = time,
    y = "fraction",
    id = "patient",
    yscale = "none",
    title = label,
    ylab = "Cell Type Percentage (%)",
    palette = c("#3399FF", "#FF6666"), 
    fill = time, 
    line.size = 0.2, 
    point.size = 1.5) +
    stat_compare_means(method="wilcox.test", 
                       paired = TRUE, size = 4, 
                       #method.args = list(alternative = "greater"),
                       aes(label = paste0("p = ", after_stat(p.format))))+
    labs(x = "") +
    theme(axis.ticks.x = element_blank(), 
          axis.text.x = element_blank(), 
          axis.title.y = element_blank(), 
          axis.text.y = element_text(size = 12, family = "sans"), 
          plot.title = element_text(size = 14, family = "sans", hjust = 0.5)) +
    scale_y_continuous(labels = custom_labels) +
    guides(fill = "none") 
  return(plot)
}

timepoint_boxplot <- function(data, label_in, time, patient, width, height, 
                              rows, cols, label_order, filename){
  label_list <- label_order 
  plot_list <- lapply(label_list, function(label) {
    fraction_table <- get_fraction_table(data, label_in, as.character(label), time, patient)
    create_paired_box_plot(fraction_table, label, time, patient)
  })
  plot_grid <- plot_list %>% patchwork::wrap_plots(nrow = rows, ncol = cols) 
  print(plot_grid)
  ggsave(plot_grid, filename=filename, device="pdf", 
         width = width, height = height, units = "cm")
}

#UMAP Plot
create_umap_plot <- function(seurat_object, metadata_category, color_mapping, #order, 
                             width, height, ptsize, alpha, filename) {
  umap_plot <- Seurat::DimPlot(
    object = seurat_object,
    group.by = metadata_category,
    cols = color_mapping, label = TRUE, 
    #order = order, 
    label.size = 3, raster = TRUE, 
    raster.dpi = c(3000,3000), pt.size = ptsize, alpha = alpha,
  ) +
    theme_void() +
    coord_fixed() +
    labs(title = NULL)
  print(umap_plot)
  ggsave(umap_plot, filename=filename, device="pdf", 
         width = width, height = height, units = "cm")
}

#FeaturePlots grid 
create_feature_plots_grid <- function(
    seurat_object,
    features,
    rows, cols, width, height, filename, ptsize,
    max.cutoff = 'q99',
    min.cutoff = 'q1', alpha
) {
  plot_list <- list()
  for (feature in features) {
    plot <- Seurat::FeaturePlot(
      object = seurat_object,
      features = feature,
      raster = TRUE,
      cols = c( "gray","red"), 
      max.cutoff = max.cutoff,
      min.cutoff = min.cutoff,
      raster.dpi = c(600,600),
      pt.size = ptsize, alpha = alpha
    ) +
      theme_void() +
      coord_fixed() +
      labs(title = bquote(bold(.(feature))))+
      theme(plot.title = element_text(hjust = 0.5, size = 22), 
            legend.text = element_text(size = 14, family = "sans"))
    plot_list <- append(plot_list, list(plot))
  }
  grid <- do.call(grid.arrange, c(plot_list, ncol = cols, nrow = rows))
  print(grid)
  ggsave(grid, filename=filename, device="pdf", 
         width = width, height = height, units = "cm")
}

create_signature_plots_grid <- function(
    seurat_object,
    features,
    rows, cols, width, height, filename, ptsize,
    max.cutoff = 'q99',
    min.cutoff = 'q1'
) {
  plot_list <- list()
  for (feature in features) {
    plot <- Seurat::FeaturePlot(
      object = seurat_object,
      features = feature,
      raster = TRUE,
      max.cutoff = max.cutoff,
      min.cutoff = min.cutoff,
      raster.dpi = c(1200,1200),
      pt.size = ptsize
    ) +
      theme_void() +
      scale_color_gradientn(colors = rev(brewer.pal(11 ,"RdBu")))+ 
      coord_fixed() +
      labs(title = bquote(bold(.(feature))))+
      theme(plot.title = element_text(hjust = 0.5, size = 22), 
            legend.text = element_text(size = 14, family = "sans"))
    plot_list <- append(plot_list, list(plot))
  }
  grid <- do.call(grid.arrange, c(plot_list, ncol = cols, nrow = rows))
  print(grid)
  ggsave(grid, filename=filename, device="pdf", 
         width = width, height = height, units = "cm")
}
#Pie chart
createSeuratPieChart <- function(seurat_obj, metadata_col, color, title = NULL, 
                                 output_file = NULL, width = 7, height = 7) {
  metadata_values <- seurat_obj@meta.data[[metadata_col]]
  metadata_freq <- table(metadata_values)
  pie_data <- data.frame(
    labels = names(metadata_freq),
    values = as.vector(metadata_freq)
  )
  p <- ggplot(pie_data, aes(x = "", y = values, fill = labels)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar(theta = "y") +
    scale_fill_manual(values = color) +
    theme(legend.text = element_text(size = 12, family = "sans")) + 
    guides(fill = guide_legend(title = NULL)) +
    theme_void() + ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5, vjust = -5.4, size = 32, family = "sans"))
  if (!is.null(output_file)) {
    ggsave(output_file, plot = p, width = width, height = height)
  }
  return(p)
}

