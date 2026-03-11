library(ggplot2)
library(Seurat)
library(ggpubr)
library(readxl)
library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(gridExtra)
setwd("G:/My Drive/University of Pennsylvania/Labs/Tan Lab/HTAN/Drug_Screening/Figures")

drug_vector <- c(rep("DMSO", 8), rep(rep(c("AZD4547", "J-104129", "Varespladib", "Danusertib"), each = 2), times = 7))
concentration_vector <- c(
  rep(0, 8),  
  4000, 4000, 1000, 1000, 4000, 4000, 8000, 8000, 
  2000, 2000, 500, 500, 2000, 2000, 4000, 4000,   
  1000, 1000, 250, 250, 1000, 1000, 2000, 2000,   
  500, 500, 125, 125, 500, 500, 1000, 1000,    
  250, 250, 62.5, 62.5, 250, 250, 500, 500,  
  125, 125, 31.25, 31.25, 125, 125, 250, 250,    
  62.5, 62.5, 15.625, 15.625, 62.5, 62.5, 125, 125)

excel_columns <- function(start, end) {
  letters <- c(LETTERS)
  # Initialize a vector to store column names
  col_names <- character(0)
  col_names <- c(col_names, letters)
  for (first in letters) {
    for (second in letters) {
      col_names <- c(col_names, paste0(first, second))
    }
  }
  return(col_names[match(start, col_names):match(end, col_names)])
}
live_index = excel_columns("C", "BN")
cytotox_index = excel_columns("BR", "EC")
caspase_index = excel_columns("EG", "GR")

#1
excel = read_excel("radiation cytotox caspase 1763.xlsx", col_names = F, sheet = "Final")
timepoints = c(0,6,12,18,24,30,36,42,48,54,60,66,72)
timepoint_starts = 22:34
cellline = "1763NR"

#2
excel = read_excel("radiation cytotox caspase 1769.xlsx", col_names = F, sheet = "Final")
timepoints = c(0,6,12,18,24,30,36,42,48,54,60,66,72)
timepoint_starts = 22:34
cellline = "1769NR"

#3
excel = read_excel("radiation cytotox caspase 3058nr.xlsx", col_names = F, sheet = "Final")
timepoints = c(0,6,12,18,24,30,36,42,48,54,60,66,72)
timepoint_starts = 22:34
cellline = "3058NR"

#4
excel = read_excel("radiation cytotox caspase 913NR.xlsx", col_names = F, sheet = "Final")
timepoints = c(0,6,12,18,24,30,36,42,48,54,60,66,72)
timepoint_starts = 22:34
cellline = "913NR"

######################

# Function to convert Excel cell references to R values
extract_indices <- function(reference, data) {
  col_letter <- gsub("\\$", "", gsub("\\d", "", reference))
  col_index <- 0
  for (i in 1:nchar(col_letter)) {
    col_index <- col_index * 26 + (utf8ToInt(toupper(substr(col_letter, i, i))) - 64)
  } 
  row_number <- as.numeric(gsub("\\D", "", reference))
  return(as.numeric(data[row_number, col_index]))
}

build_drug_concentration_df <- function(drug_vector, cell_line, concentration_vector, timepoints, timepoint_starts, column_vector, data) {
  total_rows <- length(drug_vector) * length(timepoints)
  drug_concentration_df <- data.frame(Drug = character(total_rows),
                                      CellLine = cell_line, 
                                      Concentration = numeric(total_rows),
                                      Timepoint = numeric(total_rows),
                                      Value = numeric(total_rows),
                                      stringsAsFactors = FALSE)
  row_counter <- 1
  for (i in seq_along(drug_vector)) {
    drug <- drug_vector[i]
    col <- column_vector[i]
    concentration <- concentration_vector[i]
    for (j in seq_along(timepoints)) {
      row <- timepoint_starts[j]
      reference <- paste0(col, row)
      value <- extract_indices(reference, data)
      drug_concentration_df$Drug[row_counter] <- drug
      drug_concentration_df$Concentration[row_counter] <- concentration
      drug_concentration_df$Timepoint[row_counter] <- timepoints[j]
      drug_concentration_df$Value[row_counter] <- value
      drug_concentration_df$CellLine[row_counter] <- cell_line
      row_counter <- row_counter + 1
    }
  }
  return(na.omit(drug_concentration_df))
}
concentration_filter <- list(
  "AZD4547" = c(4000, 2000, 1000,0),
  "J-104129" = c(1000, 500, 250,0),
  "Varespladib" = c(4000, 2000, 1000,0),
  "Danusertib" = c(8000, 4000, 2000,0)
)

red_concentration_df <- build_drug_concentration_df(drug_vector, cellline, 
                                                     concentration_vector, 
                                                     timepoints = timepoints,
                                                     timepoint_starts = timepoint_starts, 
                                                     column_vector = live_index,
                                                     data = excel)
# red_concentration_df <- red_concentration_df %>%
#   dplyr::filter((Drug == "AZD4547" & Concentration %in% concentration_filter[["AZD4547"]]) |
#                   (Drug == "J-104129" & Concentration %in% concentration_filter[["J-104129"]]) |
#                   (Drug == "Varespladib" & Concentration %in% concentration_filter[["Varespladib"]]) | Drug == "DMSO" | 
#                   (Drug == "Danusertib" & Concentration %in% concentration_filter[["Danusertib"]]))


cytotox_concentration_df <- build_drug_concentration_df(drug_vector, cellline, 
                                                    concentration_vector, 
                                                    timepoints = timepoints,
                                                    timepoint_starts = timepoint_starts, 
                                                    column_vector = cytotox_index,
                                                    data = excel)
# cytotox_concentration_df <- cytotox_concentration_df %>%
#   dplyr::filter((Drug == "AZD4547" & Concentration %in% concentration_filter[["AZD4547"]]) |
#                   (Drug == "J-104129" & Concentration %in% concentration_filter[["J-104129"]]) |
#                   (Drug == "Varespladib" & Concentration %in% concentration_filter[["Varespladib"]]) | Drug == "DMSO" | 
#                   (Drug == "Danusertib" & Concentration %in% concentration_filter[["Danusertib"]]))

caspase_concentration_df <- build_drug_concentration_df(drug_vector, cellline, 
                                                    concentration_vector, 
                                                    timepoints = timepoints,
                                                    timepoint_starts = timepoint_starts, 
                                                    column_vector = caspase_index,
                                                    data = excel)
# caspase_concentration_df <- caspase_concentration_df %>%
#   dplyr::filter((Drug == "AZD4547" & Concentration %in% concentration_filter[["AZD4547"]]) |
#                   (Drug == "J-104129" & Concentration %in% concentration_filter[["J-104129"]]) |
#                   (Drug == "Varespladib" & Concentration %in% concentration_filter[["Varespladib"]]) | Drug == "DMSO" | 
#                   (Drug == "Danusertib" & Concentration %in% concentration_filter[["Danusertib"]]))

drugs = c("AZD4547", "J-104129", "Varespladib", "Danusertib")
#################
plots1 = c()
for(k in 1:length(drugs)){
  print(drugs[k])
  drug_number = k
  drug = drugs[k]
  df <-  red_concentration_df[red_concentration_df$Drug %in% c("DMSO", drug), ]
  df$Drug = drug
  df$Sample = paste0(df$Drug, "_", df$Concentration, "_", df$Timepoint)
  # summary_data <- df %>%
  #   group_by(Sample) %>%
  #   summarize(Mean = mean(Value), SD = sd(Value), Concentration = Concentration, Timepoint = Timepoint, Cell_Line = CellLine)
  summary_data <- df %>%
    group_by(Sample) %>%
    summarize(Mean = mean(Value), SEM = sd(Value) / sqrt(n()), Concentration = Concentration, Timepoint = Timepoint, Cell_Line = CellLine)
  summary_data <- summary_data %>%
    distinct(Sample, .keep_all = TRUE)
  drug_curve = summary_data
  
  plots1[[k]] <-  ggplot(drug_curve, aes(x = Timepoint, y = Mean, color = as.factor(Concentration), group = as.factor(Concentration))) +
    geom_line() +  geom_point() + 
    geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), width = 2) + 
    labs(
      x = NULL,
      y = "Mean Value",
      color = "Concentration",
      title = paste0("Live Cells: ", drugs[drug_number], " (", cellline, ")")
    ) +
    theme_bw() +
    scale_color_brewer(palette = "Set1") 
}

plots2 = c()
for(k in 1:length(drugs)){
  print(drugs[k])
  drug_number = k
  drug = drugs[k]
  df <-  cytotox_concentration_df[cytotox_concentration_df$Drug %in% c("DMSO", drug), ]
  df$Drug = drug
  df$Sample = paste0(df$Drug, "_", df$Concentration, "_", df$Timepoint)
  # summary_data <- df %>%
  #   group_by(Sample) %>%
  #   summarize(Mean = mean(Value), SD = sd(Value), Concentration = Concentration, Timepoint = Timepoint, Cell_Line = CellLine)
  summary_data <- df %>%
    group_by(Sample) %>%
    summarize(Mean = mean(Value), SEM = sd(Value) / sqrt(n()), Concentration = Concentration, Timepoint = Timepoint, Cell_Line = CellLine)
  summary_data <- summary_data %>%
    distinct(Sample, .keep_all = TRUE)
  drug_curve = summary_data
  
  plots2[[k]] <-  ggplot(drug_curve, aes(x = Timepoint, y = Mean, color = as.factor(Concentration), group = as.factor(Concentration))) +
    geom_line() +  geom_point() + 
    geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), width = 2) + 
    labs(
      x = NULL,
      y = "Mean Value",
      color = "Concentration",
      title = paste0("Dead Cells: ", drugs[drug_number], " (", cellline, ")")
    ) +
    theme_bw() +
    scale_color_brewer(palette = "Set1") 
}

plots3 = c()
for(k in 1:length(drugs)){
  print(drugs[k])
  drug_number = k
  drug = drugs[k]
  df <-  caspase_concentration_df[caspase_concentration_df$Drug %in% c("DMSO", drug), ]
  df$Drug = drug
  df$Sample = paste0(df$Drug, "_", df$Concentration, "_", df$Timepoint)
  # summary_data <- df %>%
  #   group_by(Sample) %>%
  #   summarize(Mean = mean(Value), SD = sd(Value), Concentration = Concentration, Timepoint = Timepoint, Cell_Line = CellLine)
  summary_data <- df %>%
    group_by(Sample) %>%
    summarize(Mean = mean(Value), SEM = sd(Value) / sqrt(n()), Concentration = Concentration, Timepoint = Timepoint, Cell_Line = CellLine)
  summary_data <- summary_data %>%
    distinct(Sample, .keep_all = TRUE)
  drug_curve = summary_data
  
  plots3[[k]] <-  ggplot(drug_curve, aes(x = Timepoint, y = Mean, color = as.factor(Concentration), group = as.factor(Concentration))) +
    geom_line() +  geom_point() + 
    geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), width = 2) + 
    labs(
      x = NULL,
      y = "Mean Value",
      color = "Concentration",
      title = paste0("Caspase 3/7: ", drugs[drug_number], " (", cellline, ")")
    ) +
    theme_bw() +
    scale_color_brewer(palette = "Set1") 
}
p0 <- grid.arrange(
  plots1[[1]], plots1[[2]], plots1[[3]], plots1[[4]], 
  plots2[[1]], plots2[[2]], plots2[[3]], plots2[[4]],
  plots3[[1]], plots3[[2]], plots3[[3]], plots3[[4]],ncol = 4
)
ggsave(p0, filename = "Curves_Caspase_1769_SE.pdf", width = 17, height = 15/2)



