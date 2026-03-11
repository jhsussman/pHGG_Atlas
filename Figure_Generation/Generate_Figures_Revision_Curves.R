library(ggplot2)
library(Seurat)
library(ggpubr)
library(readxl)
library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(gridExtra)
setwd("G:/My Drive/University of Pennsylvania/Labs/Tan Lab/HTAN/Drug_Screening/Figures")

drugs = c("Panobinostat", "AZD4547", "Temolozolomide", "Belnacasan", "Navitoclax", "Venetoclax", 
          "J-104129", "Darifenacin", "Cevimeline", "Varespladib", "Retigabine", 
          "4-Aminopyridine", "Gabazine", "CGP52432", "Entrectinib", "Trametinib", 
          "L-741-626", "Quinipirole", "Ruxolitinib", "Danusertib", "1-Naphthyl PP1", 
          "SB366791", "Verbascoside", "Y-27632")
starting_concs = c(100,4000,200000,80000,2000,40000,100000,40000,20000,4000,
                   10000,100000,100000,4000,8000,80,4000,4000,4000,8000,20000, 
                   10000,40000,40000)

#1
excel = read_excel("1763NR.384well correct 2.xlsx", col_names = F, sheet = "Final")
timepoints = c(0,6,12,18,24,30,36,42,48,54,60,66,72)
timepoint_starts = 56:68
cellline = "1763NR"

#2
excel = read_excel("913NR 384well correct.xlsx", col_names = F, sheet = "Final")
timepoints = c(0,6,12,18,24,30,36,42,48,54,60,66,72)
timepoint_starts = 56:68
cellline = "913NR"

#3
excel = read_excel("195NR.384well correct 2.xlsx", col_names = F, sheet = "Final")
timepoints = c(0,6,12,18,24,30,36,42,48,54,60,66)
timepoint_starts = 56:67
cellline = "195NR"

#4
excel = read_excel("1769NR 384well correct 2.xlsx", col_names = F, sheet = "Final")
timepoints = c(0,6,12,18,24,30,36,42,48,54,60,66,72)
timepoint_starts = 56:68
cellline = "1769NR"

#5
excel = read_excel("3058NR.384well correct.xlsx", col_names = F, sheet = "Final")
timepoints = c(0,72)
timepoint_starts = c(56,68)
cellline = "3058NR"

#################

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

increment_excel_column <- function(column, times = 1) {
  for (i in 1:times) {
    n <- nchar(column)
    if (all(strsplit(column, NULL)[[1]] == "Z")) {
      column <- paste0(rep("A", n + 1), collapse = "")
    } else {
      column_list <- strsplit(column, NULL)[[1]]
      carry <- TRUE
      for (j in length(column_list):1) {
        if (carry) {
          if (column_list[j] == "Z") {
            column_list[j] <- "A"
          } else {
            column_list[j] <- LETTERS[which(LETTERS == column_list[j]) + 1]
            carry <- FALSE
          }
        }
      }
      if (carry) {
        column_list <- c("A", column_list)
      }
      column <- paste0(column_list, collapse = "")
    }
  }
  return(column)
}

create_data_frame <- function(refs, DMSO_refs, extract_indices, excel, timepoint, cell_line) {
  samples <- c()
  values <- c()
  dms_sample_values <- sapply(DMSO_refs, extract_indices, excel)
  samples <- c(samples, rep("DMSO", length(DMSO_refs)))
  values <- c(values, dms_sample_values)
  
  for (sample_name in names(refs)) {
    refs_list <- refs[[sample_name]]
    samples <- c(samples, rep(sample_name, length(refs_list)))
    calculated_values <- sapply(refs_list, extract_indices, excel)
    values <- c(values, calculated_values)
  }
  df <- data.frame(Sample = samples, 
                   Value = values, 
                   Cell_Line = cell_line, 
                   Timepoint = timepoint, 
                   stringsAsFactors = FALSE)
  return(na.omit(df))
}

generate_concentrations <- function(starting_conc) {
  return(c(0, starting_conc / 2^(0:6)))
}
create_drug_concentration_variables <- function(drugs, starting_concs) {
  variables <- list()
  for (i in seq_along(drugs)) {
    drug <- drugs[i]
    starting_conc <- starting_concs[i]
    concs <- generate_concentrations(starting_conc)
    for (conc in concs) {
      var_name <- paste0(drug, "_Conc_", conc, "_Refs")
      variables[[var_name]] <- conc
    }
  }
  return(variables)
}

DMSO_prefs <- c("C","D", "S","T", "AI","AJ",
                      "AY","AZ", "BO","BP", "CE",
                      "CF", "CU","CV", "DK","DL", 
                      "EA","EB", "EQ","ER", "FG","FH", 
                      "FW","FX", "GM","GN", "HC",
                      "HD", "HS","HT", "II","IJ", 
                      "IY","IZ", "JO","JP", "KE","KF", 
                      "KU","KV", "LK","LL", "MA",
                      "MB", "MQ","MR", "NG","NH")
#################
plots = c()
for(k in 1:length(drugs)){
  print(drugs[k])
  drug_number = k
  drug_curve = c()
  
  if(drugs[k] %in% c("AZD4547", "J-104129", "Varespladib", "Danusertib")){
    next
  }
  
  for(j in 1:length(timepoints)){
    start_letter <- "C"
    start_number <- timepoint_starts[j]
    timepoint = timepoints[j]
    refs <- list()
    
    for (drug in drugs) {
      for (i in 1:8) {
        var_name <- paste(drug, "Conc", i, "Refs", sep = "_")
        first_ref <- paste0(start_letter, start_number)
        second_ref <- paste0(increment_excel_column(start_letter), start_number)
        refs[[var_name]] <- c(first_ref, second_ref)
        start_letter <- increment_excel_column(start_letter, 2)
      }
    }
    
    drug_concentration_variables <- create_drug_concentration_variables(drugs, starting_concs)
    
    names(refs) = names(drug_concentration_variables)
    refs_df = t(as.data.frame(refs))
    
    DMSO_refs = paste0(DMSO_prefs, as.character(start_number))
    df <- create_data_frame(refs, DMSO_refs, extract_indices, excel, timepoint, cellline)
    
    summary_data <- df %>%
      group_by(Sample) %>%
      summarize(Mean = mean(Value), SEM = sd(Value) / sqrt(n()), Timepoint = timepoint, Cell_Line = cellline)
    # summary_data <- df %>%
    #   group_by(Sample) %>%
    #   summarize(Mean = mean(Value), SD = sd(Value), Timepoint = timepoint, Cell_Line = cellline)
    order = unique(df$Sample)
    df$Cell_Line = factor(df$Cell_Line)
    
    DMSO_df = subset(summary_data, grepl("DMSO", Sample))
    DMSO_df$Concentration = 0
    drug_df <- subset(summary_data, grepl(drugs[drug_number], Sample))
    drug_df <- drug_df[-c(1), ]
    drug_df$Concentration <- sapply(drug_df$Sample, function(x) {
      parts <- unlist(strsplit(x, "_"))
      return(parts[length(parts) - 1])
    })
    
    full_df = rbind(DMSO_df, drug_df)
    drug_curve = rbind(drug_curve, full_df)
  }
  head(drug_curve)
  drug_curve$Timepoint <- as.numeric(drug_curve$Timepoint)
  drug_curve <- drug_curve %>%
    mutate(Concentration = factor(Concentration, levels = sort(unique(as.numeric(Concentration)))))
  
  plots[[k]] <-  ggplot(drug_curve, aes(x = Timepoint, y = Mean, color = Concentration, group = Concentration)) +
    geom_line() +  geom_point() + 
    geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), width = 2) + 
    labs(
      x = NULL,
      y = "Fluorescence",
      color = "Concentration",
      title = paste0(drugs[drug_number], ": ", cellline)
    ) +
    theme_bw() +
    scale_color_brewer(palette = "Set1") 
}

# p0 = grid.arrange(
#   plots[[1]], plots[[2]], plots[[3]], plots[[4]],
#   plots[[5]], plots[[6]], plots[[7]], plots[[8]],
#   plots[[9]], plots[[10]], plots[[11]], plots[[12]],
#   plots[[13]], plots[[14]], plots[[15]], plots[[16]],
#   plots[[17]], plots[[18]], plots[[19]], plots[[20]],
#   plots[[21]], plots[[22]], plots[[23]], plots[[24]],
#   ncol = 4
# )

p0 = grid.arrange(
  plots[[1]], plots[[3]], plots[[4]],
  plots[[5]], plots[[6]], plots[[8]],
  plots[[9]],  plots[[11]], plots[[12]],
  plots[[13]], plots[[14]], plots[[15]], plots[[16]],
  plots[[17]], plots[[18]], plots[[19]], 
  plots[[21]], plots[[22]], plots[[23]], plots[[24]],
  ncol = 4
)
ggsave(p0, filename = "Curves_All_3058_SE_Supp.pdf", width = 17, height = 15 * (5/6))


