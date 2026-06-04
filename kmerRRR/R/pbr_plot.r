#!/usr/bin/env Rscript
rm(list = ls())

if (requireNamespace("tidyverse", quietly = T)){
  library(tidyverse)
} else {
  install.packages("tidyverse")
  library(tidyverse)
}
if (requireNamespace("ggplot2", quietly = T)){
  library(ggplot2)
} else {
  install.packages("ggplot2")
  library(ggplot2)
}
if (requireNamespace("data.table", quietly = T)){
  library(data.table)
} else {
  install.packages("data.table")
  library(data.table)
}

arguments <- commandArgs(trailingOnly = T)

if (length(arguments) > 0){
  argument1 <- arguments[1]
  argument2 <- arguments[2]
  argument3 <- arguments[3]
  if (!dir.exists(argument3)) {
    dir.create(argument3, recursive = T)
  }
  cat("Per base file: ", argument1, "\n")
  cat("sliding window ", argument2, "\n")
  cat("Output file prefix: ", argument3, "\n")
}

plotting_pbkr_data <- function(pbr_file_path, slide_window, output){
  
  pbr_data <- fread(pbr_file_path, sep = "\t", header = T)
  
  contig_names <- unique(pbr_data$contig)
  
  for (cname in contig_names){
    df_list <- paste0("filtered_", cname)
    filtered_data <- pbr_data %>%
      filter(contig == cname)
    filtered_data <- filtered_data %>%
      mutate(Type = ifelse(repetitive == 1 & lunique ==1, "LocalUnique", "Repetitive")) %>%
      mutate(Type = ifelse(repetitive == 0 & lunique ==1, "Unique", Type))
    filtered_data$Type <- factor(filtered_data$Type, levels = c("Unique", "Repetitive", "LocalUnique"))
    filtered_data <- filtered_data %>%
      mutate(Bins = floor(position/as.numeric(slide_window))*as.numeric(slide_window))
    filtered_data <- filtered_data %>%
      group_by(contig, Bins, Type) %>%
      summarise(Count = n(), .groups = 'drop')
    color_pal <- data.frame("Type" = c("LocalUnique", "Repetitive", "Unique"), Color = c("darkblue", "grey", "seagreen"))
    vv <- color_pal$Color
    names(vv) = color_pal$Type
    output_dir <- output
    pdf_path <- file.path(output_dir, paste0(df_list, "_kmerRRR.pdf"))
    pdf(file = pdf_path, width = 12, height = 8)
    print(
      ggplot(filtered_data, aes(x = Bins, y = Count, fill = Type)) +
        geom_bar(stat = 'identity', position = "fill") + 
        theme_minimal() + 
        scale_fill_manual(values = vv) + 
        labs(title = paste0("Distribution of Local Unique Character across the locus ", cname), x = paste0("Position ", "(", slide_window, " binned)"), y = "Counts per bin (%)", fill = "Region Type") + 
        theme(plot.title = element_text(size = 20, face = 'bold'), axis.title = element_text(size = 16), axis.text = element_text(size =12))
  )
    dev.off()
  }
}

plotting_pbkr_data(pbr_file_path = argument1, slide_window = argument2, output = argument3)