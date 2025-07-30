# ============================================================================ #
# Project: LSARP 
# Author: Thi Mui Pham, tmpham@hsph.harvard.edu
# Title: Figure 1C, Stacked Phenotypic diversity 
# Organisms: all five
# ============================================================================ #
remove(list=ls())
setwd("/Users/tm-pham/academia/hsph/lsarp/")

# Output paths
path <- "/Users/tm-pham/academia/hsph/lsarp/publications/lsarp_epi_genomics/" 
output_path <- "/Users/tm-pham/academia/hsph/lsarp/publications/epi_genomics/data_output/"

# Load functions
source(paste0(path, "code/plotting_template.R"))
source(paste0(path, "code/lsarp_function_diversity_index.R"))
source(paste0(path, "code/lsarp_function_stacked_resistance.R"))

# Load packages
library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(lubridate)

# Load data
load("data/calgary_population_estimates/lsarp_calgary_population_ageSex.RData")
load("data/LSARP_drug_bug_dictionary.RData")

load("data/figure1/lsarp_figure1C_input_df_stacked_amr.RData")

################################################################################
# PLOT
################################################################################
abx_colors <- c("#0B775E", "#FFEBAE", "#F5C600", "#DF8F44", "#D8460B", "#800000","#C51A4E")
(n_abx_perc_plot <- ggplot(df_stacked_amr, 
                           aes(x=ymd(YEAR, truncated=2), y = value, fill = variable)) + 
    facet_grid(cols = vars(organismofinterest)) + 
    geom_area(stat="identity", position = position_fill(reverse = TRUE)) + 
    labs(y="Percentage of 30-day index isolates", fill = "Antibiotic susceptibility") + 
    theme_template_white() +
    scale_fill_manual(values=abx_colors) +
    scale_x_date(date_labels = "%Y", breaks=seq(as.Date("2006-01-01"), 
                                                as.Date("2022-12-31"), by = "4 years")) + 
    scale_y_continuous(labels=scales::percent) + 
    theme(legend.position = "bottom", 
          panel.spacing = unit(2, "lines"),
          axis.title.x = element_blank(), 
          axis.text.x = element_text(size = 26),
          axis.text.y = element_text(size = 26), 
          axis.title.y = element_text(size =28), 
          legend.title = element_text(size = 30),
          legend.text = element_text(size = 28),
          strip.text = element_text(size = 30, face = "bold.italic"), 
          plot.margin = margin(10, 20, 10, 10)))
ggsave(n_abx_perc_plot, file=paste0(path, "figures/lsarp_figure1C.pdf"),
       width=42, height=8)
ggsave(n_abx_perc_plot, file=paste0(path, "figures/lsarp_figure1C.png"),
       width=42, height=8)

