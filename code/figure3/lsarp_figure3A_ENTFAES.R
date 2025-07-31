# ============================================================================ #
# Project: LSARP
# Author: Thi Mui Pham, tmpham@hsph.harvard.edu
# Title: Figure 3A
# Organism: Enterococcus faecalis
# ---------------------------------------------------------------------------- #
# Phenotype resistance profile for E faecalis 
# ============================================================================ #
remove(list=ls())
# Set working directory here
setwd("/Users/tm-pham/academia/hsph/LSARP_amr_genomics/")

# Load packges
library(dplyr)
library(stringr)
library(ggplot2)
library(lubridate)
library(cowplot)
library(ggh4x)
library(gridExtra)

# Load functions and templates
source("code/functions_templates/plotting_template.R")

# Load data
load("data/figure3/lsarp_figure3A_input_df_mlst_ab_COI_yr.RData")
load("data/figure3/lsarp_figure3A_input_df_mlst_ab_HOI_yr.RData")

################################################################################
# Figure 3A
ab_colors <- c("#1B9E77", "#D95F02" ,"#E7298A")
strip_text_size <- 34
axis_text_x <- 24
axis_title_y <- 28

################################################################################
# COMMUNITY-ONSET
(plot1.1 <- ggplot(df_mlst_ab_COI_yr %>% 
                     filter(ab_name%in%c("R-S-S", "R-R-S", "R-R-R")) %>% 
                     filter(cluster_name%in%c("1 (ST 179)", "3 (ST 16)", "2 (ST 40)","6 (ST 64)", "4 (ST 778/103)", "5 (ST 6)")),
                   aes(x=ymd(YEAR, truncated=2), y=inc, fill=ab_name, color = ab_name)) + 
    facet_grid(rows = vars(HOSPITAL_ONSET_48H), cols=vars(cluster_name)) + 
    geom_line(linetype = "dashed", linewidth = 1.0) +
    stat_smooth(method = "lm", formula = y ~ x , size = 3, expand = c(0, 0)) +
    geom_point(shape = 21, color = "black", size = 4.5, stroke = 2) + 
    labs(y="Number of index isolates\n(per 100,000 Calgary residents)", 
         fill = "", 
         color = "") + 
   scale_x_date(date_labels = "%Y", breaks = seq.Date(from = as.Date("2006-01-01"), 
                                                      to = as.Date("2022-12-31"), by = "4 years")) + 
    scale_y_continuous(limits = c(-0.2, 1.5), expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, 1.5)) + 
    scale_fill_manual(values = ab_colors) +
    scale_color_manual(values = ab_colors) +
    theme_template_white()+ 
    theme(legend.position = "top",
          legend.direction = "vertical", 
          legend.text = element_text(size = strip_text_size),
          legend.title = element_blank(),
          panel.spacing = unit(1.2, "cm"),
          strip.text.x = element_text(size = strip_text_size, face = "bold"),
          strip.text.y = element_text(size = strip_text_size-1, face = "bold"),
          axis.title.x = element_blank(), 
          axis.title.y = element_text(size = axis_title_y),
          axis.text.y = element_text(size = axis_text_x),
          axis.text.x = element_text(size = axis_text_x, angle = 0)))

################################################################################
# HOSPITAL-ONSET
(plot1.2 <- ggplot(df_mlst_ab_HOI_yr %>% 
                     filter(ab_name%in%c("R-S-S", "R-R-S", "R-R-R")) %>% 
                     filter(cluster_name%in%c("1 (ST 179)", "3 (ST 16)", "2 (ST 40)", "4 (ST 778/103)", "6 (ST 64)", "5 (ST 6)")), 
                   aes(x=ymd(YEAR, truncated=2), y=inc, fill=ab_name, color = ab_name)) + 
    facet_grid(rows = vars(HOSPITAL_ONSET_48H), cols=vars(cluster_name)) + 
    geom_line(linetype = "dashed", linewidth = 1.0) +
    stat_smooth(method = "lm", formula = y ~ x , size = 3, span = 1, fullrange=T, level = 0.95) +
    # geom_smooth(span = 2) +
    geom_point(shape = 21, color = "black", size = 4.5, stroke = 2) + 
    labs(y="Number of index isolates\n(per 100,000 patient days)", 
         fill = "", 
         color = "") + 
   scale_x_date(date_labels = "%Y", breaks = seq.Date(from = as.Date("2006-01-01"), 
                                                      to = as.Date("2022-12-31"), by = "4 years")) + 
    scale_y_continuous(limits = c(-0.4, 1.5), expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, 1.5)) + 
    scale_fill_manual(values = ab_colors) +
    scale_color_manual(values = ab_colors) +
    theme_template_white()+ 
    theme(legend.position = "none", 
          panel.spacing = unit(1.2, "cm"),
          strip.text.x = element_text(size = strip_text_size, face = "bold"),
          strip.text.y = element_text(size = strip_text_size-1, face = "bold"),
          axis.title.x = element_blank(), 
          axis.title.y = element_text(size = axis_title_y),
          axis.text.y = element_text(size = axis_text_x),
          axis.text.x = element_text(size = axis_text_x, angle = 0)))

################################################################################
# For ST 179 HOI, RRS and RSS phenotypes only 
# Segmented regression
(plot1.3 <- ggplot(df_mlst_ab_HOI_yr %>% 
                     filter(ab_name%in%c("R-S-S", "R-R-S", "R-R-R")) %>% 
                                          filter(cluster_name%in%c("1 (ST 179)")), 
                   aes(x=ymd(YEAR, truncated=2), y=inc, fill=ab_name, color = ab_name)) + 
    facet_grid(cols=vars(cluster_name)) + 
    geom_line(linetype = "dashed", linewidth = 1.0) +
    geom_smooth(data = df_mlst_ab_HOI_yr %>% 
                      filter(ab_name%in%c("R-S-S")) %>% 
                      filter(cluster_name%in%c("1 (ST 179)")), 
                method = "lm", formula = y ~ x , size = 3, span = 1, fullrange=T, level = 0.95) +
    stat_smooth(data = df_mlst_ab_HOI_yr %>% 
                      filter(ab_name%in%c("R-R-S"), cluster_name%in%c("1 (ST 179)"), YEAR<2015), 
                    method = "lm", size = 3, span = 1,level = 0.95) +
    stat_smooth(data = df_mlst_ab_HOI_yr %>% 
                  filter(ab_name%in%c("R-R-S"), cluster_name%in%c("1 (ST 179)"), YEAR>=2015), 
                method = "lm", size = 3, span = 1,  level = 0.95) +
    
    # geom_smooth(span = 2) +
    geom_point(shape = 21, color = "black", size = 4.5, stroke = 2) + 
    labs(y="Number of index isolates\n(per 100,000 patient days)", 
         fill = "", 
         color = "") + 
    scale_x_date(date_labels = "%Y", breaks = seq.Date(from = as.Date("2006-01-01"), 
                                                       to = as.Date("2022-12-31"), by = "4 years")) + 
    scale_y_continuous(limits = c(-0.4, 1.5), expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, 1.5)) + 
    scale_fill_manual(values = ab_colors) +
    scale_color_manual(values = ab_colors) +
    theme_template_white()+ 
    theme(legend.position = "none", 
          panel.spacing = unit(1.2, "cm"),
          plot.margin = margin(0, 0.5, 0, 0, "cm"),
          strip.text.x = element_text(size = strip_text_size, face = "bold"),
          strip.text.y = element_text(size = strip_text_size-1, face = "bold"),
          axis.title.x = element_blank(), 
          axis.title.y = element_text(size = 32),
          axis.text.y = element_text(size = 26),
          axis.text.x = element_text(size = 26, angle = 0)))
ggsave(plot1.3, file = "figures/figure3/lsarp_figure3A_ST179_HOI.pdf", 
       width = 7.65, height = 6.5)

################################################################################
# Combined plot
(plot <- grid.arrange(plot1.1, plot1.2, ncol = 1, heights = c(1, 0.85)))
ggsave(plot, file = "figures/figure3/lsarp_figure3A.pdf", 
       width = 36, height = 13)
ggsave(plot, file = "figures/figure3/lsarp_figure3A.png",
       width = 36, height = 13, dpi = 300)
