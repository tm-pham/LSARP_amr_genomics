# ============================================================================ #
# Project: LSARP 
# Author: Thi Mui Pham, tmpham@hsph.harvard.edu
# Title: Figure 2 A
# Organisms: S. aureus
# ---------------------------------------------------------------------------- #
# Offset: 
# - Calgary population for community-onset
# - Patient days for hospital-onset
# ============================================================================ #
remove(list=ls())
# Set working directory here
setwd("/Users/tm-pham/academia/hsph/LSARP_amr_genomics/")

# Load packages
library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(gridExtra)

# Load functions and templates
source("code/functions_templates/plotting_template.R")

# Load data
load("data/figure2/lsarp_figure2A_SA_input_df_plot_COI.RData")
load("data/figure2/lsarp_figure2A_SA_input_df_plot_HOI.RData")

################################################################################
# Figure 2A
ab_colors <- c("#1B9E77", "#B2DF8A", "#E6AB02", "#D95F02" ,"#E7298A")
strip_text_size <- 34
axis_text_x <- 24
axis_title_y <- 27
legend_text_size = 28


################################################################################
# Plot for COMMUNITY-ONSET

(plot1 <- ggplot(df_plot %>% 
                   filter(!is.na(cluster_map), HOSPITAL_ONSET_48H == "Community-onset", n_strain_ab > 30) %>% 
                   mutate(cluster_map = factor(cluster_map, levels = cluster_map_levels)),
                 aes(x=ymd(YEAR, truncated=2), y=inc_std, color=ab_name, fill = ab_name, group = ab_name)) + 
   facet_grid(cols=vars(cluster_map), rows=vars(HOSPITAL_ONSET_48H), scales = "free_y") + 
   geom_vline(xintercept=as.numeric(as.Date("2020-01-01", format="%Y-%m-%d")), 
              linetype="dashed", color="darkgrey", linewidth=1.5) + 
   # geom_line(linetype = "dashed", linewidth = 1.0) +
   geom_smooth(span = 1.5) +
   geom_point(shape = 21, color = "black", size = 4.5, stroke = 2) + 
   labs(y="Number of index isolates\n(per 100,000 Calgary population)",
        fill="Phenotypes\ncloxacillin - ciprofloxacin - erythromycin - clindamycin",
        color="Phenotypes\ncloxacillin - ciprofloxacin - erythromycin - clindamycin") + 
   
   scale_x_date(date_labels = "%Y", breaks = seq.Date(from = as.Date("2006-01-01"), 
                                                      to = as.Date("2022-01-31"), by = "4 years")) + 
   scale_y_continuous(limits = c(-1.5, 4), expand = c(0, 0)) +
   coord_cartesian(ylim = c(0, 4)) + 
   scale_color_manual(values = ab_colors) +
   scale_fill_manual(values = ab_colors) +
   theme_template_white() + 
   guides(fill = guide_legend(ncol = 1)) + 
   theme(legend.position = "top", 
         legend.direction = "vertical",
         legend.text = element_text(size = legend_text_size),
         panel.spacing = unit(1.2, "cm"),
         strip.text = element_text(size = strip_text_size, face = "bold"),
         legend.title = element_blank(), 
         axis.title.x = element_blank(), 
         axis.title.y = element_text(size = axis_title_y),
         axis.text.x = element_text(size = axis_text_x, angle = 0)))
# ggsave(plot1, file="figures/figure2/lsarp_figure2A.pdf", 
#        width = 36, height = 8)

################################################################################
# Plot for HOSPITAL-ONSET
################################################################################
(plot2 <- ggplot(df_plot_HOI %>% 
                   filter(!is.na(cluster_map), HOSPITAL_ONSET_48H == "Hospital-onset", 
                          n_strain_ab > 30), 
                 aes(x=ymd(YEAR, truncated=2), y=inc, color = ab_name, fill=ab_name, group = ab_name)) + 
   facet_grid(cols=vars(cluster_map), rows=vars(HOSPITAL_ONSET_48H), scales = "free_y") + 
   geom_vline(xintercept=as.numeric(as.Date("2020-01-01", format="%Y-%m-%d")), 
              linetype="dashed", color="darkgrey", linewidth=1.3) + 
   geom_line(linetype = "dashed", linewidth = 1.0) +
   geom_smooth(span = 1.5) +
   geom_point(shape = 21, color = "black", size = 4.5, stroke = 2) + 
   labs(y="Number of index isolates\n(per 100,000 patient days)",
        fill="StrainGST reference\n(clonal complex)",
        color="StrainGST reference\n(clonal complex)") + 
   
   scale_x_date(date_labels = "%Y", breaks = seq.Date(from = as.Date("2006-01-01"), 
                                                      to = as.Date("2022-01-31"), by = "4 years")) + 
   
   scale_y_continuous(limits = c(-1.5, 5), expand = c(0, 0)) +
   coord_cartesian(ylim = c(0, 5)) + 
   scale_color_manual(values = ab_colors) +
   scale_fill_manual(values = ab_colors) +
   theme_template_white() + 
   theme(legend.position = "none",
         legend.title = element_blank(),
         panel.spacing = unit(1.2, "cm"),
         strip.text = element_text(size = strip_text_size, face = "bold"), 
         axis.title.x = element_blank(),  
         axis.title.y = element_text(size = axis_title_y),
         axis.text.x = element_text(size = axis_text_x, angle = 0)))
# ggsave(plot2, file="figures/figure2/lsarp_figure2A_2.pdf", 
#        width = 36, height = 6)


################################################################################
# Combined plot
(plot <- grid.arrange(plot1, plot2, ncol = 1, heights = c(1, 0.8)))
ggsave(plot, file = "figures/figure2/lsarp_figure2A.pdf", 
       width = 36, height = 14)
ggsave(plot, file = "figures/figure2/lsarp_figure2A.png",
       width = 36, height = 14, dpi = 300)
