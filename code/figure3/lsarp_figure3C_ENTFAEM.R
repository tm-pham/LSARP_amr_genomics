# ============================================================================ #
# Project: LSARP
# Author: Thi Mui Pham, tmpham@hsph.harvard.edu
# Organism: Enterococcus faecium
# Title: Figure 3C, ST117 plot
# ============================================================================ #
remove(list=ls())
# Set working directory here
setwd("/Users/tm-pham/academia/hsph/LSARP_amr_genomics/")

# Load package 
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(dplyr)
library(lubridate)

# Load functions and templates
source("code/functions_templates/plotting_template.R")

# Load data
load("data/figure3/lsarp_figure3C_input_df_st117_HOI_overall_yr.RData")
load("data/figure3/lsarp_figure3C_input_df_st117_HOI_yr.RData")

################################################################################
# Variables
axis_text_x <- 26
strip_text_size <- 30

################################################################################
# Overall HOI plot with segmented Poisson regression
################################################################################
dataset <- df_st117_HOI_overall_yr %>% 
  filter(HOSPITAL_ONSET_48H=="Hospital-onset") %>% 
  mutate(ST_name = "E8328 (ST117)")

model<- glm(n ~ YEAR + offset(log(total_ptdays)), 
            data = dataset, 
            family = poisson())

davies.test(model, seg.Z = ~YEAR, k = 10)
(seg_model <- segmented::segmented(model, npsi = 1))

# Extract the breakpoints
(breaks <- seg_model$psi[, 2])

(plot1 <- ggplot(data = dataset, aes(x = ymd(YEAR, truncated=2), y = inc)) +
    facet_wrap(.~cluster_name) + 
    geom_point(color = "black", shape = 21, stroke = 1.5, size = 4, fill = "white") +
    geom_smooth(data = dataset %>% filter(YEAR<breaks), color = "black", size = 2, 
                method = "glm", method.args = list(family = "poisson")) +
    geom_smooth(data = dataset %>% filter(YEAR>=breaks), color = "black", size = 2, 
                method = "glm", method.args = list(family = "poisson")) +
    scale_y_continuous(limits = c(0, NA)) + 
    scale_x_date(date_labels = "%Y", breaks = seq.Date(from = as.Date("2006-01-01"), 
                                                       to = as.Date("2022-12-31"), by = "2 years"), 
                 limits = c(as.Date("2006-01-01"), as.Date("2022-12-31"))) + 
    labs(y="Number of index isolates\n(per 100,000 patient days)") + 
    theme_template_white() + 
    theme(
          axis.title.x = element_blank(), 
          axis.title.y = element_text(size=26, face="plain"),
          axis.text.x = element_text(size=axis_text_x, face="plain"),
          axis.text.y = element_text(size=22, face="plain"),
          strip.text = element_text(size = strip_text_size, face = "bold.italic")))

################################################################################
# Plot of ST117 by antibiogram
################################################################################
(plot2 <- ggplot(df_st117_HOI_yr %>%
                   filter(HOSPITAL_ONSET_48H=="Hospital-onset", 
                          ab_name %in% c("S-S-S", "R-S-S", "S-R-S", 
                                         "S-R-R", "R-R-S", "R-R-R")) %>% 
                   mutate(ab_name = factor(ab_name, levels =  c("S-S-S", "R-S-S", "S-R-S", 
                                                                "S-R-R", "R-R-S", "R-R-R"))), 
                 aes(x = ymd(YEAR, truncated = 2), 
                     y = inc, fill=ab_name)) + 
   facet_wrap(.~ab_name) + 
   geom_bar(stat="identity", color = "black", size = 1.2) +
   labs(y="") + 
   scale_y_continuous(limits= c(0, NA)) +
   scale_fill_manual(values= c("#1B9E77", "#E7298A", "#004B40", "#D55E00", "#EFC000FF", "#660066")) +
   scale_x_date(date_labels="%Y", breaks=seq(as.Date("2006-01-01", format = "%Y-%m-%d"), 
                                             as.Date("2022-12-31", format = "%Y-%m-%d"), 
                                             by = "4 years")) + 
   theme_template_white() + 
   theme(axis.text.x = element_text(size=18), 
         panel.spacing = unit(1, "cm"),
         strip.text = element_text(size=strip_text_size, face = "bold"),
         axis.title.y = element_text(size=22, face="plain"),
         axis.text.y = element_text(size=20, face="plain"),
         axis.title.x = element_blank(),
         legend.position = "none", 
         plot.background = element_rect(fill = "transparent", colour = NA)))

################################################################################
# Combined plot
lay <- rbind(c(1,1,1,1,1,2,2,2,2), c(1,1,1,1,1,2,2,2,2))
plot <- grid.arrange(plot1, plot2,
                     layout_matrix = lay)
ggsave(plot, file = "figures/figure3/lsarp_figure3C.pdf", 
       width = 24, height=7.5, bg = "transparent")

