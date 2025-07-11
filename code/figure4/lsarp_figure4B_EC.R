# ============================================================================ #
# Project: LSARP 
# Author: Thi Mui Pham, tmpham@hsph.harvard.edu
# Title: Figure 4B (ST1193)
# Organisms: E. coli
# ---------------------------------------------------------------------------- #
# ST1193
# Community-onset 
# Hospital-onset 
# ============================================================================ #
remove(list=ls())
# Set working directory here
setwd("/Users/tm-pham/academia/hsph/LSARP_amr_genomics/")


# Load libraries
library(dplyr)
library(reshape2)
library(MASS)
library(lubridate)
library(segmented)
library(gridExtra)

# Load functions
source("code/functions_templates/plotting_template.R")

# Load data 
load("data/figure4/lsarp_figure4A_input_df_mlst_cipro_ceftr_COI_yr.RData")
load("data/figure4/lsarp_figure4A_input_df_mlst_cipro_ceftr_HOI_yr.RData")

################################################################################
# Global variables
point_size <- 4.5
line_size <- 2.5
axis_text_x <- 24
axis_title_y <- 26
strip_text_size <- 32


################################################################################
# ST1193, COMMUNITY-ONSET
################################################################################
# R-R
dataset <- df_mlst_cipro_ceftr_COI_yr %>% filter(mlst==1193, ab_name=="R-R") %>% 
  right_join(as.data.frame(cbind(YEAR = seq(2006, 2022)))) %>% 
  mutate(ab_name = ifelse(is.na(ab_name), "R-R", ab_name), 
         mlst = 1193)
model <- glm(n ~ YEAR + offset(log(n_year)), 
             data = dataset, 
             family = poisson())
(seg_model <- segmented::segmented(model, npsi = 1))
(breaks <- seg_model$psi[, 2]) # Extract the breakpoints


(plot1 <- ggplot(data = dataset, aes(x = ymd(YEAR, truncated=2), y = 100000*n/n_year, offset = n_year)) +
    facet_wrap(.~ab_name) + 
    geom_line(linetype = "dashed", linewidth = 1.0, color = "#CD0BBC") +
    geom_point(shape = 21, color = "black", size = point_size, stroke = 2, fill = "#CD0BBC") + 
    geom_smooth(data = dataset, color = "#CD0BBC", fill = "#CD0BBC", size = line_size, 
                method = "glm", method.args = list(family = "poisson")) +
    scale_y_continuous(limits = c(0, 7)) + 
    scale_x_date(date_labels = "%Y", breaks = seq.Date(from = as.Date(paste0("2006-01-01")), 
                                                       to = as.Date(paste0("2022-12-31")), by = "4 years")) + 
    labs(x = "Calendar year", y = "Number of index isolates\n(per 100,000 Calgary residents)") +
    theme_template_white() + 
    theme(axis.title.x = element_blank(), 
          strip.text = element_text(size = strip_text_size, face = "bold"),
          axis.text.x= element_text(size = axis_text_x),
          axis.title.y= element_text(size = axis_title_y)))


# ------------------------------------------------------------------------------
# R-S
# ------------------------------------------------------------------------------
dataset <- df_mlst_cipro_ceftr_COI_yr %>% 
  filter(mlst==1193, ab_name=="S-R") %>%
  right_join(as.data.frame(cbind(YEAR = seq(2006, 2022)))) %>% 
  mutate(ab_name = ifelse(is.na(ab_name), "S-R", ab_name), 
         mlst = 1193, 
         HOSPITAL_ONSET_48H = "Community-onset")
model <- glm(n ~ YEAR + offset(log(n_year)), 
             data = dataset, 
             family = poisson())
(seg_model <- segmented::segmented(model, npsi = 1))
(seg_model <- selgmented(model, seg.Z = ~YEAR, type = "bic", Kmax = 5, check.dslope=T))
(breaks <- seg_model$psi[, 2]) # Extract the breakpoints


(plot2 <- ggplot(data = dataset, aes(x = ymd(YEAR, truncated=2), y = 100000*n/n_year, offset = n_year)) +
    facet_grid(rows = vars(HOSPITAL_ONSET_48H), cols = vars(ab_name)) + 
    geom_line(linetype = "dashed", linewidth = 1.0, color = "#003399") +
    geom_point(shape = 21, color = "black", size = point_size, stroke = 2, fill = "#003399") +
    geom_smooth(data = dataset %>% filter(YEAR<=breaks), color = "#003399", fill = "#003399", size = line_size, 
                method = "glm", method.args = list(family = "poisson")) +
    geom_smooth(data = dataset %>% filter(YEAR>breaks), color = "#003399", fill = "#003399", size = line_size, 
                method = "glm", method.args = list(family = "poisson")) +
    scale_y_continuous(limits = c(0, 7)) + 
    scale_x_date(date_labels = "%Y", breaks = seq.Date(from = as.Date(paste0("2006-01-01")), 
                                                       to = as.Date(paste0("2022-12-31")), by = "4 years")) + 
    labs(x = "Calendar year", y = "Number of index isolates\n(per 100,000 Calgary residents)") +
    theme_template_white() + 
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank(), 
          strip.text = element_text(size = strip_text_size, face = "bold"),
          axis.text.x= element_text(size = axis_text_x)))


plot_COI <- grid.arrange(plot1, plot2, nrow=1, widths = c(1, 1))
ggsave(plot_COI, file = "figures/figure4/lsarp_figure_4B_EC_COI.pdf", 
       width = 13, height=5.5)

################################################################################
# ST1193, Hospital-onset
################################################################################
# R-R
(plot2.1 <- ggplot(data = df_mlst_cipro_ceftr_HOI_yr %>% 
                   filter(mlst==1193, ab_name=="R-R") %>% 
                   right_join(as.data.frame(cbind(YEAR = seq(2006, 2022)))) %>% 
                   mutate(ab_name = ifelse(is.na(ab_name), "R-R", ab_name), 
                          mlst = 1193, 
                          ST_name = "ST1193",
                          HOSPITAL_ONSET_48H = "Hospital-onset"), 
                 aes(x = ymd(YEAR, truncated=2), y = 100000*n/n_year, offset = n_year)) +
    facet_wrap(.~ab_name) + 
    geom_bar(stat= "identity", fill = "#CD0BBC", color = "black") +

    scale_y_continuous(limits = c(0, 7)) + 
    scale_x_date(date_labels = "%Y", breaks = seq.Date(from = as.Date(paste0("2006-01-01")), 
                                                       to = as.Date(paste0("2022-12-31")), by = "4 years"), 
                 limits = c(as.Date(paste0("2006-01-01")), as.Date(paste0("2022-12-31")))) + 
    labs(x = "Calendar year", y = "Number of index isolates\n(per 100,000 patient days)") +
    theme_template_white() + 
   theme(axis.title.x = element_blank(), 
         strip.text = element_text(size = strip_text_size, face = "bold"),
         axis.text.x= element_text(size = axis_text_x),
         axis.title.y= element_text(size = axis_title_y)))

# ------------------------------------------------------------------------------
# R-S
# ------------------------------------------------------------------------------
(plot2.2 <- ggplot(data = df_mlst_cipro_ceftr_HOI_yr %>% 
                   filter(mlst==1193, ab_name=="S-R") %>% 
                   right_join(as.data.frame(cbind(YEAR = seq(2006, 2022)))) %>% 
                   mutate(ab_name = ifelse(is.na(ab_name), "S-R", ab_name), 
                          mlst = 1193, 
                          ST_name = "ST1193",
                          HOSPITAL_ONSET_48H = "Hospital-onset"), 
                 aes(x = ymd(YEAR, truncated=2), y = inc)) +
    facet_grid(rows = vars(HOSPITAL_ONSET_48H), cols = vars(ab_name)) + 
    geom_bar(stat= "identity", fill = "#003399", color = "black") +
    scale_y_continuous(limits = c(0, 8)) + 
    scale_x_date(date_labels = "%Y", breaks = seq.Date(from = as.Date(paste0("2006-01-01")), 
                                                       to = as.Date(paste0("2022-12-31")), by = "4 years"), 
                 limits = c(as.Date(paste0("2006-01-01")), as.Date(paste0("2022-12-31")))) + 
    labs(x = "Calendar year", y = "Number of 30-day index isolates\n(per 100,000 Calgary residents)") +
    theme_template_white() + 
   theme(axis.title.x = element_blank(), 
         axis.title.y = element_blank(), 
         strip.text = element_text(size = strip_text_size, face = "bold"),
         axis.text.x= element_text(size = axis_text_x)))


plot_HOI <- grid.arrange(plot2.1, plot2.2, nrow=1, widths = c(1, 1))
ggsave(plot_HOI, file = "figures/figure4/lsarp_figure_4B_EC_HOI.pdf", 
       width = 13, height=5.5)

################################################################################
# Combined plot
plot <- grid.arrange(plot_COI, plot_HOI, ncol = 1, heights = c(1, 1))
ggsave(plot, file = "figures/figure4/lsarp_figure_4B.pdf", 
       width = 13, height=11)


