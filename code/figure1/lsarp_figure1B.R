# ============================================================================ #
# Project: LSARP 
# Author: Thi Mui Pham, tmpham@hsph.harvard.edu
# Title: Incidence model using negative binomial regression for COMMUNITY-ONSET
# Organisms: Overall across 5 target species
# ---------------------------------------------------------------------------- #
# Offset: 
# - Calgary population for community-onset
# ============================================================================ #
remove(list=ls())
# Set working directory here
setwd("/Users/tm-pham/academia/hsph/LSARP_amr_genomics/")

# Load libraries
library(dplyr)
library(reshape2)
library(lubridate)
library(segmented)
library(gridExtra)

# Load functions
source("code/functions_templates/lsarp_function_aapc.R")
source("code/functions_templates/plotting_template.R")

# Load data 
load("data/LSARP_drug_bug_dictionary.RData") # drugs
load("data/figure1/lsarp_figure1B_input_df_inc_HOI_year.RData")
load("data/figure1/lsarp_figure1B_input_df_inc_org_HOI_year.RData")

################################################################################
# Colors
HOI_color <- "darkred"
point_size <- 4
line_size <- 2.5
point_shape <- 21
point_stroke <- 2
x_text_size = 26
x_title_size = 28
strip_text_size = 34

################################################################################
# Total incidence (across all five species)
################################################################################
dataset <- df_inc_HOI_year
dataset$organism <- "Total"

(plot1 <- ggplot(data = dataset, aes(x = ymd(YEAR, truncated=2), y = inc)) +
    facet_wrap(.~organism) + 
    geom_point(color = "black", fill = HOI_color, size = point_size, shape = point_shape, stroke = point_stroke) +
    geom_smooth(data = dataset, color = HOI_color, fill = HOI_color, size = line_size,
                method = "glm.nb") +
    scale_y_continuous(limits = c(0, NA)) + 
    scale_x_date(date_labels = "%Y", breaks = seq.Date(from = as.Date("2006-01-01"), 
                                                       to = as.Date("2022-12-31"), by = "4 years")) + 
    labs(x = "Calendar year", y = "Number of 30-day index isolates\n(per 100,000 patient days)") +
    theme_template_white() + 
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_text(size = x_text_size),
          axis.text.y = element_text(size = x_text_size), 
          axis.title.y = element_text(size = x_title_size), 
          strip.text = element_text(size = strip_text_size), 
          plot.margin = margin(10, 20, 10, 10)))

################################################################################
# Staphylococcus aureus
################################################################################
dataset <- df_inc_org_HOI_year %>% filter(ORG_LONG_NAME == "Staphylococcus aureus")
model_org_HOI_year <- glm(n ~ YEAR + offset(log(n_year)), 
                          data = dataset, 
                          family = poisson())

davies.test(model_org_HOI_year, seg.Z = ~YEAR, k = 5)
seg_model <- segmented::segmented(model_org_HOI_year, npsi = 1)
segmented::selgmented(model_org_HOI_year, seg.Z = ~YEAR, type = "bic", Kmax = 5, bonferroni = T)

# Extract the breakpoints
(breaks <- seg_model$psi[, 2])

(plot2 <- ggplot(data = dataset, aes(x = ymd(YEAR, truncated=2), y = inc)) +
    facet_wrap(.~ORG_LONG_NAME) + 
    geom_point(color = "black", fill = HOI_color, size = point_size, shape = point_shape, stroke = point_stroke) +
    geom_smooth(data = dataset %>% filter(YEAR<=breaks), color = HOI_color, fill = HOI_color, size = line_size, 
                method = "glm", method.args = list(family = "poisson")) +
    geom_smooth(data = dataset %>% filter(YEAR>breaks), color = HOI_color, fill = HOI_color, size = line_size, 
                method = "glm", method.args = list(family = "poisson")) +
    scale_y_continuous(limits = c(0, NA)) + 
    scale_x_date(date_labels = "%Y", breaks = seq.Date(from = as.Date("2006-01-01"), 
                                                       to = as.Date("2022-12-31"), by = "4 years")) + 
    geom_vline(xintercept = breaks, linetype = "dashed", color = HOI_color) +
    labs(y = "Number of 30-day index isolates\n(per 100,000 Calgary residents)") +
    theme_template_white() + 
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank(), 
          axis.text.x = element_text(size = x_text_size),
          axis.text.y = element_text(size = x_text_size), 
          strip.text = element_text(size = strip_text_size, face="bold.italic"), 
          plot.margin = margin(10, 20, 10, 10)))

################################################################################
# Enterococcus faecalis
################################################################################
dataset <- df_inc_org_HOI_year %>% filter(ORG_LONG_NAME == "Enterococcus faecalis")

(plot3 <- ggplot(data = dataset, aes(x = ymd(YEAR, truncated=2), y = inc)) +
    facet_wrap(.~ORG_LONG_NAME) + 
    geom_point(color = "black", fill = HOI_color, size = point_size, shape = point_shape, stroke = point_stroke) +
    geom_smooth(data = dataset, color = HOI_color, fill = HOI_color, size = line_size,
                method = "glm.nb") +
    scale_y_continuous(limits = c(0, NA)) + 
    scale_x_date(date_labels = "%Y", breaks = seq.Date(from = as.Date("2006-01-01"), 
                                                       to = as.Date("2022-12-31"), by = "4 years")) + 
    labs(y = "Number of 30-day index isolates\n(per 100,000 Calgary residents)") +
    theme_template_white() + 
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank(), 
          axis.text.x = element_text(size = x_text_size),
          axis.text.y = element_text(size = x_text_size), 
          strip.text = element_text(size = strip_text_size, face="bold.italic"), 
          plot.margin = margin(10, 20, 10, 10)))

################################################################################
# Enterococcus faecium
################################################################################
dataset <- df_inc_org_HOI_year %>% filter(ORG_LONG_NAME == "Enterococcus faecium")
model_org_HOI_year <- glm(n ~ YEAR + offset(log(n_year)), 
                          data = dataset, 
                          family = poisson())

davies.test(model_org_HOI_year, seg.Z = ~YEAR, k = 10)
seg_model <- segmented::segmented(model_org_HOI_year, npsi = 1)
segmented::selgmented(model_org_HOI_year, seg.Z = ~YEAR, type = "bic", Kmax = 3, bonferroni = T)

# Extract the breakpoints
(breaks <- seg_model$psi[, 2])

(plot4 <- ggplot(data = dataset, aes(x = ymd(YEAR, truncated=2), y = inc)) +
    facet_wrap(.~ORG_LONG_NAME) + 
    geom_point(color = "black", fill = HOI_color, size = point_size, shape = point_shape, stroke = point_stroke) +
    geom_smooth(data = dataset %>% filter(YEAR<=breaks), color = HOI_color, fill = HOI_color, size = line_size, 
                method = "glm", method.args = list(family = "poisson")) +
    geom_smooth(data = dataset %>% filter(YEAR>breaks), color = HOI_color, fill = HOI_color, size = line_size, 
                method = "glm", method.args = list(family = "poisson")) +
    scale_y_continuous(limits = c(0, NA)) + 
    scale_x_date(date_labels = "%Y", breaks = seq.Date(from = as.Date("2006-01-01"), 
                                                       to = as.Date("2022-12-31"), by = "4 years")) + 
    geom_vline(xintercept = breaks, linetype = "dashed", color = HOI_color) +
    labs(y = "Number of 30-day index isolates\n(per 100,000 Calgary residents)") +
    theme_template_white() + 
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank(), 
          axis.text.x = element_text(size = x_text_size),
          axis.text.y = element_text(size = x_text_size), 
          strip.text = element_text(size = strip_text_size, face="bold.italic"), 
          plot.margin = margin(10, 20, 10, 10)))

################################################################################
# K. pneumoniae
################################################################################
# Use negative binomial
dataset <- df_inc_org_HOI_year %>% filter(ORG_LONG_NAME == "Klebsiella pneumoniae")

(plot5 <- ggplot(data = dataset, aes(x = ymd(YEAR, truncated=2), y = inc)) +
    facet_wrap(.~ORG_LONG_NAME) + 
    geom_point(color = "black", fill = HOI_color, size = point_size, shape = point_shape, stroke = point_stroke) +
    geom_smooth(data = dataset, color = HOI_color, fill = HOI_color, size = line_size,
                method = "glm.nb") +
    scale_y_continuous(limits = c(0, NA)) + 
    scale_x_date(date_labels = "%Y", breaks = seq.Date(from = as.Date("2006-01-01"), 
                                                       to = as.Date("2022-12-31"), by = "4 years"), 
                 limits = c(as.Date("2006-01-01"), as.Date("2022-12-31"))) + 
    labs(y = "Number of 30-day index isolates\n(per 100,000 Calgary residents)") +
    theme_template_white() + 
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank(), 
          axis.text.x = element_text(size = x_text_size),
          axis.text.y = element_text(size = x_text_size), 
          strip.text = element_text(size = strip_text_size, face="bold.italic"), 
          plot.margin = margin(10, 20, 10, 10)))

################################################################################
# E. coli
################################################################################
dataset <- df_inc_org_HOI_year %>% filter(ORG_LONG_NAME == "Escherichia coli")

(plot6 <- ggplot(data = dataset, aes(x = ymd(YEAR, truncated=2), y = inc)) +
    facet_wrap(.~ORG_LONG_NAME) + 
    geom_point(color = "black", fill = HOI_color, size = point_size, shape = point_shape, stroke = point_stroke) +
    geom_smooth(data = dataset, color = HOI_color, fill = HOI_color, size = line_size,
                method = "glm.nb") +
    scale_y_continuous(limits = c(0, NA)) + 
    scale_x_date(date_labels = "%Y", breaks = seq.Date(from = as.Date("2006-01-01"), 
                                                       to = as.Date("2022-12-31"), by = "4 years"), 
                 limits = c(as.Date("2006-01-01"), as.Date("2022-12-31"))) + 
    labs(y = "Number of 30-day index isolates\n(per 100,000 Calgary residents)") +
    theme_template_white() + 
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank(), 
          axis.text.x = element_text(size = x_text_size),
          axis.text.y = element_text(size = x_text_size), 
          strip.text = element_text(size = strip_text_size, face="bold.italic"), 
          plot.margin = margin(10, 20, 10, 10)))

################################################################################
# Combine plot
################################################################################
plot <- grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, nrow=1, 
                     widths = c(1.1, 1, 1, 1, 1, 1))
ggsave(plot, file = "figures/figure1/lsarp_figure_1B.pdf", 
       width = 42, height=7)
ggsave(plot, file = "figures/figure1/lsarp_figure_1B.png", 
       width = 42, height=7)

