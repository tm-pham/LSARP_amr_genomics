# ============================================================================ #
# Project: LSARP 
# Author: Thi Mui Pham, tmpham@hsph.harvard.edu
# Title: Figure 1A Segmented regression for COMMUNITY-ONSET
# Organisms: Overall across 5 target species
# ---------------------------------------------------------------------------- #
# Offset: 
# - Calgary population for community-onset
# ============================================================================ #
remove(list=ls())
# Set working directory here
setwd("/Users/tm-pham/academia/hsph/LSARP_amr_genomics/")

# Output paths
path <- "/Users/tm-pham/academia/hsph/lsarp/publications/lsarp_epi_genomics/" 
output_path <- "results/"

# Load libraries
library(dplyr)
library(reshape2)
library(MASS)
library(lubridate)
library(segmented)
library(gridExtra)

# Load functions
source("code/functions_templates/lsarp_function_aapc.R")
source("code/functions_templates/plotting_template.R")

# Load data 
load("data/LSARP_drug_bug_dictionary.RData") # drugs
load("data/figure1/lsarp_figure1A_input_df_inc_COI_year.RData")
load("data/figure1/lsarp_figure1A_input_df_inc_org_COI_year.RData")

################################################################################
# Colors for community- and hospital-acquired infections
COI_color <- "darkblue"
HOI_color <- "darkred"
point_size <- 4
line_size <- 2.5
point_shape <- 21
point_stroke <- 2
x_text_size = 26
x_title_size = 28
strip_text_size = 30

################################################################################
# Total incidence (across all five species)
################################################################################
dataset <- df_inc_COI_year 
model_COI_year <- glm(n ~ YEAR + offset(log(n_year)), 
                      data = dataset, 
                      family = poisson())
davies.test(model_COI_year, seg.Z = ~YEAR, k = 10)
seg_model <- segmented::segmented(model_COI_year, npsi = 1)

# Extract the breakpoints
breaks <- seg_model$psi[, 2]
dataset$organism <- "Total"

(plot1 <- ggplot(data = dataset, aes(x = ymd(YEAR, truncated=2), y = 100000*n/n_year, offset = n_year)) +
    facet_wrap(.~organism) + 
    geom_point(color = "black", fill = COI_color, size = point_size, shape = point_shape, stroke = point_stroke) +
    geom_smooth(data = dataset %>% filter(YEAR<=breaks), color = COI_color, fill = COI_color, size = line_size, 
                method = "glm", method.args = list(family = "poisson")) +
    geom_smooth(data = dataset %>% filter(YEAR>breaks), color = COI_color, fill = COI_color, size = line_size, 
                method = "glm", method.args = list(family = "poisson")) +
    scale_y_continuous(limits = c(0, NA)) + 
    scale_x_date(date_labels = "%Y", breaks = seq.Date(from = as.Date(paste0("2006-01-01")), 
                                                       to = as.Date(paste0("2022-12-31")), by = "4 years")) + 
    labs(x = "Calendar year", y = "Number of 30-day index isolates\n(per 100,000 Calgary residents)") +
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
dataset <- df_inc_org_COI_year %>% filter(ORG_LONG_NAME == "Staphylococcus aureus")
model_org_COI_year <- glm(n ~ YEAR + offset(log(n_year)), 
                          data = dataset, 
                          family = poisson())

davies.test(model_org_COI_year, seg.Z = ~YEAR, k = 5)
seg_model <- segmented::segmented(model_org_COI_year, npsi = 1)

# Extract the breakpoints
breaks <- seg_model$psi[, 2]

(plot2 <- ggplot(data = dataset, aes(x = ymd(YEAR, truncated=2), y = inc)) +
    facet_wrap(.~ORG_LONG_NAME) + 
    geom_point(color = "black", fill = COI_color, size = point_size, shape = point_shape, stroke = point_stroke) +
    geom_smooth(data = dataset %>% filter(YEAR<=breaks), color = COI_color, fill = COI_color, size = line_size, 
                method = "glm", method.args = list(family = "poisson")) +
    geom_smooth(data = dataset %>% filter(YEAR>breaks), color = COI_color, fill = COI_color, size = line_size, 
                method = "glm", method.args = list(family = "poisson")) +
    scale_y_continuous(limits = c(0, NA)) + 
    scale_x_date(date_labels = "%Y", breaks = seq.Date(from = as.Date(paste0("2006-01-01")), 
                                                       to = as.Date(paste0("2022-12-31")), by = "4 years")) + 
    labs(y = "Number of 30-day index isolates\n(per 100,000 Calgary residents)") +
    theme_template_white() + 
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank(), 
          axis.text.x = element_text(size = x_text_size),
          axis.text.y = element_text(size = x_text_size), 
          strip.text = element_text(size = strip_text_size, face = "bold.italic"), 
          plot.margin = margin(10, 20, 10, 10)))

################################################################################
# Enterococcus faecalis
################################################################################
dataset <- df_inc_org_COI_year %>% filter(ORG_LONG_NAME == "Enterococcus faecalis")

(plot3 <- ggplot(data = dataset, aes(x = ymd(YEAR, truncated=2), y = inc)) +
    facet_wrap(.~ORG_LONG_NAME) + 
    geom_point(color = "black", fill = COI_color, size = point_size, shape = point_shape, stroke = point_stroke) +
    geom_smooth(color = COI_color, size = line_size, fill = COI_color,
                method = "glm.nb") +
    scale_y_continuous(limits = c(0, NA)) + 
    scale_x_date(date_labels = "%Y", breaks = seq.Date(from = as.Date(paste0("2006-01-01")), 
                                                       to = as.Date(paste0("2022-12-31")), by = "4 years")) + 
    labs(y = "Number of 30-day index isolates\n(per 100,000 Calgary residents)") +
    theme_template_white() + 
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank(), 
          axis.text.x = element_text(size = x_text_size),
          axis.text.y = element_text(size = x_text_size), 
          strip.text = element_text(size = strip_text_size, face = "bold.italic"), 
          plot.margin = margin(10, 20, 10, 10)))

################################################################################
# Enterococcus faecium
################################################################################
dataset <- df_inc_org_COI_year %>% filter(ORG_LONG_NAME == "Enterococcus faecium")
model_org_COI_year <- glm(n ~ YEAR + offset(log(n_year)), 
                          data = dataset, 
                          family = poisson())

davies.test(model_org_COI_year, seg.Z = ~YEAR, k = 10)
seg_model <- segmented::segmented(model_org_COI_year, npsi = 1)

# Extract the breakpoints
(breaks <- seg_model$psi[, 2])

(plot4 <- ggplot(data = dataset, aes(x = ymd(YEAR, truncated=2), y = inc)) +
    facet_wrap(.~ORG_LONG_NAME) + 
    geom_point(color = "black", fill = COI_color, size = point_size, shape = point_shape, stroke = point_stroke) +
    geom_smooth(data = dataset %>% filter(YEAR<=breaks), color = COI_color, fill = COI_color, size = line_size, 
                method = "glm", method.args = list(family = "poisson")) +
    geom_smooth(data = dataset %>% filter(YEAR>breaks), color = COI_color, fill = COI_color, size = line_size, 
                method = "glm", method.args = list(family = "poisson")) +
    scale_y_continuous(limits = c(0, NA)) + 
    scale_x_date(date_labels = "%Y", breaks = seq.Date(from = as.Date(paste0("2006-01-01")), 
                                                       to = as.Date(paste0("2022-12-31")), by = "4 years")) + 
    geom_vline(xintercept = breaks, linetype = "dashed", color = "blue") +
    labs(y = "Number of 30-day index isolates\n(per 100,000 Calgary residents)") +
    theme_template_white() + 
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank(), 
          axis.text.x = element_text(size = x_text_size),
          axis.text.y = element_text(size = x_text_size), 
          strip.text = element_text(size = strip_text_size, face = "bold.italic"), 
          plot.margin = margin(10, 20, 10, 10)))

################################################################################
# Klebsiella pneumoniae
################################################################################
# Use negative binomial
dataset <- df_inc_org_COI_year %>% filter(ORG_LONG_NAME == "Klebsiella pneumoniae") %>% 
  mutate(inc = ifelse(YEAR>2022, NA, inc))
model_org_COI_year <- glm.nb(n ~ YEAR + offset(log(n_year)), 
                             data = dataset)

(plot5 <- ggplot(data = dataset, aes(x = ymd(YEAR, truncated=2), y = inc)) +
    facet_wrap(.~ORG_LONG_NAME) + 
    geom_point(color = "black", fill = COI_color, size = point_size, shape = point_shape, stroke = point_stroke) +
    geom_smooth(color = COI_color, size = line_size, fill = COI_color,
                method = "glm.nb") +
    scale_y_continuous(limits = c(0, NA)) + 
    scale_x_date(date_labels = "%Y", breaks = seq.Date(from = as.Date(paste0("2006-01-01")), 
                                                       to = as.Date(paste0("2022-12-31")), by = "4 years")) + 
    labs(y = "Number of 30-day index isolates\n(per 100,000 Calgary residents)") +
    theme_template_white() + 
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank(), 
          axis.text.x = element_text(size = x_text_size),
          axis.text.y = element_text(size = x_text_size), 
          strip.text = element_text(size = strip_text_size, face = "bold.italic"), 
          plot.margin = margin(10, 20, 10, 10)))

################################################################################
# E. coli
################################################################################
dataset <- df_inc_org_COI_year %>% filter(ORG_LONG_NAME == "Escherichia coli")
model_org_COI_year <- glm(n ~ YEAR + offset(log(n_year)), 
                          data = dataset, 
                          family = poisson())

davies.test(model_org_COI_year, seg.Z = ~YEAR, k = 10)
seg_model <- segmented::segmented(model_org_COI_year, npsi = 1)
# Extract the breakpoints
(breaks <- seg_model$psi[, 2])

(plot6 <- ggplot(data = dataset, aes(x = ymd(YEAR, truncated=2), y = inc)) +
    facet_wrap(.~ORG_LONG_NAME) + 
    geom_point(color = "black", fill = COI_color, size = point_size, shape = point_shape, stroke = point_stroke) +
    geom_smooth(data = dataset %>% filter(YEAR<=breaks), color = COI_color, fill = COI_color, size = line_size, 
                method = "glm", method.args = list(family = "poisson")) +
    geom_smooth(data = dataset %>% filter(YEAR>breaks), color = COI_color, fill = COI_color, size = line_size, 
                method = "glm", method.args = list(family = "poisson")) +
    scale_y_continuous(limits = c(0, NA)) + 
    scale_x_date(date_labels = "%Y", breaks = seq.Date(from = as.Date(paste0("2006-01-01")), 
                                                       to = as.Date(paste0("2022-12-31")), by = "4 years")) + 
    labs(y = "Number of 30-day index isolates\n(per 100,000 Calgary residents)") +
    theme_template_white() + 
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank(), 
          axis.text.x = element_text(size = x_text_size),
          axis.text.y = element_text(size = x_text_size), 
          strip.text = element_text(size = strip_text_size, face = "bold.italic"), 
          plot.margin = margin(10, 20, 10, 10)))


################################################################################
# Combine plots
################################################################################
plot <- grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, nrow=1, 
                     widths = c(1.12, 1, 1, 1, 1, 1))
ggsave(plot, file = "figures/figure1/lsarp_figure_1A.pdf", 
       width = 42, height=7)
ggsave(plot, file = "figures/figure1/lsarp_figure_1A.png", 
       width = 42, height=7)

