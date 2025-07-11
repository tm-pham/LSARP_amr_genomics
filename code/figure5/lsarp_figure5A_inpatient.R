# ============================================================================ #
# Project: LSARP
# Author: Thi Mui Pham, tmpham@hsph.harvard.edu
# Title: Figure 5A (Inpatient prescribing)
# ---------------------------------------------------------------------------- #
# Data for this is created in: 221017-tmp_comm_prescription_plots.R
# ============================================================================ #
remove(list=ls())
# Set working directory here
setwd("/Users/tm-pham/academia/hsph/LSARP_amr_genomics/")

# Load packges
library(dplyr)
library(ggplot2)

# Load functions and templates
source("code/functions_templates/plotting_template.R")

################################################################################
# Load data 
# Created in: lsarp_community_abx_use_total_data_prep.R
load("/Users/tm-pham/academia/hsph/lsarp/data/abx_use/inpatient/calgary_inpatient_abx_use_processed_2003_to_2022.RData") 

# Define antibiotic classes of interest
drug_classes <- c("Fluoroquinolones", "Penicillins", "Glycopeptides", "Aminoglycoside", "Tetracyclines", "Folate pathway inhibitor", 
                  "Cephalosporins", "Lincosamides", "Macrolides")
drug_labels <- c("Fluoroquinolones", "Penicillins", "Glycopeptides", "Aminoglycosides", "Tetracyclines", "Folate pathway\ninhibitor", 
             "Cephalosporins", "Lincosamides", "Macrolides")
drug_types <- c("Fluoroquinolones", "Beta-lactamase sensitive penicillins", "Beta-lactamase resistant penicillins", "Penicillins with extended spectrum",
                "Glycopeptides", "Aminoglycoside", "Streptomycin", "Tetracyclines", "Folate pathway inhibitor", 
                "CPH-03 (cephalosporin 3rd gen)", "Lincosamides", "Macrolides")


df_inp_plot <- df_inp_abx %>% filter(drug_class%in%drug_classes, 
                                     !drug_subclass%in%c("Siderophore cephalosporin", 
                                                         "CPH-MRSA (cephalosporin with anti-mrsa activity)", 
                                                         "CPH-01 (cephalosporin 1st gen)", 
                                                         "CPH-02 (cephalosporin 2nd gen)", 
                                                         "CPH-04 (cephalosporin 4th gen)", 
                                                         "CPH-05 (cephalosporin 5th gen)"), 
                                     drug_subclass%in%drug_types, Year >=2006) %>% 
  mutate(drug_class = factor(drug_class, levels = drug_classes, labels = drug_labels),
         drug_subclass = factor(drug_subclass, levels = drug_types),
         type = "Inpatient", 
         DDD_rate = ifelse(drug_subclass == "Beta-lactamase sensitive penicillins" & Year == 2022, NA, DDD_rate))

(plot_abx_inpatient <- ggplot(df_inp_plot, 
                              aes(x = ymd(Year, truncated = 2), y= DDD_rate, 
                                  color = factor(drug_subclass), 
                                  group = drug_subclass)) + 
    facet_grid(cols = vars(drug_class), rows = vars(type), scales = "fixed") + 
    geom_point(aes(color = factor(drug_subclass), group = drug_subclass), size=3) + 
    geom_smooth(linewidth=2, span=1.1) + 
    labs(y = "Defined daily doses\nper 100 patient days", color = "Antibiotic subclasses") + 
    scale_x_date(date_labels = "%Y", breaks = seq.Date(from = as.Date(paste0("2006-01-01")), 
                                                       to = as.Date(paste0("2022-12-31")), by = "4 years")) + 
    scale_y_continuous(limits = c(0, NA)) + 
    scale_color_manual(values = c("black", "darkblue", "lightblue", "dodgerblue4",rep("black", 7))) + 
    theme_template_white() + 
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_text(size=20, face = "plain"),
          axis.text.y = element_text(size=16),
          axis.text.x = element_text(size=18, angle=0),
          strip.text = element_text(size=24, face = "bold"),
          panel.spacing = unit(1.2, "cm"),
          legend.position = "none", 
          legend.title = element_text(size=14), 
          legend.text = element_text(size=12)))
ggsave(plot_abx_inpatient, file = "figures/figure5/lsarp_figure5A_inaptient.pdf", 
       width=36, height=4.3)
