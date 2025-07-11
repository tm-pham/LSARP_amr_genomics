# ============================================================================ #
# Project: LSARP
# Author: Thi Mui Pham, tmpham@hsph.harvard.edu
# Title: Figure 5A (Community antibiotic prescribing)
# ---------------------------------------------------------------------------- #
# Data for this is created in: 221017-tmp_comm_prescription_plots.R
# ============================================================================ #
remove(list=ls())
# Set working directory here
setwd("/Users/tm-pham/academia/hsph/LSARP_amr_genomics/")

# Load functions and templates
source("code/functions_templates/plotting_template.R")

################################################################################
# Load data 
# Created in: lsarp_community_abx_use_total_data_prep.R
load("data/figure5/lsarp_community_abx_age_standardized_per_year.RData") # df_comm_abx_std

df_comm_abx_std <- df_comm_abx_std %>%
  mutate(drug.class = ifelse(drug.class == "Sulfonamides & Trimethorpim", "Folate path inhibitor", as.character(drug.class)), 
         drug.class = ifelse(drug.class == "Glycopeptide", "Glycopeptides", as.character(drug.class)),
         drug.class = ifelse(drug.class == "Other Aminoglycosides", "Aminoglycosides", as.character(drug.class)),
         Drug.Type = ifelse(Drug.Type == "Sulfonamides & Trimethorpim", "Folate path inhibitor", as.character(Drug.Type)), 
         Drug.Type = ifelse(Drug.Type == "Glycopeptide", "Glycopeptides", as.character(Drug.Type)),
         Drug.Type = ifelse(Drug.Type == "Other Aminoglycosides", "Aminoglycosides", as.character(Drug.Type)))

drug_classes <- c("Fluoroquinolones", "Penicillins", "Glycopeptides", "Aminoglycosides", "Tetracyclines", "Folate path inhibitor", 
                  "Cephalosporins", "Lincosamides", "Macrolides")
drug_labels <- c("Fluoroquinolones", "Penicillins", "Glycopeptides", "Aminoglycosides", "Tetracyclines", "Folate pathway\ninhibitor", 
                 "Cephalosporins", "Lincosamides", "Macrolides")
drug_types <- c("Fluoroquinolones", "Beta-L-Sensitive-Penicillin", "Beta-L-Resistant-Penicillin", 
                "Glycopeptides", "Aminoglycosides", "Streptomycin", "Tetracyclines", "Folate path inhibitor", 
                "Third-G-Cephalosporins", "Lincosamides", "Macrolides")

df <- df_comm_abx_std %>% 
  filter(Sex=="BOTH", drug.class%in%drug_classes, Drug.Type%in%drug_types) %>%
  mutate(drug.class = factor(drug.class, levels = drug_classes, labels = drug_labels),
         Drug.Type = factor(Drug.Type, levels = drug_types), 
         type = "Community")

df_temp <- as.data.frame(cbind(Geography = "Z2", Year = 2006:2022, Sex = "BOTH", Drug.Type = "Glycopeptides", 
                               Dispensation.Rate= 0, Total.Dispensation = 0, Unique.Dispensation = 0, Total.Population = 0, 
                               Standard.Error = 0, Standard.Score = 0, Alberta.Rate = 0, drug.class = "Glycopeptides", 
                               type = "Community")) %>% 
  mutate(across(setdiff(colnames(df), c("Geography", "Sex","Drug.Type", "drug.class", "type")), as.numeric))
df <- rbind(df, df_temp)


(plot_abx_community <- ggplot(df, 
                              aes(x=ymd(Year, truncated=2), y=Alberta.Rate, color = Drug.Type)) + 
    facet_grid(cols = vars(drug.class), rows = vars(type)) + 
    geom_point(size=3) +
    geom_smooth(linewidth=2, span = 1) +
    labs(y="Age-standardized \nprescribing rate in CHZ", color = "Antibiotic subclass") + 
    scale_x_date(date_labels = "%Y", 
                 breaks = seq(as.Date("01-01-2006", format = "%d-%m-%Y"), as.Date("31-12-2022", format = "%d-%m-%Y"), by = "4 years")) + 
    scale_y_continuous(limits = c(0, NA)) + 
    scale_color_manual(values = c("black", "darkblue", "lightblue", rep("black", 8))) + 
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
ggsave(plot_abx_community, file = "figures/figure5/lsarp_figure5A_community.pdf", 
       width = 36, height = 4.3)
