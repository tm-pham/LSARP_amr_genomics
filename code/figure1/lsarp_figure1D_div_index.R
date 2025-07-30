# ============================================================================ #
# Project: LSARP 
# Author: Thi Mui Pham, tmpham@hsph.harvard.edu
# Title: Diversity index for Figure 1D
# Organisms: All 5 species
# ---------------------------------------------------------------------------- #
# Shannon diversity index for each organism
# ============================================================================ #
remove(list=ls())
setwd("/Users/tm-pham/academia/hsph/LSARP_amr_genomics/")

# Packages
library(dplyr)
library(tidyr)

# Load data 
load("data/figure1/lsarp_figure1D_div_index_5species.RData")

# Load functions
source("code/functions_templates/plotting_template.R")

################################################################################
# Transform from list to data frame 
df_div <- do.call("rbind", df_div_list)
rownames(df_div) <- NULL

df_mlst_div_long <- data.frame(df_div) %>% melt(id=c("organismofinterest", "year")) %>% 
  mutate(organismofinterest = factor(organismofinterest, 
                                     levels = c("Staphylococcus aureus", 
                                                "Enterococcus faecalis", 
                                                "Enterococcus faecium", 
                                                "Klebsiella pneumoniae", 
                                                "Escherichia coli"), 
                                     labels = c("S aureus", 
                                                "E faecalis", 
                                                "E faecium", 
                                                "K pneumoniae", 
                                                "E coli")))

################################################################################
# Plot (Shannon diversity index)
(div_plot <- ggplot(df_mlst_div_long, aes(x=organismofinterest, y = as.numeric(value))) + 
    geom_boxplot(fill='#A4A4A4', color = "black", size = 4., outlier.size = 3.5, outlier.shape = 1, outlier.stroke = 2) + 
    geom_jitter(shape=21, position=position_jitter(0.2), size = 6, stroke = 2, fill = "white") + 
    labs(y = "Shannon diversity index (per year)") + 
    theme_template_white() + 
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_text(size = 32),
          axis.text.y = element_text(size = 30),
          axis.text.x = element_text(face="italic", size = 26, angle = 45, hjust=1),
          plot.margin = margin(10, 10, 5, 5)))
ggsave(div_plot, file = paste0(path, "figures/lsarp_fig1D_1.pdf"), 
       width = 10, height=9)
