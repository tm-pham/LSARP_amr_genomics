# ============================================================================ #
# Plotting functions
# ============================================================================ #
library(ggthemes)
theme_template <- function(base_size=14, base_family="helvetica") {
  (theme_foundation(base_size=base_size, base_family=base_family)
   + theme_bw() 
   + theme(axis.line = element_line(colour="black"),
           axis.title.y = element_text(size=20),
           axis.title.x = element_text(size=20),
           axis.text=element_text(size=18),
           axis.text.x = element_text(angle=0, size=18, hjust=0.5), 
           axis.ticks.length=unit(.25, "cm"), 
           strip.text = element_text(size=20, face="bold"),
           legend.title = element_text(size=22, face="bold"), 
           legend.text = element_text(size=20),
           panel.background = element_rect(colour = NA),
           plot.background = element_rect(colour = NA),
           panel.grid.minor=element_line(color="gray94", linetype = 'solid'),
           panel.grid.major=element_line(color="gray94", linetype = 'solid')
   ))
}


theme_template_white <- function(base_size=14, base_family="helvetica") {
  (theme_foundation(base_size=base_size, base_family=base_family)
   + theme_bw() 
   + theme(axis.line = element_line(colour="black"),
           axis.title.y = element_text(size=20),
           axis.title.x = element_text(size=20),
           axis.text=element_text(size=18),
           axis.text.x = element_text(angle=0, size=18, hjust=0.5), 
           axis.ticks.length=unit(.25, "cm"), 
           strip.text = element_text(size=20, face="bold"),
           legend.title = element_text(size=22, face="bold"), 
           legend.text = element_text(size=20),
           panel.background = element_rect(colour = NA),
           plot.background = element_rect(colour = NA), 
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank()
   ))
}
