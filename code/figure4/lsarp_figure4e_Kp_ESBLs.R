library(tidyverse)
library(cowplot)

# input table containing per isolate collection dates and hospital-onset/community-onset status
metadata <- read_csv("/bulk/LSARP/datasets/APL/versions/230924/230924-sw__APL-results-INTERP.csv")

# input table with per isolate WGS quality metrics
sequence_qc <- read_csv("/bulk/LSARP/genomics/pipeline/Klebsiella_pneumoniae/samplemetrics.csv")

# input table with per isolate straingst clusters with simplified cluster numbering
renamed_clusters <- read_tsv("/bulk/LSARP/datasets/250106-tm__rename_straingst_clusters/Kp_straingst_clusters.tsv")
renamed_clusters <- renamed_clusters %>% mutate(cluster_id = as_factor(cluster_id))

sequence_qc <- sequence_qc %>% filter(status == "PASS") %>% select(sample, straingst_ref) 
metadata <- metadata %>% select(BI_NBR, INDEX_30DAYS, COLLECT_DTM, HOSPITAL_ONSET_48H)
metadata <- metadata %>% filter(INDEX_30DAYS) %>% rename(sample = BI_NBR)
sequence_qc <- sequence_qc %>% left_join(metadata)

sequence_qc <- sequence_qc %>% drop_na()

# select isolates from K. pneumoniae specifically rather than all Kp complex isolates

kp <- sequence_qc %>% select(-INDEX_30DAYS) %>% filter(str_detect(straingst_ref, "_pneumoniae"))

# read amr data generated with AMRFinderPlus

column_types <- cols(Name = col_character(), `Protein identifier` = col_logical(), `Contig id` = col_double(), Start = col_double(), Stop = col_double(), Strand = col_character(), `Gene symbol` = col_character(), `Sequence name` = col_character(), Scope = col_character(), `Element type` = col_character(), `Element subtype` = col_character(), Class = col_character(), Method = col_character(), `Target length` = col_double(), `Reference sequence length` = col_double(), `% Coverage of reference sequence` = col_double(), `% Identity to reference sequence` = col_double(), `Alignment length` = col_double(), `Accession of closest sequence` = col_character(), `Name of closest sequence` = col_character(), `HMM id` = col_logical(), `HMM description` = col_logical())
genotype <- list.files(path="/bulk/LSARP/genomics/analyses/Klebsiella_pneumoniae/data/amrfinder", full.names=TRUE) %>% map_dfr(read_tsv, col_types = column_types)
genotype <- unique(genotype)
genotype_wider <- genotype %>% count(Name, `Gene symbol`) %>% pivot_wider(names_from = "Gene symbol", values_from = "n", values_fill = 0) %>% rename(sample = Name)

esbl_alleles <- genotype %>% filter(str_detect(Subclass, "CEPHALOSPORIN")) %>% count(`Gene symbol`) %>% pull(`Gene symbol`)
kp <- kp %>% left_join(genotype_wider)
kp_esbl <- kp %>% select(sample, straingst_ref, COLLECT_DTM, HOSPITAL_ONSET_48H, any_of(esbl_alleles))
kp_esbl %>% select_if(is.numeric) %>% map_dbl(sum)
kp_esbl <- kp_esbl %>% mutate(year = year(COLLECT_DTM))

# remove rare alleles or alleles without notable experimental evidence for association with cephalosporin resistance

kp_esbl <- kp_esbl %>% select(-`blaSHV-100`,-`blaSHV-164`,-`blaSHV-187`,-`blaSHV-38`,-`blaSHV-65`,-`blaOXA-1`,-`blaOXA-10`,-`blaDHA-1`,-`blaSHV_C-112A`)


kp_esbl <- kp_esbl %>% select(-`blaDHA`,-`blaSHV-2A`,-`blaSHV-12`,-`blaCMY-2`)
kp_esbl_long <- kp_esbl %>% pivot_longer(`blaCTX-M-14`:`blaSHV-2`, names_to = "gene", values_to = "present")
kp_esbl_long <- kp_esbl_long %>% filter(present != 0)
esbl_counts <- kp_esbl_long %>% count(straingst_ref, year, gene, HOSPITAL_ONSET_48H) %>% arrange(year)
esbl_counts <- esbl_counts %>% left_join(renamed_clusters)

# identify clusters with at least 3 ESBL+ index isolates
notable_cluster_ids <- esbl_counts %>% group_by(cluster_id, gene) %>% summarize(total = sum(n)) %>% filter(total >= 3) %>% pull(cluster_id)

# plot (only include years prior to 2021 because of incomplete data collection post-2020

p <- esbl_counts %>% filter(year < 2021) %>% 
    group_by(straingst_ref, year, gene, cluster_id) %>% summarize(total = sum(n)) %>%
    ggplot(aes(x=year,y=fct_relevel(cluster_id, rev))) + 
           geom_tile(aes(fill = as_factor(total))) + 
           facet_wrap(~gene,nrow=1) + 
           theme_cowplot(10) + 
           theme(strip.background=element_blank(),  
                 strip.text.x = element_blank(), 
                 axis.text.x =element_text(angle = 90, vjust = 0.5, hjust = 1),
                 plot.margin = margin(0,1,0,0, "cm")) +
           xlab("Year") +
           ylab("Clusters") +
           scale_fill_manual(values = c("1" = "#9e9ac8", "2" = "#756bb1", "3" = "#54278f"), name = "Count") +
           scale_y_discrete(labels = ~ ifelse(.x %in% notable_cluster_ids, .x, ""))

p_histogram <- esbl_counts %>% filter(year < 2021) %>% group_by(year,gene,HOSPITAL_ONSET_48H) %>% summarize(total = sum(n)) %>%
    ggplot(aes(x=year, y = total, fill=HOSPITAL_ONSET_48H)) +
           geom_bar(stat= "identity", position = "stack") +
           scale_fill_grey(end = 0.5, name = "Onset (48 h)", labels = c("Community", "Hospital")) + 
           facet_wrap(~gene,nrow=1) +
           theme_cowplot(10) +
           scale_y_continuous(breaks=c(2,4,6,8,10)) +
           xlab("") + 
           ylab("Isolate Count") + 
           theme(axis.ticks.x=element_blank(), axis.text.x=element_blank(), plot.margin = margin(c(0,1,0,0, "cm")))
aligned <- plot_grid(p_histogram, p, ncol = 1, align="v", axis="lr", rel_heights = c(1, 3))
ggsave("figures/figure4/lsarp_figure_4E_kp_esbls.pdf", aligned, width = 7, height = 5, units = "in")


