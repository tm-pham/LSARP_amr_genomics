library(tidyverse)
library(cowplot)

# input table containing per isolate collection dates and hospital-onset/community-onset status
metadata <- read_csv("/bulk/LSARP/datasets/APL/versions/230924/230924-sw__APL-results-INTERP.csv")

# input table with per isolate WGS quality metrics
sequence_qc <- read_csv("/bulk/LSARP/genomics/pipeline/Escherichia_coli/samplemetrics.csv")

# input table with per isolate straingst clusters with simplified cluster numbering
renamed_clusters <- read_tsv("/bulk/LSARP/datasets/250106-tm__rename_straingst_clusters/Ec_straingst_clusters.tsv")
renamed_clusters <- renamed_clusters %>% mutate(cluster_id = as_factor(cluster_id))

sequence_qc <- sequence_qc %>% filter(status == "PASS") %>% select(sample, straingst_ref, mlst) 
metadata <- metadata %>% select(BI_NBR, INDEX_30DAYS, COLLECT_DTM, HOSPITAL_ONSET_48H)
metadata <- metadata %>% filter(INDEX_30DAYS) %>% rename(sample = BI_NBR)
sequence_qc <- sequence_qc %>% left_join(metadata)

sequence_qc <- sequence_qc %>% drop_na()

ec <- sequence_qc %>% select(-INDEX_30DAYS)

# read amr data generated with AMRFinderPlus

column_types <- cols(Name = col_character(), `Protein identifier` = col_logical(), `Contig id` = col_double(), Start = col_double(), Stop = col_double(), Strand = col_character(), `Gene symbol` = col_character(), `Sequence name` = col_character(), Scope = col_character(), `Element type` = col_character(), `Element subtype` = col_character(), Class = col_character(), Method = col_character(), `Target length` = col_double(), `Reference sequence length` = col_double(), `% Coverage of reference sequence` = col_double(), `% Identity to reference sequence` = col_double(), `Alignment length` = col_double(), `Accession of closest sequence` = col_character(), `Name of closest sequence` = col_character(), `HMM id` = col_logical(), `HMM description` = col_logical())
genotype <- list.files(path="/bulk/LSARP/genomics/analyses/Escherichia_coli/data/amrfinder", full.names=TRUE) %>% map_dfr(read_tsv, col_types = column_types)
genotype <- unique(genotype)
genotype_wider <- genotype %>% count(Name, `Gene symbol`) %>% pivot_wider(names_from = "Gene symbol", values_from = "n", values_fill = 0) %>% rename(sample = Name)

esbl_alleles <- genotype %>% filter(str_detect(Subclass, "CEPHALOSPORIN")) %>% count(`Gene symbol`) %>% pull(`Gene symbol`)
ec <- ec %>% left_join(genotype_wider)
ec_esbl <- ec %>% select(sample, straingst_ref, COLLECT_DTM, HOSPITAL_ONSET_48H, any_of(esbl_alleles))
ec_esbl <- ec_esbl %>% drop_na()
ec_esbl %>% select_if(is.numeric) %>% map_dbl(sum)
ec_esbl <- ec_esbl %>% mutate(year = year(COLLECT_DTM))

ec_esbl <- ec_esbl %>% select(sample, straingst_ref, year, HOSPITAL_ONSET_48H, starts_with("blaCTX"), starts_with("blaSHV"))

# remove rare alleles or alleles without notable experimental evidence for association with cephalosporin resistance

ec_esbl <- ec_esbl %>% select(-`blaCTX-M`,-`blaCTX-M-1`,-`blaCTX-M-104`,-`blaCTX-M-189`,-`blaCTX-M-192`,-`blaCTX-M-231`,-`blaCTX-M-24`,-`blaCTX-M-3`,-`blaCTX-M-32`,-`blaCTX-M-65`,-`blaSHV-12`,-`blaSHV-2`,-`blaSHV-2A`)
ec_esbl_long <- ec_esbl %>% pivot_longer(`blaCTX-M-14`:`blaCTX-M-55`, names_to = "gene", values_to = "present")
ec_esbl_long <- ec_esbl_long %>% filter(present != 0)
esbl_counts <- ec_esbl_long %>% count(straingst_ref, year, HOSPITAL_ONSET_48H, gene) %>% arrange(year)

esbl_counts <- esbl_counts %>% left_join(renamed_clusters)

# identify clusters with at least 10 ESBL+ index isolates
notable_cluster_ids <- esbl_counts %>% group_by(cluster_id, gene) %>% summarize(total = sum(n)) %>% filter(total >= 10) %>% pull(cluster_id)


# plot

p <- esbl_counts %>%
    group_by(straingst_ref, year, gene, cluster_id) %>% summarize(total = sum(n)) %>%
    ggplot(aes(x=year,y=fct_relevel(cluster_id,rev))) + 
    geom_tile(aes(fill = total)) + 
    facet_wrap(~gene,nrow=1) + 
    theme_cowplot(10) + 
    theme(strip.background=element_blank(), strip.text.x = element_blank(), axis.ticks.y=element_blank(), axis.text.x =element_text(angle = 90, vjust = 0.5, hjust = 1), plot.margin = margin(c(0,0,0,0, "cm"))) +
    xlab("Year") +
    ylab("Clusters") +
    scale_y_discrete(labels = ~ ifelse(.x %in% notable_cluster_ids, .x, ""))

p_histogram <- esbl_counts %>% group_by(year,gene,HOSPITAL_ONSET_48H) %>% summarize(total = sum(n)) %>%
    ggplot(aes(x=year, y = total, fill = HOSPITAL_ONSET_48H)) +
           geom_bar(stat= "identity", position = "stack") +
           scale_fill_grey(end = 0.5, name = "Onset (48 h)", labels = c("Community", "Hospital")) + 
           facet_wrap(~gene,nrow=1) +
           theme_cowplot(10) +
           xlab("") + 
           ylab("Isolate Count") + 
           theme(axis.ticks.x=element_blank(), axis.text.x=element_blank(), plot.margin = margin(c(0,0,0,0, "cm")))
aligned <- plot_grid(p_histogram, p, ncol = 1, align="v", axis="lr", rel_heights = c(1, 4))
ggsave("figures/figure4/lasrp_figure_4C_ec_esbls.pdf", aligned, width = 6, height = 8, units = "in")


