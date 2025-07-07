# Thi Mui Pham, tmpham@hsph.harvard.edu
# ============================================================================ #
# Stacked resistance function
# ============================================================================ #
stacked.resistance.function <- function(df_APL, antibiotics, df_population){
  df_res_abx <- df_APL %>% 
    group_by_at(all_of(antibiotics)) %>% 
    mutate_at(antibiotics, ~ ifelse(.=="R", 1, 0)) %>%
    ungroup() %>% 
    mutate(res_sum = rowSums(.[, antibiotics]), 
           n_susc = ifelse(res_sum==0, 1, 0), 
           n_1R = ifelse(res_sum ==1, 1, 0), 
           n_2R = ifelse(res_sum ==2, 1, 0),
           n_3R = ifelse(res_sum ==3, 1, 0),
           n_4R = ifelse(res_sum ==4, 1, 0),
           n_5R = ifelse(res_sum ==5, 1, 0),
           n_6R = ifelse(res_sum ==6, 1, 0)) %>% 
    group_by(YEAR, GENDER, AGE_GRP) %>% 
    summarise(n_susc = sum(n_susc, na.rm=TRUE), 
              n_1R = sum(n_1R, na.rm=TRUE), 
              n_2R = sum(n_2R, na.rm=TRUE), 
              n_3R = sum(n_3R, na.rm=TRUE), 
              n_4R = sum(n_4R, na.rm=TRUE), 
              n_5R = sum(n_5R, na.rm=TRUE), 
              n_6R = sum(n_6R, na.rm=TRUE), 
              n_total = dplyr::n()) %>% 
    left_join(df_population %>% dplyr::select(year, sex, age_grp, n_pop, weight2006), by=c("YEAR"="year", 
                                                                                    "GENDER"="sex", 
                                                                                    "AGE_GRP"="age_grp")) %>% 
    mutate_at(vars(-YEAR, -GENDER, -AGE_GRP, -weight2006, -n_pop), ~ ./n_pop) %>% 
    ungroup() %>% 
    group_by(YEAR) %>% 
    summarise(n_susc = sum(n_susc*weight2006, na.rm=TRUE)*100000, 
              n_1R = sum(n_1R*weight2006, na.rm=TRUE)*100000, 
              n_2R = sum(n_2R*weight2006, na.rm=TRUE)*100000, 
              n_3R = sum(n_3R*weight2006, na.rm=TRUE)*100000, 
              n_4R = sum(n_4R*weight2006, na.rm=TRUE)*100000, 
              n_5R = sum(n_5R*weight2006, na.rm=TRUE)*100000, 
              n_6R = sum(n_6R*weight2006, na.rm=TRUE)*100000, 
              n_total = dplyr::n()) %>% unique()
  
  return(df_res_abx)
}
