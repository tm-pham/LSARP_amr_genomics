# ============================================================================ #
# Diverity index function
# ============================================================================ #
loadRData <- function(fileName, return_first = T){
  #loads an RData file, and returns it
  load(fileName)
  if(return_first) mget(ls()[ls() != "fileName"])[["df_APL"]]
  else mget(ls()[ls() != "fileName"])
}

diversity.index.function <- function(df_APL, 
                                     org = NULL, 
                                     time_period =seq(2006, 2022), 
                                     time = F, 
                                     variable = "mlst", 
                                     metric = "simpson"){
  df_APL <- df_APL %>% filter(YEAR%in%time_period)
  if(time){
    df_APL_div <- df_APL %>% ungroup() %>% 
      dplyr::select(all_of(c("YEAR", variable))) %>% 
      mutate(n=1) %>% 
      pivot_wider(names_from=variable, values_from = n, values_fn=sum, values_fill=0)
    
    if(is.null(org)){
      df_div <- data.frame(cbind(organismofinterest = "All species", 
                                 year = sort(unique(df_APL$YEAR)), 
                                 div.index = vegan::diversity(df_APL_div[,-1], metric)))
      
    }else{
      org <- unique(df_APL$ORG_LONG_NAME)
      df_div <- data.frame(cbind(organismofinterest = org, 
                                 year = sort(unique(df_APL$YEAR)), 
                                 div.index = vegan::diversity(df_APL_div[,-1], metric)))
    }
  }else{
    df_APL_div <- df_APL %>% ungroup() %>% 
      dplyr::select(all_of(variable)) %>% 
      mutate(n=1) %>% 
      pivot_wider(names_from=variable, values_from = n, values_fn=sum, values_fill=0)
    
    if(is.null(org)){
      df_div <- data.frame(cbind(organismofinterest = "All species", 
                                 div.index = vegan::diversity(df_APL_div[,-1], metric)))
      
    }else{
      org <- unique(df_APL$ORG_LONG_NAME)
      df_div <- data.frame(cbind(organismofinterest = org, 
                                 div.index = vegan::diversity(df_APL_div[,-1], metric)))
    }
    
  }
  return(df_div)
}
