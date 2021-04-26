read_panel_files <- function(panel_list       = list.files("datasets", pattern = "*.xlsx", full.names = T),
                             sample_size_file = list.files("datasets", pattern = "*.csv", full.names = T),
                             strata           = "clasificacion", 
                             drop_vars = c(), guess_max = 100000){
  
  #' @description Function that reads files from panel survey listed in panel_files
  #' and adds the sample sizes and the weights as contained in sample_size_file
  #' 
  #' @param panel_list List of Excel files with panel survey
  #' 
  #' @param sample_size_file File with the following structure
  #' 
  #'        strata   | Panel 1  | Panel 2 | ... | Panel K
  #'        strata 1 |  N_{1,1} | N_{1,2} | ... | N_{1,K}
  #'        strata 2 |  N_{2,1} | N_{2,2} | ... | N_{2,K}
  #'          .
  #'          .
  #'          .
  #'        strata m | N_{m,1} | N_{m,2} | ... | N_{m,K}
  #' 
  #' where N_{i,j} denotes the population size of strata i at panel j
  #' It is assumed that the sampling weights are given by N_{i,j} / N_j
  #' where N_j = sum_j N_{i,j} for each strata i.
  #' 
  #' @param strata Name of the strata both in sample_size_file and in each
  #' file of panel_list. Will be renamed to "strata".
  #' 
  #' @drop_vars Variables not to consider when reading panel data
  #' 
  #' @guess_max Maximum number of columns to read before guessing composition
  #' 
  #' @example 
  #' panel_list       <- list.files(pattern = "^resultado_*")
  #' sample_size_file <- list.files(pattern ="^Tamaños totales_*")[1]
  #' strata           <- "clasificacion"
  #' drop_vars        <- "valor_ct_2"
  #' read_panel_files(panel_list = panel_list, sample_size_file = sample_size_file,
  #'                  strata = strata, drop_vars = drop_vars)
  #'                  
  #' @return List with variables:
  #' n_dropped      .- Number of observations dropped as variable had no strata
  #' encuesta_panel .- Dataset with panel_survey and strata identified as strata
  #'                   and weights identified as sampling_weights
  
  #Read sample size file                   
  if (file.exists(sample_size_file)){
    N <- read_csv(sample_size_file)
  } else {
    stop(paste0("File ", sample_size_file, "not in current path. Set sample_size_file."))
  }
  
  #Check sample size matches number of panel files
  if (ncol(N) - 1 != length(panel_list)){
    stop(paste0("Invalid number of panels. Panels registered in ", sample_size_file, " differ in number from panel_list"))
  }
  
  #Check strata in colnames
  if (!(strata %in% colnames(N))){
    stop(paste0("strata ", strata, " not found in sample_size_file"))
  }
  
  panel_num <- 1
  for (fname in panel_list){
    
    if (file.exists(sample_size_file)){
      #temp_panel <- read_delim(fname, ";", escape_double = FALSE, trim_ws = TRUE)
      temp_panel <- read_excel(fname, guess_max = guess_max)
    } else {
      stop(paste0("Not all files in panel_list found. ", fname, " is not in current directory."))
    }
    
    #Add panel number
    temp_panel <- temp_panel %>% mutate(Ciclo = !!panel_num) %>%
      mutate(DIAS_HOSPITALIZADO = as.numeric(DIAS_HOSPITALIZADO)) %>%
      mutate(TIEMPO_NO_FUMAR_ANIOS = as.numeric(TIEMPO_NO_FUMAR_ANIOS)) %>%
      mutate(TIEMPO_NO_FUMAR_MESES = as.numeric(TIEMPO_NO_FUMAR_MESES)) %>%
      mutate(valor_ct_2 = as.character(valor_ct_2))
    panel_num  <- panel_num + 1
    
    #Delete vars in drop_vars
    if (length(drop_vars) > 0){
      for (i in 1:length(drop_vars)){
        if (drop_vars[i] %in% colnames(temp_panel)){
          temp_panel <- temp_panel %>% select(-!!sym(drop_vars[i]))
        }
      }
    }
    
    #Create panel
    if (panel_num == 2){
      encuesta_panel <- temp_panel
    } else {
      encuesta_panel <- encuesta_panel %>% bind_rows(temp_panel)
    }
  }
  
  #Drop observations without strata
  n_dropped      <- encuesta_panel %>% 
    filter(is.na(!!sym(strata))) %>% 
    tally() %>% as.numeric()
  encuesta_panel <- encuesta_panel %>% drop_na(!!sym(strata))
  
  #Bind N sizes to encuesta_panel
  for (i in 2:ncol(N)){
    
    #Get column N's
    colnum_st           <- which(colnames(N) == strata)
    sub_panel           <- N[,c(colnum_st,i)]
    colnames(sub_panel) <- c(strata, "N_muestra")
    
    #Get complete sample
    N_size <- sum(sub_panel[,"N_muestra"])
    
    #Get panel
    sub_panel           <- sub_panel %>% 
      mutate(Ciclo = !!(i - 1)) %>%
      mutate(N_size = !!N_size) %>% 
      mutate(sampling_weights = N_muestra/N_size)
    
    if (i == 2){
      temp_subpanel <- sub_panel
    } else {
      temp_subpanel <- temp_subpanel %>% bind_rows(sub_panel)
    }
    
  }
  
  #Bind columns
  encuesta_panel      <- encuesta_panel %>% left_join(temp_subpanel, by = c(strata, "Ciclo"))
  
  #Recode strata 
  encuesta_panel <- encuesta_panel %>% rename(strata  = !!sym(strata)) 
  
  return(list(n_dropped = n_dropped, encuesta_panel = encuesta_panel))
  
}