#PARA BINOMIAL 0
#---------------------------------------------------------

get_nsize_zsize <- function(nfilter, sfilter, ciclo){
  
  nsize <- panel_datasets %>% 
    filter(Ciclo == !!ciclo) %>% 
    filter(eval_tidy(parse_expr(nfilter))) %>%
    group_by(strata) %>% tally() %>% 
    full_join(stratum, by = c("strata")) %>% 
    arrange(strata) %>% 
    mutate(n = replace_na(n, 0)) %>% 
    select(n) %>% unlist() %>% as.numeric() 
  
  zsize <- panel_datasets %>% 
    filter(Ciclo == !!ciclo) %>% 
    filter(eval_tidy(parse_expr(nfilter))) %>%
    filter(eval_tidy(parse_expr(sfilter))) %>%
    group_by(strata) %>% tally() %>% 
    full_join(stratum, by = "strata") %>% 
    arrange(strata) %>% 
    mutate(n = replace_na(n, 0)) %>% 
    select(n) %>% unlist() %>% as.numeric()
  
  return(list(nsize = nsize, zsize = zsize)) 
}

genera_binomial_0 <- function(nfilter, sfilter, varname = "VARIABLE",
                              delta = 1, gamma = 1, filedir = getwd(),
                              mu = 0.5, sigma = 1/12, 
                              burnin = 5000, niter = 10000, nchains = 4){
  
  p_interes <- list()
  q_interes <- list()
  
  for (ciclo in ciclos){
    enes <- get_nsize_zsize(nfilter, sfilter, ciclo = ciclo)
    sims <- calculo_binomial_0(enes$nsize, enes$zsize, compilar_julia = T,
                               gamma = gamma, delta = delta,
                               mu = mu, sigma = sigma, 
                               burnin = burnin, niter = niter, nchains = nchains)
    p_interes[[ciclo]]  <- sims %>% select(starts_with("p")) %>% 
      mutate(Ciclo = ciclo) %>% 
      mutate(iteration = 1:n())
    q_interes[[ciclo]]  <- sims %>% select(starts_with("q"))  %>% 
      mutate(Ciclo = ciclo) %>%
      mutate(iteration = 1:n())
  }
  
  #Calculamos los pvals
  pvals           <- do.call("rbind", p_interes)
  colnames(pvals) <- c(stratum$strata, "Ciclo", "iteration")
  pvals <- pvals %>% pivot_longer(1:length(ciclos), names_to = "strata", values_to = "p") 
  pvals <- pvals %>% left_join(pondef, by = c("strata","Ciclo"))
  
  totales <- pvals %>% ungroup() %>%
    group_by(iteration, Ciclo) %>%
    dplyr::summarise(p = sum(p*sampling_weights), .groups = "keep")  %>%
    mutate(strata = "Total")
  
  pvals <- pvals %>% bind_rows(totales) %>% select(-sampling_weights)
  
  #Calculamos los qvals
  qvals <- do.call("rbind", q_interes)
  colnames(qvals) <- c(stratum$strata, "Ciclo", "iteration")
  qvals <- qvals %>% pivot_longer(1:length(ciclos), names_to = "strata", values_to = "q")
  qvals <- qvals %>% left_join(pondef, by = c("strata","Ciclo"))
  
  totales <- qvals %>% ungroup() %>%
    group_by(iteration, Ciclo) %>%
    dplyr::summarise(q = sum(q*sampling_weights), .groups = "keep")  %>%
    mutate(strata = "Total")
  
  qvals <- qvals %>% bind_rows(totales) %>% select(-sampling_weights)
  
  #Calculamos las medias
  medias_p <- pvals %>% group_by(Ciclo, strata) %>% mean_hdi(p) %>%
    select(Ciclo, strata, p, .lower, .upper, .width) %>%
    mutate(`Estimador` = 
             paste0(percent(p, 0.01)," [", 
                    percent(.lower, 0.01, suffix = ""), ", ",
                    percent(.upper, 0.01, suffix = ""), "]")) %>%
    mutate(.width = percent(.width)) %>%
    rename(Confianza = .width) %>%
    mutate(Tipo = "No ajustado")
  
  medias_q <- qvals %>% group_by(Ciclo, strata) %>% mean_hdi(q) %>%
    rename(p = q) %>%
    select(Ciclo, strata, p, .lower, .upper, .width) %>%
    mutate(`Estimador` = 
             paste0(percent(p, 0.01)," [", 
                    percent(.lower, 0.01, suffix = ""), ", ",
                    percent(.upper, 0.01, suffix = ""), "]")) %>%
    mutate(.width = percent(.width)) %>%
    rename(Confianza = .width) %>% 
    mutate(Tipo = "Ajustado")
  
  resultados <- medias_p %>% bind_rows(medias_q) %>% 
    mutate(Variable = !!varname) %>% 
    select(Variable, Tipo, Ciclo, strata, Estimador, Confianza, p, .lower, .upper) 
  
  resultados %>% write_csv(paste0(filedir,varname,".csv"))
  return(resultados)
}

#PARA BINOMIAL 1
#---------------------------------------------------------
get_nsize_zsize_tsize <- function(nfilter, sfilter, tfilter, ciclo){
  
  nsize_zsize <- get_nsize_zsize(nfilter, sfilter, ciclo)
  
  tsize <- panel_datasets %>% 
    filter(Ciclo == !!ciclo) %>% 
    filter(eval_tidy(parse_expr(nfilter))) %>%
    filter(eval_tidy(parse_expr(sfilter))) %>%
    filter(eval_tidy(parse_expr(tfilter))) %>%
    group_by(strata) %>% tally() %>% 
    full_join(stratum, by = "strata") %>%
    arrange(strata) %>% 
    mutate(n = replace_na(n, 0)) %>% 
    select(n) %>% unlist() %>% as.numeric()
  
  if (any(nsize_zsize$zsize > nsize_zsize$nsize) | any(tsize > nsize_zsize$zsize)){
    stop("Probabilidades inválidas. No se cumple que tsize <= zsize <= nsize") 
  }
  
  return(list(nsize = nsize_zsize$nsize, zsize = nsize_zsize$zsize, tsize = tsize)) 
}

genera_binomial_1_iter <- function(nfilter, sfilter, tfilter, varname = "VARIABLE",
                                   pname = "p", qname = "q",
                                   mu_p = 0.5, sigma_p = 1/12, mu_q = 0.5, sigma_q = 1/12, 
                                   delta = 1, gamma = 1, 
                                   burnin = 5000, niter = 10000, nchains = 4){
  
  p_interes      <- list()
  q_interes      <- list()
  ajus_q_interes <- list()
  
  for (ciclo in ciclos){
    enes <- get_nsize_zsize_tsize(nfilter, sfilter, tfilter, ciclo = ciclo)
    sims <- calculo_binomial_1(enes$nsize, enes$zsize, enes$tsize, compilar_julia = T,
                               gamma = gamma, delta = delta,
                               mu_p = mu_p, sigma_p = sigma_p, mu_q = mu_q, sigma_q = sigma_q, 
                               burnin = burnin, niter = niter, nchains = nchains)
    p_interes[[ciclo]]  <- sims %>% select(starts_with("p")) %>% 
      mutate(Ciclo = ciclo) %>% 
      mutate(iteration = 1:n())
    q_interes[[ciclo]]  <- sims %>% select(starts_with("q"))  %>% 
      mutate(Ciclo = ciclo) %>%
      mutate(iteration = 1:n())
    ajus_q_interes[[ciclo]]  <- sims %>% select(starts_with("ajus_q"))  %>% 
      mutate(Ciclo = ciclo) %>%
      mutate(iteration = 1:n())
  }
  
  #Calculamos los pvals
  pvals           <- do.call("rbind", p_interes)
  colnames(pvals) <- c(stratum$strata, "Ciclo", "iteration")
  pvals <- pvals %>% pivot_longer(1:length(ciclos), names_to = "strata", values_to = "p") 
  pvals <- pvals %>% left_join(pondef, by = c("strata","Ciclo"))
  
  totales <- pvals %>% ungroup() %>%
    group_by(iteration, Ciclo) %>%
    dplyr::summarise(p = sum(p*sampling_weights), .groups = "keep")  %>%
    mutate(strata = "Total")
  
  pvals <- pvals %>% bind_rows(totales) %>% select(-sampling_weights)
  
  #Calculamos los qvals
  qvals <- do.call("rbind", q_interes)
  colnames(qvals) <- c(stratum$strata, "Ciclo", "iteration")
  qvals <- qvals %>% pivot_longer(1:length(ciclos), names_to = "strata", values_to = "q")
  qvals <- qvals %>% left_join(pondef, by = c("strata","Ciclo"))
  
  totales <- qvals %>% ungroup() %>%
    group_by(iteration, Ciclo) %>%
    dplyr::summarise(q = sum(q*sampling_weights), .groups = "keep")  %>%
    mutate(strata = "Total")
  
  qvals <- qvals %>% bind_rows(totales) %>% select(-sampling_weights)
  
  #Calculamos los ajus_q
  ajus_qvals <- do.call("rbind", ajus_q_interes)
  colnames(ajus_qvals) <- c(stratum$strata, "Ciclo", "iteration")
  ajus_qvals <- ajus_qvals %>% pivot_longer(1:length(ciclos), names_to = "strata", values_to = "ajus_q")
  ajus_qvals <- ajus_qvals %>% left_join(pondef, by = c("strata","Ciclo"))
  
  totales <- ajus_qvals %>% ungroup() %>%
    group_by(iteration, Ciclo) %>%
    dplyr::summarise(ajus_q = sum(ajus_q*sampling_weights), .groups = "keep")  %>%
    mutate(strata = "Total")
  
  ajus_qvals <- ajus_qvals %>% bind_rows(totales) %>% select(-sampling_weights)
  
  #Calculamos las medias
  medias_p <- pvals %>% group_by(Ciclo, strata) %>% mean_hdi(p) %>%
    select(Ciclo, strata, p, .lower, .upper, .width) %>%
    mutate(`Estimador` = 
             paste0(percent(p, 0.01)," [", 
                    percent(.lower, 0.01, suffix = ""), ", ",
                    percent(.upper, 0.01, suffix = ""), "]")) %>%
    mutate(.width = percent(.width)) %>%
    rename(Confianza = .width) %>%
    mutate(Tipo = "No ajustado") %>% 
    mutate(Estimando = pname)
  
  medias_q <- qvals %>% group_by(Ciclo, strata) %>% mean_hdi(q) %>%
    rename(p = q) %>%
    select(Ciclo, strata, p, .lower, .upper, .width) %>%
    mutate(`Estimador` = 
             paste0(percent(p, 0.01)," [", 
                    percent(.lower, 0.01, suffix = ""), ", ",
                    percent(.upper, 0.01, suffix = ""), "]")) %>%
    mutate(.width = percent(.width)) %>%
    rename(Confianza = .width) %>% 
    mutate(Tipo = "No ajustado") %>%
    mutate(`Estimando` = qname)
  
  medias_q_ajus <- ajus_qvals %>% group_by(Ciclo, strata) %>% mean_hdi(ajus_q) %>%
    rename(p = ajus_q) %>%
    select(Ciclo, strata, p, .lower, .upper, .width) %>%
    mutate(`Estimador` = 
             paste0(percent(p, 0.01)," [", 
                    percent(.lower, 0.01, suffix = ""), ", ",
                    percent(.upper, 0.01, suffix = ""), "]")) %>%
    mutate(.width = percent(.width)) %>%
    rename(Confianza = .width) %>% 
    mutate(Tipo = "Ajustado") %>% 
    mutate(Estimando = qname)
  
  resultados <- medias_p %>% bind_rows(medias_q) %>% 
    bind_rows(medias_q_ajus) %>% 
    mutate(Variable = !!varname) %>% 
    select(Variable, Estimando, Tipo, Ciclo, strata, Estimador, Confianza, p, .lower, .upper) 
  
  #resultados %>% write_csv(paste0(filedir,varname,".csv"))
  return(resultados)
}

genera_binomial_1 <- function(nfilter, sfilter, tfilter, varname = "VARIABLE",
                              mu_p = 0.5, sigma_p = 1/12, 
                              pname = c("Probabilidad (p1)", "Probabilidad (p2)"),
                              qname = c("Probabilidad (q1)", "Probabilidad (q2)"),
                              filedir = getwd(),
                              mu_q = rep(0.5, 2), sigma_q = rep(1/12, 2), 
                              delta = 1, gamma = 1, 
                              burnin = 5000, niter = 10000, nchains = 4){
  
  resultados_1 <- genera_binomial_1_iter(nfilter = nfilter, 
                                         sfilter = sfilter, 
                                         tfilter = tfilter, 
                                         varname = varname, 
                                         pname = pname[1],
                                         qname = qname[1],
                                         mu_p = mu_p, sigma_p = sigma_p,
                                         mu_q = mu_q[1], sigma_q = sigma_q[1],
                                         delta = delta, gamma = gamma, 
                                         burnin = burnin, niter = niter, 
                                         nchains = nchains)
  
  
  # resultados_2 <- genera_binomial_1_iter(nfilter = nfilter, 
  #                                        sfilter = paste0("!(",sfilter,")"),
  #                                        tfilter = tfilter, 
  #                                        varname = varname, 
  #                                        pname = pname[2],
  #                                        qname = qname[2],
  #                                        mu_p = 1 - mu_p, sigma_p = sigma_p,
  #                                        mu_q = mu_q[2], sigma_q = sigma_q[2],
  #                                        delta = delta, gamma = gamma, 
  #                                        burnin = burnin, niter = niter, 
  #                                        nchains = nchains)
  # 
  
  resultados <- resultados_1 #%>% bind_rows(resultados_2)
  
  resultados %>% write_csv(paste0(filedir,varname,".csv"))
  return(resultados)
  
}

#PARA BINOMIAL 2
#---------------------------------------------------------
get_nsize_s1_s2_s3 <- function(nfilter, s1filter, s2filter, s3filter, ciclo){
  
  nsize_zsize_tsize <- get_nsize_zsize_tsize(nfilter, s1filter, s2filter, ciclo)
  
  s3size <- panel_datasets %>% 
    filter(Ciclo == !!ciclo) %>% 
    filter(eval_tidy(parse_expr(nfilter))) %>%
    filter(eval_tidy(parse_expr(s1filter))) %>%
    filter(eval_tidy(parse_expr(s2filter))) %>%
    filter(eval_tidy(parse_expr(s3filter))) %>%
    group_by(strata) %>% tally() %>% 
    full_join(stratum, by = "strata") %>% 
    arrange(strata) %>%  
    mutate(n = replace_na(n, 0)) %>% 
    select(n) %>% unlist() %>% as.numeric()
  
  if (any(s3size > nsize_zsize_tsize$tsize)){
    stop("Probabilidades inválidas. No se cumple que s3 <= s2 <= s1 <= nsize") 
  }
  
  return(list(nsize  = nsize_zsize_tsize$nsize, 
              s1size = nsize_zsize_tsize$zsize, 
              s2size = nsize_zsize_tsize$tsize, 
              s3size = s3size)) 
}

genera_binomial_2_iter <- function(nfilter, s1filter, s2filter, s3filter, 
                                   varname = "VARIABLE",
                                   p1name = "p1", p2name = "p2", p3name = "p3", 
                                   mu_p1 = 0.5, sigma_p1 = 1/12, 
                                   mu_p2 = 0.5, sigma_p2 = 1/12, 
                                   mu_p3 = 0.5, sigma_p3 = 1/12, 
                                   delta = 1, gamma = 1, 
                                   seeds = NULL,
                                   burnin = 5000, niter = 10000, nchains = 4){
  
  p1_interes      <- list()
  p2_interes      <- list()
  p3_interes      <- list()
  ajus_p3_interes <- list()
  
  if (length(burnin) == 1){
    burnin <- rep(burnin, length(ciclos))
  }

  if (length(niter) == 1){
    niter <- rep(niter, length(ciclos))
  }
  
  for (ciclo in ciclos){
    
    enes <- get_nsize_s1_s2_s3(nfilter, s1filter, s2filter, s3filter, ciclo = ciclo)
    
    if (!is.null(seeds)){
      julia_eval(paste0("Random.seed!(", seeds[ciclo],");")) 
    }
    
    sims <- calculo_binomial_2(enes$nsize, enes$s1size, enes$s2size, enes$s3size,
                               compilar_julia = T,
                               gamma = gamma, delta = delta,
                               mu_p1 = mu_p1, mu_p2 = mu_p2, mu_p3 = mu_p3,
                               sigma_p1 = sigma_p1, sigma_p2 = sigma_p2, sigma_p3 = sigma_p3,
                               burnin = burnin[ciclo], niter = niter[ciclo], nchains = nchains)
    
    p1_interes[[ciclo]]  <- sims %>% select(starts_with("p1")) %>% 
      mutate(Ciclo = ciclo) %>% 
      mutate(iteration = 1:n())
    p2_interes[[ciclo]]  <- sims %>% select(starts_with("p2")) %>% 
      mutate(Ciclo = ciclo) %>% 
      mutate(iteration = 1:n())
    p3_interes[[ciclo]]  <- sims %>% select(starts_with("p3")) %>% 
      mutate(Ciclo = ciclo) %>% 
      mutate(iteration = 1:n())
    ajus_p3_interes[[ciclo]]  <- sims %>% select(starts_with("ajus"))  %>% 
      mutate(Ciclo = ciclo) %>%
      mutate(iteration = 1:n())
  }
  
  #Calculamos los p1
  p1vals           <- do.call("rbind", p1_interes)
  colnames(p1vals) <- c(stratum$strata, "Ciclo", "iteration")
  p1vals <- p1vals %>% pivot_longer(1:length(ciclos), names_to = "strata", values_to = "p1") 
  p1vals <- p1vals %>% left_join(pondef, by = c("strata","Ciclo"))
  
  totales <- p1vals %>% ungroup() %>%
    group_by(iteration, Ciclo) %>%
    dplyr::summarise(p1 = sum(p1*sampling_weights), .groups = "keep")  %>%
    mutate(strata = "Total")
  
  p1vals <- p1vals %>% bind_rows(totales) %>% select(-sampling_weights)
  
  #Calculamos los p2
  p2vals <- do.call("rbind", p2_interes)
  colnames(p2vals) <- c(stratum$strata, "Ciclo", "iteration")
  p2vals <- p2vals %>% pivot_longer(1:length(ciclos), names_to = "strata", values_to = "p2")
  p2vals <- p2vals %>% left_join(pondef, by = c("strata","Ciclo"))
  
  totales <- p2vals %>% ungroup() %>%
    group_by(iteration, Ciclo) %>%
    dplyr::summarise(p2 = sum(p2*sampling_weights), .groups = "keep")  %>%
    mutate(strata = "Total")
  
  p2vals <- p2vals %>% bind_rows(totales) %>% select(-sampling_weights)
  
  #Calculamos los p3
  p3vals <- do.call("rbind", p3_interes)
  colnames(p3vals) <- c(stratum$strata, "Ciclo", "iteration")
  p3vals <- p3vals %>% pivot_longer(1:length(ciclos), names_to = "strata", values_to = "p3")
  p3vals <- p3vals %>% left_join(pondef, by = c("strata","Ciclo"))
  
  totales <- p3vals %>% ungroup() %>%
    group_by(iteration, Ciclo) %>%
    dplyr::summarise(p3 = sum(p3*sampling_weights), .groups = "keep")  %>%
    mutate(strata = "Total")
  
  p3vals <- p3vals %>% bind_rows(totales) %>% select(-sampling_weights)
  
  #Calculamos los qvals
  ajus_p3vals <- do.call("rbind", ajus_p3_interes)
  colnames(ajus_p3vals) <- c(stratum$strata, "Ciclo", "iteration")
  ajus_p3vals <- ajus_p3vals %>% pivot_longer(1:length(ciclos), names_to = "strata", values_to = "ajus_p3")
  ajus_p3vals <- ajus_p3vals %>% left_join(pondef, by = c("strata","Ciclo"))
  
  totales <- ajus_p3vals %>% ungroup() %>%
    group_by(iteration, Ciclo) %>%
    dplyr::summarise(ajus_p3 = sum(ajus_p3*sampling_weights), .groups = "keep")  %>%
    mutate(strata = "Total")
  
  ajus_p3vals <- ajus_p3vals %>% bind_rows(totales) %>% select(-sampling_weights)
  
  #Calculamos las medias
  medias_p1 <- p1vals %>% group_by(Ciclo, strata) %>% mean_hdi(p1) %>%
    rename(p = p1) %>%
    select(Ciclo, strata, p, .lower, .upper, .width) %>%
    mutate(`Estimador` = 
             paste0(percent(p, 0.01)," [", 
                    percent(.lower, 0.01, suffix = ""), ", ",
                    percent(.upper, 0.01, suffix = ""), "]")) %>%
    mutate(.width = percent(.width)) %>%
    rename(Confianza = .width) %>%
    mutate(Tipo = "No ajustado") %>% 
    mutate(Estimando = p1name)
  
  medias_p2 <- p2vals %>% group_by(Ciclo, strata) %>% mean_hdi(p2) %>%
    rename(p = p2) %>%
    select(Ciclo, strata, p, .lower, .upper, .width) %>%
    mutate(`Estimador` = 
             paste0(percent(p, 0.01)," [", 
                    percent(.lower, 0.01, suffix = ""), ", ",
                    percent(.upper, 0.01, suffix = ""), "]")) %>%
    mutate(.width = percent(.width)) %>%
    rename(Confianza = .width) %>%
    mutate(Tipo = "No ajustado") %>% 
    mutate(Estimando = p2name)
  
  medias_p3 <- p3vals %>% group_by(Ciclo, strata) %>% mean_hdi(p3) %>%
    rename(p = p3) %>%
    select(Ciclo, strata, p, .lower, .upper, .width) %>%
    mutate(`Estimador` = 
             paste0(percent(p, 0.01)," [", 
                    percent(.lower, 0.01, suffix = ""), ", ",
                    percent(.upper, 0.01, suffix = ""), "]")) %>%
    mutate(.width = percent(.width)) %>%
    rename(Confianza = .width) %>%
    mutate(Tipo = "No ajustado") %>% 
    mutate(Estimando = p3name)
  
  medias_ajus_p3 <- ajus_p3vals %>% group_by(Ciclo, strata) %>% 
    mean_hdi(ajus_p3) %>%
    rename(p = ajus_p3) %>%
    select(Ciclo, strata, p, .lower, .upper, .width) %>%
    mutate(`Estimador` = 
             paste0(percent(p, 0.01)," [", 
                    percent(.lower, 0.01, suffix = ""), ", ",
                    percent(.upper, 0.01, suffix = ""), "]")) %>%
    mutate(.width = percent(.width)) %>%
    rename(Confianza = .width) %>%
    mutate(Tipo = "Ajustado") %>% 
    mutate(Estimando = p3name)
  
  
  resultados <- medias_p1 %>% bind_rows(medias_p2) %>% 
    bind_rows(medias_p3) %>% bind_rows(medias_ajus_p3) %>% 
    mutate(Variable = !!varname) %>% 
    select(Variable, Estimando, Tipo, Ciclo, strata, Estimador, Confianza, p, .lower, .upper) 
  
  #resultados %>% write_csv(paste0(filedir,varname,".csv"))
  return(resultados)
}



genera_binomial_2 <- function(nfilter, s1filter, s2filter, s3filter, 
                              varname = "VARIABLE",
                              p1name = c("p1","p1"), 
                              p2name = c("p2","p2"), 
                              p3name = c("p3","p3"),
                              mu_p1 = 0.5, sigma_p1 = 1/12, 
                              mu_p2 = rep(0.5, 2), sigma_p2 = rep(1/12, 2),
                              mu_p3 = rep(0.5, 2), sigma_p3 = rep(1/12, 2),
                              filedir = getwd(),
                              seeds = NULL,
                              delta = 1, gamma = 1, 
                              burnin = 5000, niter = 10000, nchains = 4){
  
    
  resultados_1 <- genera_binomial_2_iter(nfilter = nfilter, 
                                         s1filter = s1filter, 
                                         s2filter = s2filter, 
                                         s3filter = s3filter,
                                         varname = varname, 
                                         p1name = p1name[1], 
                                         p2name = p2name[1],
                                         p3name = p3name[1],
                                         seeds = seeds,
                                         mu_p1 = mu_p1, sigma_p1 = sigma_p1,
                                         mu_p2 = mu_p2[1], sigma_p2 = sigma_p2[1],
                                         mu_p3 = mu_p3[1], sigma_p3 = sigma_p3[1],
                                         delta = delta, gamma = gamma, 
                                         burnin = burnin, niter = niter, 
                                         nchains = nchains)
  
  
  resultados <- resultados_1
  
  resultados %>% write_csv(paste0(filedir,varname,".csv"))
  return(resultados)
  
}

#PARA BINOMIAL 3
#---------------------------------------------------------
get_nsize_s1_s2_s3_s4 <- function(nfilter, s1filter, s2filter, s3filter, s4filter, ciclo){
  
  nsize_s1_s2_s3 <- get_nsize_s1_s2_s3(nfilter, s1filter, s2filter, s3filter, ciclo)
  
  s4size <- panel_datasets %>% 
    filter(Ciclo == !!ciclo) %>% 
    filter(eval_tidy(parse_expr(nfilter))) %>%
    filter(eval_tidy(parse_expr(s1filter))) %>%
    filter(eval_tidy(parse_expr(s2filter))) %>%
    filter(eval_tidy(parse_expr(s3filter))) %>%
    filter(eval_tidy(parse_expr(s4filter))) %>%
    group_by(strata) %>% tally() %>% 
    full_join(stratum, by = "strata") %>% 
    arrange(strata) %>%  
    mutate(n = replace_na(n, 0)) %>% 
    select(n) %>% unlist() %>% as.numeric()
  
  if (any(s4size > nsize_s1_s2_s3$s3size)){
    stop("Probabilidades inválidas. No se cumple que s4 <= s3 <= s2 <= s1 <= nsize") 
  }
  
  return(list(nsize  = nsize_s1_s2_s3$nsize, 
              s1size = nsize_s1_s2_s3$s1size, 
              s2size = nsize_s1_s2_s3$s2size, 
              s3size = nsize_s1_s2_s3$s3size,
              s4size = s4size)) 
}

genera_binomial_3_iter <- function(nfilter, s1filter, s2filter, s3filter, s4filter,
                                   varname = "VARIABLE",
                                   p1name = "p1", p2name = "p2", p3name = "p3", 
                                   p4name = "p4", 
                                   mu_p1 = 0.5, sigma_p1 = 1/12, 
                                   mu_p2 = 0.5, sigma_p2 = 1/12, 
                                   mu_p3 = 0.5, sigma_p3 = 1/12,
                                   mu_p4 = 0.5, sigma_p4 = 1/12, 
                                   delta = 1, gamma = 1, 
                                   seeds = NULL,
                                   burnin = 5000, niter = 10000, nchains = 4){

  
  p1_interes      <- list()
  p2_interes      <- list()
  p3_interes      <- list()
  p4_interes      <- list()
  ajus_p4_interes <- list()
  
  if (length(burnin) == 1){
    burnin <- rep(burnin, length(ciclos))
  }
  
  if (length(niter) == 1){
    niter <- rep(niter, length(ciclos))
  }
  
  for (ciclo in ciclos){
    
    enes <- get_nsize_s1_s2_s3_s4(nfilter, s1filter, s2filter, s3filter, s4filter, ciclo = ciclo)
    
    if (!is.null(seeds)){
      julia_eval(paste0("Random.seed!(", seeds[ciclo],");")) 
    }
    
    sims <- calculo_binomial_3(enes$nsize, enes$s1size, enes$s2size, enes$s3size, enes$s4size,
                               compilar_julia = T,
                               gamma = gamma, delta = delta,
                               mu_p1 = mu_p1, mu_p2 = mu_p2, 
                               mu_p3 = mu_p3, mu_p4 = mu_p4,
                               sigma_p1 = sigma_p1, sigma_p2 = sigma_p2, 
                               sigma_p3 = sigma_p3, sigma_p4 = sigma_p4,
                               burnin = burnin[ciclo], niter = niter[ciclo], 
                               nchains = nchains)
    
    p1_interes[[ciclo]]  <- sims %>% select(starts_with("p1")) %>% 
      mutate(Ciclo = ciclo) %>% 
      mutate(iteration = 1:n())
    p2_interes[[ciclo]]  <- sims %>% select(starts_with("p2")) %>% 
      mutate(Ciclo = ciclo) %>% 
      mutate(iteration = 1:n())
    p3_interes[[ciclo]]  <- sims %>% select(starts_with("p3")) %>% 
      mutate(Ciclo = ciclo) %>% 
      mutate(iteration = 1:n())
    p4_interes[[ciclo]]  <- sims %>% select(starts_with("p4")) %>% 
      mutate(Ciclo = ciclo) %>% 
      mutate(iteration = 1:n())
    ajus_p4_interes[[ciclo]]  <- sims %>% select(starts_with("ajus"))  %>% 
      mutate(Ciclo = ciclo) %>%
      mutate(iteration = 1:n())
  }
  
  #Calculamos los p1
  p1vals           <- do.call("rbind", p1_interes)
  colnames(p1vals) <- c(stratum$strata, "Ciclo", "iteration")
  p1vals <- p1vals %>% pivot_longer(1:length(ciclos), names_to = "strata", values_to = "p1") 
  p1vals <- p1vals %>% left_join(pondef, by = c("strata","Ciclo"))
  
  totales <- p1vals %>% ungroup() %>%
    group_by(iteration, Ciclo) %>%
    dplyr::summarise(p1 = sum(p1*sampling_weights), .groups = "keep")  %>%
    mutate(strata = "Total")
  
  p1vals <- p1vals %>% bind_rows(totales) %>% select(-sampling_weights)
  
  #Calculamos los p2
  p2vals <- do.call("rbind", p2_interes)
  colnames(p2vals) <- c(stratum$strata, "Ciclo", "iteration")
  p2vals <- p2vals %>% pivot_longer(1:length(ciclos), names_to = "strata", values_to = "p2")
  p2vals <- p2vals %>% left_join(pondef, by = c("strata","Ciclo"))
  
  totales <- p2vals %>% ungroup() %>%
    group_by(iteration, Ciclo) %>%
    dplyr::summarise(p2 = sum(p2*sampling_weights), .groups = "keep")  %>%
    mutate(strata = "Total")
  
  p2vals <- p2vals %>% bind_rows(totales) %>% select(-sampling_weights)
  
  #Calculamos los p3
  p3vals <- do.call("rbind", p3_interes)
  colnames(p3vals) <- c(stratum$strata, "Ciclo", "iteration")
  p3vals <- p3vals %>% pivot_longer(1:length(ciclos), names_to = "strata", values_to = "p3")
  p3vals <- p3vals %>% left_join(pondef, by = c("strata","Ciclo"))
  
  totales <- p3vals %>% ungroup() %>%
    group_by(iteration, Ciclo) %>%
    dplyr::summarise(p3 = sum(p3*sampling_weights), .groups = "keep")  %>%
    mutate(strata = "Total")
  
  p3vals <- p3vals %>% bind_rows(totales) %>% select(-sampling_weights)
  
  #Calculamos los p3
  p4vals <- do.call("rbind", p4_interes)
  colnames(p4vals) <- c(stratum$strata, "Ciclo", "iteration")
  p4vals <- p4vals %>% pivot_longer(1:length(ciclos), names_to = "strata", values_to = "p4")
  p4vals <- p4vals %>% left_join(pondef, by = c("strata","Ciclo"))
  
  totales <- p4vals %>% ungroup() %>%
    group_by(iteration, Ciclo) %>%
    dplyr::summarise(p4 = sum(p4*sampling_weights), .groups = "keep")  %>%
    mutate(strata = "Total")
  
  p4vals <- p4vals %>% bind_rows(totales) %>% select(-sampling_weights)
  
  #Calculamos los qvals
  ajus_p4vals <- do.call("rbind", ajus_p4_interes)
  colnames(ajus_p4vals) <- c(stratum$strata, "Ciclo", "iteration")
  ajus_p4vals <- ajus_p4vals %>% pivot_longer(1:length(ciclos), names_to = "strata", values_to = "ajus_p4")
  ajus_p4vals <- ajus_p4vals %>% left_join(pondef, by = c("strata","Ciclo"))
  
  totales <- ajus_p4vals %>% ungroup() %>%
    group_by(iteration, Ciclo) %>%
    dplyr::summarise(ajus_p4 = sum(ajus_p4*sampling_weights), .groups = "keep")  %>%
    mutate(strata = "Total")
  
  ajus_p4vals <- ajus_p4vals %>% bind_rows(totales) %>% select(-sampling_weights)
  
  #Calculamos las medias
  medias_p1 <- p1vals %>% group_by(Ciclo, strata) %>% mean_hdi(p1) %>%
    rename(p = p1) %>%
    select(Ciclo, strata, p, .lower, .upper, .width) %>%
    mutate(`Estimador` = 
             paste0(percent(p, 0.01)," [", 
                    percent(.lower, 0.01, suffix = ""), ", ",
                    percent(.upper, 0.01, suffix = ""), "]")) %>%
    mutate(.width = percent(.width)) %>%
    rename(Confianza = .width) %>%
    mutate(Tipo = "No ajustado") %>% 
    mutate(Estimando = p1name)
  
  medias_p2 <- p2vals %>% group_by(Ciclo, strata) %>% mean_hdi(p2) %>%
    rename(p = p2) %>%
    select(Ciclo, strata, p, .lower, .upper, .width) %>%
    mutate(`Estimador` = 
             paste0(percent(p, 0.01)," [", 
                    percent(.lower, 0.01, suffix = ""), ", ",
                    percent(.upper, 0.01, suffix = ""), "]")) %>%
    mutate(.width = percent(.width)) %>%
    rename(Confianza = .width) %>%
    mutate(Tipo = "No ajustado") %>% 
    mutate(Estimando = p2name)
  
  medias_p3 <- p3vals %>% group_by(Ciclo, strata) %>% mean_hdi(p3) %>%
    rename(p = p3) %>%
    select(Ciclo, strata, p, .lower, .upper, .width) %>%
    mutate(`Estimador` = 
             paste0(percent(p, 0.01)," [", 
                    percent(.lower, 0.01, suffix = ""), ", ",
                    percent(.upper, 0.01, suffix = ""), "]")) %>%
    mutate(.width = percent(.width)) %>%
    rename(Confianza = .width) %>%
    mutate(Tipo = "No ajustado") %>% 
    mutate(Estimando = p3name)
  
  medias_p4 <- p4vals %>% group_by(Ciclo, strata) %>% mean_hdi(p4) %>%
    rename(p = p4) %>%
    select(Ciclo, strata, p, .lower, .upper, .width) %>%
    mutate(`Estimador` = 
             paste0(percent(p, 0.01)," [", 
                    percent(.lower, 0.01, suffix = ""), ", ",
                    percent(.upper, 0.01, suffix = ""), "]")) %>%
    mutate(.width = percent(.width)) %>%
    rename(Confianza = .width) %>%
    mutate(Tipo = "No ajustado") %>% 
    mutate(Estimando = p4name)
  
  medias_ajus_p4 <- ajus_p4vals %>% group_by(Ciclo, strata) %>% 
    mean_hdi(ajus_p4) %>%
    rename(p = ajus_p4) %>%
    select(Ciclo, strata, p, .lower, .upper, .width) %>%
    mutate(`Estimador` = 
             paste0(percent(p, 0.01)," [", 
                    percent(.lower, 0.01, suffix = ""), ", ",
                    percent(.upper, 0.01, suffix = ""), "]")) %>%
    mutate(.width = percent(.width)) %>%
    rename(Confianza = .width) %>%
    mutate(Tipo = "Ajustado") %>% 
    mutate(Estimando = p4name)
  
  
  resultados <- medias_p1 %>% bind_rows(medias_p2) %>% 
    bind_rows(medias_p3) %>% bind_rows(medias_p4) %>% 
    bind_rows(medias_ajus_p4) %>% 
    mutate(Variable = !!varname) %>% 
    select(Variable, Estimando, Tipo, Ciclo, strata, Estimador, Confianza, p, .lower, .upper) 
  
  #resultados %>% write_csv(paste0(filedir,varname,".csv"))
  return(resultados)
}



genera_binomial_3 <- function(nfilter, s1filter, s2filter, s3filter, s4filter, 
                              varname = "VARIABLE",
                              p1name = c("p1","p1"), 
                              p2name = c("p2","p2"), 
                              p3name = c("p3","p3"),
                              p4name = c("p3","p3"),
                              mu_p1 = 0.5, sigma_p1 = 1/12, 
                              mu_p2 = rep(0.5, 2), sigma_p2 = rep(1/12, 2),
                              mu_p3 = rep(0.5, 2), sigma_p3 = rep(1/12, 2),
                              mu_p4 = rep(0.5, 2), sigma_p4 = rep(1/12, 2),
                              filedir = getwd(),
                              seeds = NULL,
                              delta = 1, gamma = 1, 
                              burnin = 5000, niter = 10000, nchains = 4){
  
  
  resultados_1 <- genera_binomial_3_iter(nfilter = nfilter, 
                                         s1filter = s1filter, 
                                         s2filter = s2filter, 
                                         s3filter = s3filter,
                                         s4filter = s4filter,
                                         varname = varname, 
                                         p1name = p1name[1], 
                                         p2name = p2name[1],
                                         p3name = p3name[1],
                                         p4name = p4name[1],
                                         seeds = seeds,
                                         mu_p1 = mu_p1, sigma_p1 = sigma_p1,
                                         mu_p2 = mu_p2[1], sigma_p2 = sigma_p2[1],
                                         mu_p3 = mu_p3[1], sigma_p3 = sigma_p3[1],
                                         mu_p4 = mu_p4[1], sigma_p4 = sigma_p4[1],
                                         delta = delta, gamma = gamma, 
                                         burnin = burnin, niter = niter, 
                                         nchains = nchains)
  
  
  resultados <- resultados_1
  
  resultados %>% write_csv(paste0(filedir,varname,".csv"))
  return(resultados)
  
}
