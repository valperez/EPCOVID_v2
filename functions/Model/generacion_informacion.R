#rm(list = ls()) #clear all
library(ggplot2)
library(ggthemes)
library(cowplot) #install.packages("cowplot")
library(combinat) #install.packages("combinat")
library(rBeta2009)
library(tidyverse)
library(dplyr)
library(readxl)
library(data.table)
library(kableExtra)
set.seed(99)

# Tenemos el problema de tipo 
# k ~ Bin(n, p)
# p ~ Beta(alpha_p, beta_p)
# t ~ Bin(k, q)
# q ~ Beta(alpha_q, beta_q)

# Función que obtiene la mu y la sigma que necesitas
# df_aprioris es el excel (o donde sea que estén los datos) donde viene las mu
# y sigma a prioris
# vars_prod son las variables de las que hay que sacarle el producto
# vars_division son las variable a la que hay que usar para divir mu final
# Ojo que depende cada caso si hay division o no
mu_y_sigma <- function(df_aprioris, vars_prod, vars_division = NULL, nombre_nueva){
  
  nombre_dada <- c()
  for (i in 1:(length(vars_prod))){           
    nombre_dada[i]           <- vars_prod[i]
  }
  
  string_parametros <- ""
  
  for (i in 1:(length(nombre_dada) - 1)){
    string_parametros <- paste0(string_parametros, " Variable == '", nombre_dada[i], "' | ")
  }
  
  string_parametros <- paste0(string_parametros, "Variable == '", 
                              nombre_dada[length(nombre_dada)],
                              "'")
  
  parametros_dados <- df_aprioris %>%
    filter(eval(parse(text = string_parametros))) %>%
    mutate(`Proporción esperada (resultado)` = 
             as.numeric(`Proporción esperada (resultado)`))
  
  aux <- prod(parametros_dados$`Proporción esperada (resultado)`)
  
  if (!is.null(vars_division)){
  
    string_parametros <- paste0("Variable == '", vars_division, "'")
    
    parametros_div <- df_aprioris %>%
      filter(eval(parse(text = string_parametros))) %>%
      mutate(`Proporción esperada (resultado)` = 
               as.numeric(`Proporción esperada (resultado)`))
    aux_div <- prod(parametros_div$`Proporción esperada (resultado)`)
    
    var_nueva <- aux/aux_div
  } else {
    var_nueva <- aux
  }
  
  auxiliar <- c(nombre_nueva, var_nueva, var_nueva*(1-var_nueva)/2, NA, NA)
  
  df_aprioris <- rbind(df_aprioris, auxiliar)
  
  return(df_aprioris)
  
}


# Función que obtiene la mu y la sigma que necesitas
# df_aprioris es el excel (o donde sea que estén los datos) donde viene las mu
# y sigma a prioris
# vars_prod son las variables de las que hay que sacarle el producto
# vars_division son las variable a la que hay que usar para divir mu final
# Ojo que depende cada caso si hay division o no
busqueda_mu_y_sigma <- function(df_aprioris, variables){

  
  string_parametros <- ""
  
  if(length(variables) > 1){
    for (i in 1:(length(variables) - 1)){
      string_parametros <- paste0(string_parametros, " Variable == '", variables[i], "' | ")
    }
    
    string_parametros <- paste0(string_parametros, "Variable == '", 
                                variables[length(variables)], "'")
  } else {
    string_parametros <- paste0(string_parametros, "Variable == '", variables, "' ")
  }
  
  parametros_dados <- df_aprioris %>%
    filter(eval(parse(text = string_parametros))) %>%
    select(Variable, `Proporción esperada (resultado)`, `Error +- de la proporción esperada (calculo uno en automático pero relléname si esto no te gusta)`) %>%
    rename(Mu = `Proporción esperada (resultado)`) %>%
    rename(Sigma = `Error +- de la proporción esperada (calculo uno en automático pero relléname si esto no te gusta)`) %>%
    mutate(Mu = as.numeric(Mu)) %>%
    mutate(Sigma = as.numeric(Sigma)) %>%
    select(-Variable) %>%
    as.matrix()
  
  return(parametros_dados)
  
}


# Funcion para obtener s y t 
# función para obtener a la gente que contestó que sí o de forma positiva al 
# aspecto k y
# obtener la gente que contestó que sí o de forma positiva al aspecto de k + el
# aspecto de t en los datos en crudo de la encuesta Panel
# Es importante que la primera variable que le damos sea la t tanto en vars_analisis
# como en igualdades_variables
k_y_t <- function(df_informacion, vars_analisis, igualdades_variables){
  variable_condicionada <- vars_analisis[1]
  igualdad_condicionada <- igualdades_variables[1]

  variable_dada <- c()
  igualdad_dada <- c()
  for (i in 2:length(vars_analisis)){            #Esto me hace que la primera entrada sea NA 
    variable_dada[i - 1]         <- vars_analisis[i]
    igualdad_dada[i - 1]         <- igualdades_variables[i]
  }
  
  string_aux <- ""
  
  if (length(variable_dada) > 1) {
    for (i in 1:(length(variable_dada)-1)){
      string_aux <- paste0(string_aux, variable_dada[i], " == '", igualdad_dada[i], "' & ")
    }
  }
    string_aux <- paste0(string_aux, 
                         variable_dada[length(variable_dada)], " == '", igualdad_dada[length(variable_dada)], "'")
    
  
  var_dada <- df_informacion %>%
    filter(!is.na(!!sym(variable_dada[1])))
  
  if (length(variable_dada) > 1){
    for (i in 2:length(variable_dada)){
      var_dada <- var_dada %>%
        filter(!is.na(!!sym(variable_dada[i])))
    }
  }
  
  # La gente que tuvo la variable dada
  # Ya está correcto
  var_dada <- var_dada %>%
    filter(eval(parse(text = string_aux))) %>%
    group_by(clasificacion) %>%
    summarise(n()) %>%
    right_join(centros, by = "clasificacion") %>%
    rename(`variables_dadasK` = `n()`) %>%
    mutate(`variables_dadasK` = replace_na(`variables_dadasK`, 0))
  
  string_aux <- paste0(string_aux, "& ", 
                       variable_condicionada, " == '", igualdad_condicionada, "'")
  
  comb_vars <- df_informacion %>%
    filter(!is.na(!!sym(variable_condicionada)))
  
  for (i in 1:length(variable_dada)){
    comb_vars <- comb_vars %>%
      filter(!is.na(!!sym(variable_dada[i])))
  }
  
  # Los que tuvieron variable_condicionada | variable_dada
  comb_vars <- comb_vars %>%
    filter(eval(parse(text = string_aux))) %>%
    group_by(clasificacion) %>%
    summarise(n()) %>%
    right_join(centros, by = "clasificacion") %>%
    rename(`k_dadot` = `n()`) %>%
    mutate(`k_dadot` = replace_na(`k_dadot`, 0))
  
  #Unimos ambas tablas y arreglamos de acuerdo a la clasificación
  s <- full_join(var_dada, comb_vars, by = "clasificacion")
  s <- s %>% arrange(clasificacion) %>% select(-clasificacion) %>% as.matrix()
  
  return(s)
}

# Funcion para obtener cualquier s que queramos
# donde vector_na es un vector para cuando necesitemos que calcule na's de algunas variables
# e string_igualdades para cuando necesitemos igualar algo
sizes <- function(df_informacion, vector_na = NULL, string_igualdades = NULL){
  
  resultado_na <- df_informacion %>%
    filter(!is.na(!!sym(vector_na[1])))
  
  if (length(vector_na) > 1){
    for (i in 2:length(vector_na)){
      resultado_na <- resultado_na %>%
        filter(!is.na(!!sym(vector_na[i])))
    }
  }
  
  # Lo guardo porque necesito que quite los NA's
  aux <- resultado_na
  
  resultado_na <- resultado_na %>%
    group_by(clasificacion) %>%
    summarise(n()) %>%
    right_join(centros, by = "clasificacion") %>%
    rename(`resultado_na` = `n()`) %>%
    mutate(`resultado_na` = replace_na(`resultado_na`, 0))
  
  
  resultado_igualdades <- aux %>%
    filter(eval(parse(text = string_igualdades))) %>%
    group_by(clasificacion) %>%
    summarise(n()) %>%
    right_join(centros, by = "clasificacion") %>%
    rename(`resultado_igualdades` = `n()`) %>%
    mutate(`resultado_igualdades` = replace_na(`resultado_igualdades`, 0))
  
  #Unimos ambas tablas y arreglamos de acuerdo a la clasificación
  s <- full_join(resultado_na, resultado_igualdades, by = "clasificacion")
  s <- s %>% arrange(clasificacion) %>% select(-clasificacion) %>% as.matrix()
  
  return(s)
}

# Función que genera la información para una sola binomial
# donde df_informacion es la base de datos del ciclo que estás haciendo
# y df_aprioris es el data frame de información (equivalente al que nos dió David)
# variable_analisis es la variable de la que estás haciendo el análisis
# nombre_apriori es como se llama la variable en el excel de David
# igualdades_variables es a lo que cada una de las variable_analisis 
# hay que igualarlas en el filter
informacion_unabin <- 
  function(df_informacion, df_aprioris, 
           variable_analisis, nombre_apriori, 
           igualdades_variables){
    
    parametros <- df_aprioris %>%
      filter(Variable == nombre_apriori)
    
    mu <- as.numeric(parametros$`Proporción esperada (resultado)`)
    
    sigma2 <- as.numeric(parametros$`Error +- de la proporción esperada (calculo uno en automático pero relléname si esto no te gusta)`)
    
    s <- df_informacion %>%
      filter(!!sym(variable_analisis) == igualdades_variables) %>%
      group_by(clasificacion) %>%
      summarise(n()) %>%
      right_join(centros, by = "clasificacion") %>%
      rename(`variable_analisis` = `n()`) %>%
      mutate(`variable_analisis` = replace_na(`variable_analisis`, 0)) %>%
      arrange(clasificacion) %>% 
      select(-clasificacion) %>% 
      as.matrix()
    
    gente_encuestada <- df_informacion %>%
      filter(!is.na(!!sym(variable_analisis))) %>%
      group_by(clasificacion) %>%
      tally() %>%
      drop_na() %>%
      arrange(clasificacion) %>%
      select(-clasificacion) %>% 
      as.matrix()
    
    return(list(mu = mu, sigma2 = sigma2, s = s, gente_encuestada = gente_encuestada))
  }



# Función para generar la media y cuantiles de las simulaciones
# lista es la lista de simulaciones 
# ciclo es el número de ciclo del que se está haciendo la simulación 
# nombres son los nombres del centro de distribución que por default está dados 
# por la clasificación de los centros
# nombre_archivo es el nombre que quieres que se llame el archivo ej. IGG_dadoNecesidad
# informativo es si se utilizaron los datos de David o no
# ajustado es si está ajustado por sensitividad y especificidad
media_cuantiles_creacionArchivos <- function(lista, 
                                             ciclo, 
                                             nombres = centros$clasificacion,
                                             nombre_archivo,
                                             informativo,
                                             ajustado){
  
  aux <- matrix(NA, nrow = length(lista) + 1, ncol = 4)
  
    for (i in 1:length(lista)){
      aux [i, 1] <- mean(lista[[i]])
      aux [i, 2] <- quantile(lista[[i]], prob = 0.0)
      aux [i, 3] <- quantile(lista[[i]], prob = 0.95)
    }
  
  
  distr_bimbo <- lista[[1]]*pesos[1] + lista[[2]]*pesos[2] + 
    lista[[3]]*pesos[3]
  
  aux[4, 1] <- mean(distr_bimbo)
  aux[4, 2] <- quantile(distr_bimbo, prob = 0.0)
  aux[4, 3] <- quantile(distr_bimbo, prob = 0.95)
  aux[ , 4] <- ciclo
  
  aux <- data.frame(aux, 
                    row.names = c(nombres, "Total"))
  
  colnames(aux) <- c("mean", "q5", "q95", "Ciclo")
  
  if (informativo == T){
    nombre_archivo <- paste0(nombre_archivo, "_Informativa")
  }
  
  if (ajustado == T){
    nombre_archivo <- paste0(nombre_archivo, "_Ajustado")
  }
  
  nombre_archivo <- paste0(nombre_archivo, "Ciclo", ciclo, ".csv")
  
  write.csv(aux, file = nombre_archivo)
  
}




