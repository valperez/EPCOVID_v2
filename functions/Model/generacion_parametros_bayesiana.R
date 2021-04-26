rm(list = ls()) #clear all
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


setwd("~/IMSS/Bimbo")

#Esto es para las preguntas con la siguiente estructura:
# PREGUNTA 1: SI | NO
# Luego analizar PCR (digamos) condicional en que PREGUNTA == SI
# PCR no tiene nada de especial, puede ser cualquier otra cosa
# este demo no está ajustado por sensibilidad / especificidad

# # # # ---- MÉTODOS NO INFORMATIVOS ---- # # # #

# Supongamos que tenemos el problema 
# sexo ~ Binomial (n, p)
# p    ~ Beta(alpha_p, beta_p)
# t    ~ Binomal(s,q)
# q    ~ Beta(alpha_q, beta_q)



# # # # ---- MÉTODOS INFORMATIVOS ---- # # # #

# Vamos a hacer las distribuciones con a priori's informativas

#Original
ciclo1 <- read_delim("~/IMSS/Bimbo/EPCOVID_c1_20210225.txt", ";", 
                    escape_double = FALSE, locale = locale(encoding = "UTF-8"), 
                    trim_ws = TRUE)

# Lo lee con los headers mal
aprioris <- read_excel("Aprioris David.xlsx")
#Entonces lo arreglamos
colnames(aprioris) <- aprioris[1, ]
aprioris <- aprioris[-1, ]

# Ya tengo los centros
centros <- ciclo1 %>%
  select(clasificacion) %>%
  distinct() %>%
  drop_na() %>%
  arrange(clasificacion)

tamanios_totales <- read_csv("Tamaños totales de muestra 29 enero 2020.csv") %>%
  select(-X1) %>%
  mutate(Total = rowSums(.)) 

# Ahora obtenemos las ponderaciones para cada centro y ciclo
# Calculamos entre el total
pesos_totales <- tamanios_totales %>%
  mutate(`Centro de ventas` = `Centro de ventas`/`Total`) %>%
  mutate(`Centro de Distribución` = `Centro de Distribución` /`Total`) %>%
  mutate(`Planta de producción` = `Planta de producción` / `Total`) %>%
  t()%>%
  as.data.frame()
#Ponemos bien los nombres de columnas y obtenemos los nombres de los centros
colnames(pesos_totales) <- c("Ciclo_1", "Ciclo_2", "Ciclo_3")
setDT(pesos_totales, keep.rownames = T)
# Arreglamos para qe todo esté de acuerdo a la clasificación (los centros)
pesos_totales <- pesos_totales %>%
  rename(clasificacion = rn) %>%
  arrange(clasificacion)


#Seleccionamos el peso que necesitamos para el Ciclo 1
pesos <- pesos_totales$Ciclo_1 #Seleccionamos el ciclo 1
pesos <- pesos[1:3]  #Quitamos la columna de total

# Valores de sensitividad y especificidad de Gelman
# en este link: http://www.stat.columbia.edu/~gelman/research/published/specificity.pdf
sensitividad <- 0.8
specificidad <- 0.995

setwd("~/IMSS/Bimbo/Resultados")

#Función que toma el archivo de David y utiliza la fórmula de Bayes 
# para cambiar de x|y a y|x
# nombres_apriori es un parámetro que te dice como se llaman las variables en 
# el excel de David. Siempre va primero la condicional que ya tenemos
# Nombres_apriori le tengo que dar en orden lo que necesito, 
# al final dándole entre lo que voy a dividir
bayes_aprioris <- function(nombres_apriori, nombre_nueva, aprioris){
  
  nombre_condicionada   <- nombres_apriori[1]
  nombre_division       <- nombres_apriori[length(nombres_apriori)]
  
  nombre_dada <- c()
  for (i in 2:(length(nombres_apriori)-1)){            #Esto me hace que la primera entrada sea NA 
    nombre_dada[i - 1]           <- nombres_apriori[i]
  }
  
  string_parametros <- paste0("Variable == '", nombre_condicionada, "'")
  
    for (i in 1:(length(nombre_dada))){
      string_parametros <- paste0(string_parametros, " | Variable == '", nombre_dada[i], "'")
    }
  
  parametros_dados <- aprioris %>%
    filter(eval(parse(text = string_parametros))) %>%
    mutate(`Proporción esperada (resultado)` = 
             as.numeric(`Proporción esperada (resultado)`))
  
  aux <- prod(parametros_dados$`Proporción esperada (resultado)`)
  
  string_parametros <- paste0("Variable == '", nombre_division, "'")
  
  parametros_div <- aprioris %>%
    filter(eval(parse(text = string_parametros))) %>%
    mutate(`Proporción esperada (resultado)` = 
             as.numeric(`Proporción esperada (resultado)`))
  aux_div <- prod(parametros_div$`Proporción esperada (resultado)`)
  
  var_nueva <- aux/aux_div
  
  auxiliar <- c(nombre_nueva, var_nueva, var_nueva*(1-var_nueva)/2, NA, NA)
  
  aprioris <- rbind(aprioris, auxiliar)
  
  return(aprioris)
  
}

aprioris <- bayes_aprioris(nombres_apriori = c("Necesidad | IGG", "IGG", "Necesidad en salud 3 meses"), 
                           nombre_nueva = "IGG | Necesidad")
aprioris <- bayes_aprioris(nombres_apriori = c("Sintomas | IGG", "IGG", "Necesidad en salud 3 meses"),
                           nombre_nueva = "IGG | Sintomas")
aprioris <- bayes_aprioris(nombres_apriori = c("Buscó | IGG Y Necesitó", "Necesidad | IGG", "IGG", "Busqueda atención"),
                           nombre_nueva = "IGG y Necesitó | Buscó")
aprioris <- bayes_aprioris(nombres_apriori = c("Buscó atención | Síntomas e IGG", "Sintomas | IGG", "IGG", "Busqueda atención"),
                           nombre_nueva = "Sintomas e IGG | Buscó")


# Función para responder una sola pregunta binomial, ej PCR+ ó Hombres
# mu y sigma2 es lo que nos dió David,
# s es la cantidad de gente que respondió por centro de trabajo
# pesos es la ponderación por centro de trabajo
# ns es la cantidad de gente encuestada por centro de trabajo

# OJO : la sigma2 en el excel de David está divida entre 2 porque si no, no funciona
una_binomial <- function(s, pesos, ns, 
                         mu = NULL, sigma2 = NULL, 
                         nsim = 5000,
                         nombres, 
                         sens = 1, spec = 1){
  if (is.null(mu)){
    mu = 0.5
  }
  
  if (is.null(sigma2)){
    sigma2 <- 1/12
  }
  if (length(s) != length(pesos) | length(pesos) != length(ns)){
    stop("Las longitudes están mal")
  }
  
  if (mu < 0 | mu > 1){
    stop("el promedio debe estar entre 0 y 1!")
  }
  
  if (sigma2 <= 0 | sigma2 >= mu*(1-mu)){
    stop("la varianza está maaaaal")
  }
  
  alpha <- -(mu*(sigma2 + mu^2 - mu))/sigma2
  beta <- (sigma2 + mu^2 - mu)*(mu - 1)/sigma2
  
  distr_beta <- list()
  
  for (i in 1:length(ns)){
    distr_beta[[i]] <- rbeta(nsim, alpha + s[i], beta + ns[i] - s[i])
  }
  
  for (i in 1:length(ns)){
    distr_beta[[i]] <- distr_beta[[i]]*sens + (1 - spec)*(1- distr_beta[[i]])
  }

  aux <- matrix(NA, nrow = length(ns) + 1, ncol = 3)
  
  for (i in 1:length(ns)){
    aux [i, 1] <- mean(distr_beta[[i]])
    aux [i, 2] <- quantile(distr_beta[[i]], prob = 0.0)
    aux [i, 3] <- quantile(distr_beta[[i]], prob = 0.95)
  }
  
  distr_bimbo <- distr_beta[[1]]*pesos[1] + distr_beta[[2]]*pesos[2] + 
    distr_beta[[3]]*pesos[3]
  
  aux[4, 1] <- mean(distr_bimbo)
  aux[4, 2] <- quantile(distr_bimbo, prob = 0.0)
  aux[4, 3] <- quantile(distr_bimbo, prob = 0.95)

  aux <- data.frame(aux, 
                    row.names = c(nombres, "Total"))
  
  colnames(aux) <- c("mean", "q5", "q95")
  
  return(aux)
  
}


# Función para responder una dos preguntas binomiales condicionales, ej Necesidad | IGG
# mu y sigma2 es lo que nos dió David (vectores ambos),
# s es la cantidad de gente que respondió a la pregunta 1 y pregunta 2 por centro de trabajo (matriz)
# osea en la columna_i de s van a venir las personas que respondieron a la pregunta_i
# pesos es la ponderación por centro de trabajo (vector)
# ns es la cantidad de gente encuestada por centro de trabajo (vector)

# OJO : la sigma2 en el excel de David está divida entre 2 porque si no, no funciona

dos_binomial <- function(s, pesos, ns, 
                         mu = NULL, sigma2 = NULL, 
                         nsim = 5000, 
                         nombres,
                         sens = 1, spec = 1){
  if (is.null(mu)){
    mu <- rep(0.5, ncol(s)) 
  }
  
  if (is.null(sigma2)){
    sigma2 <- rep(1/12, ncol(s)) 
  }
  if (nrow(s) != length(pesos) | length(pesos) != nrow(ns)){
    stop("Las longitudes están mal")
  }
  
  if (all(mu < 0) | all(mu > 1)){
    stop("el promedio debe estar entre 0 y 1!")
  }
  
  if (all(sigma2 <= 0) | all(sigma2 >= mu*(1-mu))){
    stop("la varianza está maaaaal")
  }
  
  alpha <- -(mu*(sigma2 + mu^2 - mu))/sigma2
  beta <- (sigma2 + mu^2 - mu)*(mu - 1)/sigma2
  
  distr_beta <- list()
  vec_beta <- list()

  # Quiero generar una lista de listas 
  # donde el primer renglón es la lista de simulaciones para la primera pregunta
  # y el segundo renglón es la lista de simulaciones para la segunda pregunta
  
  for (j in 1:length(alpha)){
    for (i in 1:nrow(ns)){
        distr_beta[[i]] <- rbeta(nsim, 
                                 alpha[j] + s[i, 2], 
                                 beta[j] + s[i, 1] - s[i, 2])
    }
    vec_beta[[j]] <- distr_beta
  }
  
  
  if (sens != 1 | spec != 1){
    for (i in 1:nrow(ns)){
      vec_beta[[1]][[i]] <- vec_beta[[1]][[i]]*sens + 
        (1 - spec)*(1 - vec_beta[[1]][[i]])
    }
  }
  
  
  vec_beta <- Map("*", vec_beta[[1]], vec_beta[[2]])
  
  aux <- matrix(NA, nrow = nrow(ns) + 1, ncol = 3)
  
  for (j in 1:length(alpha)){
    for (i in 1:nrow(ns)){
      aux [i, 1] <- mean(vec_beta[[i]])
      aux [i, 2] <- quantile(vec_beta[[i]], prob = 0.0)
      aux [i, 3] <- quantile(vec_beta[[i]], prob = 0.95)
    }
  }
  
  distr_bimbo <- vec_beta[[1]]*pesos[1] + vec_beta[[2]]*pesos[2] + 
    vec_beta[[3]]*pesos[3]
  
  aux[4, 1] <- mean(distr_bimbo)
  aux[4, 2] <- quantile(distr_bimbo, prob = 0.0)
  aux[4, 3] <- quantile(distr_bimbo, prob = 0.95)
  
  aux <- data.frame(aux, 
                    row.names = c(nombres, "Total"))
  
  colnames(aux) <- c("mean", "q5", "q95")
  
  return(aux)
  
  }


# Función que genera las distribuciones de manera automática #quefancy
# variable_analisis es la variable de la que estás haciendo el análisis
# nombre_apriori es como se llama la variable en el excel de David
# igualdades_variables es a lo que cada una de las variable_analisis 
# hay que igualarlas en el filter
generacion_distribucion_una_bin <- 
  function(variable_analisis, nombre_apriori, igualdades_variables){
  
  parametros <- aprioris %>%
    filter(Variable == nombre_apriori)
  
  mu <- as.numeric(parametros$`Proporción esperada (resultado)`)
  
  sigma2 <- as.numeric(parametros$`Error +- de la proporción esperada (calculo uno en automático pero relléname si esto no te gusta)`)
  
  s <- ciclo1 %>%
    filter(!!sym(variable_analisis) == igualdades_variables) %>%
    group_by(clasificacion) %>%
    summarise(n()) %>%
    right_join(centros, by = "clasificacion") %>%
    rename(`variable_analisis` = `n()`) %>%
    mutate(`variable_analisis` = replace_na(`variable_analisis`, 0)) %>%
    arrange(clasificacion) %>% 
    select(-clasificacion) %>% 
    as.matrix()
  
  gente_encuestada <- ciclo1 %>%
    filter(!is.na(!!sym(variable_analisis))) %>%
    group_by(clasificacion) %>%
    tally() %>%
    drop_na() %>%
    arrange(clasificacion) %>%
    select(-clasificacion) %>% 
    as.matrix()
  
  out_informativa <- una_binomial(s = s, 
                                  pesos = pesos, 
                                  ns = gente_encuestada, 
                                  mu = mu, 
                                  sigma2 = sigma2, 
                                  nsim = 5000, 
                                  nombres = centros$clasificacion)
  
  out_noinformativa <- una_binomial(s = s, 
                                    pesos = pesos, 
                                    ns = gente_encuestada, 
                                    mu = NULL, 
                                    sigma2 = NULL, 
                                    nsim = 5000, 
                                    nombres = centros$clasificacion)
  
  out_informativa_ajustada <- una_binomial(s = s, 
                                           pesos = pesos, 
                                           ns = gente_encuestada, 
                                           mu = mu, 
                                           sigma2 = sigma2, 
                                           nsim = 5000, 
                                           nombres = centros$clasificacion,
                                           sens = sensitividad, spec = specificidad)
  out_noinformativa_ajustada <- una_binomial(s = s, 
                                             pesos = pesos, 
                                             ns = gente_encuestada, 
                                             mu = NULL, 
                                             sigma2 = NULL, 
                                             nsim = 5000, 
                                             nombres = centros$clasificacion,
                                             sens = sensitividad, spec = specificidad)
  
  write.csv(out_informativa, file = paste0(nombre_apriori, "_Informativa.csv"))
  write.csv(out_noinformativa, file = paste0(nombre_apriori, "_NoInformativa.csv"))
  write.csv(out_informativa_ajustada, file = paste0(nombre_apriori, "_Informativa_Ajustada.csv"))
  write.csv(out_noinformativa_ajustada, file = paste0(nombre_apriori, "_NoInformativa_Ajustada.csv"))
  }

# Función que genera la distribución para dos binomiales (orita generalizamos)
# vars_analisis son las variables de las que estás haciendo el análisis según su nombre
# en el raw data
# igualdades_vars son las igualdades a las que hay que ponerle a las variables de análisis
# en el raw data
# nombres_apriori son los nombres de las variables que queremos en la tabla de David

generacion_distribucion_dos_bin<- function(vars_analisis, 
                                           igualdades_vars, 
                                           nombres_apriori, ajuste_sen = F){
  
  variable_condicionada <- vars_analisis[1]
  variable_dada         <- vars_analisis[2]
  
  igualdad_condicionada <- igualdades_vars[1]
  igualdad_dada         <- igualdades_vars[2]
  
  nombre_condicionada   <- nombres_apriori[1]
  nombre_dada           <- nombres_apriori[2]
  
  aux <- strsplit(nombre_condicionada, "[|]")
  nombres_cond <- paste0(aux[[1]][1], "_dado", aux[[1]][2])
  nombres_cond <- str_remove_all(nombres_cond, " ")
  
  # La gente que tuvo la variable dada
  var_dada <- ciclo1 %>%
    filter(!!sym(variable_dada) == igualdad_dada) %>%
    group_by(clasificacion) %>%
    summarise(n()) %>%
    right_join(centros, by = "clasificacion") %>%
    rename(variable_dada = `n()`) %>%
    mutate(variable_dada = replace_na(variable_dada, 0))
  
  # Los que tuvieron variable_condicionada | variable_dada
  comb_vars <- ciclo1 %>%
    filter(!!sym(variable_condicionada) == igualdad_condicionada 
           & !!sym(variable_dada) == igualdad_dada) %>%
    group_by(clasificacion) %>%
    summarise(n()) %>%
    right_join(centros, by = "clasificacion") %>%
    rename(nombre_condicionada = `n()`) %>%
    mutate(nombre_condicionada = replace_na(nombre_condicionada, 0))
  
  #Unimos ambas tablas y arreglamos de acuerdo a la clasificación
  s <- full_join(var_dada, comb_vars, by = "clasificacion")
  s <- s %>% arrange(clasificacion) %>% select(-clasificacion) %>% as.matrix()
  
  # Es la gente que responde la primera pregunta 
  gente_encuestada <- ciclo1 %>%
    filter(!is.na(!!sym(variable_dada))) %>%
    group_by(clasificacion) %>%
    tally() %>%
    drop_na() %>%
    arrange(clasificacion)
  
  parametros <- aprioris %>%
    filter(Variable == nombre_condicionada | Variable == nombre_dada)
  
  mu <- as.numeric(parametros$`Proporción esperada (resultado)`)
  
  sigma2 <- as.numeric(parametros$`Error +- de la proporción esperada (calculo uno en automático pero relléname si esto no te gusta)`)
  
  out_noinformativa <- dos_binomial(s = s, 
                                    pesos = pesos,
                                    ns = gente_encuestada, 
                                    mu = NULL, 
                                    sigma2 = NULL, 
                                    nsim = 5000, 
                                    nombres = centros$clasificacion)
  
  
  out_informativa <- dos_binomial(s = s, 
                                  pesos = pesos,
                                  ns = gente_encuestada, 
                                  mu = mu, 
                                  sigma2 = sigma2, 
                                  nsim = 500, 
                                  nombres = centros$clasificacion)
  
  write.csv(out_informativa, file = paste0(nombres_cond, "_Informativa.csv"))
  write.csv(out_noinformativa, file = paste0(nombres_cond, "_NoInformativa.csv"))
  
  if(ajuste_sen == T){
  
  out_informativa_ajustada <- dos_binomial(s = s, 
                                           pesos = pesos, 
                                           ns = gente_encuestada, 
                                           mu = mu, 
                                           sigma2 = sigma2, 
                                           nsim = 5000, 
                                           nombres = centros$clasificacion,
                                           sens = sensitividad, spec = specificidad)
  
  out_noinformativa_ajustada <- dos_binomial(s = s, 
                                             pesos = pesos, 
                                             ns = gente_encuestada, 
                                             mu = NULL, 
                                             sigma2 = NULL, 
                                             nsim = 5000, 
                                             nombres = centros$clasificacion,
                                             sens = sensitividad, spec = specificidad)
  
  
  write.csv(out_informativa_ajustada, file = paste0(nombres_cond, "_Informativa_Ajustada.csv"))
  write.csv(out_noinformativa_ajustada, file = paste0(nombres_cond, "_NoInformativa_Ajustada.csv"))
  }
  
}

# Función que genera la distribución para dos binomiales (orita generalizamos)
# vars_analisis son las variables de las que estás haciendo el análisis según su nombre
# en el raw data
# igualdades_vars son las igualdades a las que hay que ponerle a las variables de análisis
# en el raw data
# nombres_apriori son los nombres de las variables que queremos en la tabla de David
generacion_distribucion_varias_bin<- function(vars_analisis, igualdades_vars, nombres_apriori){
  
  variable_condicionada <- vars_analisis[1]
  igualdad_condicionada <- igualdades_vars[1]
  nombre_condicionada   <- nombres_apriori[1]
  
  variable_dada <- c()
  igualdad_dada <- c()
  nombre_dada <- c()
  for (i in 2:length(vars_analisis)){            #Esto me hace que la primera entrada sea NA 
    variable_dada[i - 1]         <- vars_analisis[i]
    igualdad_dada[i - 1]         <- igualdades_vars[i]
    nombre_dada[i - 1]           <- nombres_apriori[i]
  }
  
  nombres_cond <- str_remove_all(nombre_condicionada, " ")
  nombres_cond <- gsub("\\|", "_dado", nombres_cond)
  
  string_aux <- ""
  
  for (i in 1:(length(variable_dada)-1)){
    string_aux <- paste0(string_aux, variable_dada[i], " == '", igualdad_dada[i], "' & ")
  }
  string_aux <- paste0(string_aux, 
                       variable_dada[length(variable_dada)], " == '", igualdad_dada[length(variable_dada)], "'")
  
  var_dada <- ciclo1 %>%
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
    rename(`variables_dadas` = `n()`) %>%
    mutate(`variables_dadas` = replace_na(`variables_dadas`, 0))
  
  string_aux <- paste0(string_aux, "& ", 
                       variable_condicionada, " == '", igualdad_condicionada, "'")
  
  comb_vars <- ciclo1 %>%
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
    rename(nombre_condicionada = `n()`) %>%
    mutate(nombre_condicionada = replace_na(nombre_condicionada, 0))
  
  #Unimos ambas tablas y arreglamos de acuerdo a la clasificación
  s <- full_join(var_dada, comb_vars, by = "clasificacion")
  s <- s %>% arrange(clasificacion) %>% select(-clasificacion) %>% as.matrix()
  
  gente_encuestada <- ciclo1 
  
  # Es la gente que responde la primera pregunta 
  # Solo es saber si es NA o no
  # FIXME
  for (i in 1:length(variable_dada)){
    gente_encuestada <- gente_encuestada %>%
      filter(!is.na(!!sym(variable_dada[i])))
  }
  
  gente_encuestada <- gente_encuestada %>%
    group_by(clasificacion) %>%
    tally() %>%
    drop_na() %>%
    arrange(clasificacion)
  
  string_parametros <- ""
  
  for (i in 1:(length(nombre_dada) - 1)){
    string_parametros <- paste0(string_parametros, "Variable == '", nombre_dada[i], "' | ")
  }
  string_parametros <- paste0(string_parametros, 
                              "Variable == '", nombre_dada[length(nombre_dada)], "'")
  
  parametros <- aprioris %>%
    filter(eval(parse(text = string_parametros))) %>%
    mutate(`Proporción esperada (resultado)` = 
             as.numeric(`Proporción esperada (resultado)`))
  aux <- prod(parametros$`Proporción esperada (resultado)`)
  
  parametros <- rbind(parametros, c("Variables_dadas", aux, aux*(1 - aux)/2, NA, NA))
  parametros <- rbind(parametros, (aprioris %>% filter(Variable == nombre_condicionada)))
  parametros <- parametros[-c(1:(length(variable_dada))), ]
  
  mu <- as.numeric(parametros$`Proporción esperada (resultado)`)
  
  sigma2 <- as.numeric(parametros$`Error +- de la proporción esperada (calculo uno en automático pero relléname si esto no te gusta)`)
  
  out_Noinformativa <- dos_binomial(s = s, 
                                    pesos = pesos,
                                    ns = gente_encuestada, 
                                    mu = NULL, 
                                    sigma2 = NULL, 
                                    nsim = 5000, 
                                    nombres = centros$clasificacion)
  
  
  out_informativa <- dos_binomial(s = s, 
                                  pesos = pesos,
                                  ns = gente_encuestada, 
                                  mu = mu, 
                                  sigma2 = sigma2, 
                                  nsim = 500, 
                                  nombres = centros$clasificacion)
  
  write.csv(out_informativa, file = paste0(nombres_cond, "_InformativaP.csv"))
  write.csv(out_Noinformativa, file = paste0(nombres_cond, "_NoInformativaP.csv"))
  
}

# Ya están ajustados por sens y spec

# --- Pregunta IGG ---
generacion_distribucion_una_bin("resultado_cualit_sero", "IGG", "Positive")

# ---- Pregunta PCR+
generacion_distribucion_una_bin("resultado_sarscov2", "PCR+", "Positivo")

# ---- Pregunta Necesitó
generacion_distribucion_una_bin("necesidad_3meses", "Necesidad en salud 3 meses", "Si")

# ---------------------------- Una variable dado otra variable ----------


# --- Pregunta Necesidad | IGG <<
# Tomaremos necesidad_3meses
generacion_distribucion_dos_bin(vars_analisis = c("necesidad_3meses", "resultado_cualit_sero"),
                                igualdades_vars = c("Si", "Positive"),
                                nombres_apriori = c("Necesidad | IGG", "IGG"), 
                                ajuste_sen = T)
# --- Pregunta IGG | Necesidad


# --- Pregunta Busqueda | Necesidad <<
# Tomaremos necesidad_3meses
generacion_distribucion_dos_bin(vars_analisis = c("busqueda_salud", "necesidad_3meses"),
                                igualdades_vars = c("Si", "Si"),
                                nombres_apriori = c("Busqueda | Necesidad", "Necesidad en salud 3 meses"), 
                                ajuste_sen = F)


# --- Pregunta Sintomas | IGG <<
generacion_distribucion_dos_bin(vars_analisis = c("covid_sintoma", "resultado_cualit_sero"),
                                igualdades_vars = c("Si", "Positive"),
                                nombres_apriori = c("Sintomas | IGG", "IGG"),
                                ajuste_sen = T)

# CHECKME esta sí sale medio diferente especialmente en centro de distribución y en Bimbo
# ---- Pregunta Buscó atención | Sintomas << 
generacion_distribucion_dos_bin(vars_analisis = c("busqueda_atencion_covid", "covid_sintoma"),
                                igualdades_vars = c("Si", "Si"),
                                nombres_apriori = c("Buscó atención | Síntomas", "Sintomas"),
                                ajuste_sen = T)

# --- Pregunta Buscó | IGG y necesitó <<
generacion_distribucion_varias_bin(vars_analisis = c("busqueda_atencion_covid", "resultado_cualit_sero", "necesidad_3meses"),
                                   igualdades_vars = c("Si", "Positive", "Si"),
                                   nombres_apriori = c("Buscó | IGG Y Necesitó", "Necesidad | IGG", "IGG"))

# --- Pregunta Recibió atención | Búsqueda y necesidad <<
generacion_distribucion_varias_bin(vars_analisis = c("institucion_atencion", "busqueda_salud", "necesidad_3meses"),
                                   igualdades_vars = c("Si", "Si", "Si"),
                                   nombres_apriori = c("Recibió atención | Búsqueda y Necesidad", "Busqueda | Necesidad", "Necesidad en salud 3 meses"))


# --- Pregunta Recibió atención | Búsqueda y sintomas <<
generacion_distribucion_varias_bin(vars_analisis = c("atencion_busqueda", "busqueda_atencion_covid", "covid_sintoma"),
                                   igualdades_vars = c("Si", "Si", "Si"),
                                   nombres_apriori = c("Recibió atención | Buscó y síntomas", "Buscó atención | Síntomas", "Sintomas"))

# --- Pregunta Recibió atención IMSS | Búsqueda y sintomas <<
generacion_distribucion_varias_bin(vars_analisis = c("lugar_busqueda", "busqueda_atencion_covid", "covid_sintoma"),
                                   igualdades_vars = c("IMSS", "Si", "Si"),
                                   nombres_apriori = c("Recibió atención en IMSS | Buscó y síntomas", "Buscó atención | Síntomas", "Sintomas"))

# Pregunta Buscó atención | Sintomas e IGG+ <<
generacion_distribucion_varias_bin(vars_analisis = c("busqueda_atencion_covid", "covid_sintoma", "resultado_cualit_sero"),
                                   igualdades_vars = c("Si", "Si", "Positive"),
                                   nombres_apriori = c("Buscó atención | Síntomas e IGG", "Sintomas | IGG", "IGG"))



# Pregunta Recibió atención | Buscó, necesitó e IGG <<
generacion_distribucion_varias_bin(vars_analisis = c("atencion_busqueda", "busqueda_atencion_covid", "necesidad_3meses", "resultado_cualit_sero"),
                                igualdades_vars = c("Si", "Si", "Si", "Positive"),
                                nombres_apriori = c("Recibió | Buscó, Necesitó e IGG", "Buscó | IGG Y Necesitó", "Necesidad | IGG", "IGG"))


# Pregunta Recibió atención | Buscó, sintomas e IGG <<
generacion_distribucion_varias_bin(vars_analisis = c("atencion_busqueda", "busqueda_atencion_covid", "covid_sintoma", "resultado_cualit_sero"),
                                   igualdades_vars = c("Si", "Si", "Si", "Positive"),
                                   nombres_apriori = c("Recibió atención | Buscó y síntomas e IGG+", "Buscó atención | Síntomas e IGG", "Sintomas | IGG", "IGG"))


# Pregunta Recibió atención IMSS | Buscó, sintomas e IGG <<
generacion_distribucion_varias_bin(vars_analisis = c("lugar_busqueda", "busqueda_atencion_covid", "covid_sintoma", "resultado_cualit_sero"),
                                   igualdades_vars = c("IMSS", "Si", "Si", "Positive"),
                                   nombres_apriori = c("Recibió atención | Buscó y síntomas e IGG+", "Buscó atención | Síntomas e IGG", "Sintomas | IGG", "IGG"))


# Vamos a hacer P(IGG | Necesitó, Busco)
variable_condicionada <- vars_analisis[1]
igualdad_condicionada <- igualdades_vars[1]
nombre_condicionada   <- nombres_apriori[1]

variable_dada <- c()
igualdad_dada <- c()
nombre_dada <- c()
for (i in 2:length(vars_analisis)){            #Esto me hace que la primera entrada sea NA 
  variable_dada[i - 1]         <- vars_analisis[i]
  igualdad_dada[i - 1]         <- igualdades_vars[i]
  nombre_dada[i - 1]           <- nombres_apriori[i]
}

nombres_cond <- str_remove_all(nombre_condicionada, " ")
nombres_cond <- gsub("\\|", "_dado", nombres_cond)

string_aux <- ""

for ( i in 1:(length(variable_dada)-1)){
  string_aux <- paste0(string_aux, variable_dada[i], " == '", igualdad_dada[i], "' & ")
}
string_aux <- paste0(string_aux, 
                     variable_dada[length(variable_dada)], " == '", igualdad_dada[length(variable_dada)], "'")

var_dada <- ciclo1 %>%
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
  rename(`variables_dadas` = `n()`) %>%
  mutate(`variables_dadas` = replace_na(`variables_dadas`, 0))

string_aux <- paste0(string_aux, "& ", 
                     variable_condicionada, " == '", igualdad_condicionada, "'")

comb_vars <- ciclo1 %>%
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
  rename(nombre_condicionada = `n()`) %>%
  mutate(nombre_condicionada = replace_na(nombre_condicionada, 0))

#Unimos ambas tablas y arreglamos de acuerdo a la clasificación
s <- full_join(var_dada, comb_vars, by = "clasificacion")
s <- s %>% arrange(clasificacion) %>% select(-clasificacion) %>% as.matrix()

gente_encuestada <- ciclo1 

# Es la gente que responde la primera pregunta 
# Solo es saber si es NA o no
# FIXME
for (i in 1:length(variable_dada)){
  gente_encuestada <- gente_encuestada %>%
    filter(!is.na(!!sym(variable_dada[i])))
}

gente_encuestada <- gente_encuestada %>%
  group_by(clasificacion) %>%
  tally() %>%
  drop_na() %>%
  arrange(clasificacion)





## --- Para el caso multinomial ---
# Tengo el problema de tipo
# t ~ Binomial(n, q)
# q ~ Beta(alpha_q, beta_q)
# k ~ Multinomial(t, p1,..., pm)
# p ~ Dir(alpha)

# Para la beta
mu_q_t <- runif(n,min = 0, max= 1)
sigma2_q_t <- runif(n, min = 0, max = mu_q_t*(1-mu_q_t)) 

alpha_q_t <- -(mu_q_t*(sigma2_q_t + mu_q_t^2 - mu_q_t))/sigma2_q_t
beta_q_t <- (sigma2_q_t + mu_q_t^2 - mu_q_t)*(mu_q_t - 1)/sigma2_q_t

#Para la Dirichlet
k <- 7 #Numero de categorias
t <- 2000 #cantidad de positivos
alpha <- runif(k, 0, 1)
k_i <- c(1:7) #vector de los contados en cada categoria

# Entonces la distribucion de q es Beta(alpha_q + t, beta_q + n - t) 
# y p es Dirichlet(alpha + k )
distr_mult_betaq <- list()
distr_mult_dir   <- list()

for (i in 1:n){
  distr_mult_betaq[[i]]<- rbeta(n, alpha_q_t[i] + t, beta_q_t[i] + n - t)
  distr_mult_dir[[i]] <- rdirichlet(n, alpha + k_i) 
}

for (i in 1:n){
  mean_mult_betaq <- mean(mean(distr_mult_betaq[[i]]))
  mean_mult_dir <- mean(mean(distr_mult_dir[[i]]))#tienen que salir las proporciones de k_i/t
}

quantiles_mult_beta_q <- sapply(distr_mult_betaq, quantile, prob=c(0.025, .975), na.rm=TRUE)
quantiles_mult_dir <- sapply(distr_mult_dir, quantile, prob=c(0.025, .975), na.rm=TRUE)
