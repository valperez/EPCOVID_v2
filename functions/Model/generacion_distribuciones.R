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

# Supongamos que tenemos el problema 
# s ~ Binomial (n, p)
# p    ~ Beta(alpha_p, beta_p)
# t    ~ Binomal(s,q)
# q    ~ Beta(alpha_q, beta_q)


# Ahora haré una función que me genere las simulaciones de la distribución beta
# mu y sigma2 es lo que nos dió David (vectores ambos),
# s es la cantidad de gente que respondió a la pregunta 1 y pregunta 2 por centro de trabajo (matriz)
# osea en la columna_i de s van a venir las personas que respondieron a la pregunta_i
# pesos es la ponderación por centro de trabajo (vector)
# n es la cantidad de gente encuestada por centro de trabajo (vector) que respondió logicamente
# Es importante que la primera variable que le damos sea la t tanto en vars_analisis
# como en igualdades_variables

# OJO : la sigma2 en el excel de David está divida entre 2 porque si no, no funciona
generacion_unabeta <- function(s, pesos, ns, 
                               mu = NULL, sigma2 = NULL, 
                               nsim = 5000){
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
  return(distr_beta)
}

#función que genera una lista de listas
# La primera lista es la informativa mientras que la segunda es la no informativa 
# los paremetros son los mismos que para la función anterior

generacion_p <- function(s, pesos, ns, 
                             mu, sigma2, 
                             nsim = 5000){
  
  #Verificacion de longitudes de las variables
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
  beta_informativa <- list()
  
  for (i in 1:nrow(ns)){
    distr_beta[[i]] <- rbeta(nsim, 
                             alpha + s[i, 1], 
                             beta + ns[i, 1] - s[i, 1])
  }
  beta_informativa <- distr_beta
  
  mu <- NULL
  sigma2 <- NULL
  
   if (is.null(mu)){
    mu <- rep(0.5, ncol(s)) 
  }
  
  if (is.null(sigma2)){
    sigma2 <- rep(1/12, ncol(s)) 
  }
  
  alpha <- -(mu*(sigma2 + mu^2 - mu))/sigma2
  beta <- (sigma2 + mu^2 - mu)*(mu - 1)/sigma2
  
  distr_beta <- list()
  beta_noInformativa <- list()
  
  for (i in 1:nrow(ns)){
    distr_beta[[i]] <- rbeta(nsim, 
                             alpha + s[i, 1], 
                             beta + ns[i, 1] - s[i, 1])
  }
  beta_noInformativa <- distr_beta
  
  
  return(list(beta_informativa = beta_informativa, 
              beta_noInformativa = beta_noInformativa))
}

#Función que genera la distribucion posteriori de q. Los parámetros son los mismos
# que las funciones anteriores
#

generacion_q <- function(s, pesos, ns, 
                         mu, sigma2, 
                         nsim = 5000){
  
  #Verificacion de longitudes
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
  beta_informativa <- list()

  for (i in 1:nrow(ns)){
    distr_beta[[i]] <- rbeta(nsim, 
                             alpha + s[i, 2], 
                             beta + s[i, 1] - s[i, 2])
  }
  beta_informativa <- distr_beta
  
  mu <- NULL
  sigma2 <- NULL
  
  if (is.null(mu)){
    mu <- rep(0.5, ncol(s)) 
  }
  
  if (is.null(sigma2)){
    sigma2 <- rep(1/12, ncol(s)) 
  }
  
  alpha <- -(mu*(sigma2 + mu^2 - mu))/sigma2
  beta <- (sigma2 + mu^2 - mu)*(mu - 1)/sigma2
  
  distr_beta <- list()
  beta_noInformativa <- list()
  
  for (i in 1:nrow(ns)){
    distr_beta[[i]] <- rbeta(nsim, 
                             alpha + s[i, 2], 
                             beta + s[i, 1] - s[i, 2])
  }
  beta_noInformativa <- distr_beta

  return(list(beta_informativa = beta_informativa, 
              beta_noInformativa = beta_noInformativa))
}

# Función que me ajusta por sensitividad y especificidad
# lista_betas es una lista de betas que tengo generadas por 
# la función generacion_betas
# indice es al elemento de la lista que le quiero aplicar el ajuste
# sens y spec son los valores que le daré a sensitividad y especificidad
ajusta_sens_spec <- function(lista_betas, indice, sens = 1, spec = 1){
  
  if (sens != 1 | spec != 1){
    for (i in 1:length(lista_betas[[indice]])){
      lista_betas[[indice]][[i]] <- lista_betas[[indice]][[i]]*sens + 
        (1 - spec)*(1 - lista_betas[[indice]][[i]])
    }
  }
  return (lista_betas)
}


