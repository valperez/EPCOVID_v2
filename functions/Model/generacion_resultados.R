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


setwd("~/IMSS/EPCOVID/Model/Model")

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

source("generacion_distribuciones.R", encoding = "UTF-8")
source("generacion_informacion.R", encoding = "UTF-8")
source("funciones_finales.R", encoding = "UTF-8")

# # # # ---- MÉTODOS INFORMATIVOS ---- # # # #

# Vamos a hacer las distribuciones con a priori's informativas

#Original
ciclo1 <- read_delim("~/IMSS/EPCOVID/Cuestionario/EPCOVID_c1_20210225.txt", ";", 
                     escape_double = FALSE, locale = locale(encoding = "UTF-8"), 
                     trim_ws = TRUE)
ciclo2 <- read_delim("~/IMSS/EPCOVID/Cuestionario/EPCOVID_c2_20210225.txt", ";", 
                     escape_double = FALSE, locale = locale(encoding = "UTF-8"), 
                     trim_ws = TRUE)
ciclo3 <- read_delim("~/IMSS/EPCOVID/Cuestionario/EPCOVID_c3_20210225.txt", ";", 
                     escape_double = FALSE, locale = locale(encoding = "UTF-8"), 
                     trim_ws = TRUE)


# Lo lee con los headers mal
aprioris <- read_excel("~/IMSS/EPCOVID/Cuestionario/aprioris_conAgregados.xlsx")
#Entonces lo arreglamos
aprioris <- aprioris[, -1 ]

# Ya tengo los centros
centros <- ciclo1 %>%
  select(clasificacion) %>%
  distinct() %>%
  drop_na() %>%
  arrange(clasificacion)

tamanios_totales <- read_csv("~/IMSS/EPCOVID/Cuestionario/Tamaños totales de muestra 29 enero 2020.csv") %>%
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

setwd("~/IMSS/EPCOVID/Model/Model")

# IGG
enes <- sizes(ciclo1, 
                vector_na = c("resultado_cualit_sero"), 
                string_igualdades = "resultado_cualit_sero == 'Positive' ")
parametros <- busqueda_mu_y_sigma(df_aprioris = aprioris,
                                  variables = "IGG")


IGG <- calculo_binomial_0(nsize = enes[, 1], zsize = enes[, 2], 
                   mu = parametros$`Proporción esperada (resultado)`, 
                   sigma = parametros$`Error +- de la proporción esperada (calculo uno en automático pero relléname si esto no te gusta)`, 
                   delta = 0.8, gamma = 0.95, 
                               compilar_julia = T,
                               burnin = 100000, niter = 200000, nchains = 4, 
                               starting_point = c(1, 1, 0, 0, 0), 
                               proba_nom = "IGG")

# IGG | Necesidad
enes <- sizes(ciclo1, 
              vector_na = c("resultado_cualit_sero", "necesidad_3meses"), 
              string_igualdades = "resultado_cualit_sero == 'Positive' & necesidad_3meses == 'Si'")
z <- sizes(ciclo1, 
           vector_na = c("necesidad_3meses"), 
           string_igualdades = "necesidad_3meses == 'Si'")

parametros <- busqueda_mu_y_sigma(df_aprioris = aprioris,
                                  variables = c("IGG", "IGG | Necesidad"))

IGG_dadoNecesidad <- calculo_binomial_1(nsize = enes[,1], zsize = z[, 2], tsize = enes[, 2], 
                                        mu_p = parametros[1, 1], sigma_p = parametros[1, 2],  
                                        mu_q = parametros[2, 1], sigma_q = parametros[2, 2], 
                                        delta = 0.8, gamma = 0.95, compilar_julia = T,
                                        burnin = 100000, niter = 200000, nchains = 4, 
                                        starting_point = c(1, 1, 1, 1, 0, 0, 0, 0, 0, 0), 
                                        proba_nom = "IGGdadoNecesidad")
