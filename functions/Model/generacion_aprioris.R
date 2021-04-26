# Script para generar las probabilidad a priori que 
# necesitamos para que el modelo funcione
# Autores: Rodrigo Zepeda, Valeria Perez

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
source("generacion_informacion.R", encoding = "UTF-8")

setwd("~/IMSS/EPCOVID/Cuestionario")


# Lo lee con los headers mal
aprioris <- read_excel("~/IMSS/EPCOVID/Cuestionario/Aprioris David.xlsx")
#Entonces lo arreglamos
colnames(aprioris) <- aprioris[1, ]
aprioris <- aprioris[-1, ]

# Necesitamos *IGG | Necesidad* pero para eso, necesitamos (IGG, Necesidad)

# IGG, Necesidad = (Necesidad | IGG)(IGG)
aprioris <- mu_y_sigma(aprioris, 
                       vars_prod = c("Necesidad | IGG", "IGG"), vars_division = NULL, 
                       "Necesidad, IGG")
# IGG | Necesidad = Necesidad, IGG / Necesidad
aprioris <- mu_y_sigma(aprioris, 
                       vars_prod =c("Necesidad, IGG"), 
                       vars_division = "Necesidad en salud 3 meses", 
                       nombre_nueva = "IGG | Necesidad")

# De seguro necesitamos Necesidad | Busqueda 
# Necesidad, Busqueda = (Busqueda | Necesidad) * (Necesidad)
aprioris <- mu_y_sigma(aprioris, 
                       vars_prod = c("Busqueda | Necesidad", "Necesidad en salud 3 meses"), 
                       vars_division = NULL, 
                       "Necesidad, Busqueda")

# Necesidad | Busqueda = Necesidad, Busqueda / Busqueda
aprioris <- mu_y_sigma(aprioris, 
                       vars_prod =c("Necesidad, Busqueda"), 
                       vars_division = "Busqueda atención", 
                       nombre_nueva = "Necesidad | Busqueda")

# IGG, Busco, Necesito = (Busco | IGG, Necesito)(IGG, Necesito)
aprioris <- mu_y_sigma(aprioris, 
                       vars_prod =c("Buscó | IGG Y Necesitó", "Necesidad, IGG"), 
                       vars_division = NULL, 
                       nombre_nueva = "IGG, Busco, Necesito")

# IGG | Busco y Necesitó = (IGG, Busco, Necesitó) /(Buscó, Necesitó)
aprioris <- mu_y_sigma(aprioris, 
                       vars_prod = "IGG, Busco, Necesito",
                       vars_division = "Necesidad, Busqueda",
                       nombre_nueva = "IGG | Busco y Necesitó")

# Necesitó, Busco, Recibio = (Recibio | Busqueda y Necesidad) (Busqueda y Necesidad)
aprioris <- mu_y_sigma(aprioris, 
                       vars_prod = c("Recibió atención | Búsqueda y Necesidad", "Necesidad, Busqueda"),
                       vars_division = NULL, 
                       nombre_nueva = "Necesito, Busco, Recibio")

# IGG, Necesitó, Buscó, Recibió = (Recibio | Busco, Necesitó, IGG)(IGG, Busco, Necesito)
aprioris <- mu_y_sigma(aprioris, 
                       vars_prod = c("Recibió | Buscó, Necesitó e IGG", "IGG, Busco, Necesito"),
                       vars_division = NULL, 
                       nombre_nueva = "IGG, Necesitó, Buscó, Recibió")

# IGG | Necesitó, Buscó, Recibió = (IGG, Necesitó, Busco, Recibió) / (Necesitó, Buscó, Recibio)
aprioris <- mu_y_sigma(aprioris, 
                       vars_prod = c("IGG, Necesitó, Buscó, Recibió"),
                       vars_division = c("Necesito, Busco, Recibio"), 
                       nombre_nueva = "IGG | Necesitó, Buscó, Recibió")

# Sintomas, IGG = (Sintomas | IGG) (IGG)
aprioris <- mu_y_sigma(aprioris, 
                       vars_prod = c("Sintomas | IGG", "IGG"),
                       vars_division = NULL, 
                       nombre_nueva = "Sintomas, IGG")

# IGG | Sintomas = (Sintomas, IGG) / Sintomas
aprioris <- mu_y_sigma(aprioris, 
                       vars_prod = "Sintomas, IGG", 
                       vars_division = "Sintomas", 
                       nombre_nueva = "IGG | Sintomas")

# IGG, Busco atencion, Sintomas = (Busco | Sintomas e IGG) * (Sintomas, IGG)
aprioris <- mu_y_sigma(aprioris, 
                       vars_prod = c("Buscó atención | Síntomas e IGG", "Sintomas, IGG"), 
                       vars_division = NULL,
                       nombre_nueva = "IGG, Busco atencion, Sintomas")

# Busco, Sintomas = (Busco atención | Sintomas) * (Sintomas)
aprioris <- mu_y_sigma(aprioris, 
                       vars_prod = c("Buscó atención | Síntomas", "Sintomas"), 
                       vars_division = NULL,
                       nombre_nueva = "Busco, Sintomas")

# IGG | Busco atencion, Sintomas = (IGG, Busco, Sintomas) / (Busco, Sintomas)
aprioris <- mu_y_sigma(aprioris, 
                       vars_prod = "IGG, Busco atencion, Sintomas", 
                       vars_division = "Busco, Sintomas", 
                       nombre_nueva = "IGG | Busco atencion, Sintomas")

# IGG, Sintomas, Busco, Recibio = (Recibio | Busco, Sintomas, IGG) * (Busco, Sintomas, IGG)
aprioris <- mu_y_sigma(aprioris, 
                       vars_prod = c("Recibió atención | Buscó y síntomas e IGG+", "IGG, Busco atencion, Sintomas"), 
                       vars_division = NULL,
                       nombre_nueva = "IGG, Sintomas, Busco, Recibio")

# Sintomas, Busco, Recibio = (Recibio | Busco y sintomas) * (Busco, Sintomas)
aprioris <- mu_y_sigma(aprioris, 
                       vars_prod = c("Recibió atención | Buscó y síntomas", "Busco, Sintomas"), 
                       vars_division = NULL, 
                       nombre_nueva = "Sintomas, Busco, Recibio")

# IGG | Sintomas, Busco, Recibio = (IGG, Sintomas, Busco, Recibio) / (Sintomas, Busco, Recibio)
aprioris <- mu_y_sigma(aprioris, 
                       vars_prod = "IGG, Sintomas, Busco, Recibio",  
                       vars_division = "Sintomas, Busco, Recibio", 
                       nombre_nueva = "IGG | Sintomas, Busco, Recibio")

# IGG, Recibio IMSS, Busco, Sintomas = (Recibio IMSS | Busco, Sintomas, IGG) * (Busco, Sintomas, IGG)
aprioris <- mu_y_sigma(aprioris, 
                       vars_prod = c("Recibió atención en IMSS | Buscó y síntomas e IGG+", "IGG, Busco atencion, Sintomas"), 
                       vars_division = NULL,
                       nombre_nueva = "IGG, Recibio IMSS, Busco, Sintomas")

# Recibio IMSS, Busco, Sintomas = (Recibió atención en IMSS | Buscó y síntomas) * (Busco, Sintomas)
aprioris <- mu_y_sigma(aprioris, 
                       vars_prod = c("Recibió atención en IMSS | Buscó y síntomas", "Busco, Sintomas"), 
                       vars_division = NULL,
                       nombre_nueva = "Recibio IMSS, Busco, Sintomas")

# IGG | Recibió IMSS, Busco, Sintomas = (IGG, Recibio IMSS, Busco, Sintomas) / (Recibio IMSS, Busco, Sintomas)
aprioris <- mu_y_sigma(aprioris, 
                       vars_prod = "IGG, Recibio IMSS, Busco, Sintomas", 
                       vars_division = "Recibio IMSS, Busco, Sintomas", 
                       nombre_nueva = "IGG | Recibió IMSS, Busco, Sintomas")


write.csv(aprioris, "aprioris_conAgregados.csv")
