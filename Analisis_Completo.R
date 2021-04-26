rm(list = ls())

library(tidyverse)
library(readxl)
library(tidybayes)
library(JuliaCall)
library(rlang)
library(scales)
library(purrr)

setwd("/home")

#Número de ciclos de la encuesta
nciclos <- 3

#Sourcing functions
source("functions/julia_setup_calls.R")
source("functions/read_panel_files.R")
source("functions/raux.R")
source("functions/Model/funciones_finales.R")

mantener.convergidos <- T #SI lo pones a FALSE recalcula todo;
                             #en true sólo calcula los q faltan de converger

#https://www.cochrane.org/CD013652/INFECTN_what-diagnostic-accuracy-antibody-tests-detection-infection-covid-19-virus
gamma_igg <- 0.907 #especificidad
delta_igg <- 0.98 #sensibilidad

#PCR
gamma_pcr <- 1.00 #especificidad
delta_pcr <- 0.98 #sensibilidad


#Leemos todos los archivos que están en datasets y los procesamos
panel_datasets <- read_panel_files()$encuesta_panel

#Obtenemos los estratos de la muestra
stratum <- panel_datasets %>% select(strata) %>% 
  arrange(strata) %>% distinct()
ciclos  <- unique(panel_datasets$Ciclo)

pondef  <- panel_datasets %>% select(Ciclo, sampling_weights, strata) %>% 
  distinct() %>%
  arrange(strata)
#hats aqui siempre

#----------------------------------------------------------------------------
#PROPORCIONES 
#----------------------------------------------------------------------------

#Proporción de PCR+ en la muestra
#----------------------------------------------------------------------------
#CONVERGENCE = TRUE
if(!mantener.convergidos){
  message("Calculando PCR+")
  julia_eval("Random.seed!(237);")
  genera_binomial_0("!is.na(resultado_sarscov2) & resultado_sarscov2 != 'Inconcluso'",
                    "str_detect(resultado_sarscov2,'Positivo')", 
                    delta = delta_pcr, gamma = gamma_pcr,
                    burnin = 200000, niter = 400000,
                    varname = "Proporción de PCR+", filedir = "resultados/")
}

#Proporción de IGG+ en la muestra
#----------------------------------------------------------------------------
if(!mantener.convergidos){
  message("Calculando IGG+")
  julia_eval("Random.seed!(23854);")
  genera_binomial_0("!is.na(resultado_cualit_sero)",
                    "resultado_cualit_sero == 'Positive'", 
                    delta = delta_igg, gamma = gamma_igg,
                    burnin = 40000, niter = 80000,
                    varname = "Proporción de IGG+", filedir = "resultados/")
}
#Proporción de Necesidad 3 meses en la muestra
#----------------------------------------------------------------------------
if(!mantener.convergidos){
  message("Calculando Necesidad")
  julia_eval("Random.seed!(23854);")
  genera_binomial_0("!is.na(NECESIDAD_3MESES)",
                    "NECESIDAD_3MESES == 'Si'", 
                    burnin = 40000, niter = 80000,
                    varname = "Proporción que tuvo necesidad en los últimos 3 meses", filedir = "resultados/")
}
#Proporción de ALCOHOL CONSUMO 
#----------------------------------------------------------------------------
if(!mantener.convergidos){
  message("Calculando Alcohol")
  julia_eval("Random.seed!(93);")
  genera_binomial_0("!is.na(ALCOHOL_CONSUMO) & ALCOHOL_CONSUMO != 'No sabe' & ALCOHOL_CONSUMO != 'No deseo responder' & ALCOHOL_CONSUMO != '2'",
                    "ALCOHOL_CONSUMO == 'Si'", 
                    burnin = 10000, niter = 50000,
                    varname = "Proporción que consumió alcohol durante el confinamiento", filedir = "resultados/")
}
#Proporción de COVID_SINTOMA
#----------------------------------------------------------------------------
if(!mantener.convergidos){
  message("Calculando Covid síntoma")
  julia_eval("Random.seed!(274247);")
  genera_binomial_0("!is.na(COVID_SINTOMA) & COVID_SINTOMA != 'No sabe' & COVID_SINTOMA != 'No deseo responder'",
                    "COVID_SINTOMA == 'Si'", 
                    burnin = 100000, niter = 200000,
                    varname = "Proporción que ha tenido enfermedad respiratoria", filedir = "resultados/")
}

#Proporción de COVID_SINTOMA
#----------------------------------------------------------------------------
if(!mantener.convergidos){
  message("Calculando Covid síntoma")
  julia_eval("Random.seed!(274247);")
  genera_binomial_0("!is.na(COVID_SINTOMA) & COVID_SINTOMA != 'No sabe' & COVID_SINTOMA != 'No deseo responder'",
                    "COVID_SINTOMA == 'Si'", 
                    burnin = 100000, niter = 200000,
                    varname = "Proporción que ha tenido enfermedad respiratoria", filedir = "resultados/")
}

#Proporción de COVID_SINTOMAM CONFIRMA
#----------------------------------------------------------------------------
if(!mantener.convergidos){
  nfilter <- "!is.na(COVID_SINTOMA_CONFIRMACION)"
  sfilter <- "str_detect(COVID_SINTOMA_CONFIRMACION,'Dolor|Pérdida|nasal|cabeza|Ojos|garganta|articulaciones|Tos|Diarrea|Falta|Vómito')"
  message("Calculando Covid síntoma")
  julia_eval("Random.seed!(18111);")
  genera_binomial_0(nfilter,
                    sfilter, 
                    burnin = 10000, niter = 20000,
                    varname = "Proporción que tuvo síntoma confirmatorio de COVID-19", filedir = "resultados/")
}

#Proporción de SINTOMA_DISNEA
#----------------------------------------------------------------------------
if(!mantener.convergidos){
  message("Calculando disnea síntoma")
  julia_eval("Random.seed!(274247);")
  genera_binomial_0("!is.na(SINTOMA_DISNEA) & SINTOMA_DISNEA != 'No sabe' & SINTOMA_DISNEA != 'No deseo responder'",
                    "SINTOMA_DISNEA == 'Si'", 
                    burnin = 100000, niter = 200000,
                    varname = "Proporción que presentó disnea (dificultad respiratoria)", filedir = "resultados/")
}

#Proporción de itt 
#----------------------------------------------------------------------------
if(!mantener.convergidos){
  message("Calculando itt")
  julia_eval("Random.seed!(26);")
  genera_binomial_0("TRUE",
                    "!is.na(fh_inicio_itt)", 
                    burnin = 200000, niter = 400000,
                    varname = "Proporción que presentó al menos una incapacidad temporal en el trabajo", 
                    filedir = "resultados/")
}

#Proporción de hombres en la muestra
#----------------------------------------------------------------------------
if(!mantener.convergidos){
  message("Calculando % hombres")
  julia_eval("Random.seed!(1999);")
  genera_binomial_0("!is.na(SEXO)","SEXO == 'Hombre'", 
                    burnin = 100000, niter = 200000,
                    varname = "Proporción de hombres", filedir = "resultados/")
}

#----------------------------------------------------------------------------
#PROPORCIONES CON UNA CONDICIÓN
#----------------------------------------------------------------------------

#Hombre e IGG
#----------------------------------------------------------------------------
if(!mantener.convergidos){
  message("Calculando hombres e igg")
  julia_eval("Random.seed!(26);")
  genera_binomial_1("!is.na(SEXO) & !is.na(resultado_cualit_sero)", 
                    "SEXO == 'Hombre'", "resultado_cualit_sero == 'Positive'",
                    pname = c("Proporción de hombres", "Proporción de mujeres"),
                    qname = c("IGG+ en el grupo de hombres", "IGG+ en el grupo de mujeres"),
                    delta = delta_igg, gamma = gamma_igg,
                    burnin = 10000, niter = 20000,
                    filedir = "resultados/", varname = "IGG+ por sexo hombre")
}

#Mujer e IGG
#----------------------------------------------------------------------------
if(!mantener.convergidos){
  message("Calculando mujeres e igg")
  julia_eval("Random.seed!(26);")
  genera_binomial_1("!is.na(SEXO) & !is.na(resultado_cualit_sero)", 
                    "SEXO != 'Hombre'", 
                    "resultado_cualit_sero == 'Positive'",
                    pname = c("Proporción de mujeres"),
                    qname = c("IGG+ en el grupo de mujeres"),
                    delta = delta_igg, gamma = gamma_igg,
                    burnin = 10000, niter = 20000,
                    filedir = "resultados/", varname = "IGG+ por sexo mujer")
}

#Sexo por PCR
#----------------------------------------------------------------------------
if(!mantener.convergidos){
  message("Calculando hombres pcr")
  julia_eval("Random.seed!(26);") #71717
  genera_binomial_1("!is.na(SEXO) & !is.na(resultado_sarscov2) & resultado_sarscov2 != 'Inconcluso'", 
                    "SEXO == 'Hombre'", "str_detect(resultado_sarscov2,'Positivo')",
                    pname = c("Proporción de hombres", "Proporción de mujeres"),
                    qname = c("PCR+ en el grupo de hombres", "PCR+ en el grupo de mujeres"),
                    delta = delta_pcr, gamma = gamma_pcr,
                    burnin = 10000, niter = 20000,
                    filedir = "resultados/", varname = "PCR+ por sexo hombre")
}

if(!mantener.convergidos){
  message("Calculando pcr mujer")
  julia_eval("Random.seed!(26);") #71717
  genera_binomial_1("!is.na(SEXO) & !is.na(resultado_sarscov2) & resultado_sarscov2 != 'Inconcluso'", 
                    "SEXO != 'Hombre'", "str_detect(resultado_sarscov2,'Positivo')",
                    pname = c("Proporción de mujeres", "Proporción de homres"),
                    qname = c("PCR+ en el grupo de mujeres", "PCR+ en el grupo de hombres"),
                    delta = delta_pcr, gamma = gamma_pcr,
                    burnin = 10000, niter = 20000,
                    filedir = "resultados/", varname = "PCR+ por sexo mujer")
}

#NECESITO e IGG+
#----------------------------------------------------------------------------
if(!mantener.convergidos){
  message("Calculando necesitó e igg")
  julia_eval("Random.seed!(2747);")
  genera_binomial_1("!is.na(NECESIDAD_3MESES) & !is.na(resultado_cualit_sero)", 
                    "NECESIDAD_3MESES == 'Si'", "resultado_cualit_sero == 'Positive'",
                    pname = c("Proporción con necesidad en salud en últimos 3 meses", 
                              "Proporción sin necesidad en salud en últimos 3 meses"),
                    qname = c("IGG+ en el grupo que tuvo una necesidad en salud", 
                              "IGG+ en el grupo que no tuvo una necesidad en salud"),
                    delta = delta_igg, gamma = gamma_igg,
                    burnin = 10000, niter = 20000,
                    filedir = "resultados/", varname = "IGG+ por NECESIDAD_3MESES 1")
}

#No NECESITO | IGG+
#----------------------------------------------------------------------------
if(!mantener.convergidos){
  message("Calculando no necesitó e igg")
  julia_eval("Random.seed!(82);")
  genera_binomial_1("!is.na(NECESIDAD_3MESES) & !is.na(resultado_cualit_sero)", 
                    "NECESIDAD_3MESES == 'No'", "resultado_cualit_sero == 'Positive'",
                    pname = c("Proporción sin necesidad en salud en últimos 3 meses", 
                              "Proporción con necesidad en salud en últimos 3 meses"),
                    qname = c("IGG+ en el grupo que no tuvo una necesidad en salud", 
                              "IGG+ en el grupo que tuvo una necesidad en salud"),
                    delta = delta_igg, gamma = gamma_igg,
                    burnin = 10000, niter = 20000,
                    filedir = "resultados/", varname = "IGG+ por NECESIDAD_3MESES 2")
}
  
#NECESITO SEXO
#----------------------------------------------------------------------------
if(!mantener.convergidos){
  message("Calculando hombres que necesitaron")
  julia_eval("Random.seed!(82);")
  genera_binomial_1("!is.na(SEXO) & !is.na(NECESIDAD_3MESES)", 
                    "SEXO == 'Hombre'", "NECESIDAD_3MESES == 'Si'",
                    pname = c("Proporción de hombres", "Proporción de mujeres"),
                    qname = c("Tuvo necesidad en salud en ultimos 3 meses en el grupo de hombres", "Tuvo necesidad en salud en ultimos 3 meses en el grupo de mujeres"),
                    filedir = "resultados/", 
                    varname = "Tuvo necesidad en salud en ultimos 3 meses por sexo hombre")
}

if(!mantener.convergidos){
  message("Calculando mujeres que necesitaron")
  julia_eval("Random.seed!(82);")
  genera_binomial_1("!is.na(SEXO) & !is.na(NECESIDAD_3MESES)", 
                    "SEXO != 'Hombre'", "NECESIDAD_3MESES == 'Si'",
                    pname = c("Proporción de mujeres", "Proporción de mujeres"),
                    qname = c("Tuvo necesidad en salud en ultimos 3 meses en el grupo de mujeres", "Tuvo necesidad en salud en ultimos 3 meses en el grupo de mujeres"),
                    filedir = "resultados/", 
                    burnin = 10000, niter = 20000,
                    varname = "Tuvo necesidad en salud en ultimos 3 meses por sexo mujer")
}

#NECESITO POR ITT
#----------------------------------------------------------------------------
if(!mantener.convergidos){
  julia_eval("Random.seed!(82);")
  message("Calculando itt que necesitaron")
  genera_binomial_1("!is.na(NECESIDAD_3MESES)", 
                    "is.na(fh_inicio_itt)", 
                    "NECESIDAD_3MESES == 'Si'",
                    pname = c("Proporción sin ITT", "Proporción con ITT"),
                    qname = c("Tuvo necesidad en salud en ultimos 3 meses pero no registró ITT", 
                              "Tuvo necesidad en salud en ultimos 3 meses y registró ITT"),
                    filedir = "resultados/", 
                    varname = "Tuvo necesidad en salud en ultimos 3 meses e ITT")
}

#BUSCO | NECESITO
#----------------------------------------------------------------------------
if(!mantener.convergidos){
  message("Calculando busco necesitaro")
  julia_eval("Random.seed!(82);")
  genera_binomial_1("(NECESIDAD_3MESES == 'Si' & !is.na(BUSQUEDA_SALUD)) | (NECESIDAD_3MESES == 'No' & is.na(BUSQUEDA_SALUD))", 
                    "NECESIDAD_3MESES == 'Si'", 
                    "BUSQUEDA_SALUD == 'Si'",
                    burnin = 10000, niter = 20000,
                    pname = c("Proporción con necesidad en salud en últimos 3 meses"),
                    qname = c("Buscó atender la necesidad en salud que tuvo los últimos 3 meses"),
                    filedir = "resultados/", 
                    varname = "Busco dado que tuvo necesidad")
}

#ER | IGG #FIXME
# Así como está, a mí me marcó bien todo atte: Val :)
#----------------------------------------------------------------------------
if(!mantener.convergidos){
  message("Calculando Covid síntoma igg+")
  julia_eval("Random.seed!(8181);")
  genera_binomial_1("!is.na(resultado_cualit_sero) & !is.na(COVID_SINTOMA) & COVID_SINTOMA != 'No sabe' & COVID_SINTOMA != 'No deseo responder'",
                    "COVID_SINTOMA == 'Si'", 
                    "resultado_cualit_sero == 'Positive'",
                    burnin = 10000, niter = 20000,
                    pname = c("Tuvo una enfermedad respiratoria"),
                    qname = c("IGG+ dado que tuvo una enfermedad respiratoria"),
                    delta = delta_igg, gamma = gamma_igg,
                    varname = "Proporción que tuvo enfermedad respiratoria e igg+", filedir = "resultados/")
}

#SINTOMA | IGG  #FIXME
#----------------------------------------------------------------------------
if(!mantener.convergidos){
  message("Calculando Covid síntoma igg+")
  julia_eval("Random.seed!(8181);")
  sfilter  <- "str_detect(COVID_SINTOMA_CONFIRMACION,'Dolor|Pérdida|nasal|cabeza|Ojos|garganta|articulaciones|Tos|Diarrea|Falta|Vómito')"
  nfilter  <- paste0("(str_detect(COVID_SINTOMA_CONFIRMACION,'Dolor|Pérdida|nasal|cabeza|Ojos|garganta|articulaciones|Tos|Diarrea|Falta|Vómito') & !is.na(BUSQUEDA_ATENCION_COVID)) | str_detect(COVID_SINTOMA_CONFIRMACION,'NINGUNO')")
  genera_binomial_1(nfilter,
                    sfilter, 
                    "resultado_cualit_sero == 'Positive'",
                    burnin = 10000, niter = 20000,
                    pname = c("Tuvo una enfermedad respiratoria"),
                    qname = c("IGG+ dado que tuvo una enfermedad respiratoria"),
                    delta = delta_igg, gamma = gamma_igg,
                    varname = "Proporción que tuvo enfermedad respiratoria e igg+", filedir = "resultados/")
}

#COVID SINTOMAS | busco #FIXME
#----------------------------------------------------------------------------
if(!mantener.convergidos){
  message("Calculando Covid BUSCO sintoma")
  julia_eval("Random.seed!(8181);")
  s1filter <- "str_detect(COVID_SINTOMA_CONFIRMACION,'Dolor|Pérdida|nasal|cabeza|Ojos|garganta|articulaciones|Tos|Diarrea|Falta|Vómito')"
  nfilter  <- paste0("(str_detect(COVID_SINTOMA_CONFIRMACION,'Dolor|Pérdida|nasal|cabeza|Ojos|garganta|articulaciones|Tos|Diarrea|Falta|Vómito') & !is.na(BUSQUEDA_ATENCION_COVID)) | str_detect(COVID_SINTOMA_CONFIRMACION,'NINGUNO')")
  genera_binomial_1(nfilter,
                    s1filter, 
                    "BUSQUEDA_ATENCION_COVID == 'Si'",
                    burnin = 10000, niter = 20000,
                    delta = delta_igg, gamma = gamma_igg,
                    pname = c("Tuvo sintoma confirmacion covid"),
                    qname = c("Busco atencion dado que tuvo sintoma confirmatorio"),
                    varname = "Proporción que tuvo sintoma y busco atencion", filedir = "resultados/")
}


#----------------------------------------------------------------------------
#PROPORCIONES CON DOS CONDICIONES
#----------------------------------------------------------------------------

#Buscó dado que necesitó e IGG+
# a pesar de los cambios en la sigma parece que todo bien :)
#----------------------------------------------------------------------------
nfilter <- "!is.na(resultado_cualit_sero) & ((NECESIDAD_3MESES == 'Si' & !is.na(BUSQUEDA_SALUD)) | (NECESIDAD_3MESES == 'No' & is.na(BUSQUEDA_SALUD)))"
if(!mantener.convergidos){
  message("Busco dado necesito e igg")
  julia_eval("Random.seed!(3771);")
  genera_binomial_2(nfilter,
                    "NECESIDAD_3MESES == 'Si'",
                    "BUSQUEDA_SALUD == 'Si'",
                    "resultado_cualit_sero == 'Positive'",
                    p1name = c("Tuvo necesidad"),
                    burnin = 10000, niter = 20000,
                    seeds = c(7373, 3472, 22),
                    p2name = c("Buscó dado que tuvo necesidad"),
                    p3name = c("IGG+ en grupo que buscó y necesitó"),
                    delta = delta_igg, gamma = gamma_igg,
                    filedir = "resultados/", varname = "IGG+ con busqueda y necesidad")
}

if(!mantener.convergidos){
  message("No Busco aunque necesito e igg")
  julia_eval("Random.seed!(3771);")
  genera_binomial_2(nfilter,
                    "NECESIDAD_3MESES == 'Si'",
                    "BUSQUEDA_SALUD == 'No'",
                    "resultado_cualit_sero == 'Positive'",
                    p1name = c("Tuvo necesidad"),
                    p2name = c("No Buscó aunque tuvo necesidad"),
                    p3name = c("IGG+ en grupo que no buscó pero necesitó"),
                    burnin = c(10000,10000,10000), niter = c(20000,20000,50000),
                    seeds = c(7373, 7373, 3771),
                    delta = delta_igg, gamma = gamma_igg,
                    filedir = "resultados/", 
                    varname = "IGG+ sin busqueda con necesidad")
}


# #----------------------------------------------------------------------------
# #PROPORCIONES CON TRES CONDICIONES
# #----------------------------------------------------------------------------
# 
# #Recibió dado que Buscó dado que necesitó e IGG+
# #----------------------------------------------------------------------------
if(!mantener.convergidos){
  nfilter <- paste("!is.na(resultado_cualit_sero) &", 
                   "((NECESIDAD_3MESES == 'Si' & BUSQUEDA_SALUD == 'Si' & !is.na(ATENCION_BUSQUEDA)) | ", 
                   "(NECESIDAD_3MESES == 'Si' & BUSQUEDA_SALUD == 'No'  & is.na(ATENCION_BUSQUEDA)) | ", 
                   "(NECESIDAD_3MESES == 'No' & is.na(BUSQUEDA_SALUD) & is.na(ATENCION_BUSQUEDA)))")
  s1filter <- "NECESIDAD_3MESES == 'Si'"
  s2filter <- "BUSQUEDA_SALUD == 'Si'"
  s3filter <- "ATENCION_BUSQUEDA == 'Si'"
  s4filter <- "resultado_cualit_sero == 'Positive'"
  julia_eval("Random.seed!(3771);") #CADENA 3 PEDOS
  message("Busco necesito recibio e igg+")
  genera_binomial_3(nfilter, s1filter, s2filter, s3filter, s4filter,
                    p1name = "Tuvo necesidad",
                    p2name = "Buscó atención dado que tuvo necesidad",
                    p3name = "Recibió atención dado que buscó atención y tuvo necesidad",
                    p4name = "IGG+ en grupo que recibió atención que buscó y necesitó",
                    seeds = c(9919, 27189, 22),
                    burnin = 50000, niter = c(200000, 100000, 100000),
                    delta = delta_igg, gamma = gamma_igg,
                    filedir = "resultados/", varname = "IGG+ dado que recibió buscó y necesitó")
}

