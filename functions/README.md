# EPCOVID
Este repositorio contiene los archivos necesarios para generar distribuciones a posteriori de un modelo Beta-Binomial. 

El archivo más importante es el de "resultados_finales.R". En el se encuentran las funciones para calcular las probabilidades 
posteriores de los modelos beta binomial así como el chequeo de convergencia (check_convergencia.jl). 
Los archivos que comienzan con "betabinomial" o "multinomial" son los códigos en Julia que calculan las distribuciones necesarias. 

Los otros scripts de R son intentos fallidos de resolver el problema. Favor de no hacerles mucho caso. 

----- IGNORE THIS

El primer archivo que se debe correr es el de generación_resultados.R hasta la línea 86 (setwd("~/IMSS/Bimbo/Resultados")) 
para que cargue los archivos necesarios para que funcione el modelo. Este archivo también contiene un procedimiento (fallido) para
calcular una probabilidad multicondicional de la forma x | y1, y2, y3

El segundo archivo que se debe correr es el R2Julia.R que contiene los comandos que funcionan en Julia para correr un modelo (no fallido)
que calcula las probailidades a posteriori de una distribución con un solo condicionamiento. 

El tercer archivo que hay que correr es resultados_turing.jl que contiene los modelos en Julia así como los test de convergencia que 
estamos tomando en cuenta para evaluar los resultados. 

--- Los siguientes archivos son intentos fallidos de responder las preguntas de interés. Por lo tanto, solo se utilizan algunas partes
--- o como backup para comparar resultados. 

El archivo generacion_distribuciones.R contiene las funciones que generan las distribuciones posteriores de una distribución beta binomial.
Es decir, contiene el calculo de una Beta(alpha + s, beta + n - s) y Beta(alpha + t, beta + s -t).

El archivo de generacion_informacion.R contiene las funciones para obtener la información necesaria para cada probabilidad sin tener que 
calcular caso por caso. Además, contiene una función que te calcula la media, los cuantiles y genera el archivo de salida con los resultados. 

Por último (y, en este caso el menos importante) fue el primer acercamiento a la solución del problema. Se llama generacion_parametros_bayesiana.R
y contiene todas funciones que vienen en los otros archivos pero mucho menos accesibles. Se sugiere ampliamente no utilizarlo (pero lo dejare hasta
que los resultados salgan).
