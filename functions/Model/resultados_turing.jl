using Random
using Plots
using Turing
using Distributions
using StatsPlots
using DataFrames
using MCMCChains
using CSV
#using TableView

Random.seed!(99);
nchains = 4
niter = 10000
nburnini = 5000

cd("/Users/rod/Dropbox/Coronavirus/EPCOVID_V2/functions/Model/")

include("betabinomial.jl")
include("betabinomialdos.jl")
include("betabinomialtres.jl")
include("betabinomialcuatro_v2.jl")
include("multinomial_cero_uncentro.jl")
include("multinomial_cero.jl")
include("multinomial_uno.jl")

#---Pruebas de convergencia
#Ponemos el working directory
#cd("C:/Users/Valeria/Documents/IMSS/EPCOVID/Resultados/Convergencia")
# Buen link donde explica los tests
#https://www2.math.su.se/matstat/reports/master/2011/rep2/report.pdf

include("check_convergencia.jl")
include("calculo_xy.jl")

#IGG
#delta = 0.8
#gamma = 0.995

#PCR
delta = 0.86
gamma = 0.9999

x = calculo_xy(0.05, 0.02375)[1]
y = calculo_xy(0.05, 0.02375)[2]

n = [165, 879, 473]
z = [0, 12, 3]

sim = betabinomial(n, z, x, y, delta, gamma)
hmcsample = sample(sim, HMC(0.01, 5), MCMCThreads(),
            burnin = nburnini, niter, nchains,
            init_theta = [1, 1, 0, 0, 0]) #starting point
check_convergencia(hmcsample, niter, burnin, proba_num)

# Vamos a hacer IGG | Necesidad
s = [34, 272, 132]
t = [15, 114, 65]
n = [162, 737, 450]

delta = 0.8       #sensitividad
gamma = 0.995     #especificidad

x_p = calculo_xy(0.15, 0.06375)[1]
y_p = calculo_xy(0.15, 0.06375)[2]

x_q = calculo_xy(0.4, 0.12)[1]
y_q = calculo_xy(0.4, 0.12)[2]

sim = betabinomialdos(n, s, t, x_p, y_p, x_q, y_q, delta, gamma)
hmcsample = sample(sim, HMC(0.01, 5), MCMCThreads(),
            burnin = 100000, 200000, nchains,
            init_theta = [x_p, x_q, y_p, y_q, s./n, t./s])
check_convergencia(hmcsample, 200000, 100000, 4, "prueba_subirVarianza")
df = DataFrame(hmcsample)b
CSV.write("C:/Users/Valeria/Documents/IMSS/EPCOVID/Resultados/Convergencia\\IGGdadoNecesidad.csv", df)


#---Vamos a hacer PCR | Hombres
s = [152, 824, 338]
t = [0, 9, 2]
n = [165, 879, 474]

delta = 0.86       #sensitividad
gamma = 0.9999     #especificidad

#Hombres
x_p = calculo_xy(0.8, 0.08)[1]
y_p = calculo_xy(0.8, 0.08)[2]

#PCR | Hombres
x_q = calculo_xy(0.82, 0.0738)[1]
y_q = calculo_xy(0.82, 0.0738)[2]

sim = betabinomialdos(n, s, x_p, y_p, x_q, y_q, delta, gamma, t)
hmcsample = sample(sim, HMC(0.01, 5), MCMCThreads(),
            burnin = nburnini, niter, nchains,
            init_theta = [1, 1, 1, 1, 0, 0, 0, 0, 0, 0]) #starting point =#


#----Vamos a probarlo con IGG | Necesitó y Buscó

# n es la cantidad de gente que no dio NA a necesito + no dio NA a IGG + no dio NA a Busco dado que Necesito = Si
# FIXME
n = [162, 737, 450]
# s1 es la cantidad de gente que necesitó
s1 = [34, 272, 132]
#s2 es la gente que Busco dado que necesitó
s2 = [31, 213, 109]
# s3 es la gente que tuvo IGG dado que Busco
# FIXME que el filtro que sea de los tres
s3 = [13, 99, 56]

#Necesitó
x_p1 = calculo_xy(0.15, 0.06375)[1]
y_p1 = calculo_xy(0.15, 0.06375)[2]

#Busco | Necesito
x_p2 = calculo_xy(0.7, 0.105)[1]
y_p2 = calculo_xy(0.7, 0.105)[2]

#IGG | Busco, Necesito
#esto no está correcto, solo quiero ver que funcione
x_p3 = calculo_xy(0.5, 0.0833)[1]
y_p3 = calculo_xy(0.5, 0.0833)[2]

delta = 0.8
gamma = 0.95

sim = betabinomialtres(n, s1, s2, s3, x_p1, y_p1, x_p2, y_p2, x_p3, y_p3, delta, gamma)

#OJO QUE AUMENTÉ EL NUMERO DE ITERACIONES
# TAL VEZ sea por la proba que estaba mal(?)
hmcsample = sample(sim, HMC(0.01, 5), MCMCThreads(),
            burnin = 6000, 12000, nchains,
            init_theta = [1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0])
check_convergencia(hmcsample, 12000, 6000, "IGG_dadoNecesitoyBusco")

#----Vamos a probarlo con IGG | Necesitó, busco y recibio

# n es la cantidad de gente que no dio NA a necesito + no dio NA a IGG + no dio NA a Busco dado que Necesito = Si
# FIXME
n = [162, 737, 450]
# s1 es la cantidad de gente que necesitó
s1 = [34, 272, 132]
#s2 es la gente que Busco dado que necesitó
s2 = [22, 177, 93]
# s3 es la gente que recibio dado que necesitó y buscó
# FIXME que el filtro que sea de los tres
s3 = [21, 166, 91]
s4 = [10, 83, 51]

#k_y_t(ciclo1, c("busqueda_atencion_covid", "necesidad_3meses"), c("Si", "Si"))
#k_y_t(ciclo1, c("atencion_busqueda","busqueda_atencion_covid", "necesidad_3meses"), c("Si", "Si", "Si"))
#k_y_t(ciclo1, c("resultado_cualit_sero", "atencion_busqueda","busqueda_atencion_covid", "necesidad_3meses"), c("Positive","Si", "Si", "Si"))

#Necesitó
x_p1 = calculo_xy(0.15, 0.06375)[1]
y_p1 = calculo_xy(0.15, 0.06375)[2]

#Busco | Necesito
x_p2 = calculo_xy(0.7, 0.105)[1]
y_p2 = calculo_xy(0.7, 0.105)[2]

#IGG | Busco, Necesito
#esto no está correcto, solo quiero ver que funcione
x_p3 = calculo_xy(0.4571, 0.1241)[1]
y_p3 = calculo_xy(0.4571, 0.1241)[2]

#IGG | Busco, Necesito, Recibio
#esto no está correcto, solo quiero ver que funcione
x_p4 = calculo_xy(0.4668, 0.1244)[1]
y_p4 = calculo_xy(0.4668, 0.1244)[2]

delta = 0.8
gamma = 0.95

sim = betabinomialcuatro_v2(n, s1, s2, s3, s4, x_p1, y_p1, x_p2, y_p2, x_p3, y_p3, x_p4, y_p4, delta, gamma)

nchains = 6
#OJO QUE AUMENTÉ EL NUMERO DE ITERACIONES
# TAL VEZ sea por la proba que estaba mal(?)
@time hmcsample = sample(sim, HMC(0.01, 5), MCMCThreads(),
            burnin = 50000, 100000, 6,
            init_theta = [1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
check_convergencia(hmcsample, 100000, 50000, "IGG_dadoNecesitoyBuscoyRecibio")
# bueno nunca paso la prueba de estacionariedad pero segun las imagenes,
# sí convirgio

# ---- Prueba para con un solo centro
# Ya saliooooo
n = 162
alpha_vec = [0.1, 0.9] #OJO QUE ESTO TIENE QUE SUMAR 1
k = [19, 140] # y esto tiene que sumar la n!!!!!!!!!!!!!!
delta = 0.8
gamma = 0.95
nchains = 4
prueba_uno = multinomial_cero_uncentro(n, alpha_vec, k, delta, gamma)
hmcsample = sample(prueba_uno, HMC(0.01, 5), MCMCThreads(),
            burnin = 1000, 2000, 4)

# -------------- Pruebas para la multinomial_cero
# QUEEEEEEEEEEEN YA SALIO
 #DUDA: n es gente que no puso NA en los centros o en la ocupacion?? AMBOS
ns = [162, 737, 450] #ocupacion
alpha_vec = [0.1, 0.15, 0.2, 0.05, 0.01, 0.015, 0.015, 0.08, 0.2, 0.08, 0.1]
# que pedo con jefe de departamento y jefe, supervisor
# y con trabajadores en servicios en lo de David
k = [0 31 1 52 0 17 8 6 23 16 8
    19 35 372 79 2 49 57 1 101 18 4
     4 191 1 22 0 20 39 104 44 22 3]
delta = 0.8
gamma = 0.95
nchains = 4

prueba = multinomial_cero(ns, alpha_vec, k)
@time hmcsample = sample(prueba, HMC(0.01, 5), MCMCThreads(),
            burnin = 350000, 700000, 4)
check_convergencia(hmcsample, 700000, 350000, 4, "prueba_ocupacion")




# ---- Vamos a hacer la prueba para multinomial_uno
# Gente que sale diferente de NA al test de IGG
n = [165, 902, 473]
#Gente que sale positiva a IGG
t = [40, 264, 181]
# Gente de las edades
# que raro, el maximo es 60 osea sobra una categoria ????
k = [17 13 10
     65 67 132
     42 69 70]
alpha_vec = [0.3, 0.3, 0.4]
x_q = calculo_xy(0.2, 0.08)[1]
y_q = calculo_xy(0.2, 0.08)[2]
delta = 0.8
gamma = 0.95
prueba_dos = multinomial_uno(n, t, k, alpha_vec, x_q, y_q, delta, gamma)
hmcsample = sample(prueba_dos, HMC(0.01, 5), MCMCThreads(),
            burnin = 1000, 2000, 4)


# Vamos a responder la pregunta "Cuantas personas viven en su casa?"
# donde s es la cantidad de gente que respondió a esa pregunta por centro
# x_p, y_p vienen de los aprioris de David, maybe??????
#No está bien, solo es para saber si funciona
s = [100, 200, 300]
x_p = calculo_xy(0.5, 0.0833)[1]
y_p = calculo_xy(0.5, 0.0833)[2]

sim = gammapoissoncero(s, x_p, y_p)
@time hmcsample = sample(sim, HMC(0.01, 5), MCMCThreads(),
            burnin = 5000, 10000, 4)
            #Hacer el chequeo de convergencia
