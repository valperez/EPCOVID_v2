##### SCRIPT PARA OBTENER LAS PROBABILIDADES (FINALLY)
## de Valeria Pérez para Rodrigo Zepeda

# delta es la sensitividad y 
# gamma es la especificidad para todos los métodos
# adjusted es si queremos que sea proba ajustada o no


delta <- 0.8       #sensitividad
gamma <- 0.995     #especificida


#función para hacer el ajuste de sensitividad y especificidad
ajuste_sens <- function(x, delta = 1, gamma = 1){
  return(delta*x + (1 - gamma)*(1 - x))
}


check_sum <- function(x){
  suma = sum(x)
  if (suma != 1){
    x <- x / suma
  }
  return(x)
}

# NOTA: los resultados de cuando lo corres te lo guarda en este directorio a menos que lo cambies

# Función para calcular las n's para una probabilidad sin nada condicionado
# ej P(IGG)
# IGG = s ~ Binomial(n, p) con p ~ Beta(alpha, beta) 
calculo_binomial_0 <- function(nsize, zsize, mu = 0.5, sigma = 1/12, 
                               delta = 1, gamma = 1, 
                               compilar_julia = F,
                               burnin = 5000, niter = 10000, nchains = 4, 
                               starting_point = NULL, 
                               proba_nom = str_remove(as.character(runif(1)),"\\.")){
  
  # Carga los metodos que tiene que cargar
  if (compilar_julia){ 
    julia_eval('include("functions/Model/betabinomial.jl")')
    julia_eval('include("functions/Model/calculo_xy.jl")')
    julia_eval('include("functions/Model/check_convergencia.jl")')
  }
  
  tiene_ceros <- is.element(0, nsize)
  
  nsize <- JuliaObject(nsize) #OJO AQUI, tiene que decir "Julia Object of type Array{Float64,1}"
  zsize <- JuliaObject(zsize) # para los 3. Si no, ponle el transpuesto  
  
  julia_assign("nsize", nsize)
  julia_assign("zsize", zsize)
  
  # Calculo de x, y
  julia_assign("mu", mu)
  julia_assign("sigma", sigma)
  x <- julia_eval("x = calculo_xy(mu, sigma)[1]")
  y <- julia_eval("y = calculo_xy(mu, sigma)[2]")
  
  if (is.null(starting_point) & tiene_ceros == F){
    starting_point <- julia_eval("starting_point = [x, y, zsize ./ nsize]")
  } else { 
    if (is.null(starting_point) & tiene_ceros == T){
      starting_point <- julia_eval("starting_point = [x, y, 0, 0, 0]")
    }
  }
  
  # Convertimos el starting point a JuliaObject
  starting_point <- JuliaObject(starting_point)
  julia_assign("starting_point", starting_point)
  
  julia_assign("proba_nom", proba_nom)
  julia_eval("proba_nom = string(proba_nom)")
  burnin <- format(burnin, scientific = F)
  niter  <- format(niter, scientific = F)
  julia_assign("niter", niter)
  julia_assign("burnin", burnin)
  julia_assign("nchains", nchains)
  
  # Aquí solo le estás diciendo "oye guardame en sim el método betabinomial con estas variables"
  julia_eval(paste0("sim       = betabinomial(",nsize,",",zsize,",", 
                    x,",", y,",",
                    delta,",", gamma,")"))
  
  # Aquí es cuando ya corres las simulaciones como lo hice yo
  julia_eval(paste0("hmcsample = sample(sim, HMC(0.01, 5), MCMCThreads(), burnin = ", 
                    burnin, ", ", niter, ", ", nchains,  
                    ", init_theta = starting_point)"))
  
  # En caso de no convergencia: 
  # Paso 1. Checar que tus n's estén bien (neta, checalo)
  # Paso 2. Aumentar el número de iteraciones
  # Paso 3. Si solamente falla la estacionariedad, checa las gráficas y ve donde está fallando
  julia_eval(paste0("check_convergencia(hmcsample, ", niter, ", ", burnin,
                    ", ", nchains, ", proba_nom)"))
  
  # guardamos muestras como DataFrame
  muestras <- julia_eval("DataFrame(hmcsample)") 
  
  # creamos un tibble auxiliar y lo modificamos
  aux <- as_tibble(muestras) %>%
    select(starts_with("p")) %>%
    purrr::modify(., ~ ajuste_sens(., delta, gamma))
  
  # guardamos los nombres y modificamos
  nombres <- colnames(aux)
  nombres <- str_replace(nombres, "p", "q")
  colnames(aux) <- nombres
  
  #Bindeamos todo
  muestras <- cbind(muestras, aux)
  
  return(muestras)
}
# Función para calcular las n's para una probabilidad con una sola variable condicionada
# ej P(IGG | Necesidad)
# donde s = Necesidad ~ Binomial(n, p) con p ~ Beta(alpha_p, beta_p)
# y     t = IGG       ~ Binomial(s, q) con q ~ Beta(alpha_q, beta_q)
# proba_nom es la probabilidad que estás haciendo. Es el nombre con el cual vas a guardar los plots
calculo_binomial_1 <- function(nsize, zsize, tsize, 
                               mu_p = 0.5, sigma_p = 1/12, mu_q = 0.5, sigma_q = 1/12, 
                               delta = 1, gamma = 1, compilar_julia = T,
                               burnin = 50000, niter = 100000, nchains = 4, 
                               starting_point = NULL, 
                               proba_nom = str_remove(as.character(runif(1)),"\\.")){
  
  tiene_ceros <- (is.element(0, nsize) | is.element(0, zsize))
  
  # Carga los metodos que tiene que cargar
  if (compilar_julia){ 
    julia_eval('include("functions/Model/betabinomialdos.jl")')
    julia_eval('include("functions/Model/calculo_xy.jl")')
    julia_eval('include("functions/Model/check_convergencia.jl")')
  }
  
  nsize <- JuliaObject(nsize) #OJO AQUI, tiene que decir "Julia Object of type Array{Float64,1}"
  zsize <- JuliaObject(zsize) # para los 3. Si no, ponle el transpuesto osea 
  tsize <- JuliaObject(tsize) #  nsize <- JuliaObject(t(nsize))
  
  julia_assign("nsize", nsize)
  julia_assign("zsize", zsize)
  julia_assign("tsize", tsize)
  
  # Calculo de x_p, y_p
  julia_assign("mu_p", mu_p)
  julia_assign("sigma_p", sigma_p)
  x_p <- julia_eval("x_p = calculo_xy(mu_p, sigma_p)[1]")
  y_p <- julia_eval("y_p = calculo_xy(mu_p, sigma_p)[2]")
  
  # Calculo de x_q, y_q
  julia_assign("mu_q", mu_q)
  julia_assign("sigma_q", sigma_q)
  x_q <- julia_eval("x_q = calculo_xy(mu_q, sigma_q)[1]")
  y_q <- julia_eval("y_q = calculo_xy(mu_q, sigma_q)[2]")
  
  if (is.null(starting_point) & tiene_ceros == F){
    starting_point <- julia_eval("starting_point = [x_p, x_q, y_p, y_q, 
                                 zsize ./ nsize, tsize ./zsize]")
  } else { 
    if (is.null(starting_point) & tiene_ceros == T){
      starting_point <- julia_eval("starting_point = [x_p, x_q, y_p, y_q, 
                                 0, 0, 0, 0, 0, 0]")
    }
  }
  
  # Convertimos el starting point a JuliaObject
  starting_point <- JuliaObject(starting_point)
  julia_assign("starting_point", starting_point)
  
  
  julia_assign("proba_nom", proba_nom)
  julia_eval("proba_nom = string(proba_nom)")
  burnin <- format(burnin, scientific = F)
  niter  <- format(niter, scientific = F)
  julia_assign("niter", niter)
  julia_assign("burnin", burnin)
  julia_assign("nchains", nchains)
  
  
  # Aquí solo le estás diciendo "oye guardame en sim el método betabinomialdos con estas variables"
  julia_eval(paste0("sim       = betabinomialdos(",nsize,",",zsize,",", tsize, ",",
                    x_p,",", y_p,",",x_q,",", y_q,",",
                    delta,",", gamma,")"))
  
  # Aquí es cuando ya corres las simulaciones como lo hice yo
  julia_eval(paste0("hmcsample = sample(sim, HMC(0.01, 5), MCMCThreads(), burnin = ", 
                    burnin, ", ", niter, ", ", nchains,  
                    ", init_theta = starting_point)"))
  
  
  # En caso de no convergencia: 
  # Paso 1. Checar que tus n's estén bien (neta, checalo)
  # Paso 2. Aumentar el número de iteraciones
  # Paso 3. Si solamente falla la estacionariedad, checa las gráficas y ve donde está fallando
  julia_eval(paste0("check_convergencia(hmcsample, ", niter, ", ", burnin, ", ", nchains, ", proba_nom)"))
  
  # guardamos muestras como DataFrame
  muestras <- julia_eval("DataFrame(hmcsample)") 
  
  # creamos un tibble auxiliar y lo modificamos
  aux <- as_tibble(muestras) %>%
    select(starts_with("q")) %>%
    purrr::modify(., ~ ajuste_sens(., delta, gamma))
  
  # guardamos los nombres y modificamos
  nombres <- colnames(aux)
  nombres <- str_replace(nombres, "q", "ajus_q")
  colnames(aux) <- nombres
  
  #Bindeamos todo
  muestras <- cbind(muestras, aux)
  
  return(muestras)
}

# Función para calcular las n's para una probabilidad con dos variables condicionadas
# ej P(IGG | Necesidad y Buscó)
# donde s1 = Necesidad ~ Binomial(s0, p1) con p1 ~ Beta(alpha_p1, beta_p1)
#       s2 = Buscó     ~ Binomail(s1, p2) con p2 ~ Beta(alpha_p2, beta_p2)
# y     s3 = IGG       ~ Binomial(s2, p3) con p3 ~ Beta(alpha_p3, beta_p3)

calculo_binomial_2 <- function(nsize, s1, s2, s3, 
                               mu_p1 = 0.5, sigma_p1 = 1/12, 
                               mu_p2 = 0.5, sigma_p2 = 1/12, 
                               mu_p3 = 0.5, sigma_p3 = 1/12, 
                               delta = 1, gamma = 1, compilar_julia = T,
                               burnin = 50000, niter = 100000, nchains = 4, 
                               starting_point = NULL, 
                               proba_nom = str_remove(as.character(runif(1)),"\\.")){
  # s0 es el total de gente que sí contesto algo en Necesidad
  # s1 es el total de gente que tuvo necesidad
  # s2 es el total de gente que busco atencion medica y tuvo necesidad
  # s3 es la gente que dió positivo a covid, tuvo necesidad y buscó atención médica
  
  tiene_ceros <- (is.element(0, nsize) | is.element(0, s1) | is.element(0, s2))
  
  # Carga los metodos que tiene que cargar
  if (compilar_julia){ 
    julia_eval('include("functions/Model/betabinomialtres_v2.jl")')
    julia_eval('include("functions/Model/calculo_xy.jl")')
    julia_eval('include("functions/Model/check_convergencia.jl")')
  }
  
  nsize <- JuliaObject(nsize) #OJO AQUI, tiene que decir "Julia Object of type Array{Float64,1}"
  s1 <- JuliaObject(s1) # para los 3. Si no, ponle el transpuesto osea 
  s2 <- JuliaObject(s2) #  nsize <- JuliaObject(t(nsize))
  s3 <- JuliaObject(s3)
  
  julia_assign("nsize", nsize)
  julia_assign("s1", s1)
  julia_assign("s2", s2)
  julia_assign("s3", s3)
  
  # Calculo de x_p, y_p
  # p1
  julia_assign("mu_p1", mu_p1)
  julia_assign("sigma_p1", sigma_p1)
  x_p1 <- julia_eval("x_p1 = calculo_xy(mu_p1, sigma_p1)[1]")
  y_p1 <- julia_eval("y_p1 = calculo_xy(mu_p1, sigma_p1)[2]")
  
  # p2
  julia_assign("mu_p2", mu_p2)
  julia_assign("sigma_p2", sigma_p2)
  x_p2 <- julia_eval("x_p2 = calculo_xy(mu_p2, sigma_p2)[1]")
  y_p2 <- julia_eval("y_p2 = calculo_xy(mu_p2, sigma_p2)[2]")
  
  # p3
  julia_assign("mu_p3", mu_p3)
  julia_assign("sigma_p3", sigma_p3)
  x_p3 <- julia_eval("x_p3 = calculo_xy(mu_p3, sigma_p3)[1]")
  y_p3 <- julia_eval("y_p3 = calculo_xy(mu_p3, sigma_p3)[2]")
  
  if (is.null(starting_point) & tiene_ceros == F){
    starting_point <- julia_eval("starting_point = [x_p1, x_p2, x_p3, y_p1, y_p2, y_p3, 
                                 s1 ./ nsize, s2 ./ s1, s3 ./s2]")
  } else { 
    if (is.null(starting_point) & tiene_ceros == T){
      starting_point <- julia_eval("starting_point = [x_p1, x_p2, x_p3, y_p1, y_p2, y_p3, 
                                 0, 0, 0, 0, 0, 0, 0, 0, 0]")
    }
  }
  
  # Convertimos el starting point a JuliaObject
  starting_point <- JuliaObject(starting_point)
  julia_assign("starting_point", starting_point)
  
  
  julia_assign("proba_nom", proba_nom)
  julia_eval("proba_nom = string(proba_nom)")
  burnin <- format(burnin, scientific = F)
  niter <- format(niter, scientific = F)
  julia_assign("niter", niter)
  julia_assign("burnin", burnin)
  julia_assign("nchains", nchains)
  
  # Aquí solo le estás diciendo "oye guardame en sim el método betabinomialtres con estas variables"
  julia_eval(paste0("sim       = betabinomialtres(",nsize,",",s1,",", s2, ",", s3, ",",
                    x_p1,",", y_p1,",", x_p2,",", y_p2, ",", x_p3,",", y_p3, ",",
                    delta,",", gamma,")"))
  
  # Aquí es cuando ya corres las simulaciones como lo hice yo
  julia_eval(paste0("hmcsample = sample(sim, HMC(0.01, 5), MCMCThreads(), burnin = ", 
                    burnin, ", ", niter, ", ", nchains,  
                    ", init_theta = starting_point)"))
  
  # En caso de no convergencia: 
  # Paso 1. Checar que tus n's estén bien (neta, checalo)
  # Paso 2. Aumentar el número de iteraciones
  # Paso 3. Si solamente falla la estacionariedad, checa las gráficas y ve donde está fallando
  julia_eval(paste0("check_convergencia(hmcsample, ", niter, ", ", burnin, ", ", nchains, ", proba_nom)"))
  
  # guardamos muestras como DataFrame
  muestras <- julia_eval("DataFrame(hmcsample)") 
  
  # creamos un tibble auxiliar y lo modificamos
  aux <- as_tibble(muestras) %>%
    select(starts_with("p3")) %>%
    purrr::modify(., ~ ajuste_sens(., delta, gamma))
  
  # guardamos los nombres y modificamos
  nombres <- colnames(aux)
  nombres <- str_replace(nombres, "p3", "ajus_p3")
  colnames(aux) <- nombres
  
  #Bindeamos todo
  muestras <- cbind(muestras, aux)
  
  return(muestras)
  
}

# Función para calcular las n's para una probabilidad con tres variables condicionadas
# ej P(IGG | Necesidad, Buscó y Recibió)
# donde s1 = Necesidad ~ Binomial(s0, p1) con p1 ~ Beta(alpha_p1, beta_p1)
#       s2 = Buscó     ~ Binomail(s1, p2) con p2 ~ Beta(alpha_p2, beta_p2)
#       s3 = Recibió   ~ Binomial(s2, p3) con p3 ~ Beta(alpha_p3, beta_p3)
#       s4 = IGG       ~ Binomial(s3, p4) con p4 ~ Beta(alpha_p4, beta_p4)
# Nota: le subi el default de iteraciones por experiencia pasada
calculo_binomial_3 <- function(nsize, s1, s2, s3, s4, 
                               mu_p1 = 0.5, sigma_p1 = 1/12,
                               mu_p2 = 0.5, sigma_p2 = 1/12,
                               mu_p3 = 0.5, sigma_p3 = 1/12, 
                               mu_p4 = 0.5, sigma_p4 = 1/12,
                               delta = 1, gamma = 1, compilar_julia = T,
                               burnin = 50000, niter = 100000, nchains = 4, 
                               starting_point = NULL, 
                               proba_nom = str_remove(as.character(runif(1)),"\\.")){
  # s0 es el total de gente que sí contesto algo en Necesidad
  # s1 es el total de gente que tuvo necesidad
  # s2 es el total de gente que busco atencion medica y tuvo necesidad
  # s3 es la gente que tuvo necesidad, buscó atención médica y la recibió
  # s4 es la gente que tuvo necesidad, buscó atención médica, la recibió y dio + a IGG
  
  tiene_ceros <- (is.element(0, nsize) | is.element(0, s1) | is.element(0, s2) | is.element(0, s3))
  
  # Carga los metodos que tiene que cargar
  if (compilar_julia){ 
    julia_eval('include("functions/Model/betabinomialcuatro_v2.jl")')
    julia_eval('include("functions/Model/calculo_xy.jl")')
    julia_eval('include("functions/Model/check_convergencia.jl")')
  }
  
  nsize <- JuliaObject(nsize) #OJO AQUI, tiene que decir "Julia Object of type Array{Float64,1}"
  s1 <- JuliaObject(s1) # para los 3. Si no, ponle el transpuesto osea 
  s2 <- JuliaObject(s2) #  nsize <- JuliaObject(t(nsize))
  s3 <- JuliaObject(s3)
  s4 <- JuliaObject(s4)
  
  julia_assign("nsize", nsize)
  julia_assign("s1", s1)
  julia_assign("s2", s2)
  julia_assign("s3", s3)
  julia_assign("s4", s4)
  
  # Calculo de x_p, y_p
  # p1
  julia_assign("mu_p1", mu_p1)
  julia_assign("sigma_p1", sigma_p1)
  x_p1 <- julia_eval("x_p1 = calculo_xy(mu_p1, sigma_p1)[1]")
  y_p1 <- julia_eval("y_p1 = calculo_xy(mu_p1, sigma_p1)[2]")
  
  # p2
  julia_assign("mu_p2", mu_p2)
  julia_assign("sigma_p2", sigma_p2)
  x_p2 <- julia_eval("x_p2 = calculo_xy(mu_p2, sigma_p2)[1]")
  y_p2 <- julia_eval("y_p2 = calculo_xy(mu_p2, sigma_p2)[2]")
  
  # p3
  julia_assign("mu_p3", mu_p3)
  julia_assign("sigma_p3", sigma_p3)
  x_p3 <- julia_eval("x_p3 = calculo_xy(mu_p3, sigma_p3)[1]")
  y_p3 <- julia_eval("y_p3 = calculo_xy(mu_p3, sigma_p3)[2]")
  
  # p4
  julia_assign("mu_p4", mu_p4)
  julia_assign("sigma_p4", sigma_p4)
  x_p4 <- julia_eval("x_p4 = calculo_xy(mu_p4, sigma_p4)[1]")
  y_p4 <- julia_eval("y_p4 = calculo_xy(mu_p4, sigma_p4)[2]")
  
  if (is.null(starting_point) & tiene_ceros == F){
    starting_point <- julia_eval("starting_point = [x_p1, x_p2, x_p3, x_p4,
                                 y_p1, y_p2, y_p3, y_p4, 
                                 s1 ./ nsize, s2 ./ s1, s3 ./ s2, s4 ./ s3]")
  } else { 
    if (is.null(starting_point) & tiene_ceros == T){
      starting_point <- julia_eval("starting_point = [x_p1, x_p2, x_p3, x_p4,
                                 y_p1, y_p2, y_p3, y_p4, 
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]")
    }
  }
  
  # Convertimos el starting point a JuliaObject
  starting_point <- JuliaObject(starting_point)
  julia_assign("starting_point", starting_point)
  
  
  julia_assign("proba_nom", proba_nom)
  julia_eval("proba_nom = string(proba_nom)")
  burnin <- format(burnin, scientific = F)
  niter <- format(niter, scientific = F)
  julia_assign("niter", niter)
  julia_assign("burnin", burnin)
  julia_assign("nchains", nchains)
  
  # Aquí solo le estás diciendo "oye guardame en sim el método betabinomialcuatro con estas variables"
  julia_eval(paste0("sim       = betabinomialcuatro_v2(",nsize,",",s1,",", s2, ",", s3, ",", s4, ",",
                    x_p1,",", y_p1,",", x_p2,",", y_p2, ",", x_p3,",", y_p3, ",", x_p4,",", y_p4, ",",
                    delta,",", gamma,")"))
  
  # Aquí es cuando ya corres las simulaciones como lo hice yo
  julia_eval(paste0("hmcsample = sample(sim, HMC(0.01, 5), MCMCThreads(), burnin = ", 
                    burnin, ", ", niter, ", ", nchains,  
                    ", init_theta = starting_point)"))
  
  # En caso de no convergencia: 
  # Paso 1. Checar que tus n's estén bien (neta, checalo)
  # Paso 2. Aumentar el número de iteraciones
  # Paso 3. Si solamente falla la estacionariedad, checa las gráficas y ve donde está fallando
  julia_eval(paste0("check_convergencia(hmcsample, ", niter, ", ", burnin, ", ", nchains, ", proba_nom)"))
  
  # guardamos muestras como DataFrame
  muestras <- julia_eval("DataFrame(hmcsample)") 
  
  # creamos un tibble auxiliar y lo modificamos
  aux <- as_tibble(muestras) %>%
    select(starts_with("p4")) %>%
    purrr::modify(., ~ ajuste_sens(., delta, gamma))
  
  # guardamos los nombres y modificamos
  nombres <- colnames(aux)
  nombres <- str_replace(nombres, "p4", "ajus_p4")
  colnames(aux) <- nombres
  
  #Bindeamos todo
  muestras <- cbind(muestras, aux)
  
  return(muestras)
  
}


# Función para calcular las n's para una probabilidad con tres variables condicionadas
# ej P(IGG | Necesidad, Buscó, Recibió, Recibió IMSS)
# donde s1 = Necesidad ~ Binomial(s0, p1) con p1 ~ Beta(alpha_p1, beta_p1)
#       s2 = Buscó     ~ Binomail(s1, p2) con p2 ~ Beta(alpha_p2, beta_p2)
#       s3 = Recibió   ~ Binomial(s2, p3) con p3 ~ Beta(alpha_p3, beta_p3)
#       s4 = Rec IMSS  ~ Binomial(s3, p4) con p4 ~ Beta(alpha_p4, beta_p4)
#       s5 = IGG       ~ Binomial(s4, p5) con p5 ~ Beta(alpha_p5, beta_p5)
# Nota: le subi el default de iteraciones por experiencia pasada
calculo_binomial_4 <- function(nsize, s1, s2, s3, s4, s5, 
                               mu_p1 = 0.5, sigma_p1 = 1/12,
                               mu_p2 = 0.5, sigma_p2 = 1/12,
                               mu_p3 = 0.5, sigma_p3 = 1/12, 
                               mu_p4 = 0.5, sigma_p4 = 1/12,
                               mu_p5 = 0.5, sigma_p5 = 1/12,
                               delta = 1, gamma = 1, compilar_julia = T,
                               burnin = 50000, niter = 100000, nchains = 4, 
                               starting_point = NULL, 
                               proba_nom = str_remove(as.character(runif(1)),"\\.")){
  
  tiene_ceros <- (is.element(0, nsize) | is.element(0, s1) 
                  | is.element(0, s2) | is.element(0, s3) |
                    is.element(0, s4))
  
  # Carga los metodos que tiene que cargar
  if (compilar_julia){ 
    julia_eval('include("functions/Model/betabinomialcinco.jl")')
    julia_eval('include("functions/Model/calculo_xy.jl")')
    julia_eval('include("functions/Model/check_convergencia.jl")')
  }
  
  nsize <- JuliaObject(nsize) #OJO AQUI, tiene que decir "Julia Object of type Array{Float64,1}"
  s1 <- JuliaObject(s1) # para los 3. Si no, ponle el transpuesto osea 
  s2 <- JuliaObject(s2) #  nsize <- JuliaObject(t(nsize))
  s3 <- JuliaObject(s3)
  s4 <- JuliaObject(s4)
  s5 <- JuliaObject(s5)
  
  julia_assign("nsize", nsize)
  julia_assign("s1", s1)
  julia_assign("s2", s2)
  julia_assign("s3", s3)
  julia_assign("s4", s4)
  julia_assign("s5", s5)
  
  # Calculo de x_p, y_p
  # p1
  julia_assign("mu_p1", mu_p1)
  julia_assign("sigma_p1", sigma_p1)
  x_p1 <- julia_eval("x_p1 = calculo_xy(mu_p1, sigma_p1)[1]")
  y_p1 <- julia_eval("y_p1 = calculo_xy(mu_p1, sigma_p1)[2]")
  
  # p2
  julia_assign("mu_p2", mu_p2)
  julia_assign("sigma_p2", sigma_p2)
  x_p2 <- julia_eval("x_p2 = calculo_xy(mu_p2, sigma_p2)[1]")
  y_p2 <- julia_eval("y_p2 = calculo_xy(mu_p2, sigma_p2)[2]")
  
  # p3
  julia_assign("mu_p3", mu_p3)
  julia_assign("sigma_p3", sigma_p3)
  x_p3 <- julia_eval("x_p3 = calculo_xy(mu_p3, sigma_p3)[1]")
  y_p3 <- julia_eval("y_p3 = calculo_xy(mu_p3, sigma_p3)[2]")
  
  # p4
  julia_assign("mu_p4", mu_p4)
  julia_assign("sigma_p4", sigma_p4)
  x_p4 <- julia_eval("x_p4 = calculo_xy(mu_p4, sigma_p4)[1]")
  y_p4 <- julia_eval("y_p4 = calculo_xy(mu_p4, sigma_p4)[2]")
  
  # p4
  julia_assign("mu_p5", mu_p5)
  julia_assign("sigma_p5", sigma_p5)
  x_p5 <- julia_eval("x_p5 = calculo_xy(mu_p5, sigma_p5)[1]")
  y_p5 <- julia_eval("y_p5 = calculo_xy(mu_p5, sigma_p5)[2]")
  
  if (is.null(starting_point) & tiene_ceros == F){
    starting_point <- julia_eval("starting_point = [x_p1, x_p2, x_p3, x_p4,
                                 y_p1, y_p2, y_p3, y_p4, 
                                 s1 ./ nsize, s2 ./ s1, s3 ./ s2, s4 ./ s3, s5 ./ s4]")
  } else { 
    if (is.null(starting_point) & tiene_ceros == T){
      starting_point <- julia_eval("starting_point = [x_p1, x_p2, x_p3, x_p4, x_p5,
                                 y_p1, y_p2, y_p3, y_p4, y_p5, 
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]")
    }
  }
  
  # Convertimos el starting point a JuliaObject
  starting_point <- JuliaObject(starting_point)
  julia_assign("starting_point", starting_point)
  
  
  julia_assign("proba_nom", proba_nom)
  julia_eval("proba_nom = string(proba_nom)")
  burnin <- format(burnin, scientific = F)
  niter <- format(niter, scientific = F)
  julia_assign("niter", niter)
  julia_assign("burnin", burnin)
  julia_assign("nchains", nchains)
  
  # Aquí solo le estás diciendo "oye guardame en sim el método betabinomialcuatro con estas variables"
  julia_eval(paste0("sim       = betabinomialcinco(",nsize,",",s1,",", s2, ",", s3, ",", s4, ",", s5, ",",
                    x_p1,",", y_p1,",", x_p2,",", y_p2, ",", x_p3,",", y_p3, ",", x_p4,",", y_p4, ",", x_p5,",", y_p5, ",",
                    delta,",", gamma,")"))
  
  # Aquí es cuando ya corres las simulaciones como lo hice yo
  julia_eval(paste0("hmcsample = sample(sim, HMC(0.01, 5), MCMCThreads(), burnin = ", 
                    burnin, ", ", niter, ", ", nchains,  
                    ", init_theta = starting_point)"))
  
  # En caso de no convergencia: 
  # Paso 1. Checar que tus n's estén bien (neta, checalo)
  # Paso 2. Aumentar el número de iteraciones
  # Paso 3. Si solamente falla la estacionariedad, checa las gráficas y ve donde está fallando
  julia_eval(paste0("check_convergencia(hmcsample, ", niter, ", ", burnin, ", ", nchains, ", proba_nom)"))
  
  # guardamos muestras como DataFrame
  muestras <- julia_eval("DataFrame(hmcsample)") 
  
  # creamos un tibble auxiliar y lo modificamos
  aux <- as_tibble(muestras) %>%
    select(starts_with("p5")) %>%
    purrr::modify(., ~ ajuste_sens(., delta, gamma))
  
  # guardamos los nombres y modificamos
  nombres <- colnames(aux)
  nombres <- str_replace(nombres, "p5", "ajus_p5")
  colnames(aux) <- nombres
  
  #Bindeamos todo
  muestras <- cbind(muestras, aux)
  
  return(muestras)
  
}


# Función para calcular las n's para una probabilidad de ocupacion
# ej P(Ocupacion)
# donde k ~ Multinomial(nsize, p) y p ~ Dirichlet(alpha_vec) 
# Nota: le subi el default de iteraciones por experiencia pasada
calculo_multinomial_0 <- function(nsize, alpha_vec = NULL, k_matrix, 
                                  compilar_julia = T,
                                  burnin = 50000, niter = 100000, nchains = 4, 
                                  proba_nom = str_remove(as.character(runif(1)),"\\.")){
  
  # Carga los metodos que tiene que cargar
  if (compilar_julia){ 
    julia_eval('include("functions/Model/multinomial_cero.jl")')
    julia_eval('include("functions/Model/check_convergencia.jl")')
  }
  
  if (is.null(alpha_vec)){
    alpha_vec <- rep(1, ncol(k_matrix))
  } else {
    alpha_vec <- check_sum(alpha_vec)
  }
  
  nsize <- JuliaObject(nsize) #OJO AQUI, tiene que decir "Julia Object of type Array{Float64,1}"
  # para los 3. Si no, ponle el transpuesto
  alpha_vec <- JuliaObject(alpha_vec)
  k_matrix <- JuliaObject(k_matrix)
  
  julia_assign("nsize", nsize)
  julia_assign("alpha_vec", alpha_vec)
  julia_assign("k_matrix", k_matrix)
  
  julia_assign("proba_nom", proba_nom)
  julia_eval("proba_nom = string(proba_nom)")
  burnin <- format(burnin, scientific = F)
  niter <- format(niter, scientific = F)
  julia_assign("niter", niter)
  julia_assign("burnin", burnin)
  julia_assign("nchains", nchains)
  
  # Aquí solo le estás diciendo "oye guardame en sim el método multinomial_cero con estas variables"
  julia_eval(paste0("sim       = multinomial_cero(",nsize,",",alpha_vec,",", k_matrix, ")"))
  
  # Aquí es cuando ya corres las simulaciones como lo hice yo
  julia_eval(paste0("hmcsample = sample(sim, HMC(0.01, 5), MCMCThreads(), burnin = ", 
                    burnin, ", ", niter, ", ", nchains, ")"))
  
  # En caso de no convergencia: 
  # Paso 1. Checar que tus n's estén bien (neta, checalo)
  # Paso 2. Aumentar el número de iteraciones
  # Paso 3. Si solamente falla la estacionariedad, checa las gráficas y ve donde está fallando
  julia_eval(paste0("check_convergencia(hmcsample, ", niter, ", ", burnin, ", ", nchains, ", proba_nom)"))
  
  # guardamos muestras como DataFrame
  muestras <- julia_eval("DataFrame(hmcsample)") 
  
  #Recuerda que para este tipo de casos no tenemos que ajustar con sens y espec
  
  muestras <- as_tibble(muestras)
  
  return(muestras)
  
}

# Función para calcular las n's para una probabilidad multinomial condicionada a una binomial
# ej P(Ocupacion)
# donde t ~ Binomial(nsize, q) con q ~ Beta(alpha_q, beta_q)
# y     k ~ Multinomial(nsize, p) y p ~ Dirichlet(alpha_vec) 
# Nota: le subi el default de iteraciones por experiencia pasada
calculo_multinomial_1 <- function(nsize, tsize, alpha_vec = NULL, k_matrix,
                                  mu_q = 0.5, sigma_q = 1/12,
                                  delta = 1, gamma = 1,
                                  compilar_julia = T,
                                  burnin = 50000, niter = 100000, nchains = 4, 
                                  proba_nom = str_remove(as.character(runif(1)),"\\.")){
  
  # Carga los metodos que tiene que cargar
  if (compilar_julia){ 
    julia_eval('include("functions/Model/multinomial_uno.jl")')
    julia_eval('include("functions/Model/calculo_xy.jl")')
    julia_eval('include("functions/Model/check_convergencia.jl")')
  }
  
  if (is.null(alpha_vec)){
    alpha_vec <- rep(1, ncol(k_matrix))
  } else {
    alpha_vec <- check_sum(alpha_vec)
  }
  
  nsize <- as.integer(nsize)
  tsize <- as.integer(tsize)
  
  alpha_vec <- JuliaObject(alpha_vec)
  k_matrix <- JuliaObject(k_matrix)
  
  julia_assign("nsize", nsize)
  julia_assign("tsize", tsize)
  julia_assign("alpha_vec", alpha_vec)
  julia_assign("k_matrix", k_matrix)
  
  # Calculo de x_q, y_q
  julia_assign("mu_q", mu_q)
  julia_assign("sigma_q", sigma_q)
  x_q <- julia_eval("x_q = calculo_xy(mu_q, sigma_q)[1]")
  y_q <- julia_eval("y_q = calculo_xy(mu_q, sigma_q)[2]")
  
  julia_assign("proba_nom", proba_nom)
  julia_eval("proba_nom = string(proba_nom)")
  burnin <- format(burnin, scientific = F)
  niter <- format(niter, scientific = F)
  julia_assign("niter", niter)
  julia_assign("burnin", burnin)
  julia_assign("nchains", nchains)
  
  # Aquí solo le estás diciendo "oye guardame en sim el método multinomial_cero con estas variables"
  julia_eval(paste0("sim       = multinomial_uno(",nsize,",", tsize,"," ,k_matrix,",",
                    alpha_vec,",", x_q, ",", y_q, ", ", delta, ", ", gamma,  ")"))
  
  # Aquí es cuando ya corres las simulaciones como lo hice yo
  julia_eval(paste0("hmcsample = sample(sim, HMC(0.01, 5), MCMCThreads(), burnin = ", 
                    burnin, ", ", niter, ", ", nchains, ")"))
  
  # En caso de no convergencia: 
  # Paso 1. Checar que tus n's estén bien (neta, checalo)
  # Paso 2. Aumentar el número de iteraciones
  # Paso 3. Si solamente falla la estacionariedad, checa las gráficas y ve donde está fallando
  julia_eval(paste0("check_convergencia(hmcsample, ", niter, ", ", burnin, ", ", nchains, ", proba_nom)"))
  
  # guardamos muestras como DataFrame
  muestras <- julia_eval("DataFrame(hmcsample)") 
  
  # creamos un tibble auxiliar y lo modificamos
  aux <- as_tibble(muestras) %>%
    select(starts_with("q")) %>%
    purrr::modify(., ~ ajuste_sens(., delta, gamma))
  
  # guardamos los nombres y modificamos
  nombres <- colnames(aux)
  nombres <- str_replace(nombres, "q", "ajus_q")
  colnames(aux) <- nombres
  
  #Bindeamos todo
  muestras <- cbind(muestras, aux)
  
  return(muestras)
}

# Función para calcular las n's para la gente que vive en su casa 
# donde s ~ Poisson(lambda) y lambda ~ Gamma(x_p_gamma, y_p_gamma) 
calculo_poisson_0 <- function(zsize, x_p_gamma = 1/3, y_p_gamma = 0, 
                                  compilar_julia = T,
                                  burnin = 50000, niter = 100000, nchains = 4, 
                                  proba_nom = str_remove(as.character(runif(1)),"\\.")){
  
  # Carga los metodos que tiene que cargar
  if (compilar_julia){ 
    julia_eval('include("functions/Model/gammapoissoncero.jl")')
    julia_eval('include("functions/Model/check_convergencia.jl")')
  }
  
  
  zsize <- JuliaObject(zsize) #OJO AQUI, tiene que decir "Julia Object of type Array{Float64,1}"
  # para los 3. Si no, ponle el transpuesto
  
  julia_assign("zsize", zsize)
  
  # Calculo de x_p_gamma, y_p_gamma
  julia_assign("x_p_gamma", x_p_gamma)
  julia_assign("y_p_gamma", y_p_gamma)
  
  julia_assign("proba_nom", proba_nom)
  julia_eval("proba_nom = string(proba_nom)")
  burnin <- format(burnin, scientific = F)
  niter <- format(niter, scientific = F)
  julia_assign("niter", niter)
  julia_assign("burnin", burnin)
  julia_assign("nchains", nchains)
  
  # Aquí solo le estás diciendo "oye guardame en sim el método multinomial_cero con estas variables"
  julia_eval(paste0("sim       = gammapoissoncero(",zsize,",",x_p_gamma,",", y_p_gamma, ")"))
  
  # Aquí es cuando ya corres las simulaciones como lo hice yo
  julia_eval(paste0("hmcsample = sample(sim, HMC(0.01, 5), MCMCThreads(), burnin = ", 
                    burnin, ", ", niter, ", ", nchains, ")"))
  
  # En caso de no convergencia: 
  # Paso 1. Checar que tus n's estén bien (neta, checalo)
  # Paso 2. Aumentar el número de iteraciones
  # Paso 3. Si solamente falla la estacionariedad, checa las gráficas y ve donde está fallando
  julia_eval(paste0("check_convergencia(hmcsample, ", niter, ", ", burnin, ", ", nchains, ", proba_nom)"))
  
  # guardamos muestras como DataFrame
  muestras <- julia_eval("DataFrame(hmcsample)") 
  
  #Recuerda que para este tipo de casos no tenemos que ajustar con sens y espec
  
  muestras <- as_tibble(muestras)
  
  return(muestras)
}

# Función para calcular como se distribuyen la edad de los trabajadores 
# donde k ~ LogNormal(m, p) con
# m ~ Normal(mu, tau)
# p ~ Gamma(alpha, beta)
# OJO: k es matriz donde las columnas deben ser 3 (1 por centro)
calculo_lognormal <- function(k, alpha = 1, beta = 1, 
                              tau = 1, mu = 0,
                              compilar_julia = T,
                              burnin = 5000, niter = 10000, nchains = 4, 
                              proba_nom = str_remove(as.character(runif(1)),"\\.")){
  
  # Carga los metodos que tiene que cargar
  if (compilar_julia){ 
    julia_eval('include("functions/Model/lognormal.jl")')
    julia_eval('include("functions/Model/check_convergencia.jl")')
  }
  
  k <- JuliaObject(k) 
  
  julia_assign("k", k)
  
  # Asigacion de parametros
  julia_assign("alpha", alpha)
  julia_assign("beta", beta)
  julia_assign("tau", tau)
  julia_assign("mu", mu)
  
  julia_assign("proba_nom", proba_nom)
  julia_eval("proba_nom = string(proba_nom)")
  burnin <- format(burnin, scientific = F)
  niter <- format(niter, scientific = F)
  julia_assign("niter", niter)
  julia_assign("burnin", burnin)
  julia_assign("nchains", nchains)
  
  # Aquí solo le estás diciendo "oye guardame en sim el método multinomial_cero con estas variables"
  julia_eval(paste0("sim       = lognormal(",k,",",alpha,",", beta, ",", 
                    tau,",", mu,")"))
  
  # Aquí es cuando ya corres las simulaciones como lo hice yo
  julia_eval(paste0("hmcsample = sample(sim, HMC(0.01, 5), MCMCThreads(), burnin = ", 
                    burnin, ", ", niter, ", ", nchains, ")"))
  
  # En caso de no convergencia: 
  # Paso 1. Checar que tus n's estén bien (neta, checalo)
  # Paso 2. Aumentar el número de iteraciones
  # Paso 3. Si solamente falla la estacionariedad, checa las gráficas y ve donde está fallando
  julia_eval(paste0("check_convergencia(hmcsample, ", niter, ", ", burnin, ", ", nchains, ", proba_nom)"))
  
  # guardamos muestras como DataFrame
  muestras <- julia_eval("DataFrame(hmcsample)") 
  
  #Recuerda que para este tipo de casos no tenemos que ajustar con sens y espec
  
  muestras <- as_tibble(muestras)
  
  return(muestras)
}