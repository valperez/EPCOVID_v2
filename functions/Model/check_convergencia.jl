# Función que checa si se cumplen o no las pruebas de convergencia
# que estamos aplicando
# hmcsample es la sample que corrimos con el modelo debido
# proba_nom es el nombre de la probabilidad que estamos calculando (solo para
# guardar los plots con nombres bonitos)
function check_convergencia(hmcsample, niter, nburnini, nchains, proba_nom)
    println("Empieza el chequeo de convergencia")
    # Plot
    #plot(hmcsample)
    #savefig("plot"*proba_nom)

    # Autocorrelacion
    #High autocorrelations within chains indicate slow mixing and
    #slow convergence
    #autocorplot(hmcsample)
    #savefig("autocor"*proba_nom)

    #Gelman
    #Usualmente necesitamos que sea menor a 1.1
    gelman = gelmandiag(hmcsample, transform = true)
    psrf = gelman[:, :psrf]
    renglones = size(gelman)[1] #obtiene el numero de renglones
    for i in 1:renglones
        if psrf[i] > 1.1
            println("Error, un parámetro no está cumpliendo con Gelman")
        end
    end

    # Heidel
    #Te dice si la cadena de Markov viene de una distribucion estacionaria
    #statitonarity test dice convergencia con un 1 y no convergencia con un 0
    #la hip nula es que la cadena viene de una distribucion estacionaria
    heidel = heideldiag(hmcsample)
    for j in 1:nchains
        prueba_heidel = heidel[j]
        for k in 1:(size(heidel[1])[1])
            estacionariedad = prueba_heidel[:, :stationarity]
            for u in 1:length(estacionariedad)
                if estacionariedad[u] != 1
                    println("Hay un problema con la estacionariedad en la cadena ", j)
                end
            end
        end
    end

    # Ess = Effective sample size
    # is the number of effectively independent draws from the posterior
    # distribution that the Markov chain is equivalent to
    #tiene que ser menor que el numero de iteraciones menos el warmup
    ess_sample = ess(hmcsample)
    ess_num = ess_sample[:, :ess]
    for v in 1:length(ess_num)
        if (ess_num[v] > (niter - nburnini))
            println("ess es mayor de lo que debería")
        end
    end

    # No estoy segura que hacer con el summary
    summary = summarize(hmcsample)

    println("Termina el chequeo de convergencia")
end
