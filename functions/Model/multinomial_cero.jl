# Vamos a hacer un modelo tipo P(Ocupacion)
# donde n es el numero de gente que tiene una ocupacion
# alpha_vec es un vector de probabilidades de que la gente esté en cierta ocupacion
# y k es matriz de gentes en las ocupaciones
# solo son las mu's para la alpha
# verificar que sumen 1 y si no, las sumo y las normalizo pi/sum(pi)
# alpha es el vector de probas de David y está normalizado para que sume 1
@model multinomial_cero(n, alpha_vec, k, ::Type{T} = Float64) where {T} = begin
    # Para p (con todos los centros)
    #OJO que p es vector de length(alpha) para cada uno de los centros
    # p matriz de 3 x length(alpha)
    p         = Matrix{T}(undef, length(n), length(alpha_vec))

    for i in 1:length(n)
        # Cada renglón se distribuye Dir(alpha_vec)
        p[i, :] ~ Dirichlet(alpha_vec) # es un vector de longitud el numero de ocupaciones

        k[i, :] ~ Multinomial(n[i], p[i, :])
        # OJO que segun yo aqui no hay porque ajustar nada
    end
end;
