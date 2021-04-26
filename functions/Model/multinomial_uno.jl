# Funcion que calcula probabilidades de tipo x | y (ambos unitarios)
# Por ejemplo, IGG | Edad < 35 años
# El problema es de tipo
# t ~ Bin(n, q) con q ~ Beta(alpha_q, beta_q)
# k ~ Multinomial(t, p1, ... , pm) con p ~ Dirichlet(alpha_vec)
# x, y se calculan con la mu y sigma que nos dió David
@model multinomial_uno(n, t, k, alpha_vec, x_q, y_q, delta, gamma, ::Type{T} = Float64) where {T} = begin

    # Para q
    alpha_q ~ truncated(Normal(x_q, 0.1), 0, Inf)
    beta_q ~ truncated(Normal(y_q, 0.1), 0, Inf)
    q_ajustada = Vector{T}(undef, length(n)) #este ajuste sí se queda asi, no?
    q          = Vector{T}(undef, length(n))
    p          = Matrix{T}(undef, length(n), length(alpha_vec))
    for i = 1:length(n)
         # Para q
         q[i] ~ Beta(alpha_q, beta_q)
         q_ajustada[i] = delta*q[i] + (1-gamma)*(1-q[i])
         t[i] ~ Binomial(n[i], q_ajustada[i])
         # Cada renglón se distribuye Dir(alpha_vec)
         p[i, :] ~ Dirichlet(alpha_vec) # es un vector de longitud el numero de ocupaciones
         k[i, :] ~ Multinomial(t[i], p[i, :])
    end
end;
