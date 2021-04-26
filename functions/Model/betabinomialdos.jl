# Funcion que calcula probabilidades de tipo x | y (ambos unitarios)
# Por ejemplo, IGG | Necesidad
# El problema es de tipo
# s ~ Bin(n, p) con p ~ Beta(alpha_p, beta_p)
# t ~ Bin(s, q) con q ~ Beta(alpha_q, beta_q)
# x, y se calculan con la mu y sigma que nos dió David
@model betabinomialdos(n, s, t, x_p, y_p, x_q, y_q, delta, gamma, ::Type{T} = Float64) where {T} = begin
    # Para p
    alpha_p ~ truncated(Normal(x_p, 0.1), 0, Inf)
    beta_p ~ truncated(Normal(y_p, 0.1), 0, Inf)
    p          = Vector{T}(undef, 3)
    # Para q
    alpha_q ~ truncated(Normal(x_q, 0.1), 0, Inf)
    beta_q ~ truncated(Normal(y_q, 0.1), 0, Inf)
    q_ajustada = Vector{T}(undef, 3)
    q          = Vector{T}(undef, 3)
    for i = 1:length(n)
        #Para p
         p[i] ~ Beta(alpha_p, beta_p)
         s[i] ~ Binomial(n[i], p[i])
         # Para q
         q[i] ~ Beta(alpha_q, beta_q)
         q_ajustada[i] = delta*q[i] + (1-gamma)*(1-q[i])
         t[i] ~ Binomial(s[i], q_ajustada[i])
    end
end;
