# Vamos a hacer un modelo tipo IGG | Necesitó y Buscó (dos variables condicionadas)
# El ejemplo es
# s1 = Necesitó ~ Bin(s0, p1)
# s2 = Buscó ~ Bin(s1, p2)
# s3 = IGG ~ Bin(s2, p3)
# OJO con el orden
@model betabinomialtres(n, s1, s2, s3, x_p1, y_p1, x_p2, y_p2, x_p3, y_p3, delta, gamma, ::Type{T} = Float64) where {T} = begin
    # Para p
    alpha_p1 ~ truncated(Normal(x_p1, 0.1), 0, Inf)
    beta_p1 ~ truncated(Normal(y_p1, 0.1), 0, Inf)
    p1         = Vector{T}(undef, 3)
    # Para p2
    alpha_p2 ~ truncated(Normal(x_p2, 0.1), 0, Inf)
    beta_p2 ~ truncated(Normal(y_p2, 0.1), 0, Inf)
    p2        = Vector{T}(undef, 3)
    # Para p3
    alpha_p3 ~ truncated(Normal(x_p3, 0.1), 0, Inf)
    beta_p3 ~ truncated(Normal(y_p3, 0.1), 0, Inf)
    p3 = Vector{T}(undef, 3)
    # p3_ajustada es el que corresponde a IGG
    p3_ajustada = Vector{T}(undef, 3)

    for i = 1:length(n)
        #Para p1
         p1[i] ~ Beta(alpha_p1, beta_p1)
         s1[i] ~ Binomial(n[i], p1[i])
         # Para p2
         p2[i] ~ Beta(alpha_p2, beta_p2)
         s2[i] ~ Binomial(s1[i], p2[i])
         #Para p3
         p3[i] ~ Beta(alpha_p3, beta_p3)
         s3[i] ~ Binomial(s2[i], p3[i])
         p3_ajustada[i] = delta*p3[i] + (1-gamma)*(1-p3[i])
         s3[i] ~ Binomial(s2[i], p3_ajustada[i])
    end
end;
