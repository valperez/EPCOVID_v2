#
# Vamos a hacer un modelo tipo IGG | Necesitó, Buscó y recibió (betabinomialtres variables condicionadas)
# El ejemplo es
# s1 = Necesitó ~ Bin(s0, p1)
# s2 = Buscó ~ Bin(s1, p2)
# s3 = Recibió ~ Bin(s2, p3)
# s4 = IGG ~ Bin(s3, p4)
# OJO con el orden
@model betabinomialcuatro_v2(n, s1, s2, s3, s4, x_p1, y_p1, x_p2, y_p2, x_p3, y_p3, x_p4, y_p4, delta, gamma, ::Type{T} = Float64) where {T} = begin
    # Para p1
    alpha_p1 ~ truncated(Normal(x_p1, 0.1), 0, Inf)
    beta_p1 ~ truncated(Normal(y_p1, 0.1), 0, Inf)
    p1      ~ filldist(Beta(alpha_p1, beta_p1), length(n))
    # Para p2
    alpha_p2 ~ truncated(Normal(x_p2, 0.1), 0, Inf)
    beta_p2 ~ truncated(Normal(y_p2, 0.1), 0, Inf)
    p2      ~ filldist(Beta(alpha_p2, beta_p2), length(n))
    # Para p3
    alpha_p3 ~ truncated(Normal(x_p3, 0.1), 0, Inf)
    beta_p3 ~ truncated(Normal(y_p3, 0.1), 0, Inf)
    p3      ~ filldist(Beta(alpha_p3, beta_p3), length(n))
    # Para p4
    alpha_p4 ~ truncated(Normal(x_p4, 0.1), 0, Inf)
    beta_p4 ~ truncated(Normal(y_p4, 0.1), 0, Inf)
    p4      ~ filldist(Beta(alpha_p4, beta_p4), length(n))
    # p4_ajustada es el que corresponde a IGG
    p4_ajustada = Vector{T}(undef, 3)
    unos = ones(length(n))
    p4_ajustada = delta*p4 + (1 - gamma)*(unos - p4)

    @. s1 ~ Binomial(n, p1)
    @. s2 ~ Binomial(s1, p2)
    @. s3 ~ Binomial(s2, p3)
    @. s4 ~ Binomial(s3, p4)
    @. s4 ~ Binomial(s3, p4_ajustada)
end;
