# Vamos a hacer un modelo tipo IGG | Necesitó, Buscó, recibió y recibió IMSS (cuatro variables condicionadas)
# El ejemplo es
# s1 = Necesitó ~ Bin(s0, p1)
# s2 = Buscó ~ Bin(s1, p2)
# s3 = Recibió ~ Bin(s2, p3)
# s4 = Recibió IMSS ~ Bin(s3, p4) (Binomial porque recibió en el imss o no)
# s5 = IGG ~ Bin(s4, p5)
# OJO con el orden
@model betabinomialcinco(n, s1, s2, s3, s4, s5, x_p1, y_p1, x_p2, y_p2, x_p3, y_p3, x_p4, y_p4, x_p5, y_p5, delta, gamma, ::Type{T} = Float64) where {T} = begin
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
    # Para p5
    alpha_p5 ~ truncated(Normal(x_p5, 0.1), 0, Inf)
    beta_p5 ~ truncated(Normal(y_p5, 0.1), 0, Inf)
    p5      ~ filldist(Beta(alpha_p5, beta_p5), length(n))

    # p5_ajustada es el que corresponde a IGG
    p5_ajustada = Vector{T}(undef, 3)
    unos = ones(length(n))
    p5_ajustada = delta*p5 + (1 - gamma)*(unos - p5)

    @. s1 ~ Binomial(n, p1)
    @. s2 ~ Binomial(s1, p2)
    @. s3 ~ Binomial(s2, p3)
    @. s4 ~ Binomial(s3, p4)
    @. s5 ~ Binomial(s4, p5)
    @. s5 ~ Binomial(s4, p5_ajustada)
end;
