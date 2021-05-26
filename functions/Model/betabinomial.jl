# Funcion para calcular la probabilidad posterior de algo tipo P(x)
# por ejemplo, P(IGG) o P(PCR)
@model betabinomial(n, z, x, y, delta, gamma, ::Type{T} = Float64) where {T} = begin
    alpha_p ~ truncated(Normal(x, 0.1), 0, Inf)
    beta_p ~ truncated(Normal(y, 0.1), 0, Inf)
    q_ajustada = Vector{T}(undef, 3)
    p          = Vector{T}(undef, 3)
    for i = 1:length(n)
         p[i] ~ Beta(alpha_p, beta_p)
         q_ajustada[i] = delta*p[i] + (1-gamma)*(1-p[i])
         z[i] ~ Binomial(n[i], q_ajustada[i])
    end
end;
