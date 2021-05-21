# Función para calcular una pregunta estilo "Cuál es tu edad?"
# k es una matriz donde cada columna representa un centro de dist
@model lognormal(k, alpha, beta, lambda, ::Type{T} = Float64) where {T} = begin
    p          = Vector{T}(undef, size(k, 2))
    m          = Vector{T}(undef, size(k, 2))
    for i = 1:size(k, 2)
         p[i] ~ Gamma(alpha, beta)
         m[i] ~ Exponential(lambda)
         k[i] ~ LogNormal(m[i], p[i])
    end
end;
