# Función para calcular una pregunta estilo "Cuantas personas viven en su casa?"
@model gammapoissoncero(s, x_p_gamma, y_p_gamma, ::Type{T} = Float64) where {T} = begin

    alpha_p ~ truncated(Normal(x_p_gamma, 0.1), 0, Inf)
    beta_p ~ truncated(Normal(y_p_gamma, 0.1), 0, Inf)
    lambda          = Vector{T}(undef, 3)
    for i = 1:length(n)
         lambda[i] ~ Gamma(alpha_p, beta_p)
         s[i] ~ Poisson(lambda[i])
    end
end;
