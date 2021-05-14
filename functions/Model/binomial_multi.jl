# Funcion para calcular la probabilidad posterior de una pregunta que ddetermina si la
# persona tomo medidas en caso de haber tenido sintomas de covid
# una pregunta sintomas = s ~ Binomial(n, p)
#              medidas  = m ~ Multinomial(s, alpha_vec)
#  donde s es la gente que dijo que sí tuvo algun sintoma (el que sea)
#  y m es la matriz que dice cuanta gente tuvo cada uno de los sintomas
@model binomial_multi(n, s, m, x, y, alpha_vec, ::Type{T} = Float64) where {T} = begin
    alpha_p ~ truncated(Normal(x, 0.1), 0, Inf)
    beta_p ~ truncated(Normal(y, 0.1), 0, Inf)
    p          = Vector{T}(undef, length(n))
    q          = Matrix{T}(undef, length(n), length(alpha_vec))
    for i = 1:length(n)
         p[i] ~ Beta(alpha_p, beta_p)
         s[i] ~ Binomial(n[i], p[i])
         q[i, :] ~ Dirichlet(alpha_vec) # es un vector de longitud el numero de ocupaciones
         m[i, :] ~ Multinomial(s[i], q[i, :])
    end
end;
