@model betabinomialdos_upgraded(n, s, t, w, x_p, y_p, x_q, y_q, delta, gamma, ::Type{T} = Float64) where {T} = begin
    # Para p
    alpha_p ~ truncated(Normal(x_p, 0.1), 0, Inf)
    beta_p  ~ truncated(Normal(y_p, 0.1), 0, Inf)
    p          = Vector{T}(undef, 3)
    # Para q
    alpha_q ~ truncated(Normal(x_q, 0.1), 0, Inf)
    beta_q ~ truncated(Normal(y_q, 0.1), 0, Inf)
    q_ajustada = Vector{T}(undef, 3)
    q          = Vector{T}(undef, 3)

    q_noajustada = Vector{T}(undef, 3)
    qno          = Vector{T}(undef, 3)
    for i = 1:length(n)
        #Para p
        # ITT
         p[i] ~ Beta(alpha_p, beta_p)
         s[i] ~ Binomial(n[i], p[i]) #GENTE CON ITT

         #q = IGG | ITT
         q[i] ~ Beta(alpha_q, beta_q)
         q_ajustada[i] = delta*q[i] + (1-gamma)*(1-q[i])
         t[i] ~ Binomial(s[i], q_ajustada[i])

         #q_no = IGG | NO ITT
         qno[i] ~ Beta(alpha_q, beta_q)
         q_noajustada[i] = delta*qno[i] + (1-gamma)*(1-qno[i])
         w[i] ~ Binomial(n[i] - s[i], q_noajustada[i]) #Gente que no tiene ITT
    end
end;
