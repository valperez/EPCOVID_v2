function calculo_xy(mu, sigma2)
    x = -(mu*(sigma2 + mu^2 - mu))/sigma2
    y = (sigma2 + mu^2 - mu)*(mu - 1)/sigma2
    return x, y
end
