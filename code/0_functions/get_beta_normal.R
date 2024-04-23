# mean to beta params
get_beta_normal <- function(mu, sd) {
  var <- sd^2
  alpha <- ((1 - mu) / var - 1 / mu) * mu^2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}
