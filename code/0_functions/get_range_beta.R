get_range_beta <- function(prop, n) {
  var <- prop * (1 - prop) / n
  sd <- sqrt(var)
  lower <- prop - qnorm(0.975) * sd
  upper <- prop + qnorm(0.975) * sd
  return(list(sd = sd, lower = lower, upper = upper))
}
