get_mean_beta <- function(lower, upper) {
  alpha <-
    ((1 - mean(c(lower, upper))) / var(c(lower, upper)) - 1 / mean(c(lower, upper))) * mean(c(lower, upper)) ^ 2
  beta <- alpha * (1 / mean(c(lower, upper)) - 1)
  mean_beta <- alpha / (alpha + beta)
  return(mean_beta)
}