`%+%` <- function(x,y) paste(x,y,sep="")

# Task 2 ------------------------------------------------------------------

simYgivenTheta = function(theta, w, n) {
  J = length(theta)
  Y = matrix(nrow = J, ncol = n)
  for (j in 1:J) {
    Y[j,] = rpois(n, w[j]*theta[j])
  }
  return(Y)
}

# Tasks 3, 4, 5 -----------------------------------------------------------

DrawTheta = function(J, mu, sd) {
  return(exp(rnorm(J, mu, sd)))
}

CalculateCoverage = function(Theta.true, logTheta.matrix) {
  # Computes the coverage using 68% and 95% posterior intervals for each
  # Theta_j using MCMC output
  #
  # Args:
  #   logTheta.true: True value of logTheta
  #   logTheta.matrix: J x ndraw logTheta matrix output by poisson.logn.mcmc
  #
  # Returns:
  #   J x 3 vector
  #   First column contains logTheta values
  #   Second column contains 1 if true value is in 68% interval, else 0
  #   Third column contains 1 if true value is in 95% interval, else 0

  logTheta.true = log(Theta.true)
  J = length(logTheta.true)
  if (J != nrow(logTheta.matrix)) {
    error("Error -- need nrow(logTheta.matrix) == J")
  }
  prob = c(0.68, 0.95)
  res = matrix(nrow = J, ncol = length(prob) + 1)
  for (i in 1:length(prob)) {
    for (j in 1:J) {
      lower = quantile(logTheta.matrix[j, ], (1-prob[i])/2)
      upper = quantile(logTheta.matrix[j, ], (1+prob[i])/2)
      res[j, i] = as.numeric(logTheta.true[j] > lower && logTheta.true[j] < upper)
    }
  }
  res = cbind(logTheta.true, res)
  return(res)
}