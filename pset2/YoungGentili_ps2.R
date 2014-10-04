source('poissonLogN_MCMC.R')
source('rASL.R')

# Task 2 ------------------------------------------------------------------

simYgivenTheta = function(theta, w, n) {
  J = length(theta)
  Y = matrix(nrow = J, ncol = n)
  for (j in 1:J) {
    Y[j,] = rpois(n, w[j]*theta[j])
  }
  return(Y)
}

# Task 3 ------------------------------------------------------------------

J = 1000
w = rep(1, J)
N = 2

params = list(c(mu=1.6,sd=0.7),c(mu=2.5,sd=1.3),c(mu=5.2,sd=1.3),c(mu=4.9,sd=1.6))

# First, you will draw Theta from the given generative distribution
DrawTheta = function(J, mu, sd) {
  return(exp(rnorm(J, mu, sd)))
}
thetas = lapply(params, function(param){
  mu = param[1]
  sd = param[2]
  drawTheta(J, mu, sd)
})

# Then you will repeatedly draw Y | Theta for a set of simulations
Ys = lapply(thetas, function(theta){
  simYgivenTheta(theta, w, N)
})

system.time(mcmc <- lapply(Ys, function(Y){
  poisson.logn.mcmc(Y, w)
}))
# user  system elapsed
# 38.732   3.729  42.477
# It takes roughly 10.5 seconds to run MCMC

CalculateCoverage = function(Theta.true, logTheta.matrix) {
  # Computes the coverage using 68% and 95% posterior intervals for each
  # Theta_j using MCMC output
  #
  # Args:
  #   logTheta.true: True value of logTheta
  #   logTheta.matrix: J x ndraw logTheta matrix output by poisson.logn.mcmc
  #
  # Returns:
  #   J x 2 vector \in {0, 1} indicating coverage for 68% and 95% intervals
  logTheta.true = log(Theta.true)
  J = length(logTheta.true)
  if (J != nrow(logTheta.matrix)) {
    error("Error -- need nrow(logTheta.matrix) == J")
  }
  prob = c(0.68, 0.95)
  res = matrix(nrow = J, ncol = length(prob))
  for (i in 1:length(prob)) {
    for (j in 1:J) {
      lower = quantile(logTheta.matrix[j, ], (1-prob[i])/2)
      upper = quantile(logTheta.matrix[j, ], (1+prob[i])/2)
      res[j, i] = as.numeric(logTheta.true[j] > lower && logTheta.true[j] < upper)
    }
  }
  return(res)
}

coverage = list()
for (i in 1:length(mcmc)) {
  coverage[[i]] = CalculateCoverage(thetas[[i]], mcmc[[i]]$logTheta)
}

