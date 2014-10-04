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
thetas = lapply(params, function(param){
  mu = param[1]
  sd = param[2]
  exp(rnorm(J, mu, sd))
})

# Then you will repeatedly draw Y | Theta for a set of simulations
Ys = lapply(thetas, function(theta){
  simYgivenTheta(theta, w, N)
})

mcmc = lapply(Ys, function(Y){
  poisson.logn.mcmc(Y, w)
})
