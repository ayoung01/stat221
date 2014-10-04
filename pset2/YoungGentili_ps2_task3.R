source('poissonLogN_MCMC.R')
source('YoungGentili_ps2_functions.R')

# Task 3 ------------------------------------------------------------------

J = 1000
w = rep(1, J)
N = 2

params = list(c(mu=1.6,sd=0.7),c(mu=2.5,sd=1.3),c(mu=5.2,sd=1.3),c(mu=4.9,sd=1.6))

# First, you will draw Theta from the given generative distribution
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

coverage = list()
for (i in 1:length(mcmc)) {
  coverage[[i]] = CalculateCoverage(thetas[[i]], mcmc[[i]]$logTheta)
}

