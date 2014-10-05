source('../poissonLogN_MCMC.R')
source('../YoungGentili_ps2_functions.R')

args <- as.numeric(commandArgs(trailingOnly = TRUE))
if(length(args) != 2) {
  stop("Not correct no. of args")
}
M = args[1] # number of Ys to draw
job.id = args[2]

J = 1000
w = scan('../weights.txt')
N = 2

params = list(c(mu=1.6,sd=0.7),c(mu=2.5,sd=1.3),c(mu=5.2,sd=1.3),c(mu=4.9,sd=1.6))

# First, you will draw Theta from the given generative distribution
thetas = lapply(params, function(param){
  mu = param[1]
  sd = param[2]
  DrawTheta(J, mu, sd)
})

# create list of 4 to store coverage for each set of parameters
coverage = as.list(rep(NA, 4))
for (m in 1:M) {
  print('Rep: '%+%m%+%', job.id: '%+%job.id)
  # Then you will repeatedly draw Y | Theta for a set of simulations
  Ys = lapply(thetas, function(theta){
    simYgivenTheta(theta, w, N)
  })

  mcmc = lapply(Ys, function(Y){
    poisson.logn.mcmc(Y, w)
  })

  for (i in 1:length(params)) {
    coverage[[i]] = rbind(coverage[[i]], CalculateCoverage(thetas[[i]], mcmc[[i]]$logTheta))
  }
}

coverage = lapply(coverage, function(x){
  x = x[-1,-ncol(x)] # Remove first row of NAs
  x = cbind(x, rep(w, M)) # Add info about weights
  colnames(x) = c('logTheta.true', 'sd1','sd2','w') # rename columns
  return(x)
})
names(coverage) = as.character(params) # add parameter info
save(coverage, file=sprintf("odyssey/coverage_%d.rda", job.id))

