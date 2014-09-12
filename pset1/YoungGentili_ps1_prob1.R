library(mvtnorm)
`%+%` <- function(x,y) paste(x,y,sep="")

createCovMat <- function(alpha, beta, dim) {
  mat <- matrix( data = -1*beta, nrow = dim, ncol = dim)
  for(i in 1:dim) {
    mat[i,i] <- alpha
  }
  return( mat )
}

dlogisticnorm <- function(u,mu,alpha,beta) {
  cov <- createCovMat(alpha, beta, length(mu))
  u.trans <- log(u/(1-sum(u)))
  res <- 1/(prod(u)*(1-sum(u)))*dmvnorm( u.trans, mean = mu, sigma = cov )
  return( res )
}

loglikelihood <- function(par, data) {
  mu <- par[1:3]
#   print(mu)
  alpha <- par[4]
  beta <- par[5]
  return(prod(apply(data, 1, function(x) dlogisticnorm(x, mu, alpha, beta))))
}

logisticnorm.mle <- function(U) {
  #mu.init <- numeric(length(U))
  mu.init <- c(0,0,0)
  alpha.init <- 1
  beta.init <- .1
  par.init <- c(mu.init, alpha.init, beta.init)
  result <- optim(par.init, loglikelihood, data = U, control=list(fnscale=-1))
  return(result)
}

data <- read.table("dataLogisticNorm3D.txt", header=TRUE)
x <- logisticnorm.mle(data)

print('mu.mle: '%+%x$par[1:3])
print('alpha.mle: '%+%x$par[4])
print('beta.mle: '%+%x$par[5])
