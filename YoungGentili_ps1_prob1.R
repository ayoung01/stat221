library(mvtnorm)

createCovMat <- function(alpha, beta, dim) {
  mat <- matrix( data = -1*beta, nrow = dim, ncol = dim)
  for(i in 1:dim) {
    mat[i,i] <- alpha
  }
  return( mat )
}

dlogisticnorm <- function(u,mu,alpha,beta) {
  cov <- createCovMat( alpha, beta, length(mu))
  u.trans <- log(u/(1-sum(u)))
  print( dmvnorm( u.trans, mean = mu, sigma = cov ) )
  res <- 1/(prod(u)*(1-sum(u)))*dmvnorm( u.trans, mean = mu, sigma = cov )
  return( res )
}

loglikelihood <- function(par, data) {
  mu <- par[0:(length(par)-2)]
  print(mu)
  alpha <- par[length(par)-1]
  beta <- par[length(par)]
  return( sum(apply( data, 1, function (x) dlogisticnorm(x, mu, alpha, beta))) )
}

logisticnorm.mle <- function(U) {
  #mu.init <- numeric(length(U))
  mu.init <- c(0,0,0)
  alpha.init <- 1
  beta.init <- 1
  par.init <- c(mu.init, alpha.init, beta.init)
  result <- optim( par.init, loglikelihood, data = U )
  return( list( "mu.hat"=result$par[1], "alpha.hat"=result$par[2], "beta.hat"=result$par[3]) )
}

data <- read.table("dataLogisticNorm3D.txt",header=TRUE)
x <- logisticnorm.mle( data )
