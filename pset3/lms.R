# Stochastic gradient descent. 
# Panos Toulis, ptoulis@fas.harvard.edu
# Example on a simple linear model.
# n = #samples, p=#covariates
#
# SYNOPSIS.
#   d = sample.data(1e5, 1e2)
#   test.lm(d)
#   test.sgd(d, alpha=0.5)
#   test.sgd(d, alpha=0.1)
#
rm(list=ls())

# gen* functions to generate the data taken from Friedman and Tibshirani (JOSS)
# 
genjerry = function(x, snr){
  # generate data according to Friedman's setup
  n=nrow(x)
  p=ncol(x)
  b=((-1)^(1:p))*exp(-2*((1:p)-1)/20)
  # b=sample(c(-0.8, -0.45, 0.45, 0.9, 1.4), size=p, replace=T)
  # ((-1)^(1:p))*(1:p)^{-0.65}#exp(-2*((1:p)-1)/20)
  f=x%*%b
  e=rnorm(n)
  k=sqrt(var(f)/(snr*var(e)))
  y=f+k*e
  return(list(y=y, beta=b))
}

genx2 = function(n,p,rho){
  #    generate x's multivariate normal with equal corr rho
  # Xi = b Z + Wi, and Z, Wi are independent normal.
  # Then Var(Xi) = b^2 + 1
  #  Cov(Xi, Xj) = b^2  and so cor(Xi, Xj) = b^2 / (1+b^2) = rho
  z=rnorm(n)
  if(abs(rho)<1){
    beta=sqrt(rho/(1-rho))
    x0=matrix(rnorm(n*p),ncol=p)
    A = matrix(z, nrow=n, ncol=p, byrow=F)
    x= beta * A + x0
  }
  if(abs(rho)==1){ x=matrix(z,nrow=n,ncol=p,byrow=F)}
  
  return(x)
}

sample.data <- function(dim.n, dim.p, rho=0.0, 
                        model="gaussian") {
  # Samples the covariates as normal with the specific correlation
  X = genx2(dim.n, dim.p, rho)
  d = genjerry(X, 3)
  beta = d$beta
  Y = d$y
  if(model=="logistic") {
    Y = as.vector(1* (runif(dim.n) < 1/(1 + exp(-Y))))
  }
  return(list(Y=Y, X=X, beta=beta))
}

sgd <- function(data, alpha=0.5) {
  # Runs SGD and returns the estimate.
  n = nrow(data$X)
  p = ncol(data$X)
  theta.sgd = matrix(0, nrow=p, ncol=1)
  I = diag(p)
  for(i in 1:n) {
    xi = matrix(data$X[i, ], ncol=1)
    yi = data$Y[i]
    yi.pred = sum(xi * theta.sgd)
    ai = alpha / (alpha + i)
    theta.sgd = theta.sgd + ai * (yi - yi.pred) * xi
  }
  return(theta.sgd)
}

test.lm <- function(data) {
  # Runs simple lm()
  tick = system.time({ fit <- lm(data$Y ~ data$X + 0)})
  print(sprintf("Time elapsed for lm() %.3f secs.", tick[["elapsed"]]))
  b.hat = as.numeric(coefficients(fit))
  plot(d$beta, d$beta, type="l", col="red", main="Fitted vs. true params")
  points(b.hat, d$beta, pch=6)
}

test.sgd <- function(data, alpha=1, add.plot=F) {
  # Plot/test sgd
  tick = system.time({ theta.sgd <- sgd(data,alpha=alpha) })
  print(sprintf("Time elapsed for sgd() %.3f secs.", tick[["elapsed"]]))
  if(add.plot) {
    points(theta.sgd, data$beta, pch=6)
  } else {
    plot(data$beta, data$beta, type="l", col="red")
    points(theta.sgd, data$beta, pch=13)
  }  
}



