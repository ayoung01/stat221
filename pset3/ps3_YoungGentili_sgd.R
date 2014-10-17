# Copyright (c) 2014
# Panos Toulis, ptoulis@fas.harvard.edu
# Helper code for pset3 for Stat 221, Fall 2014
# rm(list=ls())
#
# EXAMPLE run:
#
# d = sample.data(1e4, 50)
# sgd.2b(d)
library(mvtnorm)

random.orthogonal <- function(p) {
  # Get an orthogonal matrix.
  B = matrix(runif(p^2), nrow=p)
  qr.Q(qr(B))
}

sample.data <- function(dim.n, dim.p,
                        model="gaussian") {
  # Samples the covariates as normal with the specific correlation

  # Create A matrix (variance of the covariates xn)
  Q = random.orthogonal(p=dim.p)
  lambdas = seq(0.01, 1, length.out=dim.p)
  A = Q %*% diag(lambdas) %*% t(Q)

  X = rmvnorm(dim.n, mean=rep(0, dim.p), sigma=A)
  theta = matrix(1, ncol=1, nrow=dim.p)
  epsilon = rnorm(dim.n, mean=0, sd=1)
  # Data generation
  y = X %*% theta  + epsilon

  return(list(Y=y, X=X, A=A, theta=theta))
}

check.data <- function(data) {
  # Do this to check the data object.
  #
  nx = nrow(data$X)
  ny = length(data$Y)
  p = ncol(data$X)
  stopifnot(nx==ny, p==length(data$theta))
  lambdas = eigen(cov(data$X))$values
  print(lambdas)
  print(mean(data$Y))
  print(var(data$Y))
  print(1 + sum(cov(data$X)))
}

plot.risk <- function(data, est) {
  # est = p x niters
  niters = ncol(est)
  est.bias.sgd = apply(est, 2, function(colum)
    log(t(colum-data$theta) %*% data$A %*% (colum-data$theta)))

  est.bias.asgd = numeric(niters)
  for (i in 2:niters) {
    colum = rowMeans(est[, 1:i])
    est.bias.asgd[i] = log(t(colum-data$theta) %*% data$A %*% (colum-data$theta))
  }
  plot(est.bias.sgd, type="l", lty=3)
  lines(est.bias.asgd, type="l", lty=3)
}

plot.risk.asgd <- function(data, est) {

}

sgd.2a <- function(data) {
  A = data$A
  n = nrow(data$X)
  theta.sgd = matrix(0, nrow=p, ncol=1)
  alpha = 0.01
  tr = sum(diag(A))
  for (i in 1:n) {
    a = alpha/(alpha/tr + i)
    xi = data$X[i, ]
    theta.old = theta.sgd[, i]
    theta.new = theta.old + a * A %*% (xi - theta.old)
    theta.sgd = cbind(theta.sgd, theta.new)
  }
  return(theta.sgd)
}

sgd.2b <- function(data, plot=T) {
  # check.data(data)
  n = nrow(data$X)
  p = ncol(data$X)
  # matrix of estimates of SGD (p x iters)
  theta.sgd = matrix(0, nrow=p, ncol=1)
  # params for the learning rate seq.
  gamma0 = 1 / (sum(seq(0.01, 1, length.out=p)))
  lambda0 = 0.01

  for(i in 1:n) {
    xi = data$X[i, ]
    theta.old = theta.sgd[, i]
    ai = gamma0 / (1 + gamma0 * lambda0 * i)
    # make computations easier.
    lpred = sum(theta.old * xi)
    theta.new = theta.old - ai * (lpred * xi - data$Y[i] * xi)
    theta.sgd = cbind(theta.sgd, theta.new)
  }
  if(plot) {
    plot.risk(data, theta.sgd)
  } else {
    return(theta.sgd)
  }
}

batch.2b <- function(data) {
  n = nrow(data$X)
  p = ncol(data$X)
  mat = matrix(0, ncol=p, nrow=p)
  res = numeric(p)

  for (t in 1:n) {
    xi = data$X[t, ]
    mat = mat + (xi)%*%t(xi)
    res = res + data$Y[t] * xi
  }
  return(solve(mat) %*% res)
}