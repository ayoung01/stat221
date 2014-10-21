# Copyright (c) 2013
# Panos Toulis, ptoulis@fas.harvard.edu
#
# Using the simulation setup as in glmnet JoSS paper(Friedman, Hastie, Tibshirani)
# http://www.jstatsoft.org/v33/i01/
#
# EXAMPLE run:
#  run.glmnet(1e4, 10)
#
library(mvtnorm)
library(glmnet)

# genjerry, genx2 are functions taken from the above paper.
# These functions generate the simulation data.
# NOTE: use function sample.data(..) instead.
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

sample.data <- function(dim.n, dim.p, rho=0.0, snr=1) {
  # Samples the dataset according to Friedman et. al.
  #
  # 1. sample covariates
  X = genx2(dim.n, dim.p, rho)
  # 2. ground truth params.
  theta = ((-1)^(1:dim.p))*exp(-2*((1:dim.p)-1)/20)

  f= X %*% theta
  e = rnorm(dim.n)
  k= sqrt(var(f)/(snr*var(e)))
  y=f+k*e
  return(list(y=y, X=X, theta=theta))
}

dist <- function(x, y) {
  if(length(x) != length(y))
    stop("MSE should compare vectors of same length")
  sqrt(mean((x-y)^2))
}

# 3c. Implement the estimation method for the same problem using standard SGD
#     and create a new Table 1 for this method. Compare with the results from
#     glmnet and comment on your findings.

sgd <- function(data, implicit=F) {
  n = nrow(data$X)
  p = ncol(data$X)
  # params for the learning rate seq.
  gamma0 = 1e-6 # this just works for some reason
  lambda0 = 0.01
  theta.old = rep(0, p)

  for(i in 1:n) {
    xi = data$X[i, ]
    yi = data$y[i]
    ai = gamma0 / (1 + gamma0 * lambda0 * i)
    lpred = sum(theta.old * xi)
    if (implicit) {
      # TODO: set ai for implicit
      xi.norm = sum(xi^2)
      fi = 1 / (1 + ai * sum(xi^2))
      theta.new = (theta.old  - ai * fi * lpred * xi) +
        (ai * yi * xi - ai^2 * fi * yi * xi.norm * xi)
    }
    # make computations easier.
    else {
      # regular sgd
      theta.new = theta.old + ai * (yi - lpred) * xi
    }
    theta.old = theta.new
  }
  return(theta.new)
}

# Wrapper function to run the experiment on local machine
run.timings <- function(rho.values = c(0.0, 0.1, 0.2, 0.5, 0.9, 0.95),
                    nreps = 1,
                    ns = c(1e4, 1e5, 1e6),
                    ps = c(1e1, 1e2, 1e3),
                    glmnet = T, sgd = T, sgd.implicit = T) {
  ## Runs glmnet() for various param values.
  ##
  niters = 0
  cols = c("glmnet", "n", "p", "rho", "rep", "time", "mse")
  timings = matrix(nrow=0, ncol=length(cols))
  colnames(timings) <- cols
  rownames(timings) = NULL
  total.iters = nreps * length(rho.values) * length(ns) * length(ps) * (glmnet + sgd)
  pb = txtProgressBar(style=3)

  for (n in ns) {
    for (p in ps) {
        for (rho in rho.values) {
          timings = rbind(timings, run.timing(rho, nreps, n, p, glmnet, sgd, sgd.implicit))
          niters = niters + glmnet + sgd
          setTxtProgressBar(pb, niters/total.iters)
        }
      }

  }
  return(timings)
}

# Will submit this to Odyssey
run.timing <- function(rho, nreps, n, p, glmnet=T, sgd=T, sgd.implicit=T) {
  ## Runs glmnet() for various param values.
  ## method: 0 = glmnet, 1 = sgd, 2 = implicit sgd
  cols = c("method", "n", "p", "rho", "rep", "time", "mse")
  timings = matrix(nrow=0, ncol=length(cols))
  colnames(timings) <- cols
  rownames(timings) = NULL

  for (i in 1:nreps) {
    # 1. (For every repetition) Sample the dataset
    dataset = sample.data(dim.n=n, dim.p=p, rho=rho, snr=3.0)
    true.theta = dataset$theta
    x = dataset$X
    y = dataset$y
    stopifnot(nrow(x) == n, ncol(x) == p)
    # 1b. Define metrics:
    #    = time for the method to finish
    #   mse = Distance (e.g. RMSE) of the estimates to the ground truth.
    #         (q1, q2, q3) representing the quartiles (since glmnet returns grid of estimates)
    #         Implicit has (x, x, x) i.e., the same value in all places.

    # 2. Run the glmnet method and tabulate timing
    if (glmnet) {
      glmnet.time = system.time({ glmnet.fit = glmnet(x, y,alpha=1, standardize=FALSE, type.gaussian="naive")})[1]
      glmnet.mse = median(apply(glmnet.fit$beta, 2, function(est) dist(est, true.theta)))
      timings = rbind(timings, c(0, n, p, rho, i,
                                 glmnet.time,
                                 glmnet.mse))
    }
    # 3. Run the sgd method and tabulate timing
    if (sgd) {
      sgd.time = system.time({ sgd.est = sgd(dataset)})[1]
      sgd.mse = dist(sgd.est, true.theta)
      timings = rbind(timings, c(1, n, p, rho, i,
                                 sgd.time,
                                 sgd.mse))
    }
    # 4. Run the implicit sgd method and tabulate timing
    if (sgd.implicit) {
      sgd.implicit.time = system.time({ sgd.est = sgd(dataset, implicit=T)})[1]
      sgd.implicit.mse = dist(sgd.est, true.theta)
      timings = rbind(timings, c(2, n, p, rho, i,
                                 sgd.implicit.time,
                                 sgd.implicit.mse))
    }
  }
  return(timings)
}
#
# timings = run.timings(ns=1e2, nrep=2, ps=1e4, rho.values=.95)
# View(timings)
