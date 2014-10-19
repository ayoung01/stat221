library(mvtnorm)

#creates the data as per section 6.1 of Xu paper
sample.data.2a <- function(t.size, A.size = 100) {
  #create eigenvalues [1,1,1,0.02,...]
  lambdas = numeric(A.size)
  for( i in seq(1, A.size) ) {
    if( i <= 3) {
      lambdas[i] = 1
    } else {
      lambdas[i] = 0.02
    }
  }

  A = diag(lambdas)
  X = rmvnorm(t.size, mean = rep(0, A.size), sigma = diag(A.size))
  return(list(X=X, A=A))
}

#plots the risk
plot.risk.sgd.2a <- function(data, est) {
  A = data$A
  est.bias.sgd = apply(est, 2, function(column)
    t(column) %*% A %*% column )
  plot(est.bias.sgd, type='l', log="y", lty=3)
}

plot.risk.all.2a <- function(est) {
  A = est$param.data$A
  est.risk = lapply(est[-1], function(estimate) {
               apply(estimate, 2, function(col)
                 t(col) %*% A %*% col)})
  cl = rainbow(length(est.risk))
  plot(est.risk[[1]], type='l', log="xy", col=cl[1], xlab = "training size t", ylab = "excess risk")
  for( i in seq(2,length(est.risk))) {
    lines(est.risk[[i]], col=cl[i])
  }
  legend("bottomleft", legend=names(est.risk), col=cl, pch=15)
}

get.data.all.2a <- function(t.size) {
  data = sample.data.2a(t.size)
  batch.data = batch.sgd.2a(data)
  implicit.data = implicit.sgd.2a(data)
  average.data = average.sgd.2a(data)
  sgd.data = sgd.2a(data)
  bad.average.data = bad.average.sgd.2a(data)
  return(list(param.data=data,
              batch=batch.data,
              implicit=implicit.data,
              ASGD=average.data,
              SGD=sgd.data,
              ASGD_BAD=bad.average.data))
}

batch.sgd.2a <- function(data) {
  A = data$A
  n = nrow(data$X)
  p = ncol(data$X)
  theta.sgd = matrix(runif(p), nrow=p, ncol=1)
  for (i in 1:n) {
    xi = data$X[i, ]
    theta.old = theta.sgd[, i]
    theta.new = (1 - 1/i)*theta.old + 1/i * xi
    theta.sgd = cbind(theta.sgd, theta.new)
  }
  return(theta.sgd)
}

average.sgd.2a <- function(data) {
  A = data$A
  lambda_0 = min(eigen(A)$values)
  n = nrow(data$X)
  p = ncol(data$X)
  theta.sgd = matrix(runif(p), nrow=p, ncol=1)
  avg.theta.sgd = matrix(runif(p), nrow=p, ncol=1)
  tr = sum(diag(A))
  for (i in 1:n) {
    a = 1/(1 + lambda_0*i)^(2/3)
    xi = data$X[i, ]
    theta.old = theta.sgd[, i]
    theta.new = theta.old + a * A %*% (xi - theta.old)
    avg.theta.new = (1 - 1/i)*avg.theta.sgd[,i] + 1/i * theta.new
    theta.sgd = cbind(theta.sgd, theta.new)
    avg.theta.sgd = cbind(avg.theta.sgd, avg.theta.new)
  }
  return(avg.theta.sgd)
}

bad.average.sgd.2a <- function(data) {
  A = data$A
  n = nrow(data$X)
  p = ncol(data$X)
  theta.sgd = matrix(runif(p), nrow=p, ncol=1)
  avg.theta.sgd = matrix(runif(p), nrow=p, ncol=1)
  tr = sum(diag(A))
  for (i in 1:n) {
    a = 1/(1 + i)^(1/2)
    xi = data$X[i, ]
    theta.old = theta.sgd[, i]
    theta.new = theta.old + a * A %*% (xi - theta.old)
    avg.theta.new = (1 - 1/i)*avg.theta.sgd[,i] + 1/i * theta.new
    theta.sgd = cbind(theta.sgd, theta.new)
    avg.theta.sgd = cbind(avg.theta.sgd, avg.theta.new)
  }
  return(avg.theta.sgd)
}

implicit.sgd.2a <- function(data) {
  A = data$A
  n = nrow(data$X)
  p = ncol(data$X)
  theta.sgd = matrix(runif(p), nrow=p, ncol=1)
  alpha = 1
  tr = sum(diag(A))
  for (i in 1:n) {
    #a = alpha/(alpha/tr + i)
    a = 1/(1+0.02*i)^(2/3)
    xi = data$X[i, ]
    theta.old = theta.sgd[, i]
    inv.term = solve( diag(p) + a * A )
    theta.new = inv.term %*% (theta.old + a * A %*% xi )
    theta.sgd = cbind(theta.sgd, theta.new)
  }
  return(theta.sgd)
}

sgd.2a <- function(data) {
  A = data$A
  lambda_0 = min(eigen(A)$values)
  n = nrow(data$X)
  p = ncol(data$X)
  theta.sgd = matrix(runif(p), nrow=p, ncol=1)
  tr = sum(diag(A))
  for (i in 1:n) {
    a = 1/(1 + lambda_0*i)
    xi = data$X[i, ]
    theta.old = theta.sgd[, i]
    theta.new = theta.old + a * A %*% (xi - theta.old)
    theta.sgd = cbind(theta.sgd, theta.new)
  }
  return(theta.sgd)
}
