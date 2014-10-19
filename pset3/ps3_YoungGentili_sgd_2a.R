library(mvtnorm)


sample.data.2a <- function(t.size, A.size = 100) {
  #create eigenvalues as in 6.1
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


plot.risk.sgd.2a <- function(data, est) {
  A = data$A
  est.bias.sgd = apply(est, 2, function(column)
    t(column) %*% A %*% column )
  plot(est.bias.sgd, type='l', lty=3)
}

batch.2a <- function(data) {
  A = data$A
  n = nrow(data$X)
  p = ncol(data$X)
  theta.sgd = matrix(0, nrow=p, ncol=1)
  for (i in 1:n) {
    xi = data$X[i, ]
    theta.old = theta.sgd[, i]
    theta.new = (1 - 1/i)*theta.old + 1/i * xi
    theta.sgd = cbind(theta.sgd, theta.new)
  }
  return(theta.sgd)
}

asgd.2a <- function(data) {
  A = data$A
  n = nrow(data$X)
  p = ncol(data$X)
  theta.sgd = matrix(0, nrow=p, ncol=1)
  alpha = 0.01
  tr = sum(diag(A))
  for (i in 1:n) {
    a = alpha/(alpha/tr + i)
    xi = data$X[i, ]
    theta.old = theta.sgd[, i]
    inv.term = solve( diag(p) - a * A )
    theta.new = inv.term %*% (theta.old + a * A %*% xi )
    theta.sgd = cbind(theta.sgd, theta.new)
  }
  return(theta.sgd)
}

sgd.2a <- function(data) {
  A = data$A
  n = nrow(data$X)
  p = ncol(data$X)
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
