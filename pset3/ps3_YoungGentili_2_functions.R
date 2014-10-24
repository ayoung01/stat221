library(mvtnorm)

# Helper functions --------------------------------------------------------

random.orthogonal <- function(p) {
  # Get an orthogonal matrix.
  B = matrix(runif(p^2), nrow=p)
  qr.Q(qr(B))
}
generate.A <- function(p) {
  # Create A matrix (variance of the covariates xn)
  Q = random.orthogonal(p)
  lambdas = seq(0.01, 1, length.out=p)
  A = Q %*% diag(lambdas) %*% t(Q)
  return(A)
}
sample.data <- function(dim.n, A,
                        model="gaussian") {
  # Samples the dataset. Returns a list with (Y, X, A ,true theta)
  dim.p = nrow(A)
  # This call will make the appropriate checks on A.
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
get.data.all.2b <- function(t.size, alpha) {
  A = generate.A(100)
  data = sample.data(t.size, A)
  batch.data = batch.sgd.2b(data)
  implicit.data = implicit.sgd.2b(data, alpha, plot=F)
  average.data = average.sgd.2b(data, alpha, plot=F)
  sgd.data = sgd.2b(data, alpha, plot=F)
  return(list(param.data=data,
              batch=batch.data,
              implicit=implicit.data,
              ASGD=average.data,
              SGD=sgd.data))
}
plot.risk.all <- function(est) {
  A = est$param.data$A
  theta = est$param.data$theta
  estimates = est[-1]
  est.risk = lapply(estimates, function(estimate) {
    apply(estimate$theta, 2, function(col)
      log(t(col-theta) %*% A %*% (col-theta)))})
  cl = rainbow(length(est.risk))
  plot(estimates[[1]]$indices, est.risk[[1]], type='l', log="xy", col=cl[1], xlab = "training size t", ylab = "excess risk")
  for( i in seq(2,length(est.risk))) {
    lines(estimates[[i]]$indices, est.risk[[i]], col=cl[i])
  }
  legend("bottomleft", legend=names(est.risk), col=cl, pch=15)
}
plot.risk <- function(data, est) {

  # est = p x niters
  #est.bias = apply(est$theta, 2, function(colum)
  #  log(t(colum-data$theta) %*% data$A %*% (colum-data$theta)))

  est.bias = apply(est$theta, 2, function(colum)
    t(colum-data$theta) %*% data$A %*% (colum-data$theta))

  print(sprintf("Risk of first estimate = %.3f Risk of last estimate %.3f",
                head(est.bias, 1), tail(est.bias, 1)))
  plot(est$indices, est.bias, type="l", log='xy', lty=3)
}
lr <- function(alpha, n) {
  ## learning rate
  alpha / (alpha + n)
}
sqrt.norm <- function(X) {
  sqrt(mean(X^2)  )
}

# 2b ----------------------------------------------------------------------


sgd.2b <- function(data, alpha, plot=T) {
  n = nrow(data$X)
  p = ncol(data$X)
  # matrix of estimates of SGD (p x iters)
  theta.sgd = matrix(0, nrow=p, ncol=1)
  # params for the learning rate seq.
  # gamma0 = 1 / (sum(seq(0.01, 1, length.out=p)))
  trace = sum(diag(data$A))  # NOTE: data snooping.

  for(i in 1:n) {
    xi = data$X[i, ]
    theta.old = theta.sgd[, i]
    # NOTE: This is assuming implicit. Use rate
    #  an = alpha / (alpha * trace(A) + n)  if standard.
    ai = alpha / (alpha * trace + i)

    lpred = sum(theta.old * xi)
    yi = data$Y[i]
    theta.new = (theta.old - ai * lpred * xi) + ai * yi * xi
    theta.sgd = cbind(theta.sgd, theta.new)
  }

  if(plot) {
    print("Last estimate")
    print(theta.sgd[, ncol(theta.sgd)])
    plot.risk(data, theta.sgd)
  } else {
    return(list(theta=theta.sgd, indices=1:(n+1)))
  }
}
average.sgd.2b <- function(data, alpha, plot=T) {
  # check.data(data)
  # Implements implicit
  n = nrow(data$X)
  p = ncol(data$X)
  # matrix of estimates of SGD (p x iters)
  theta.sgd = matrix(0, nrow=p, ncol=1)
  avg.theta.sgd = matrix(0, nrow=p, ncol=1)
  # params for the learning rate seq.
  # gamma0 = 1 / (sum(seq(0.01, 1, length.out=p)))
  trace = sum(diag(data$A))  # NOTE: data snooping.

  for(i in 1:n) {
    xi = data$X[i, ]
    theta.old = theta.sgd[, i]
    ai = alpha / (alpha * trace + i)

    # make computations easier.
    lpred = sum(theta.old * xi)
    yi = data$Y[i]
    # Standard SGD
    # Uncomment this line for standard SGD.
    theta.new = (theta.old - ai * lpred * xi) + ai * yi * xi
    avg.theta.new = (1 - 1/i)*avg.theta.sgd[,i] + 1/i * theta.new
    theta.sgd = cbind(theta.sgd, theta.new)
    avg.theta.sgd = cbind(avg.theta.sgd, avg.theta.new)
  }

  if(plot) {
    print("Last estimate")
    print(avg.theta.sgd[, ncol(avg.theta.sgd)])
    plot.risk(data, avg.theta.sgd)
  } else {
    return(list(theta=avg.theta.sgd, indices=1:(n+1)))
  }
}
implicit.sgd.2b <- function(data, alpha, plot=T) {
  # check.data(data)
  # Implements implicit
  n = nrow(data$X)
  p = ncol(data$X)
  # matrix of estimates of SGD (p x iters)
  theta.sgd = matrix(0, nrow=p, ncol=1)
  # params for the learning rate seq.
  # gamma0 = 1 / (sum(seq(0.01, 1, length.out=p)))
  trace = sum(diag(data$A))  # NOTE: data snooping.
  I = diag(p)

  for(i in 1:n) {
    xi = data$X[i, ]
    theta.old = theta.sgd[, i]
    # NOTE: This is assuming implicit. Use rate
    #  an = alpha / (alpha * trace(A) + n)  if standard.
    ai = lr(alpha, i)

    # make computations easier.
    xi.norm = sum(xi^2)
    lpred = sum(theta.old * xi)
    fi = 1 / (1 + ai * sum(xi^2))
    yi = data$Y[i]
    # Implicit SGD
    theta.new = (theta.old  - ai * fi * lpred * xi) +
      (ai * yi * xi - ai^2 * fi * yi * xi.norm * xi)
    # Standard SGD
    # Uncomment this line for standard SGD.
    # theta.new = (theta.old - ai * lpred * xi) + ai * yi * xi

    theta.sgd = cbind(theta.sgd, theta.new)
  }

  if(plot) {
    print("Last estimate")
    print(theta.sgd[, ncol(theta.sgd)])
    plot.risk(data, theta.sgd)
  } else {
    return(list(theta=theta.sgd, indices=1:(n+1)))
  }
}


batch.sgd.2b <- function(data) {
  increment = 100
  n = nrow(data$X)
  p = ncol(data$X)
  theta.sgd = matrix(0, nrow=p, ncol=1)
  indices = numeric(n/increment + 1)
  #for (x in 1:n) {
  for (x in seq(1,n,increment)) {
    indices[x/increment + 1] = x
    mat = matrix(0, ncol=p, nrow=p)
    res = numeric(p)
    if(x <= p) {
      theta.sgd = cbind(theta.sgd, matrix(0, nrow=p, ncol=1))
    } else {
      for (t in 1:x) {
        xi = data$X[t, ]
        mat = mat + (xi)%*%t(xi)
        res = res + data$Y[t] * xi
      }
      theta.sgd = cbind(theta.sgd, (solve(mat) %*% res))
    }
  }
  return( list(theta=theta.sgd, indices=indices ) )
}

# 2cd ---------------------------------------------------------------------


run.sgd.2cd <- function(alpha, n=1e4, p=100, asgd=F, implicit=F, verbose=T, A=NULL) {
  if (is.null(A)) {
    A = generate.A(p)
  }
  data = sample.data(n, A)

  trace = sum(diag(data$A))
  theta.sgd = matrix(0, nrow=p, ncol=n)
  for (i in 1:(n-1)) {
    theta.old = theta.sgd[, i]
    xi = data$X[i, ]
    ai = alpha / (alpha * trace + i)
    lpred = sum(theta.old * xi)
    yi = data$Y[i]
    theta.new = theta.old + ai * (yi - lpred) * xi
    if (asgd) { # ASGD
      theta.new = (1 - 1/i)*theta.old + 1/i * ai * (yi - lpred) * xi
    }
    if (implicit) { # Implicit
      ai = lr(alpha, i)
      xi.norm = sum(xi^2)
      fi = 1 / (1 + ai * sum(xi^2))
      theta.new = (theta.old  - ai * fi * lpred * xi) +
        (ai * yi * xi - ai^2 * fi * yi * xi.norm * xi)
    }
    theta.sgd[, i+1] = theta.new
  }
  return(theta.sgd)
}

# 2e ----------------------------------------------------------------------

run.sgd.2e <- function(alpha, n, implicit=F, asgd=F,
                         p=100, verbose=F) {
  # Calculates || Empirical variance - Theoretical ||
  #
  A = generate.A(p)

  # Compute theoretical variance
  I = diag(p)
  Sigma.theoretical <- alpha * solve(2 * alpha * A - I) %*% A
  stopifnot(all(eigen(Sigma.theoretical)$values > 0))

  thetas = run.sgd.2cd(n=n, p=p, alpha=alpha, A=A, asgd=asgd, implicit=implicit)

  # Calculate empirical variance
  empirical.var = (1 / lr(alpha, n)) * cov(thetas)
  dist = sqrt.norm(empirical.var - Sigma.theoretical)
  return(dist)
}
