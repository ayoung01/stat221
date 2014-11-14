tr = function(x) sum(diag(x))
sigma2 = function(lambda, c=2) {
  return(lambda**c)
}
getSigma = function(phi, lambdas, c=2) {
  I = length(lambdas)
  sigma = matrix(0, nrow=I, ncol=I)
  for (i in 1:I) {
    sigma[i, i] = sigma2(lambdas[i], c)
  }
  return(phi * sigma)
}
getQ = function(theta, A, y, w=11, debug=0) {
  if (debug==2) {
    browser()
  }
  phi = theta[1]
  lambdas = theta[2:length(theta)]
  sigma11 = getSigma(phi, lambdas)
  sigma12 = sigma11 %*% t(A)
  sigma21 = A %*% sigma11
  sigma22 = A %*% sigma11 %*% t(A)
  mu1 = lambdas
  mu2 = A %*% lambdas
  s = 0
  for (t in 1:nrow(y)) {
    a = y[t, ]
    m_t = mu1 + sigma12 %*% solve(sigma22) %*% (a - mu2)
    s = s + t(m_t - lambdas) %*% solve(sigma11) %*% (m_t - lambdas)
  }
  R = sigma11 - sigma12 %*% solve(sigma22) %*% sigma21
  return(-w/2*(log(det(sigma11)) + tr(solve(sigma11) %*% R)) - s/2)
}

## Run EM for every window
runEM_1.4 <- function(Y, verbose=1, debug=0) {
  # debug 1: break in EM loop
  # debug 2: break in getQ
  thetas_t = list()
  for (t in 6:282) {
    if (verbose > 0) {
      print('t: '%+%t)
    }
    # initialize phi & lambdas
    ### TODO: initialize this in a smarter way
    phi.init = 1
    lambdas.init = rep(1e4, I)
    theta = c(phi.init, lambdas.init)
    # get window of 11 data points
    y = Y[(t-(w-1)/2):t+(w-1)/2, ]

    epsilon = 1e-4
    counter = 0
    # repeat while EM has not converged
    repeat {
      # E-step
      Q = getQ(theta, A, y, w=w, debug=debug)
      if (debug == 1) {
        browser()
      }

      # M-step
      fit <- optim(par=theta,
                   fn=getQ,
                   method="L-BFGS-B",
                   lower=rep(1e-4, I+1), # phi > 0 and lambda > 0
                   control=list(fnscale=-1),
                   A=A, y=y, w=w)
      theta = fit$par
      counter = counter + 1

      if (verbose > 1) print('Q at iteration '%+%counter%+%': '%+%fit$value)

      if (fit$value - Q < epsilon) {
        if (verbose) {
          print('EM exited after: '%+%counter%+%' steps')
        }
        break
      }
    }
    thetas_t$t = theta
    if (verbose > 0) {
      print(theta)
      if (verbose > 1) print(fit$value)
    }
  }
  return(thetas_t)
}