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

get.m_t = function(phi, lambdas, y_t, A) {
  sigma11 = getSigma(phi, lambdas)
  sigma12 = sigma11 %*% t(A)
  sigma21 = A %*% sigma11
  sigma22 = A %*% sigma11 %*% t(A)
  mu1 = lambdas
  mu2 = A %*% lambdas
  a = y_t
  m_t = mu1 + sigma12 %*% solve(sigma22) %*% (a - mu2)
  return(m_t)
}

get.m = function(phi, lambdas, Y, A) {
  I = length(lambdas)
  m = matrix(0, nrow=nrow(Y), ncol=I)
  for (i in 1:nrow(Y)) {
    m[i, ] = get.m_t(phi, lambdas, Y[i, ], A)
  }
  return(m)
}

get.R = function(phi, lambdas, A) {
  sigma11 = getSigma(phi, lambdas)
  sigma12 = sigma11 %*% t(A)
  sigma21 = A %*% sigma11
  sigma22 = A %*% sigma11 %*% t(A)
  R = sigma11 - sigma12 %*% solve(sigma22) %*% sigma21
  return(R)
}

get.a = function(R, m_t) {
  a = numeric(I)
  T. = length(m_t)
  for (i in 1:length(A)) {
    a[i] = R[i, i] + 1/T. * sum(m_t^2)
  }
  return(a)
}

get.b = function(m) {
  return(colMeans(m))
}

get.J = function(phi, lambdas, Y, A, c=2, debug=FALSE) {
  if(debug) browser()
  I = length(lambdas)
  J = matrix(0, nrow=I+1, ncol=I+1)
  m = get.m(phi, lambdas, Y, A)
  b = get.b(m)
  for (i in 1:I) {
    J[i, i] = phi*c^2*lambdas[i]^(c-1) + 2*(2-c)*lambdas[i] - 2*(1-c)*b[i] # populate diagonal
    J[i, I+1] = c*lambdas[i]^c # populate last column
    J[I+1, i] = c*lambdas[i]^c # populate last column
#     J[I+1, i] = (2-c)*lambdas[i]^(1-c)-(1-c)*lambdas[i]^(-c)*b[i] # populate last row
    # TODO: figure out coefficients due to taking log
  }
  return(J)
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
    lambdas.init = rep(1, I)
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
                   method="Nelder-Mead",
                   lower=c(1e-4, rep(0.1, I)), # phi > 0 and lambda > 0
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
    thetas_t[[t]] = theta
    if (verbose > 0) {
      print(theta)
      if (verbose > 1) print(fit$value)
    }
  }
  return(thetas_t)
}