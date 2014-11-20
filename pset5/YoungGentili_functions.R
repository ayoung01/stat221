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
getQ = function(theta, A, y, c=2, w=11, debug=0, verbose=0) {
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
  Q = -w/2*(sum(log(diag(sigma11))) + tr(solve(sigma11) %*% R)) - s/2
  if (is.infinite(Q)){
    browser()
  }
  if (verbose>2) {
    print('Q: '%+%Q)
  }
  return(Q)
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

# sigmas = list of sigmas indexed by t
# etas = list of etas indexed by t
# Y = list of y_t indexed by t
g = function(eta_t, etas, sigmas_cond, t, Y, h=5, debug=FALSE) {
  if (debug) {
    browser()
  }
  logprior = dmvnorm(eta_t, mean=etas[[t-1]], sigma=sigmas_cond[[t]], log=TRUE)
  etas[[t]] = eta_t
  lambdas = lapply(etas, function(x) {
    exp(x[2:17])
  })
  sigmas = lapply(etas, function(eta){
    phi = exp(eta[1])
    lambda = exp(eta[2:17])
    return(getSigma(phi, lambda))
  })
  loglik = 0
  for (i in seq(t-h, t+h)) {
    loglik = loglik + dmvnorm(Y[[i]], A%*%lambdas[[t]], A%*%sigmas[[t]]%*%t(A), log=TRUE)
#     loglik = loglik -1/2*(det(A%*%sigmas[[i]]%*%t(A))+
#                             t(Y[[i]]-A%*%lambdas[[i]])%*%solve(A%*%sigmas[[i]]%*%t(A))%*%(Y[[i]]-A%*%lambdas[[i]]))
  }
  return(logprior + loglik)
}

getQ.prior = function(theta_t, thetas, t, sigmas_cond, A, y, w=11, debug=0) {
  if (debug==2) {
    browser()
  }
  thetas[[t]] = theta_t
  logprior = dmvnorm(theta_t, mean=thetas[[t-1]], sigmas_cond[[t]], log=TRUE)
  if (is.infinite(logprior)) {
    logprior = 0 # hackish thing
    print('oops infinite logprior')
  }
  if (debug==2) {
    print('logprior'%+%logprior)
  }
  phi = theta_t[1]
  lambdas = theta_t[2:length(theta_t)]
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
  return(logprior-w/2*(log(det(sigma11)) + tr(solve(sigma11) %*% R)) - s/2)
}

## Run EM for every window
locally_iid_EM <- function(data, c=2, A, I=16, verbose=1, debug=0, reltol=1e-8) {
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
    y = data[(t-(w-1)/2):t+(w-1)/2, ]

    epsilon = 1e-4
    counter = 0
    Q = Inf

    # repeat while EM has not converged
    repeat {
      # E-step
      Q = getQ(theta, A, y, w=w, debug=debug, verbose=verbose)
      if (debug == 1) {
        browser()
      }

      # M-step
      fit <- optim(par=theta,
                   fn=getQ,
                   method="L-BFGS-B",
                   lower=c(1e-4, rep(0.1, I)), # phi > 0 and lambda > 0
#                    upper=c(1e4, rep(1e6, I)),
                   control=list(fnscale=-1, reltol=reltol),
                   A=A, y=y, w=w, c=c, verbose=verbose)


      if (verbose > 1) {
        print('Q at iteration '%+%counter%+%': '%+%fit$value)
        print('old Q: '%+%Q)
        if (verbose > 2) {
          print('old theta')
          print(theta)
          print('new theta')
          print(fit$par)
        }
      }
      theta = fit$par
      counter = counter + 1

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