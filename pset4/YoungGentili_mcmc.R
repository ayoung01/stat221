# Albert Young & Marco Gentili
# Problem Set 4
# Fall 2014, Stat221, Harvard

library(MASS)
`%+%` <- function(x,y) paste(x,y,sep="")


impala = read.table('impala.txt', header=T, quote="\"", stringsAsFactors=F)$impala
waterbuck = read.table('waterbuck.txt', header=T, quote="\"", stringsAsFactors=F)$waterbuck

NBound = 150
thetaBound = 1

log.lik <- function(N, theta, Y) {
  # Log-likelihood of the data
  sum(dbinom(Y, size=N, prob=theta, log=T))
}

log.prior <- function(N, theta) {
  return(log(1/(N*(N-1)*theta)))
}

log.posterior <- function(a, b, Y) {
  log.lik(a, b, Y) + log.prior(a, b)
}

mle <- function(y) {
  # Compute the MLE given the data
  # Note: The model y ~ Binom(N, theta) is misspecified so there are
  #   multiple MLEs.
  #
  # Args:
  #   y = vector of values

  f <- function(par) {
    log.lik(par[1], par[2], y)
  }
  x = optim(par=c(50,0.5), fn=f, control=list(fnscale=-1),
            method="L-BFGS-B", lower=c(0, 0))
  print(sprintf("Log-likelihood=%.3f", x$value))
  return(x$par)
}


plot.posterior.density <- function(data) {
  # Plots the density of the posterior density with contour lines.
  #
  # Args:
  #   data = vector of values

  a = seq(.1, kBound, length.out=400)
  b = a
  gr = expand.grid(x=a, y=b)
  z = apply(gr, 1, function(par) log.lik(par[1], par[2], data))
  z = matrix(z, nrow=length(a), ncol=length(b))
  print("Summary of likelihood values.")
  print(summary(as.vector(z)))
  z0 = quantile(z, probs = c(0.1))
  z1 = max(z)
  ncols = 10000
  cols = topo.colors(ncols)
  i = apply(gr, 1, function(par) {
    ll = log.posterior(par[1], par[2], data)
    if(ll < z0) return(0)
    return((ll - z0) / (z1 - z0))
  })
  par(mfrow=c(1, 1))
  # hist(i)
  j = ncols * i
  j[j==0] <- 1
  plot(gr[,1], gr[,2], xlab="a", ylab="b", col=cols[j], pch=".", cex=4)
  contour(a, b, z, nlevels = 50, add=T, cex=2, col="black")
}

rgamma.trunc <- function(upper.bound, s, r) {
  # Sample from truncated gamma.
  # TODO: Find a better implementation.
  x <- upper.bound + 10
  while(x > upper.bound) {
    x = rgamma(1, shape=s, rate=r)
  }
  return(x)
}

plot.chain <- function(mcmc.chain) {
  mcmc.niters = nrow(mcmc.chain)
  burnin = 0.1 * mcmc.niters
  x = mcmc.chain[, 1]
  y = mcmc.chain[, 2]
  mcmc.chain = mcmc.chain[burnin:mcmc.niters, ]
  f = kde2d(x, y, n=100)
  image(f, xlim=c(min(x), max(x)), ylim=c(min(y), max(y)))
  contour(f$x, f$y, f$z, nlevels=5, add=T, cex=2, col='black')
  #   points(x, y, pch=21, col=rgb(0,0,0,0.01))
}

mcmc.gibbs <- function(y, mcmc.niters=1e4) {
  # Runs Gibbs on the data y
  S = sum(y)  # sufficient statistic
  n = length(y)

  mcmc.chain <- matrix(3, nrow=mcmc.niters, ncol=2)

  for(i in 2:mcmc.niters) {
    beta.last <- mcmc.chain[i-1, 2]
    # Conditionals are truncated gammas.
    alpha <- rgamma.trunc(kBound, s=S+1, r=n * beta.last)
    beta <- rgamma.trunc(kBound, s=S+1,  r=n * alpha)
    mcmc.chain[i, ] <- c(alpha, beta)
  }
  # Plot empirical density
  plot.chain(mcmc.chain)
  return(mcmc.chain)
}

mcmc.mh <- function(y, mcmc.niters=1e4, init=c(90, 0.7)) {
  mcmc.chain <- matrix(0.1, nrow=mcmc.niters, ncol=2)
  mcmc.chain[1, ] = init
  nacc <- 0
  for(i in 2:mcmc.niters) {
    # 1. Current state
    N.old = mcmc.chain[i-1, 1]
    theta.old = mcmc.chain[i-1, 2]
    # 2. Propose new state
    N.new = trunc(rnorm(1, N.old, 3))
    theta.new = rnorm(1, theta.old, 0.05)

    if (N.new < 0) {
      N.new = 1
    }
    if (theta.new < 0) {
      theta.new = 0
    }
    if (theta.new > 1) {
      theta.new = 1
    }

    # 3. Ratio
    mh.ratio = min(0, log.posterior(N.new, theta.new, y) -
                     log.posterior(N.old, theta.old, y))
    if(runif(1) < exp(mh.ratio)) {
      # Accept
      mcmc.chain[i, ] <- c(N.new, theta.new)
      nacc <- nacc + 1
    } else {
      mcmc.chain[i, ] <- c(N.old, theta.old)
    }
  }
  # Cut the burn-in period.
  plot.chain(mcmc.chain)
  print(sprintf("Acceptance ratio %.2f%%", 100 * nacc / mcmc.niters))
  return(mcmc.chain)
}

generate_1.4_plots = function(y, niter=1e4, init=c(90, 0.7), title) {
  name <- function(v1) {
    deparse(substitute(v1))
  }
  par(mfrow=c(5,2),
      oma = c(5,4,0,0) + 0.6,
      mar = c(0,0,1,1) + 0.6
  )
  for (i in 1:10) {
    mcmc.mh(y, niter, init)
  }
  mtext(title%+%' MCMC chain ('%+%niter%+%' Iterations)', side=3, outer=TRUE, line=-1)
  title(xlab='N',
        ylab='theta',
        outer=TRUE,
        line=3)
}

generate_1.4_plots(impala, init=c(50, 0.6), title='impala')
generate_1.4_plots(waterbuck, init=c(90, 0.7), title='waterbuck')