# Albert Young & Marco Gentili
# Problem Set 4
# Fall 2014, Stat221, Harvard

library(MASS)

impala = read.table('impala.txt', header=T, quote="\"", stringsAsFactors=F)$impala
waterbuck = read.table('waterbuck.txt', header=T, quote="\"", stringsAsFactors=F)$waterbuck

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
  mcmc.chain = mcmc.chain[burnin:mcmc.niters, ]
  f = kde2d(x=mcmc.chain[, 1], y=mcmc.chain[, 2], n=100)
  image(f, xlim=c(0, kBound), ylim=c(0, kBound))
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

mcmc.mh <- function(y, mcmc.niters=1e4) {
  # Complete with MH.
  S = sum(y)
  n = length(y)
  mcmc.chain <- matrix(0.1, nrow=mcmc.niters, ncol=2)
  nacc <- 0
  for(i in 2:mcmc.niters) {
    mu <- rgamma.trunc(kBound**2, s=S+1, r=n)
    # 1. Current state
    alpha.old = mcmc.chain[i-1, 1]
    beta.old = mcmc.chain[i-1, 2]
    # 2. Propose new state
    #   Respect symmetry in (a,b)
    alpha.new = runif(1, min=mu/kBound, max=kBound)
    beta.new = mu / alpha.new
    if(runif(1) < 0.5) {
      beta.new = alpha.new
      alpha.new = mu / beta.new
    }
    # 3. Ratio
    mh.ratio = min(0, log.posterior(alpha.new, beta.new, y) -
                      log.posterior(alpha.old, beta.old, y))
    if(runif(1) < exp(mh.ratio)) {
      # Accept
      mcmc.chain[i, ] <- c(alpha.new, beta.new)
      nacc <- nacc + 1
    } else {
      mcmc.chain[i, ] <- c(alpha.old, beta.old)
    }
  }
  # Cut the burnin period.
  print(sprintf("Acceptance ratio %.2f%%", 100 * nacc / mcmc.niters))
  plot.chain(mcmc.chain)
  return(mcmc.chain)
}
