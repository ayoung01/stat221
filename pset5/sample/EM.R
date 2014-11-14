# Panos Toulis, ptoulis@fas.harvard.edu
# Section on EM
# Fall 2014, Stat221, Harvard

#
# Problems 
#   mm.median -> shows the Majorization-Maximization algorithm.
#                 Calculates the median.
# 
rm(list=ls())

mm.median <- function() {
  y = round(rnorm(10, mean=20, sd=5), 6)
  print("Data")
  print(y)
  print(sort(y))
  print(sprintf("Median is %.2f", median(y)))
  theta.old = 0
  inloop = T
  theta <- c()
  fobj <- function(par) sum(abs(y-par))
  
  while(inloop) {
    theta <- c(theta, theta.old)
    w = 1 / abs(y - theta.old)
    theta.next = sum(y * w) / sum(w)
    # print(sprintf("Next iteration theta=%.2f", theta.next))
    # print(round(w, 3))
    print(sprintf("f=%.3f", sum(abs(y-theta.old))))
    if(abs(theta.old-theta.next) < 1e-5) {
      inloop <- F
    }
    theta.old = theta.next
  }
  
  fvals = sapply(theta, fobj)
  plot(fvals, type="l")
  print(sprintf("MM algo=%.3f value=%.7f", theta.old, fobj(theta.old)))
  abline(h=fobj(median(y)), col="red")
  print(sprintf("median=%.3f  Value from median %.7f", median(y), fobj(median(y))))
  fit <- optim(par=runif(1, min=min(y), max=max(y)),
               fn=function(x) sum(abs(x-y)), method="L-BFGS-B",
               lower=min(y), upper=max(y), control=list(fnscale=1))
  print(sprintf("optim=%.3f Value from optim=%.7f", fit$par, fit$value))
  abline(h=fobj(fit$par), col="green")
}

EM.rao <- function() {
  y = c(125, 18, 20, 34)
  get.piVector <- function(pi) {
    c(1/2 + 1/4 * pi, 1/4 * (1-pi), 1/4 * (1-pi), 1/4 * pi)  
  }
  get.complete.piVector <- function(pi) {
    c(1/2, 1/4 * pi, 1/4 * (1-pi), 1/4 * (1-pi), 1/4 * pi)
  }
  
  pi.old = 0.0
  inloop <- T
  # We will consider 5 classes 
  while(inloop) {
    # E-step: Given pi.old, impute the sufficient statistic
    multinom.probs = get.complete.piVector(pi.old)
    x1 = y[1] * multinom.probs[1] / (multinom.probs[1] + multinom.probs[2])
    x2 = y[1] * multinom.probs[2] / (multinom.probs[1] + multinom.probs[2])
    
    x.complete = c(x1, x2, y[2], y[3], y[4])
    # M-step: Find the MLE assuming complete data.
    fobj <- function(x) {
      z = dmultinom(x.complete, size=197, prob=get.complete.piVector(x), log=T)
      if(z == -Inf) return(-1e10)
      return(z)
    }
    fit <- optim(par=pi.old, 
                 fn=fobj,
                 method="L-BFGS-B", 
                 lower=c(1e-5), upper=c(1), control=list(fnscale=-1))
    pi.next = fit$par
    # pi.next = (x2 + y[4]) / (x2 + y[4] + y[2] + y[3])
    if(abs(pi.old - pi.next) < 1e-6) {
      inloop <- F
    }
    pi.old = pi.next
    print(sprintf("New iteration pi=%.4f", pi.old))
  }
}

EM.mixture <- function(nsamples=100) {
  # Compute the MLE values from a mixture of Poissons.
  lambdas = c(1, 5, 12)
  # Sample memberships.
  a = sample(1:3, size=nsamples, prob=c(3, 1, 1), replace=T)
  # Sample data
  y = rep(0, nsamples)
  for(i in 1:length(lambdas)) {
    y <- y + as.numeric(a==i) * rpois(nsamples, lambda=lambdas[i])
  }  
  par(mfrow=c(1, 2))
  plot(a, y, main="Data from mixture of Poissons.")
  inloop = T
  lam.old <- rep(3, length(lambdas))
  membership = c()
  
  while(inloop) {
    # E-step: Impute the missing indicators.
    a.mis = sapply(1:nsamples, function(i) {
      d = dpois(y[i], lambda=lam.old)
      sample(1:3, size=1, prob=d)
    })
    # M-step: Treat imputed indicators as correct
    y1 = y[which(a.mis==1)]
    y2 = y[which(a.mis==2)]
    y3 = y[which(a.mis==3)]
    ll <- function(lam) {
      z = sum(dpois(y1, lambda=lam[1], log=T)) + 
        sum(dpois(y2, lambda=lam[2], log=T)) + 
        sum(dpois(y3, lambda=lam[3], log=T))
      if(z == -Inf) return(-1e20)
      return(z)
    }
    
    fit <- optim(par=lam.old, fn=ll, method="L-BFGS-B", 
                 lower=rep(1e-6, 3), upper=rep(max(y)^3, 3),
                 control=list(fnscale=-1))
    # print(fit$par)
    lam.next = fit$par
    if(all(abs(lam.next - lam.old) < 1e-5)) {
      inloop <- F
    }
    lam.old = lam.next
    membership <- a.mis
  }
  print("EM estimate for lambas")
  # Need to reorder
  z = match(1:length(lambdas), order(lam.old))
  lam.old = sort(lam.old)
  membership <- z[membership]
  print(head(membership))
  plot(jitter(a, amount = .1), jitter(membership, amount=.1), pch=".", cex=2)
  print(lam.old)
  print(sprintf("Hit rate %.1f%%", 100 * sum(diag(as.matrix(table(membership, a)))) / nsamples))
}