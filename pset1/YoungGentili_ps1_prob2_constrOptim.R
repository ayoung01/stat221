`%+%` <- function(x,y) paste(x,y,sep="")

# optimize ll function over g
ll_g = function(par, theta_H, theta_L, X) {
  G = matrix(par, ncol=2, byrow=TRUE)
  N = dim(G)[1]
  P = length(theta_H)

  res = 0
  for (n in 1:N) {
    for (p in 1:P) {
      theta_H_p = theta_H[[p]]
      theta_L_p = theta_L[[p]]
      y = log(G[n,1] * theta_L_p[X[n,p]+1] + G[n,2] * theta_H_p[X[n,p]+1])
      res = res + y
    }
  }
  return(res)
}

# optimize ll function over theta_H
ll_H = function(par, G, theta_L, X) {
  N = dim(G)[1]
  theta_H = relist(par, sk)
  P = length(theta_L)
  res = 0
  for (n in 1:N) {
    for (p in 1:P) {
      theta_H_p = theta_H[[p]]
      theta_L_p = theta_L[[p]]
      y = log(G[n,1] * theta_L_p[X[n,p]+1] + G[n,2] * theta_H_p[X[n,p]+1])
      res = res + y
    }
  }
  return(res)
}

# optimize ll function over theta_L
ll_L = function(par, G, theta_H, X) {
  N = dim(G)[1]
  theta_L = relist(par, sk)
  P = length(theta_H)
  res = 0
  for (n in 1:N) {
    for (p in 1:P) {
      theta_H_p = theta_H[[p]]
      theta_L_p = theta_L[[p]]
      y = log(G[n,1] * theta_L_p[X[n,p]+1] + G[n,2] * theta_H_p[X[n,p]+1])
      res = res + y
    }
  }
  return(res)
}

# Find MLE for (G, Theta) using a coordinate ascent approach, maximizing
# over parameters in the folLing sequence:
# g_L,n for n=1,...,N
# Theta_L,j for j=1,...,J
# Theta_H,j for j=1,...,J
gomMLE = function(X, G0, theta0) {
  curr.value = -Inf
  theta_H = theta0$H
  theta_L = theta0$L
  counter = 1
  epsilon = 1e-3

  # create ui matrix for G
  ui.G = matrix(0, nrow=nrow(G0), ncol=2*nrow(G0))
  j = 1
  for (i in 1:nrow(G0)) {
    ui.G[i, j] = -1
    ui.G[i, j + 1] = -1
    j = j + 2
  }
  ci.G = rep(-1, nrow(G0))

  # create ui matrix for Theta
  n.levels = length(theta0$L)
  n.theta = length(unlist(theta0$L))
  len.levels = sapply(theta0$L, length)
  ui.theta = matrix(0, nrow=n.levels, ncol=n.theta)

  j = 1
  i = 1
  for (len.level in len.levels) {
    for (k in 1:len.level) {
      ui.theta[i, j] = -1
      j = j + 1
    }
    i = i + 1
  }
  ci.theta = rep(-1, n.levels)


  repeat  {
    print('Iteration: '%+%counter)
    # Need to flatten all parameters before passing to optim
    print('Optimizing over G...')
    res = constrOptim(theta=as.vector(t(G0)), f=ll_g, grad=NULL, ui=ui.G, ci=ci.G,
                control=list(fnscale=-1,trace=1), X=X, theta_H=theta_H, theta_L=theta_L)
    G0 = matrix(res$par, ncol=2, byrow=TRUE)
    print(res$value)

    print('Optimizing over theta_L...')
    res = constrOptim(theta=unlist(theta_L), f=ll_L, grad=NULL, ui=ui.theta, ci=ci.theta,
                control=list(fnscale=-1,trace=1), X=X, G=G0, theta_H=theta_H)
    theta_L = relist(res$par, sk)
    print(res$value)

    print('Optimizing over theta_H...')
    res = constrOptim(theta=unlist(theta_H), f=ll_H, grad=NULL, ui=ui.theta, ci=ci.theta,
                control=list(fnscale=-1,trace=1), X=X, G=G0, theta_L=theta_L)
    theta_H = relist(res$par, sk)

    # Stop coordinate ascent if our log-likelihood function stops increasing by a factor of > epsilon
    if (res$value - curr.value < epsilon) {
      print('And we are done!')
      print(c(curr.value, res$value))
      break
    }
    curr.value = res$value
    print("Current value: "%+%curr.value)
    counter = counter + 1
  }
  return(list(G.hat=G0, theta.hat=list(H=theta_H, L=theta_L), maxlik=curr.value))
}

data1985 = read.delim("data1985_area2.csv", stringsAsFactors=FALSE)
X = as.matrix(data1985[,2:length(data1985)])
load('theta0list.Rdata')

# Correct gleba column to be zero-indexed
X[,'gleba'] = X[,'gleba'] - 1

# G0 = matrix(data = 1/2, nrow = nrow(X), ncol = 2)
# We can't start with the initial point in the feasibility region
# https://stat.ethz.ch/pipermail/r-devel/2010-June/057730.html
G0 = matrix(data = 0.499, nrow = nrow(X), ncol = 2)

# Reformat theta0List into list of list of vector
theta0 = list()
theta0$H = lapply(theta0List, function(x){
  return(x$high)
})
theta0$L = lapply(theta0List, function(x){
  return(x$low)
})

# Manually replace some initial probabilities so we don't start with zeroes
theta0$H[[6]] = c(.7,.28,.01,.01)
theta0$H[[7]] = c(.7,.27,.01,.01,.01)
theta0$H[[46]] = c(.01,.35,.35,.05,.08,.06,.09,.01)
theta0$H[[49]] = c(.44,.01,.11,.43,.01)
theta0$L[[7]] = c(.2,.39,.39,.01,.01)

# Decrease first element of each probability by .001 so we don't start on the boundary
theta0$H = lapply(theta0$H, function(x){
  x[[1]] = x[[1]] - 0.001
  return(x)
})
theta0$L = lapply(theta0$L, function(x){
  x[[1]] = x[[1]] - 0.001
  return(x)
})

# skeleton structure of relisting
sk = theta0$L

gomMLE.res = gomMLE(X, G0, theta0)