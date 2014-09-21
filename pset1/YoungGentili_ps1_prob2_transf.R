`%+%` <- function(x,y) paste(x,y,sep="")

# transform G to simplex
# x: N x 2 matrix
x_to_u_g = function(x) {
  sums = apply(x, 1, function(x){
    1 + exp(x[1]) + exp(x[2])
  })
  for (i in 1:nrow(x)) {
    x[i,1] = exp(x[i,1])/sums[i]
    x[i,2] = exp(x[i,2])/sums[i]
  }
  return(x)
}
# transform simplex to G
# x: N x 2 matrix
u_to_x_g = function(x) {
  sums = apply(x, 1, function(x){
    1 - (x[1] + x[2])
  })
  for (i in 1:nrow(x)) {
    x[i,1] = log(x[i,1]/sums[i])
    x[i,2] = log(x[i,2]/sums[i])
  }
  return(x)
}

# transform theta to simplex
# x: list
x_to_u_theta = function(x) {
  sums = 1 + sapply(x, function(y){sum(exp(y))})
  res = list()
  x = lapply(seq_along(x), function(x.idx) {
    return(sapply(x[[x.idx]], function(z){
      exp(z)/sums[x.idx]
    }))
  })
  return(x)
}
# transform simplex to theta
# x: list
u_to_x_theta = function(x) {
  sums = 1 - sapply(x, sum)
  lapply(seq_along(x), function(x.idx){
    sapply(x[[x.idx]], function(z){
      log(z/sums[x.idx])
    })
  })
}

# optimize ll function over g
ll_g = function(par, theta_H, theta_L, X) {
  G = x_to_u_g(matrix(par, ncol=2, byrow=TRUE))
  theta_H = x_to_u_theta(theta_H)
  theta_L = x_to_u_theta(theta_L)
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
  G = x_to_u_g(G)
  theta_H = x_to_u_theta(relist(par, sk))
  theta_L = x_to_u_theta(theta_L)
  N = dim(G)[1]
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
  G = x_to_u_g(G)
  theta_H = x_to_u_theta(theta_H)
  theta_L = x_to_u_theta(relist(par, sk))
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

# Find MLE for (G, Theta) using a coordinate ascent approach, maximizing
# over parameters in the folLing sequence:
# g_L,n for n=1,...,N
# Theta_L,j for j=1,...,P
# Theta_H,j for j=1,...,P
gomMLE = function(X, G0, theta0) {
  curr.value = -Inf
  NUM_ITER = 1
  EPSILON = 0.1
  MAX_ITER = 20
  theta_H = theta0$H
  theta_L = theta0$L

  # transform initial G and theta into the real space to run optim

  G0 = u_to_x_g(G0)
  theta_H = u_to_x_theta(theta0$H)
  theta_L = u_to_x_theta(theta0$L)

  repeat  {
    trace = 3
    lower = 0
    upper = 100
    print('ITERATION: '%+%NUM_ITER)
    # Need to flatten all parameters before passing to optim
    print('Optimizing over G...')
    res = optim(par=as.vector(t(G0)), fn=ll_g, method='L-BFGS-B', lower=lower, upper=upper,
                control=list(fnscale=-1,trace=trace), X=X, theta_H=theta_H, theta_L=theta_L)
    G0 = matrix(res$par, ncol=2, byrow=TRUE)
    print(res$value)

    print('Optimizing over theta_L...')
    res = optim(par=unlist(theta_L), fn=ll_L, method='L-BFGS-B', lower=lower, upper=upper,
                control=list(fnscale=-1,trace=trace), X=X, G=G0, theta_H=theta_H)
    theta_L = relist(res$par, sk)
    print(res$value)

    print('Optimizing over theta_H...')
    res = optim(par=unlist(theta_H), fn=ll_H, method='L-BFGS-B', lower=lower, upper=upper,
                control=list(fnscale=-1,trace=trace), X=X, G=G0, theta_L=theta_L)
    theta_H = relist(res$par, sk)

    # Stop coordinate ascent if our log-likelihood function stops increasing by a factor of > epsilon
    # or if we exceed a maximum number of iterations
    if (res$value - curr.value < EPSILON || NUM_ITER > MAX_ITER) {
      print('And we are done!')
      print(c(curr.value, res$value))
      break
    }
    curr.value = res$value
    print("Current value: "%+%curr.value)
    NUM_ITER = NUM_ITER + 1
  }
  G0 = x_to_u_g(G0)
  theta_H = x_to_u_theta(theta_H)
  theta_L = x_to_u_theta(theta_L)
  mles = list(G.hat=G0, theta.hat=list(H=theta_H, L=theta_L), maxlik=curr.value)
  save(mles, files='mles.RData')
  return(mles)
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

# subtract 0.001 from first element of each vector to satisfy condition sum over ui < 1
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