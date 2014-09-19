`%+%` <- function(x,y) paste(x,y,sep="")

# ll = function(G, theta, X) {
#   N = dim(G)[1]
#   P = dim(theta)[2]
#   res = 1
#   for (n in 1:N) {
#     for (p in 1:P) {
#       theta_p_high = theta[[p]]$high
#       theta_p_low = theta[[p]]$low
#       res = res + log(G[n,1] * theta_p_low[X[n,p]+1] + G[n,2] * theta_p_high[X[n,p]+1])
#     }
#   }
#   return(res)
# }

# optimize ll function over g
ll_g = function(par, theta_H, theta_L, X) {
  G = matrix(par, ncol=2)
  N = dim(G)[1]
  P = length(theta_H)

  res = 1
  for (n in 1:N) {
    for (p in 1:P) {
      theta_p_high = theta_H[[p]]
      theta_p_low = theta_L[[p]]
      y = log(G[n,1] * theta_p_low[X[n,p]+1] + G[n,2] * theta_p_high[X[n,p]+1])
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
      theta_p_high = theta_H[[p]]
      theta_p_low = theta_L[[p]]
      y = log(G[n,1] * theta_p_low[X[n,p]+1] + G[n,2] * theta_p_high[X[n,p]+1])
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
  res = 1
  for (n in 1:N) {
    for (p in 1:P) {
      theta_p_high = theta_H[[p]]
      theta_p_low = theta_L[[p]]
      y = log(G[n,1] * theta_p_low[X[n,p]+1] + G[n,2] * theta_p_high[X[n,p]+1])
      res = res + y
    }
  }
  return(res)
}

# Find MLE for (G, Theta) using a coordinate ascent approach, maximizing
# over parameters in the following sequence:
# g_L,n for n=1,...,N
# Theta_L,j for j=1,...,J
# Theta_H,j for j=1,...,J
gomMLE = function(X, G0, theta0) {
  curr.value = -Inf
  theta_H = theta0$high
  theta_L = theta0$low
  counter = 1

  epsilon = 1e-3

  repeat  {
    print('Iteration: '%+%counter)
    # Need to flatten all parameters before passing to optim
    dim(G0) = NULL
    print('Optimizing over G...')
    res = optim(par=G0, fn=ll_g, lower=0, upper=1, method='L-BFGS-B',
                control=list(fnscale=-1), X=X, theta_H=theta_H, theta_L=theta_L)
    G0 = matrix(res$par, ncol=2)
    print(res$value)

    print('Optimizing over theta_L...')
    res = optim(par=unlist(theta_L), fn=ll_L, lower=0, upper=1, method='L-BFGS-B',
                control=list(fnscale=-1), X=X, G=G0, theta_H=theta_H)
    theta_H = relist(res$par, sk)
    print(res$value)

    print('Optimizing over theta_H...')
    res = optim(par=unlist(theta_H), fn=ll_H, lower=0, upper=1, method='L-BFGS-B',
                control=list(fnscale=-1), X=X, G=G0, theta_L=theta_L)
    theta_L = relist(res$par, sk)

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
  return(list(G.hat=G0, theta.hat=list(high=theta_H, low=theta_L), maxlik=curr.value))
}

data1985 = read.delim("data1985_area2.csv", stringsAsFactors=FALSE)
X = as.matrix(data1985[,2:length(data1985)])
load('theta0list.Rdata')

# Correct gleba column to be zero-indexed
X[,'gleba'] = X[,'gleba'] - 1

N = nrow(X)
G0 = matrix(data = 1/2, nrow = N, ncol = 2)

# Reformat theta0List into list of list of vector
theta0 = list()
theta0$high = lapply(theta0List, function(x){
  y = x$high
  # Replace 0.00 with 0.01
  y = replace(y, which(y==0), 0.01)
  return(y)
})
theta0$low = lapply(theta0List, function(x){
  y = x$low
  # Replace 0.00 with 0.01
  y = replace(y, which(y==0), 0.01)
  return(y)
})

# skeleton structure of relisting
sk = theta0$low

# create ui matrix for constrOptim
n.levels = length(theta0$high)
n.theta = length(unlist(theta0$high))
len.levels = sapply(theta0$high, length)
ui = matrix(0, nrow=n.levels, ncol=n.theta)

col.idx = 1
row.idx = 1
for (len.level in len.levels) {
  for (i in 1:len.level) {
    ui[row.idx, col.idx] = -1
    col.idx = col.idx + 1
  }
  row.idx = row.idx + 1
}
ci = rep(-1, n.levels)

gomMLE.res = gomMLE(X, G0, theta0)