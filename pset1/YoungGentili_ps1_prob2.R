
ll = function(G, theta, X) {
  N = dim(G)[1]
  P = dim(theta)[2]
  res = 1
  for (n in 1:N) {
    for (p in 1:P) {
      theta_p = theta[[p]]
      theta_p_high = theta[[p]]$high
      theta_p_low = theta[[p]]$low
      res = res*(G[n,1] * theta_p_low[X[n,p]] + G[n,2] * theta_p_high[X[n,p]])
    }
  }
  return(res)
}

gomMLE = function(X, G0, theta0) {
  
}
# setwd('pset1')
data1985 = read.delim("data1985_area2.csv", stringsAsFactors=FALSE)
X = as.matrix(data1985[,2:length(data1985)])
load('theta0list.Rdata')

N = nrow(X)
G0 = matrix(data = 1/2, nrow = N, ncol = 2)