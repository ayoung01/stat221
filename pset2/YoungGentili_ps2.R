# theta and w are of length J
# N is a scalar
# should return matrix Y
simYgivenTheta = function(theta, w, n) {
  J = length(theta)
  Y = matrix(nrow = J, ncol = n)
  for (j in 1:J) {
    Y[j,] = rpois(n, w[j]*theta[j])
  }
  return(Y)
}