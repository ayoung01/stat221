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
getQ = function(theta, A, w=11) {
  theta[1] = phi
  theta[2:length(theta)] = lambdas
  sigma = getSigma(phi, lambdas)
  R = sigma - sigma %*% t(A) %*% solve(A%*%sigma%*%t(A)) %*% A %*% sigma
  return(-w/2*(log(det(sigma)) + tr(solve(sigma)%*%R)))
}
tr = function(x) sum(diag(x))