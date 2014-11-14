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
getQ = function(theta, A, y, w=11) {
  theta[1] = phi
  theta[2:length(theta)] = lambdas
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
    s = s + m_t
  }
  R = sigma11 - sigma12 %*% solve(sigma22) %*% sigma21
  return(-w/2*(log(det(sigma)) + tr(solve(sigma11) %*% R)) - s/2)
}