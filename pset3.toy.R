genx2 = function(n,p,rho){
  #    generate x's multivariate normal with equal corr rho
  # Xi = b Z + Wi, and Z, Wi are independent normal.
  # Then Var(Xi) = b^2 + 1
  #  Cov(Xi, Xj) = b^2  and so cor(Xi, Xj) = b^2 / (1+b^2) = rho
  z=rnorm(n)
  if(abs(rho)<1){
    beta=sqrt(rho/(1-rho))
    x0=matrix(rnorm(n*p),ncol=p)
    A = matrix(z, nrow=n, ncol=p, byrow=F)
    x= beta * A + x0
  }
  if(abs(rho)==1){ x=matrix(z,nrow=n,ncol=p,byrow=F)}

  return(x)
}

sample.data <- function(X, snr=1) {
  # Samples the dataset according to Friedman et. al.
  #
  dim.n = nrow(X)
  dim.p = ncol(X)
  # ground truth params
  theta = ((-1)^(1:dim.p))*exp(-2*((1:dim.p)-1)/20)

  f= X %*% theta
  e = rnorm(dim.n)
  k= sqrt(var(f)/(snr*var(e)))
  y=f+k*e
  return(list(y=y, X=X, theta=theta))
}

n = 1000
p = 100
rho = 0.95

X = genx2(n, p, rho)
dataset = sample.data(X, snr=3.0)

mat = dataset$X

var.sum = 0
S = matrix(0, nrow=p, ncol=p)

for (i in 1:n) {
  xi = mat[i, ]


  S = S + xi%*%t(xi)
}
for (i in 1:p) {
  var.sum = var.sum + var(mat[,i])
}
S.avg = S/n

eigenvalue.sum = sum(eigen(S.avg, only.values=T)$values)

print(var.sum)
print(eigenvalue.sum)