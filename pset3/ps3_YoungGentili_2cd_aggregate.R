library(ggplot2)

cols = c('method', 'n', 'alpha', 'bias', 'variance')
d = matrix(nrow=0, ncol=length(cols))
colnames(d) <- cols
rownames(d) = NULL
dir = '2cd'

for (dir in dir) {
  output.files = list.files(dir, full.names=T)
  for(file in output.files) {
    tmp = read.table(file, header=TRUE, quote="\"", stringsAsFactors=FALSE)
    d = rbind(d, tmp)
  }
}

sgd = d[d$method=='sgd',]
asgd = d[d$method=='sgd',]
implicit = d[d$method=='implicit',]

# compute variance
empirical.var = (1 / lr(alpha, n)) * cov(thetas)
var.trace = sum(diag(empirical.var))

# compute bias
true.theta = matrix(1, ncol=1, nrow=p)
theta.avg = rowMeans(thetas)
bias = sqrt(t(theta.avg-true.theta) %*% (theta.avg-true.theta))

p <- ggplot(sgd, aes(x=n, y=bias, group=alpha, colour=group))
p + geom_line()
