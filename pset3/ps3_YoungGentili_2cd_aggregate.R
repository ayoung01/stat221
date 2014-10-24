library(ggplot2)
`%+%` <- function(x,y) paste(x,y,sep="")
lr <- function(alpha, n) {
  ## learning rate
  alpha / (alpha + n)
}

dirs = c('2cd/sgd/', '2cd/asgd/', '2cd/implicit/')
alphas = c('0.01', '0.10', '0.50', '1.00', '2.00', '10.00', '20.00', '50.00', '100.00', '200.00')
reps = c('1','2','3','4','5')
methods = c('sgd', 'asgd', 'implicit')

p = 100
n = 1e4
true.theta = matrix(1, ncol=1, nrow=p)

res = list()
# average over reps
for (method in methods) {
  dir = '2cd/'%+%method%+%'/'
  for (alpha in alphas) {
    thetas = matrix(0, nrow=p, ncol=n)
    for (rep in reps) {
      file = dir%+%alpha%+%'_'%+%rep%+%'.txt'
      tmp = read.table(file, header=F, quote="\"", stringsAsFactors=FALSE)
      thetas = thetas + tmp
    }
    thetas = thetas/length(reps)

    var.trace = numeric(n)
    bias = numeric(n)

    for (i in 1:n) {
      thetas_n = thetas[, 1:i]

      # compute bias
      theta.avg = rowMeans(thetas_n)
      bias = sqrt(t(theta.avg-true.theta) %*% (theta.avg-true.theta))

      # compute variance
      empirical.var = (1 / lr(alpha, i)) * cov(thetas_n)
      var.trace[i] = sum(diag(empirical.var))
    }

    res[[method]][[alpha]][[var.trace]] = var.trace
    res[[method]][[alpha]][[bias]] = bias
  }
}

p <- ggplot(sgd, aes(x=n, y=bias, group=alpha, colour=group))
p + geom_line()
