library(ggplot2)
library(data.table)
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
  res[[method]] = list()
  print(method)
  dir = '2cd/'%+%method%+%'/'
  for (alpha in alphas) {
    res[[method]][[alpha]] = list()
    print(alpha)
    thetas = matrix(0, nrow=p, ncol=n)
    for (rep in reps) {
      file = dir%+%alpha%+%'_'%+%rep%+%'.txt'
      tmp = fread(file)
      thetas = thetas + tmp
    }
    thetas = as.matrix(thetas/length(reps))

    div = 1000
    var.trace = numeric(n/div)
    bias = numeric(n/div)


    for (i in 1:n) {
      if (i %% div == 0){
        print(i)
        thetas_n = thetas[, 1:i]
        # compute bias
        theta.avg = rowMeans(thetas_n)
        bias[i/div] = sqrt(t(theta.avg-true.theta) %*% (theta.avg-true.theta))

        # compute variance
        empirical.var = (1 / lr(as.numeric(alpha), i)) * cov(thetas_n)
        var.trace[i/div] = sum(diag(empirical.var))
      }
    }
    res[[method]][[alpha]][['var.trace']] = var.trace
    res[[method]][[alpha]][['bias']] = bias
  }
}

dfs = list()
# plot graphs
idx = seq(1000, 10000, 1000)
for (method in methods) {
  df = data.frame()
  for (alpha in alphas) {
    df = rbind(df, cbind(alpha=as.numeric(alpha), n=idx, log.var.trace=log(res[[method]][[alpha]][['var.trace']])))
  }
  dfs[[method]] = df
}

p <- ggplot(dfs$sgd, aes(x=n, y=log.var.trace, group=alpha, colour=factor(alpha)))
p + geom_line() + ggtitle('Log Empirical Variance Trace vs. N (Standard SGD)')


p <- ggplot(dfs$asgd, aes(x=n, y=log.var.trace, group=alpha, colour=factor(alpha)))
p + geom_line() + ggtitle('Log Empirical Variance Trace vs. N (ASGD)')


p <- ggplot(dfs$implicit, aes(x=n, y=log.var.trace, group=alpha, colour=factor(alpha)))
p + geom_line() + ggtitle('Log Empirical Variance Trace vs. N (Implicit SGD)')


for (method in methods) {
  df = data.frame()
  for (alpha in alphas) {
    df = rbind(df, cbind(alpha=as.numeric(alpha), n=idx, log.bias=log(res[[method]][[alpha]][['bias']])))
  }
  dfs[[method]] = df
}

p <- ggplot(dfs$sgd, aes(x=n, y=log.bias, group=alpha, colour=factor(alpha)))
p + geom_line() + ggtitle('Log(Bias) vs. N (Standard SGD)')


p <- ggplot(dfs$asgd, aes(x=n, y=log.bias, group=alpha, colour=factor(alpha)))
p + geom_line() + ggtitle('Log(Bias) vs. N (ASGD)')


p <- ggplot(dfs$implicit, aes(x=n, y=log.bias, group=alpha, colour=factor(alpha)))
p + geom_line() + ggtitle('Log(Bias) vs. N (Implicit SGD)')



