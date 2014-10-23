library(tables)

cols = c("method", "n", "p", "rho", "rep", "time", "mse")
timings = matrix(nrow=0, ncol=length(cols))
colnames(timings) <- cols
rownames(timings) = NULL
dir = 'out'

for (dir in dir) {
  output.files = list.files(dir, full.names=T)
  for(file in output.files) {
    d = read.table(file, header=TRUE, quote="\"", stringsAsFactors=FALSE)
    timings = rbind(timings, d)
  }
}

ns = c(100, 1000, 5000)
rhos = c(0.0, 0.1, 0.2, 0.5, 0.9, 0.95)
# average over reps
avg = with(timings, aggregate(timings, list(method, n, rho, p), mean))
drops <- c('Group.1','Group.2','Group.3','Group.4', 'rep')
avg = avg[,!names(avg)%in%drops]
avg$method[avg$method==0] = 'glmnet (naive)'
avg$method[avg$method==1] = 'glmnet (cov)'
avg$method[avg$method==2] = 'sgd (standard)'
avg$method[avg$method==3] = 'sgd (implicit)'
avg$method = as.factor(avg$method)
avg$rho = as.factor(avg$rho)

View(timings)

# 3b+c --------------------------------------------------------------------


# table1 = with(avg, subset(avg, n==1000 & p==100))
# table2 = with(avg, subset(avg, n==5000 & p==100))
# table3 = with(avg, subset(avg, n==100 & p==1000))
# table4 = with(avg, subset(avg, n==100 & p==5000))
# table5 = with(avg, subset(avg, n==100 & p==20000))
# table6 = with(avg, subset(avg, n==100 & p==50000))
#
# latex(tabular( method ~ rho*time*identity, data=table1))
# latex(tabular( method ~ rho*time*identity, data=table2))
# latex(tabular( method ~ rho*time*identity, data=table3))
# latex(tabular( method ~ rho*time*identity, data=table4))
# latex(tabular( method ~ rho*time*identity, data=table5))
# latex(tabular( method ~ rho*time*identity, data=table6))


# 3d ----------------------------------------------------------------------

table1 = with(avg, subset(avg, n==1e3 & p==50000))
table2 = with(avg, subset(avg, n==1e4 & p==5000))
table3 = with(avg, subset(avg, n==1e5 & p==1000))
table4 = with(avg, subset(avg, n==1e6 & p==100))

latex(tabular(  rho*(time*identity + mse*identity) ~ method, data=table1))
latex(tabular(  rho*(time*identity + mse*identity) ~ method, data=table2))
latex(tabular(  rho*(time*identity + mse*identity) ~ method, data=table3))
latex(tabular(  rho*(time*identity + mse*identity) ~ method, data=table4))
