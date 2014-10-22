cols = c("method", "n", "p", "rho", "rep", "time", "mse")
timings = matrix(nrow=0, ncol=length(cols))
colnames(timings) <- cols
rownames(timings) = NULL

for (dir in 'out') {
  output.files = list.files(dir, full.names=T)
  for(file in output.files) {
    d = read.table(file, header=TRUE, quote="\"", stringsAsFactors=FALSE)
    timings = rbind(timings, d)
  }
}
View(timings)

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
table1 = with(avg, subset(avg, n==1000 & p==100))
table2 = with(avg, subset(avg, n==5000 & p==100))
table3 = with(avg, subset(avg, n==100 & p==1000))
table4 = with(avg, subset(avg, n==100 & p==5000))
table5 = with(avg, subset(avg, n==100 & p==20000))
table6 = with(avg, subset(avg, n==100 & p==50000))
