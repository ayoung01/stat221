source('ps3_YoungGentili_2_functions.R')

args <- as.numeric(commandArgs(trailingOnly = TRUE))
if(length(args) != 2) {
  stop("Not correct no. of args")
}
n = args[1]
alpha = args[2]

n=100
alpha=10

res.sgd = c('sgd', run.sgd.2cd(nreps=100, alpha=alpha, n=n, p=100))
res.asgd = c('asgd', run.sgd.2cd(nreps=100, alpha=alpha, n=n, p=100))
res.implicit = c('implicit', run.sgd.2cd(nreps=100, alpha=alpha, n=n, p=100))
res = rbind(res.sgd, res.asgd, res.implicit)
colnames(res) = c('method', 'n', 'alpha', 'bias', 'variance')
write.table(res, file=sprintf("2cd/%.0f_%.0f.txt", n, alpha),
            row.names=FALSE, quote=FALSE)



