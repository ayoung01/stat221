source('ps3_YoungGentili_2_functions.R')

args <- as.numeric(commandArgs(trailingOnly = TRUE))
if(length(args) != 2) {
  stop("Not correct no. of args")
}
n = args[1]
alpha = args[2]

res.sgd = run.sgd.2cd(nreps=100, alpha=alpha, n=n, p=100)
res.asgd = run.sgd.2cd(nreps=100, alpha=alpha, n=n, p=100)
res.implicit = run.sgd.2cd(nreps=100, alpha=alpha, n=n, p=100)
res = rbind(res.sgd, res.asgd, res.implicit)
colnames(res) = c('bias', 'variance')
write.table(res, file=sprintf("2cd/%.0f_%.0f.txt", n, alpha),
            row.names=FALSE, quote=FALSE)



