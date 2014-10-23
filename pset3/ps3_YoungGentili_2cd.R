source('ps3_YoungGentili_2_functions.R')

args <- as.numeric(commandArgs(trailingOnly = TRUE))
if(length(args) != 2) {
  stop("Not correct no. of args")
}
alpha = args[1]
rep.no = args[2]

res.sgd = c('sgd', rep.no, run.sgd.2cd(alpha=alpha, n=n, p=100))
res.asgd = c('asgd', rep.no, run.sgd.2cd(alpha=alpha, n=n, p=100, asgd=T))
res.implicit = c('implicit', rep.no, run.sgd.2cd(alpha=alpha, n=n, p=100, implicit=T))
res = rbind(res.sgd, res.asgd, res.implicit)
colnames(res)[1] = 'method'
colnames(res)[2] = 'rep.no'

write.table(res, file=sprintf("2cd/%.0f_%.0f.txt", alpha, rep.no),
            row.names=FALSE, quote=FALSE)



