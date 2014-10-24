source('ps3_YoungGentili_2_functions.R')

args <- as.numeric(commandArgs(trailingOnly = TRUE))
if(length(args) != 2) {
  stop("Not correct no. of args")
}
n = args[1]
alpha = args[2]

res = run.sgd.2e(alpha=alpha, n=n)
write.table(res, file=sprintf("2e/sgd/%.0f_%.0f.txt", n, alpha),
            row.names=FALSE, quote=FALSE)

res = run.sgd.2cd(alpha=alpha, n=n, implicit=T)
write.table(res, file=sprintf("2e/implicit/%.0f_%.0f.txt", n, alpha),
            row.names=FALSE, quote=FALSE)