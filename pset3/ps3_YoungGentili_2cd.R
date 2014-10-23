source('ps3_YoungGentili_2_functions.R')

args <- as.numeric(commandArgs(trailingOnly = TRUE))
if(length(args) != 3) {
  stop("Not correct no. of args")
}
alpha = args[1]
rep.no = args[2]
n = args[3]

res = run.sgd.2cd(alpha=alpha, n=n)
write.table(res, file=sprintf("2cd/sgd/%.2f_%.0f.txt", alpha, rep.no),
            row.names=FALSE, quote=FALSE)

res = run.sgd.2cd(alpha=alpha, n=n, asgd=T)
write.table(res, file=sprintf("2cd/asgd/%.2f_%.0f.txt", alpha, rep.no),
            row.names=FALSE, quote=FALSE)

res = run.sgd.2cd(alpha=alpha, n=n, implicit=T)
write.table(res, file=sprintf("2cd/implicit/%.2f_%.0f.txt", alpha, rep.no),
            row.names=FALSE, quote=FALSE)