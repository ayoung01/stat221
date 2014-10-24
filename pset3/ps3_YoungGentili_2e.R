source('ps3_YoungGentili_2_functions.R')

args <- as.numeric(commandArgs(trailingOnly = TRUE))
if(length(args) != 2) {
  stop("Not correct no. of args")
}
n = args[1]
alpha = args[2]

res = run.sgd.2e(alpha=alpha, n=n)
write(res, file=sprintf("2e/sgd/%.0f_%.0f.txt", n, alpha))

res = run.sgd.2e(alpha=alpha, n=n, implicit=T)
write(res, file=sprintf("2e/implicit/%.0f_%.0f.txt", n, alpha))