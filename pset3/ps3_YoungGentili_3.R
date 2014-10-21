source('ps3_YoungGentili_3_functions.R')

args <- as.numeric(commandArgs(trailingOnly = TRUE))
if(length(args) != 4) {
  stop("Not correct no. of args")
}
n = args[1]
p = args[2]
rho = args[3]
nreps = args[4]

timing = run.timing(rho, nreps, n, p)

write.table(timing,
            file=sprintf("out/timing_%.0f_%.0f_%.2f_%.0f.txt",
                         n, p, rho, nreps),
            row.names=FALSE, quote=FALSE)


