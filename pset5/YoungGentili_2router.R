library(scales)
library(mvtnorm)
library(numDeriv)
library(lattice)

# 1.9 ---------------------------------------------------------------------
source('YoungGentili_functions.R')
router2 = read.csv("2router_linkcount.dat")

I2 = 64
# J x I incidence matrix
A2 = matrix(c(rep(1, 8), rep(0, 56),
              rep(0, 8), rep(1, 8), rep(0, 48),
              rep(0, 16), rep(1, 8), rep(0, 40),
              rep(0, 24), rep(1, 8), rep(0, 32),
              rep(0, 32), rep(1, 8), rep(0, 24),
              rep(0, 40), rep(1, 8), rep(0, 16),
              rep(0, 48), rep(1, 8), rep(0, 8),
              rep(0, 56), rep(1, 8),
              rep(c(1,0,0,0,0,0,0,0), 8),
              rep(c(0,1,0,0,0,0,0,0), 8),
              rep(c(0,0,1,0,0,0,0,0), 8),
              rep(c(0,0,0,1,0,0,0,0), 8),
              rep(c(0,0,0,0,1,0,0,0), 8),
              rep(c(0,0,0,0,0,1,0,0), 8),
              rep(c(0,0,0,0,0,0,1,0), 8)), ncol=64, byrow=T)

Y2 = do.call(rbind, tapply(router2$value, router2$time, function(x) x))
colnames(Y2) = router2$nme[1:16]
Y2 = Y2[, -ncol(Y2)]

thetas2_t = locally_iid_EM(data=Y2, A=A2, c=2, I=I2, debug=0, verbose=2)


df = do.call(rbind, thetas2_t)
lambdas = data.frame(df[ ,2:ncol(df)])
lambdas = stack(lambdas)

names = do.call(rbind, strsplit(router2$nme[1:8], ' '))[, 2]
names2 = character()
names2.1 = character()
names2.2 = character()
for (name1 in names) {
  for (name2 in names) {
    names2 = c(names2, paste(name1, name2, sep='->'))
    names2.1 = c(names2.1, name1)
    names2.2 = c(names2.2, name2)
  }
}

names2 = factor(names2, levels = names2)

lambdas$names2 = rep(names2, each=277)
lambdas$names2.1 = rep(names2.1, each=277)
lambdas$names2.2 = rep(names2.2, each=277)

p = ggplot(lambdas, aes(x=rep(6:282*5/60, I2), y=lambdas$values)) + geom_line() + xlab('hour of day') + ylab('bytes/sec') +
  ylim(0, 1e5) + scale_x_continuous(breaks=seq(0, 24, 4)) + ggtitle('Mean Traffic Estimates lambda_t for all OD Pairs' )
p + facet_wrap(~ names2)


