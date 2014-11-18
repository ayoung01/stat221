library(scales)
library(mvtnorm)
source('YoungGentili_functions.R')

router1 = read.csv("1router_allcount.dat")
router2 = read.csv("2router_linkcount.dat")
links = router1[grep("src|dst", router1$nme),]

# 1.1 ---------------------------------------------------------------------

graph_1.1 <- function () {
  tmp = do.call(rbind, strsplit(links$nme, ' '))
  links$dir = tmp[, 1]
  links$client = tmp[, 2]
  links$client = factor(links$client, levels=c('corp','local','switch','fddi')) # control order of panels

  for (i in 1:len(links$time)) {
    links$time[i] = strftime(strptime(links$time[i], format="(%m/%d/%y %H:%M:%S)"), format="%H:%M:%S") # convert to time string
  }

  # http://stackoverflow.com/questions/14938278/plotting-timeseries-in-ggplot2
  p = ggplot(links, aes(x=strptime(time, "%H:%M:%S"), y=value, colour=as.factor(dir))) + geom_line() +
              scale_x_datetime(limits=c(as.POSIXct("00:02:43",format="%H:%M:%S"),as.POSIXct("23:52:43",format="%H:%M:%S")),
                               breaks=date_breaks("2 hour"), labels = date_format("%H:%M")) +
              xlab("hour of day") + ylab("bytes/sec") + ggtitle('Link Measurements for the Router1 Subnetwork')
  p + facet_grid(client ~ .)
}
# graph_1.1()

# 1.2 ---------------------------------------------------------------------

graph_1.2 <- function () {
  # Router 1
  time.1130 = links[1065:1152,]
  mean.1130 = tapply(time.1130$value, time.1130$nme, function(x) {log10(mean(x))})
  var.1130 = tapply(time.1130$value, time.1130$nme, function(x) {log10(var(x))})

  time.1530 = links[1449:1536,]
  mean.1530 = tapply(time.1530$value, time.1530$nme, function(x) {log10(mean(x))})
  var.1530 = tapply(time.1530$value, time.1530$nme, function(x) {log10(var(x))})

  df = data.frame(x = c(mean.1130, mean.1530), y = c(var.1130, var.1530), panel = c(rep("time 11:30", 8), rep("time 15:30", 8)))
  p = ggplot(df, aes(x, y)) + geom_smooth(method="lm",data=df,colour="gray",se=FALSE, fullrange=TRUE)
  p + geom_point() + facet_grid(. ~ panel) + xlab('log10(mean)') + ylab('log10(var)') +
    ggtitle('Local Variances Versus Local Means for Router1 Link Measurements')

  # Router 2
  time.1130 = router2[2129:2304,]
  mean.1130 = tapply(time.1130$value, time.1130$nme, function(x) {log10(mean(x))})
  var.1130 = tapply(time.1130$value, time.1130$nme, function(x) {log10(var(x))})

  time.1530 = links[2897:3072,]
  mean.1530 = tapply(time.1530$value, time.1530$nme, function(x) {log10(mean(x))})
  var.1530 = tapply(time.1530$value, time.1530$nme, function(x) {log10(var(x))})

  df = data.frame(x = c(mean.1130, mean.1530), y = c(var.1130, var.1530), panel = c(rep("time 11:30", 8), rep("time 15:30", 8)))
  p = ggplot(df, aes(x, y)) + geom_smooth(method="lm",data=df,colour="gray",se=FALSE, fullrange=TRUE)
  p + geom_point() + facet_grid(. ~ panel) + xlab('log10(mean)') + ylab('log10(var)') +
    ggtitle('Local Variances Versus Local Means for Router2 Link Measurements')
}
# graph_1.2()

# 1.4 ---------------------------------------------------------------------

### Prepare constants

# J x I incidence matrix
A = matrix(c(rep(1, 4), rep(0, 12),
             rep(0, 4), rep(1, 4), rep(0, 8),
             rep(0, 8), rep(1, 4), rep(0, 4),
             rep(0, 12), rep(1, 4),
             rep(c(1,0,0,0), 4),
             rep(c(0,1,0,0), 4),
             rep(c(0,0,1,0), 4)), ncol=16, byrow=T)
I = 16
w = 11

# Y is a T x J matrix
Y = do.call(rbind, tapply(links$value, links$time, function(x) x))
colnames(Y) = links$nme[1:8]
Y = Y[, -ncol(Y)] # drop last column so we have linearly independent link measurements


thetas_t = runEM_1.4(Y, debug=0, verbose=2)

### not the most elegant way to plot this but it works

df = do.call(rbind, thetas_t)
lambdas = data.frame(df[ ,2:ncol(df)])
lambdas = stack(lambdas)

names = unique(router1$nme)
names = names[grep("->", names)]
names = factor(names, levels = names)

lambdas$names = rep(names, each=277)

p = ggplot(lambdas, aes(x=rep(6:282*5/60, 16), y=lambdas$values)) + geom_line() + xlab('hour of day') + ylab('bytes/sec') +
    ylim(0, 1e6) + scale_x_continuous(breaks=seq(0, 24, 4)) + ggtitle('Mean Traffic Estimates lambda_t for all OD Pairs' )
p + facet_wrap(~ names)


# 1.6 ---------------------------------------------------------------------

# sigmas = list of sigmas indexed by t
# etas = list of etas indexed by t
# Y = list of y_t indexed by t
g = function(eta_t, etas, sigmas_hat, t, Y, debug=FALSE) {
  if (debug) {
    browser()
  }
  logprior = dmvnorm(eta_t, mean=etas[[t-1]], sigma=sigmas_hat[[t-1]], log=TRUE)
  lambdas = lapply(etas, function(x) {
    exp(x[2:17])
  })
  sigmas = lapply(etas, function(eta){
    phi = exp(eta[1])
    lambda = exp(eta[2:17])
    return(getSigma(phi, lambda))
  })
  loglik = 0
  for (i in t-5:t+5) {
    loglik = loglik -1/2*(log(det(A%*%sigmas[[i]]%*%t(A)))+
                            t(log(Y[[i]])-A%*%log(lambdas[[i]]))%*%solve(A%*%sigmas[[i]]%*%t(A))%*%(log(Y[[i]])-A%*%log(lambdas[[i]])))
  }
  return(logprior + loglik)
}
# g(etas[[8]], etas, sigmas.cond, 8, Y.list, debug=TRUE)
source('YoungGentili_functions.R')
### initialize eta_0 and sigma_0
V = diag(I + 1)

etas = lapply(thetas_t, function(x) {
  if (!is.null(x)) {
    return(log(x))
  }
  else {
    return(runif(I+1))
  }
})
sigmas = list()
sigmas.cond = list()
for (i in 1:287) {
  sigmas[[i]] = 10*diag(I + 1)
  sigmas.cond[[i]] = 10*diag(I + 1)
}

Y.list = split(Y, 1:nrow(Y))

### iterate through all time points
for (t in 6:282) {
  print(t)
  sigmas.cond[[t]] = sigmas[[t-1]] + V
  fit <- optim(par=etas[[t-1]],
               fn=g,
               method="Nelder-Mead",
               control=list(fnscale=-1),
               etas=etas, sigmas_hat=sigmas.cond, t=t, Y=Y.list)
  print(fit$par)
  print(fit$value)
  etas[[t]] = fit$par
  phi = exp(etas[[t]][1])
  lambdas = exp(etas[[t]][2:length(etas[[t]])])
  J = get.J(phi=phi, lambdas=lambdas, Y=Y[(t-h):(t+h), ], A=A, c=2, debug=FALSE)
  sigmas[[t]] = solve(-solve(sigmas.cond[[t]]) + J)
}
