output.files = list.files("odyssey", full.names=T, pattern="rda")
A <- c()
numTests = 4
totalCoverage = list()
created = FALSE
for(file in output.files) {
  load(file)
  if( !created ) {
    totalCoverage = coverage[[i]]
    created = TRUE
  }
  for(i in numTests) {
    totalCoverage[[i]] = rbind(totalCoverage[[i]], coverage[[i]])
  }
}

aggdata = list()
for(i in numTests) {
  x = data.frame( totalCoverage[[i]] )
  x$bucket = round(x$logTheta.true, 2)

  aggdata[[i]] = aggregate(x, by=list(key=x$bucket),
                      FUN="mean", na.rm=TRUE)
  pdf(sprintf("SampleGraph%d.pdf", i),width=7,height=5)
  plot(aggdata[[i]]$bucket, aggdata[[i]]$sd1, main="68% coverage", xlab="log(theta)", ylab="coverage")
  plot(aggdate[[i]]$bucket, aggdata[[i]]$sd2, main="95% coverage", xlab="log(theta)", ylab="coverage")
  dev.off()
}
