#aggregate data from cluster machines to get coverages

output.files = list.files("odyssey", full.names=T, pattern="rda")

totalCoverage = list()
created = FALSE
for(file in output.files) {
  load(file)
  if( !created ) {
    totalCoverage = coverage
    created = TRUE
  }
  for(i in length(totalCoverage)) {
    totalCoverage[[i]] = rbind(totalCoverage[[i]], coverage[[i]])
  }
}

aggdata = list()
pdf("SampleGraphs.pdf",width=7,height=5)
for(i in 1:length(totalCoverage)) {
  dat = data.frame( totalCoverage[[i]] )
  dat$bucket = round(dat$logTheta.true, 1)
  aggdata[[i]] = aggregate(dat, by=list(key=dat$bucket),
                      FUN="mean", na.rm=TRUE)
  x = aggdata[[i]]$bucket
  y1 = aggdata[[i]]$sd1
  y2 = aggdata[[i]]$sd2
  plot( x, y1, main=sprintf("68pct coverage%d", i), xlab="log(theta)", ylab="coverage")
  lo = loess(y1~x)
  lines(x, predict(lo), col='red', lwd=2)
  plot( x, y2, main=sprintf("95pct coverage%d", i), xlab="log(theta)", ylab="coverage")
  lo = loess(y2~x)
  lines(x, predict(lo), col='red', lwd=2)
}
dev.off()