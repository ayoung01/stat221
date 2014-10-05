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

types = c("logTheta.true")
for(t in types) {
  for(i in 1:length(totalCoverage)) {
    dat = data.frame( totalCoverage[[i]] )
    dat$bucket = round(dat[[t]], 1)
    aggdata = aggregate(dat, by=list(key=dat$bucket),
                        FUN="mean", na.rm=TRUE)
    x = aggdata$bucket
    y1 = aggdata$sd1
    y2 = aggdata$sd2

    png(sprintf("YoungGentili_ps2_task3_plot%d_%s.png", i, t))
    par(mfrow=c(2, 1))
    plot( x, y1, main=sprintf("68pct coverage%d", i), xlab=t, ylab="coverage")
    lo = loess(y1~x)
    lines(x, predict(lo), col='red', lwd=2)
    plot( x, y2, main=sprintf("95pct coverage%d", i), xlab=t, ylab="coverage")
    lo = loess(y2~x)
    lines(x, predict(lo), col='red', lwd=2)
    dev.off()
    keeps = c("key","sd1","sd2")
    aggdata = subset(aggdata, select = keeps)
    save(aggdata, file=sprintf("YoungGentili_ps2_task3_par%d_%s.dat", i, t))
  }
}
