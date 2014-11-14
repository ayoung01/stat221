library(scales)
router1 = read.csv("1router_allcount.dat")
router2 = read.csv("2router_linkcount.dat")

# 1.1 ---------------------------------------------------------------------

links = router1[grep("src|dst", links$nme),]

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


# 1.2 ---------------------------------------------------------------------

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
