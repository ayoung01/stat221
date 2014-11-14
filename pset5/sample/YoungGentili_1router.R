library(scales)
# 1.1 ---------------------------------------------------------------------

links <- read.csv("1router_allcount.dat")
links = links[grep("src|dst", links$nme),]

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



