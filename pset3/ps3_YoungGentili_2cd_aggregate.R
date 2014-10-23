library(ggplot2)

cols = c('method', 'n', 'alpha', 'bias', 'variance')
d = matrix(nrow=0, ncol=length(cols))
colnames(d) <- cols
rownames(d) = NULL
dir = '2cd'

for (dir in dir) {
  output.files = list.files(dir, full.names=T)
  for(file in output.files) {
    tmp = read.table(file, header=TRUE, quote="\"", stringsAsFactors=FALSE)
    d = rbind(d, tmp)
  }
}

sgd = d[d$method=='sgd',]
asgd = d[d$method=='sgd',]
implicit = d[d$method=='implicit',]

p <- ggplot(sgd, aes(x=n, y=bias, group=alpha, colour=group))
p + geom_line()
