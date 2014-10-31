library(ggplot2)
library(data.table)
`%+%` <- function(x,y) paste(x,y,sep="")

alphas = c('60','70','80','90','100','110','120','130','140','150')
ns = as.character(seq(1e3,1e4,1e3))
methods = c('sgd', 'implicit')

df = data.frame()
for (method in methods) {
  dir = '2e/'%+%method%+%'/'
  for (alpha in alphas) {
    for (n in ns) {
      file = dir%+%n%+%'_'%+%alpha%+%'.txt'
      print(file)
      tmp = scan(file)
      print(tmp)
      row = cbind(method, as.numeric(alpha), as.numeric(n), as.numeric(tmp))
      df = rbind(df, row)
    }
  }
}
colnames(df) = c('method', 'alpha', 'n', 'dist')
df$dist=as.numeric(as.character(df$dist))
df.sgd = df[df$method=='sgd',]
df.imp = df[df$method=='implicit',]

p <- ggplot(df.sgd, aes(x=n, y=dist, group=alpha, colour=factor(alpha)))
p + geom_line() + ggtitle('||Empirical - Theoretical Var|| vs. N (Standard SGD)')

p <- ggplot(df.imp, aes(x=n, y=dist, group=alpha, colour=factor(alpha)))
p + geom_line() + ggtitle('||Empirical - Theoretical Var|| vs. N (Implicit SGD)')
