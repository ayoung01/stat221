setwd('/Users/ayoung/Dropbox/Fall 2014/STAT 221/psets')
for (dir in c('dat3','dat4','dat4_w','dat5','dat5_w')) {
  browser()
  output.files = list.files(dir, full.names=T, pattern="rda")
  for(file in output.files) {
    load(file)
    df.list = lapply(ls(), function(df.name) {
      eval(parse(text=df.name))
    })
    df.combined = do.call(rbind, df.list)
    save(df.combined, file=sprintf("YoungGentili_par%s.dat", dir))
    rm(list=ls())
  }
}


