source('ERGM.MCMC.R')

# florentine --------------------------------------------------------------
data(florentine)
g = as.matrix(flomarriage)
ergm.res = ergm(g~ edges + triangles)

n.iter = 1000
theta0.ET = as.matrix(c(-1,-1))

init.ET.G = vector('list', n.iter)
for (i in 1:n.iter) {
  init.ET.G[[i]] = g
}

ET.res = SGD.Monte.Carlo(init.ET.G, init.ET.G[[1]], theta0.ET, ERGM.ET.ss,
                         function(i) (1/(i)), n.draws=50, ERGM.ET.ss.diff, use.pkg=T, debug=F, verbose=T, model='ET')
print('Predicted theta: '%+%ET.res[[n.samples]])

plot(1:length(unlist(ET.res)), unlist(ET.res))

# karate data from http://www-personal.umich.edu/~mejn/netdata/
library(igraph)

g = read.graph('karate/karate.gml', format='gml')
plot(g)
g = get.adjacency(g, sparse=F)

n.iter = 1000
theta0.ET = as.matrix(c(-1,-1))

init.ET.G = vector('list', n.iter)
for (i in 1:n.iter) {
  init.ET.G[[i]] = g
}

ET.res = SGD.Monte.Carlo(init.ET.G, init.ET.G[[1]], theta0.ET, ERGM.ET.ss,
                         function(i) (1/(i)), n.draws=50, ERGM.ET.ss.diff, use.pkg=T, debug=F, verbose=T, model='ET')
print('Predicted theta: '%+%ET.res[[n.samples]])

plot(1:length(unlist(ET.res)), unlist(ET.res))