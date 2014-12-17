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

# karate ------------------------------------------------------------------

# karate data from http://www-personal.umich.edu/~mejn/netdata/
# 34 nodes
# 78 edges, 45 triangles, 528 two-stars
library(igraph)

karate = read.graph('karate/karate.gml', format='gml')
plot(karate)
karate = get.adjacency(karate, sparse=F)
karate.ET.ergm = ergm(karate~ edges + triangles)
karate.triad.ergm = ergm(karate~ edges + triangles + twopath)


# Run experiments ---------------------------------------------------------

n.nodes = 20
n.samples = 1000

# extract subnetwork of 20 random nodes
karate.samples = sample.subnetworks(karate, n.samples = 1000, n.nodes = 20)

# fit ET model
karate.theta0.ET = as.matrix(c(1, -1))
karate.ET.res = SGD.Monte.Carlo(karate.samples, karate.samples[[1]], karate.theta0.ET, ERGM.ET.ss,
                         function(i) (1/(5+i)), n.draws=10, ERGM.ET.ss.diff, use.pkg=T, debug=F, verbose=T, model='ET')
print('Predicted theta: '%+%karate.ET.res[[n.samples]])

karate.ET.res = t(do.call(cbind, karate.ET.res))
colnames(karate.ET.res) = c('edges', 'triangles')
karate.ET.res = as.data.frame(karate.ET.res)
karate.ET.res = stack(karate.ET.res)
karate.ET.res$iteration = rep(1:n.samples, 2)
colnames(karate.ET.res) = c('theta', 'stat', 'iteration')

p = ggplot(karate.ET.res, aes(x=iteration, y=theta, colour=stat)) + geom_line()
p + ggtitle('Predicted theta by iteration for ET model')


# fit triad model
karate.theta0.triad = as.matrix(c(-1, 1, -1))
karate.triad.res = SGD.Monte.Carlo(karate.samples, karate.samples[[1]], karate.theta0.triad, ERGM.triad.ss,
                                function(i) (1/(50+i)), n.draws=10, ERGM.triad.ss.diff, use.pkg=T, debug=F, verbose=T, model='triad')
print('Predicted theta: '%+%karate.triad.res[[n.samples]])

karate.triad.res = t(do.call(cbind, karate.triad.res))
colnames(karate.triad.res) = c('edges', 'triangles', 'twostar')
karate.triad.res = as.data.frame(karate.triad.res)
karate.triad.res = stack(karate.triad.res)
karate.triad.res$iteration = rep(1:n.samples, 3)
colnames(karate.triad.res) = c('theta', 'stat', 'iteration')

p = ggplot(karate.triad.res, aes(x=iteration, y=theta, colour=stat)) + geom_line()
p + ggtitle('Predicted theta by iteration for triad model')

