source('ERGM.MCMC.R')
library(ergm)

ERGM.triad.ss = function( G) {
  res = numeric(3)
  res[1] = ERGM.edges(G)
  res[2] = ERGM.triangles(G)
  res[3] = ERGM.twostars(G)
  return(as.matrix(res))
}

ERGM.triad.ss.diff = function( G, edge ) {
  res = numeric(3)
  res[1] = ERGM.edges.diff(G, edge)
  res[2] = ERGM.triangles.diff(G, edge)
  res[3] = ERGM.twostars.diff(G, edge)

  return(as.matrix(res))
}

ERGM.ET.ss = function( G) {
  res = numeric(2)
  res[1] = ERGM.edges(G)
  res[2] = ERGM.triangles(G)
  return(as.matrix(res))
}

ERGM.ET.ss.diff = function( G, edge ) {
  res = numeric(2)
  res[1] = ERGM.edges.diff(G, edge)
  res[2] = ERGM.triangles.diff(G, edge)
  return(as.matrix(res))
}

ERGM.triad.filter.degenerates = function(G.samples) {
  min.edges = 0
  max.edges = ncol(G.samples[[1]])*(ncol(G.samples[[1]] - 1))
  num.edges = lapply(G.samples, ERGM.edges)
  return( G.samples[ num.edges > min.edges && num.edges < max.edges] )
}

#how exactly to define degenerates
ERGM.ET.filter.degenerates = function(G.samples) {
  min.edges = 0
  min.triangles = 0
  n.nodes = ncol(G.samples[[1]])
  max.edges = n.nodes * (n.nodes-1)
  max.triangles = n.nodes*(n.nodes-1)*(n.nodes-2)/2
  num.edges = lapply(G.samples, ERGM.edges)
  num.triangles = lapply(G.samples,ERGM.triangles)
  return( G.samples[ num.edges > min.edges && num.edges < max.edges &&
                       num.triangles > min.triangles && num.triangles < max.triangles])
}

Bartz.experiment = function(n.nodes=100, n.samples=1000, n.draws=10, G.samples, theta_0=as.matrix(c(-0.1,-0.1)),
                               theta.actual, ss, ss.diff, lr, use.pkg=T, debug=F, verbose=F, model) {
  res = SGD.Monte.Carlo(G.samples, G.samples[[1]], theta_0, ss,
                        lr, n.draws, ss.diff, use.pkg=use.pkg, debug=debug, verbose, model)
  print('Actual theta: '%+%theta.actual)
  print('Predicted theta: '%+%res[[n.samples]])
  return(res)
}

n.nodes = 20
n.samples = 1000

# ET ----------------------------------------------------------------------

theta.ET.actual = as.matrix(rnorm(2, -1, 1))
print(theta.ET.actual)
# theta.ET.working = theta.ET.actual
theta0.ET = as.matrix(c(0,0))
init.ET.G = ERGM.generate.samples(n.nodes, n.samples, theta.ET.actual, use.pkg=T, model='ET')
init.ET.G.filtered = ERGM.ET.filter.degenerates(init.ET.G)
ET.res = Bartz.experiment(n.nodes, n.samples, n.draws=10, init.ET.G.filtered, theta0.ET, theta.ET.actual,
                       ss=ERGM.ET.ss, ss.diff=ERGM.ET.ss.diff, lr=simple.lr, use.pkg=F, debug=F, verbose=T, model='ET')

# look at graph statistics
ET.edges = sapply(init.ET.G, ERGM.edges)
ET.triangles = sapply(init.ET.G, ERGM.triangles)
ET.twostars = sapply(init.ET.G, ERGM.twostars)

print(table(ET.edges))
print(table(ET.triangles))
print(table(ET.twostars))

plot(1:length(unlist(ET.res)), unlist(ET.res))

# triad -------------------------------------------------------------------

theta.triad.actual = as.matrix(rnorm(3, 0, 1))
theta0.triad = as.matrix(c(0,0,0))
init.triad.G = ERGM.generate.samples(n.nodes, n.samples, theta.triad.actual, use.pkg=T, model='triad')

# look at graph statistics
triad.edges = sapply(init.triad.G, ERGM.edges)
triad.triangles = sapply(init.triad.G, ERGM.triangles)
triad.twostars = sapply(init.triad.G, ERGM.twostars)

# filter out graphs with zero triangles
# notriangles = which(triad.triangles==0)
# init.triad.G.filtered = init.triad.G[-notriangles]

# print(sprintf('Edges: %s, Triangles: %s, Two-stars: %s', triad.edges))
#triad.log.probs = lapply( init.triad.G, function (x) ERGM.log.prob(x, theta.triad.actual, ERGM.triad.ss))
#sum(unlist(triad.log.probs))

triad.res = Bartz.experiment(n.nodes, n.samples, n.draws=10, init.triad.G, theta0.triad, theta.triad.actual,
                             ss=ERGM.triad.ss, ss.diff=ERGM.triad.ss.diff, lr=simple.lr, use.pkg=F, debug=F, verbose=T, model='triad')
plot(1:length(unlist(triad.res)), unlist(triad.res))

