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
                               theta.actual, ss, ss.diff, lr, use.pkg=T, debug=F) {
  res = SGD.Monte.Carlo(G.samples, G.samples[[1]], theta_0, ss,
                        lr, n.draws, ss.diff, use.pkg=use.pkg, debug=debug)
  print('Actual theta: '%+%theta.actual)
  print('Predicted theta: '%+%res[[n.samples]])
  return(res)
}

n.nodes = 50
n.samples = 100





# ET ----------------------------------------------------------------------

# theta.ET.actual = as.matrix(rnorm(2, -1, 1))
theta0.ET = as.matrix(c(0,0))
init.ET.G = ERGM.ET.generate.samples(n.nodes, n.samples, theta.ET.actual, use.pkg=T)
init.ET.G.filtered = ERGM.ET.filter.degenerates(init.ET.G)
ET.res = Bartz.experiment(n.nodes, n.samples, n.draws=10, init.ET.G, theta_0, theta.ET.actual,
                       ss=ERGM.ET.ss, ss.diff=ERGM.ET.ss.diff, lr=simple.lr, use.pkg=F, debug=F)



plot(1:length(unlist(res)), unlist(res))

# triad -------------------------------------------------------------------

theta.triad.actual = as.matrix(rnorm(3, -1, 1))
theta0.triad = as.matrix(c(0,0,0))
init.triad.G = ERGM.ET.generate.samples(n.nodes, n.samples, theta.ET.actual, use.pkg=T)

triad.res = Bartz.experiment(n.nodes, n.samples, n.draws=10, init.ET.G, theta_0, theta.ET.actual,
                             ss=ERGM.triad.ss, ss.diff=ERGM.triad.ss.diff, lr=simple.lr, use.pkg=F, debug=F)


