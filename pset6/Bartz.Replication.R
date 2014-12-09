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
  res[2] = ERGM.twostars(G)
  return(as.matrix(res))
}

ERGM.ET.ss.diff = function( G, edge ) {
  res = numeric(2)
  res[1] = ERGM.edges.diff(G, edge)
  res[2] = ERGM.twostars.diff(G, edge)
  return(as.matrix(res))
}

ERGM.triad.generate.samples = function(n.nodes, n.samples, theta.actual) {
  G_0 = generate.random.graph(n.nodes, 0.5)
  G.samples = vector("list", n.samples)

  #let markov chain mix a lot for first sample
  G.samples[[1]] = ERGM.MCMC.fast( G_0, theta.actual, ERGM.triad.ss, ERGM.triad.ss.diff, n.nodes**3 )
  for( i in 2:n.samples) {
    G.samples[[i]] = ERGM.MCMC.fast( G.samples[[i-1]], theta.actual, ERGM.triad.ss, ERGM.triad.ss.diff, n.nodes**2 )
  }
  return(G.samples)
}

ERGM.ET.generate.samples = function(n.nodes, n.samples, theta.actual) {
  G_0 = generate.random.graph(n.nodes, 0.5)
  G.samples = vector("list", n.samples)

  #let markov chain mix a lot for first sample
  G.samples[[1]] = ERGM.MCMC.fast( G_0, theta.actual, ERGM.ET.ss, ERGM.ET.ss.diff, n.nodes**3 )
  for( i in 2:n.samples) {
    G.samples[[i]] = ERGM.MCMC.fast( G.samples[[i-1]], theta.actual, ERGM.ET.ss, ERGM.ET.ss.diff, n.nodes**2 )
  }
  return(G.samples)
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
  min.twostars = 0
  n.nodes = ncol(G.samples[[1]])
  max.edges = n.nodes * (n.nodes-1)
  max.twostars = n.nodes*(n.nodes-1)*(n.nodes-2)/2
  num.edges = lapply(G.samples, ERGM.edges)
  num.twostars = lapply(G.samples,ERGM.twostars)
  return( G.samples[ num.edges > min.edges && num.edges < max.edges &&
                       num.twostars > min.twostars && num.twostars < max.twostars])
}

Bartz.ET.experiment = function(n.nodes = 100, n.samples = 1000, init.G=NULL,
                               init.theta.actual=NULL, debug=F) {
  theta.actual = as.matrix(rnorm(2, -1, 1))
  if (!is.null(init.theta.actual)) {
    theta.actual = init.theta.actual
  }
  if (is.null(init.G)) {
    print('Generating samples...')
    G.samples = ERGM.ET.generate.samples(n.nodes, n.samples, theta.actual)
    print(sprintf("Num samples before: %s", length(G.samples)))
    G.samples = ERGM.ET.filter.degenerates(G.samples)
    print(sprintf("Num samples after: %s", length(G.samples)))
  } else {
    print('Initial samples supplied!')
    G.samples = init.G
  }

  theta_0 = as.matrix(c(-0.1,-0.1))
  n.draws = 10
  res = SGD.Monte.Carlo(G.samples, G.samples[[1]], theta_0, ERGM.ET.ss,
                        simple.learning.rate, n.draws, ERGM.ET.ss.diff, debug=debug)
  print('Actual theta: '%+%theta.actual)
  print('Predicted theta: '%+%res[[n.samples]])
  return(res)
}

Bartz.triad.experiment = function(n.nodes = 100, n.samples = 1000) {
  theta.actual = as.matrix(rnorm(3, -1, 1))
  print('Generating samples...')
  G.samples = ERGM.triad.generate.samples(n.nodes, n.samples, theta.actual)
  #x = avg.over.list(lapply(G.samples, ERGM.triad.ss))
  #check for degeneracy of ERGM samples (meaning that won't be able to use for model)
  print('Filtering degenerates...')
  G.samples = ERGM.triad.filter.degenerates(G.samples)
  print('Done!')
  theta_0 = as.matrix(c(-0.1,-0.1, -0.1))
  n.draws = 1
  res = SGD.Monte.Carlo(G.samples, G.samples[[1]], theta_0, ERGM.triad.ss,
                        simple.learning.rate, n.draws, ERGM.triad.ss.diff)
  print('Actual theta: '%+%theta.actual)
  print('Predicted theta: '%+%res[[n.samples]])
  return(res)
}

n.nodes = 50
n.samples = 1000
# res = Bartz.triad.experiment(n.nodes, n.samples)

theta.ET.actual = as.matrix(rnorm(2, -1, 1))
# init.ET.G = ERGM.ET.generate.samples(n.nodes, n.samples, theta.ET.actual)
# res = Bartz.ET.experiment(n.nodes, n.samples, init.ET.G, theta.ET.actual, debug=F)

# res = Bartz.ET.experiment(n.nodes, 100, init.ET.G[1:100], theta.ET.actual, debug=F)




g.sim = simulate(network(n.nodes, directed=F) ~ edges + twopath, nsim= 10, coef=theta.ET.actual)
g.sim.mat = lapply(g.sim, as.matrix)
res = Bartz.ET.experiment(n.nodes, 10, g.sim.mat, theta.ET.actual, debug=F)




