source('ERGM.MCMC.R')

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
  if(G[edge[1],edge[2]] == 1) {
    res = -res
  }
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
  if(G[edge[1],edge[2]] == 1) {
    res = -res
  }
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

Bartz.ET.experiment = function(n.nodes = 100, n.samples = 1000) {
  theta.actual = as.matrix(rnorm(2, -1, 1))
  G.samples = ERGM.ET.generate.samples(n.nodes, n.samples, theta.actual)
  print(sprintf("Num samples before: ", length(G.samples)))
  G.samples = ERGM.ET.filter.degenerates(G.samples)
  print(sprintf("Num samples after: ", length(G.samples)))
  
  theta_0 = as.matrix(c(-0.1,-0.1))
  n.samples.per.iter = 1
  res = SGD.Monte.Carlo(G.samples, G.samples[[1]], theta_0, ERGM.ET.ss, 
                        simple.learning.rate, n.samples.per.iter, ERGM.ET.ss.diff)
}

Bartz.triad.experiment = function(n.nodes = 100, n.samples = 1000) {
  theta.actual = as.matrix(rnorm(3, -1, 1))
  G.samples = ERGM.triad.generate.samples(n.nodes, n.samples, theta.actual)
  #x = avg.over.list(lapply(G.samples, ERGM.triad.ss))
  #check for degeneracy of ERGM samples (meaning that won't be able to use for model)
  G.samples = ERGM.triad.filter.degenerates(G.samples)
  
  theta_0 = as.matrix(c(-0.1,-0.1, -0.1))
  n.samples.per.iter = 1
  res = SGD.Monte.Carlo(G.samples, G.samples[[1]], theta_0, ERGM.triad.ss, 
                        simple.learning.rate, n.samples.per.iter, ERGM.triad.ss.diff)
} 

n.nodes = 50
n.samples = 1000
#Bartz.triad.experiment(n.nodes, n.samples)
res = Bartz.ET.experiment(n.nodes, n.samples)





