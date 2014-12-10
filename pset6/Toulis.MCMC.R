source('ERGM.MCMC.R')

ERGM.edges.ss = function( G ) {
  return(as.matrix(ERGM.edges(G)))
}

ERGM.edges.ss.diff = function( G, edge ) {
  res = ERGM.edges.diff(G, edge)
  if(G[edge[1],edge[2]] == 1) {
    res = -1*res
  }
  return(as.matrix(res))
}

#Bernoulli ERGM model where every edge is present at fixed probability
#e^theta/(1+e^theta)
Toulis.generate.graph = function(n, theta) {
  p = exp(theta)/(1+exp(theta))
  return(generate.random.graph(n,p))
}

graph.generating.sanity = function( theoretical.samples, my.samples ) {
  theoretical.avg.edges = avg.over.list(lapply(theoretical.samples, ERGM.edges))[1,1]
  my.avg.edges = avg.over.list(lapply(my.samples, ERGM.edges))[1,1]
  print(sprintf("Theoretical # edges %s vs actual # edges %s", theoretical.avg.edges, my.avg.edges))
  if( abs(theoretical.avg.edges - my.avg.edges)/theoretical.avg.edges > 0.02) {
    print(">2% difference in generation!")
  }
}

MCMC.test = function(n.iters = 1000, theta = -0.9, n.nodes = 18, n.samples = 3) {
  G.data = vector("list", n.iters)
  for( i in 1:n.iters) {
    G.data[[i]] = Toulis.generate.graph(n, theta)
  }

  G.samples = vector("list", n.iters)
  G.samples[[1]] = ERGM.MCMC.fast( G_0, as.matrix(theta), ERGM.edges.ss, ERGM.edges.ss.diff, n**3 )
  for( i in 2:n.iters ) {
    G.samples[[i]] = ERGM.MCMC.fast( G.samples[[i-1]], as.matrix(theta), ERGM.edges.ss, ERGM.edges.ss.diff, n*(n-1) )
  }
  graph.generating.sanity(G.data, G.samples)
}

#MCMC.test()

Toulis.Bernoulli.Experiment = function(n.iters = 10000, theta = -0.9, n.nodes = 18, n.samples = 1) {
  G.data = vector("list", n.iters)
  for(i in 1:n.iters) {
    G.data[[i]] = Toulis.generate.graph(n.nodes, theta)
  }

  theta_0 = as.matrix(-0.1)
  start = proc.time()
  G_0 = generate.random.graph(n.nodes, 0.5)
  res = SGD.Monte.Carlo(G.data, G_0, theta_0, ERGM.edges.ss, function(x) parameter.lr( x, 5), n.samples, ERGM.edges.ss.diff)
  end = proc.time()
  print( end - start )
  return(res)
}

res = Toulis.Bernoulli.Experiment()
