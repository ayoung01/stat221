#graph is stored as adjaceny matrix

ERGM.prob = function( G, theta, ss ) {
  return(exp(t(theta) %*% ss(G))[1,1])
}

ERGM.log.prob = function( G, theta, ss) {
  return( (t(theta) %*% ss(G))[1,1])
}

#takes in ss.diff (the difference in sufficient statistics given that edge(i,j) is toggled)
ERGM.log.prob.diff = function( G, edge, theta, ss.diff) {
    return( (t(theta) %*% ss.diff(G, edge))[1,1])
}

toggle.random.edge = function( G ) {
  n_vertices = ncol(G)
  edge = sample(1:n_vertices, 2) #guaranteed x[0] != x[1]
  return(toggle.edge(G, edge))
}

toggle.edge = function( G, edge) {
  G[edge[1],edge[2]] = 1 - G[edge[1], edge[2]]
  G[edge[2],edge[1]] = G[edge[1],edge[2]]
  return(G)
}

propose.edge = function( G ) {
  n_vertices = ncol(G)
  return( sample(1:n_vertices, 2))
}

#based off of Bartz paper
#todo: are we considering directed edges?
#using MCMC
#using gibbs sampling: given current graph, what is conditional probability of graph with one edge changed?
ERGM.MCMC = function( G_0, theta_0, ss, n_iters = ncol(G_0)*(ncol(G_0)-1)) {
  #n = ncol(G_0)
  graphs = vector("list", n_iters)
  graph.old = G_0
  #graphs[[1]] = G_0
  for (i in 2:n_iters) {
    graph.new = toggle.random.edge(graph.old)
    acceptance.prob = exp(ERGM.log.prob(graph.new, theta_0, ss) - ERGM.log.prob(graph.old, theta_0, ss))
    #print(sprintf("%s %s %s %s", theta_0, ERGM.log.prob(graph.new, theta_0, ss), ERGM.log.prob(graphs[[i-1]], theta_0, ss), acceptance.prob))
    if( runif(1) < acceptance.prob) {
     graph.old = graph.new
    }
  }
  return(graph.old)
}

ERGM.MCMC.fast = function( G_0, theta_0, ss, ss.diff, n_iters = ncol(G_0)*(ncol(G_0)-1) ) {
  graphs = vector("list", n_iters)
  graph.old = G_0
  for( i in 2:n_iters) {
    edge.proposal = propose.edge(graph.old)
    acceptance.prob = exp(ERGM.log.prob.diff(graph.old, edge.proposal, theta_0, ss.diff))
    if( runif(1) < acceptance.prob) {
      graph.old = toggle.edge(graph.old, edge.proposal)
    }
  }
  
  return(graph.old)
}

#takes in graph and outputs sufficient statistics
ERGM.sufficient.statistic = function( G ) {
  return(as.matrix(sum(G)/2)) #number of edges for now
}

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
ERGM.triad.ss = function( G) {
  res = numeric(3)
  res[1] = ERGM.edges(G)
  res[2] = ERGM.triangles(G)
  res[3] = ERGM.twostars(G)
  return(as.matrix(res))
}

ERGM.edges = function( G ) {
  return(sum(G)/2)
}

#defined as three edges that are mutually adjacent
ERGM.triangles = function( G ) {
  res = 0
  n = ncol(G)
  for( i in 3:n) {
    for( j in 1:(i-1)) {
      for( k in 1:(j-1))
        res = res + G[i,j]*G[j,k]*G[k,i]
    }
  }
  return(res)
}

#two star is defined by two edges being connected at a vertex
ERGM.twostars = function( G ) {
  res = 0
  n = ncol(G)
  for( i in 3:n ) {
    for( j in 1:(i-1)) {
      for( k in 1:(j-1)) {
        res = res + G[i,j]*G[j,k] + G[j,k]*G[k,i] + G[k,i]*G[i,j]
      }
    }
  }
  return(res)
}

ERGM.triad.ss.diff = function( G, edge ) {
  res = numeric(3)
  res[1] = ERGM.edges.diff(G, edge)
  res[2] = ERGM.triangles.diff(G, edge)
  res[3] = ERGM.twostars.diff(G, edge)
  if(G[edge[1],edge[2]] == 1) {
    res = -1*res
  }
  return(as.matrix(res))
}

ERGM.edges.diff = function(G, edge) {
  res = 1
  return(res)
}

ERGM.triangles.diff = function(G, edge) {
  res = 0
  for(k in 1:ncol(G)) {
    if( k != edge[1] && k != edge[2]) {
      res = res + G[edge[2],k]*G[k,edge[1]]
    }
  }
}

ERGM.twostars.diff = function(G, edge) {
  res = 0
  for(k in 1:ncol(G)) {
    if( k != edge[1] && k != edge[2]) {
      res = res + G[edge[1],k] + G[k,edge[2]]
    }
  }
  return(res)
}

avg.over.list = function( l ) {
  s = 0
  for( i in l) {
    s = s + i
  }
  return( s/length(l))
}

#todo: How to choose G_0?
SGD.Monte.Carlo = function(G.data, G_0, theta_0, ss, learning.rate, n.samples = 1, ss.diff = NULL) {
  n.iters = length(G.data)
  thetas = vector("list", n.iters)
  thetas[[1]] = theta_0
  for( i in 2:n.iters ) {
    a = learning.rate(i)
    G.samples = vector("list", n.samples)
    for( j in 1:n.samples) {
      if( typeof(ss.diff) == "closure" ) {
        G.samples[[j]] = ERGM.MCMC.fast(G_0, thetas[[i-1]], ss, ss.diff)
      } else {
        G.samples[[j]] = ERGM.MCMC(G_0, thetas[[i-1]], ss)  
      }
    }
    #G.old = ERGM.MCMC(G_0, thetas[[i-1]], ss) # only obtain one sample each time
    #s.avg = ss(G.old)
    s.avg = avg.over.list(lapply(G.samples, ss))
    #print(s.avg)
    C = diag(ncol(s.avg))
    thetas[[i]] = thetas[[i-1]] + a*C%*%( ss(G.data[[i]])- s.avg )
  }
  return(thetas)
}

#Bernoulli ERGM model where every edge is present at fixed probability
#e^theta/(1+e^theta)
Toulis.generate.graph = function(n, theta) {
  x = matrix(0, n, n)
  for(i in 2:n) {
    for(j in 1:(i-1)) {
      #print(sprintf("%s %s", i, j))
      if( runif(1) < exp(theta)/(1+exp(theta))) {
        x[j,i] = 1
        x[i,j] = 1
      }
    }
  }
  return(x)
}

simple.learning.rate = function( i ) {
  return(1/i)
}

n.iters = 100
theta = -0.9
n = 18
n.samples = 3
G.data = vector("list", n.iters)
for(i in 1:n.iters) {
  G.data[[i]] = Toulis.generate.graph(n, theta)
}

theta_0 = as.matrix(-0.1)
start = proc.time()
#res = SGD.Monte.Carlo(G.data, G.data[[1]], theta_0, ERGM.edges.ss, simple.learning.rate, n.samples)
res = SGD.Monte.Carlo(G.data, G.data[[1]], theta_0, ERGM.edges.ss, simple.learning.rate, n.samples, ERGM.edges.ss.diff)
end = proc.time()
print( end - start )
