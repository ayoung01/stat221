#graph is stored as adjaceny matrix?

ERGM.prob = function( G, theta, ss ) {
  if(t(theta) %*% ss(G) < -1000){
    browser()
  }
  return(exp(t(theta) %*% ss(G))[1,1])
}

ERGM.log.prob = function( G, theta, ss) {
  return( (t(theta) %*% ss(G))[1,1])
}

toggle.edge = function( G ) {
  n_vertices = ncol(G)
  x = sample(1:n_vertices, 2) #guaranteed x[0] != x[1]
  G[x[1],x[2]] = 1-G[x[1], x[2]]
  G[x[2],x[1]] = G[x[1],x[2]]
  return(G)
}

#based off of Bartz paper
#todo: are we considering directed edges?
#using MCMC
#using gibbs sampling: given current graph, what is conditional probability of graph with one edge changed?
ERGM.MCMC = function( G_0, theta_0, ss ) {
  n = ncol(G_0)
  graphs = vector("list", n*(n-1))
  graph.old = G_0
  #graphs[[1]] = G_0
  for (i in 2:(n*(n-1))) {
    graph.new = toggle.edge(graph.old)
    #acceptance.prob = ERGM.prob(graph.new, theta_0, ss)/ERGM.prob(graphs[[i-1]], theta_0, ss)
    acceptance.prob = exp(ERGM.log.prob(graph.new, theta_0, ss) - ERGM.log.prob(graph.old, theta_0, ss))
    #print( graph.new - graphs[[i-1]] )
    #print(sprintf("%s %s %s %s", theta_0, ERGM.log.prob(graph.new, theta_0, ss), ERGM.log.prob(graphs[[i-1]], theta_0, ss), acceptance.prob))
    if( runif(1) < acceptance.prob) {
     graph.old = graph.new
    }
    #else {
    #  graphs[[i]] = graphs[[i-1]]
    #}
  }
  return(graph.old)
}

#takes in graph and outputs sufficient statistics
ERGM.sufficient.statistic = function( G ) {
  return(as.matrix(sum(G)/2)) #number of edges for now
}

#todo: How to choose G_0?
SGD.Monte.Carlo = function(G.data, G_0, theta_0, ss, n.samples = 1) {
  n.iters = length(G.data)
  thetas = vector("list", n.iters)
  thetas[[1]] = theta_0
  for( i in 2:n.iters ) {
    a = 1/i
    G.old = ERGM.MCMC(G_0, thetas[[i-1]], ss) # only obtain one sample each time
    s.avg = ss(G.old)
    #s.avg = ss(G.old[[length(G.old)]])
    C = diag(ncol(s.avg))
    thetas[[i]] = thetas[[i-1]] + a*C%*%( ss(G.data[[i]])- s.avg )
  }
  return(thetas)
}

#Bernoulli ERGM model where every edge is present at fixed probability
#e^theta/(1+e^theta)
Toulis.generate.graph = function(n, theta) {
  x = matrix(0, n, n)
  for(i in 1:n) {
    for(j in 1:(i-1)) {
      if( runif(1) < exp(theta)/(1+exp(theta))) {
        x[j,i] = 1
        x[i,j] = 1
      }
    }
  }
  return(x)
}

n.iters = 1000
theta = -0.9
n = 18
G.data = vector("list", n.iters)
for(i in 1:n.iters) {
  G.data[[i]] = Toulis.generate.graph(n, theta)
}

theta_0 = as.matrix(-0.1)
res = SGD.Monte.Carlo(G.data, G.data[[1]], theta_0, ERGM.sufficient.statistic)
