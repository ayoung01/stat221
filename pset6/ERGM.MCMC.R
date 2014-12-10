`%+%` <- function (x, y) {
  paste(x, y, sep = "")
}
#graph is stored as an adjacency matrix

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

ERGM.edges = function( G ) {
  return(sum(G)/2)
}

#defined as three edges that are mutually adjacent
ERGM.triangles = function( G ) {
  r1 = 0
  n = ncol(G)
  for( i in 3:n) {
    for( j in 2:(i-1)) {
      for( k in 1:(j-1)) {
        r1 = r1 + G[i,j]*G[j,k]*G[k,i]
      }
    }
  }
  return(r1)
}

#two star is defined by two edges being connected at a vertex
ERGM.twostars = function( G ) {
  res = 0
  n = ncol(G)
  for( i in 3:n ) {
    for( j in 2:(i-1)) {
      for( k in 1:(j-1)) {
        res = res + G[i,j]*G[j,k] + G[j,k]*G[k,i] + G[k,i]*G[i,j]
      }
    }
  }
  return(res)
}

ERGM.edges.diff = function(G, edge) {
  res = 1
  if(G[edge[1],edge[2]] == 1) {
    res = -res
  }
  return(res)
}

ERGM.triangles.diff = function(G, edge) {
  res = 0
  for(k in 1:ncol(G)) {
    if( k != edge[1] && k != edge[2]) {
      res = res + G[edge[2],k]*G[k,edge[1]]
    }
  }
  if(G[edge[1],edge[2]] == 1) {
    res = -res
  }
  return( res )
}

ERGM.twostars.diff = function(G, edge) {
  res = 0
  for(k in 1:ncol(G)) {
    if( k != edge[1] && k != edge[2]) {
      res = res + G[edge[1],k] + G[k,edge[2]]
    }
  }
  if(G[edge[1],edge[2]] == 1) {
    res = -res
  }
  return(res)
}

#based off of Bartz paper
#using MCMC
ERGM.MCMC = function( G_0, theta_0, ss, n_iters = ncol(G_0)*(ncol(G_0)-1)) {
  graphs = vector("list", n_iters)
  graph.old = G_0
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

ERGM.MCMC.fast = function( G_0, theta_0, ss, ss.diff, n_iters = ncol(G_0)*(ncol(G_0)-1), verbose=F ) {
  graphs = vector("list", n_iters)
  graph.old = G_0
  acceptance.percentage = 0
  for( i in 2:n_iters) {
    edge.proposal = propose.edge(graph.old)
    acceptance.prob = exp(ERGM.log.prob.diff(graph.old, edge.proposal, theta_0, ss.diff))
    if( runif(1) < acceptance.prob) {
      graph.old = toggle.edge(graph.old, edge.proposal)
      acceptance.percentage = acceptance.percentage + 1
    }
  }
  if (verbose) {
    print('Acceptance percentage: '%+%(acceptance.percentage/n_iters*100)%+%'%')
  }
  return(graph.old)
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

ERGM.ET.generate.samples = function(n.nodes, n.samples, theta, use.pkg=T) {
  if (use.pkg) {
    G.samples = simulate(network(n.nodes, directed=F) ~ edges + triangles, nsim=n.samples, coef=theta, sequential=F)
    G.samples = lapply(G.samples, as.matrix)
  }
  else {
    G_0 = generate.random.graph(n.nodes, 0.5)
    G.samples = vector("list", n.samples)

    #let markov chain mix a lot for first sample
    G.samples[[1]] = ERGM.MCMC.fast( G_0, theta, ERGM.ET.ss, ERGM.ET.ss.diff, n.nodes**3 )
    for( i in 2:n.samples) {
      G.samples[[i]] = ERGM.MCMC.fast( G.samples[[i-1]], theta, ERGM.ET.ss, ERGM.ET.ss.diff, n.nodes**2 )
    }
  }
  return(G.samples)
}

avg.over.list = function( l ) {
  s = 0
  for( i in l) {
    s = s + i
  }
  return( as.matrix(s/length(l)))
}

#TODO: How to choose G_0 and C?
SGD.Monte.Carlo = function(G.data, G_0, theta_0, ss, lr, n.draws = 1,
                           ss.diff = NULL, use.pkg=T, debug=F) {
  n.nodes = ncol(G_0)
  n.iters = length(G.data)
  thetas = vector("list", n.iters)
  thetas[[1]] = theta_0
  pb <- txtProgressBar(min = 0, max = n.iters, style = 3)
  for( i in 2:n.iters ) {
    a = lr(i)
    if (use.pkg) {
      G.samples = ERGM.ET.generate.samples(n.nodes, n.draws, thetas[[i-1]], use.pkg=T)
    } else {
      G.samples = vector("list", n.draws)
      for( j in 1:n.draws) {
        if( typeof(ss.diff) == "closure" ) {
          G.samples[[j]] = ERGM.MCMC.fast(G_0, thetas[[i-1]], ss, ss.diff)
        } else {
          G.samples[[j]] = ERGM.MCMC(G_0, thetas[[i-1]], ss)
        }
      }
    }

    ss.list = lapply(G.samples, ss)
    ss.mat = do.call(cbind, ss.list)
    if (debug){
      browser()
    }
    # calculate condition matrix
    Fisher.hat = cov(t(ss.mat))
    C = solve(Fisher.hat)

    s.avg = avg.over.list(ss.list)

    # condition matrix should be inverse of Fisher information
    # Fisher information = variance of ss
    # can run separate SGD on Fisher information
    thetas[[i]] = thetas[[i-1]] + a*C%*%( (ss(G.data[[i]])- s.avg) )
    setTxtProgressBar(pb, i)
  }
  close(pb)
  return(thetas)
}

#n = number of nodes, p= probability of particular edge existing
generate.random.graph = function(n, p) {
  x = matrix(0, n, n)
  for(i in 2:n) {
    for(j in 1:(i-1)) {
      if( runif(1) < p) {
        x[j,i] = 1
        x[i,j] = 1
      }
    }
  }
  return(x)
}

simple.lr = function( i ) {
  return(1/i)
}

parameter.lr = function( i, a = 1) {
  return(a/i)
}