library(ergm)


# Graph statistics --------------------------------------------------------

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
  r1 = 1
  if(G[edge[1],edge[2]] == 1) {
    r1 = -r1
  }
  return(r1)
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


# Bartz -------------------------------------------------------------------

Bartz.experiment = function(G.samples, theta.actual, theta0, verbose=F, n.reps=25, model) {
  n.samples = length(G.samples)
  if (verbose) {
    # look at graph statistics
    edges = sapply(G.samples, ERGM.edges)
    triangles = sapply(G.samples, ERGM.triangles)
    twostars = sapply(G.samples, ERGM.twostars)

    print(table(edges))
    print(table(triangles))
    print(table(twostars))
  }
  if (model=='ET') {
    ss = ERGM.ET.ss
    ss.diff = ERGM.ET.ss.diff
  }
  if (model=='triad') {
    ss = ERGM.triad.ss
    ss.diff = ERGM.triad.ss.diff
  }

  res = list()
  for (i in 1:n.reps) {
    res[[i]] = SGD.Monte.Carlo(G.samples, G.samples[[1]], theta0, ss,
                               function(i) (1/(5+i)), n.draws=10, ss.diff, use.pkg=T,
                               debug=F, verbose=verbose, model=model)
    print('Actual theta: '%+%theta.actual)
    print('Predicted theta: '%+%res[[i]][[n.samples]])
    plot(1:length(unlist(res[[i]])), unlist(res[[i]]))
  }
  return(res)
}

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


# Monte Carlo SGD ---------------------------------------------------------------

SGD.Monte.Carlo = function(G.data, G_0, theta_0, ss, lr, n.draws = 1,
                           ss.diff = NULL, use.pkg=F, debug=F, verbose=F, model) {
  n.nodes = ncol(G_0)
  n.iters = length(G.data)
  thetas = vector("list", n.iters)
  thetas[[1]] = theta_0
  pb <- txtProgressBar(min = 0, max = n.iters, style = 3)
  # Fisher information = variance of ss
  Fisher.hat = diag(nrow(theta_0))
  for( i in 2:n.iters ) {
    a = lr(i)
    if (use.pkg) {
      G.samples = ERGM.generate.samples(n.nodes, n.draws, thetas[[i-1]], use.pkg=T, model=model)
    } else {
      G.samples = vector("list", n.draws)
      for( j in 1:n.draws) {
        if( typeof(ss.diff) == "closure" ) {
          G.samples[[j]] = ERGM.MCMC(G_0, thetas[[i-1]], ss, ss.diff)
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
    if(nrow(ss.mat) == 1) {
      Fisher.hat = diag(1)
    } else {
      Fisher.hat = 0.4 * Fisher.hat + 0.6 * cov(t(ss.mat))
    }
    C = solve(Fisher.hat)
    s.avg = avg.over.list(ss.list)
    thetas[[i]] = thetas[[i-1]] + a*C%*%( (ss(G.data[[i]])- s.avg) )

    if (verbose) {
      print('theta_'%+%(i-1))
      print(thetas[[i-1]])
      print('theta_'%+%i)
      print(thetas[[i]])
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  return(thetas)
}


# Graph utility functions --------------------------------------------------------

#based off of Bartz paper
ERGM.MCMC = function( G_0, theta_0, ss, ss.diff, n_iters = ncol(G_0)*(ncol(G_0)-1), verbose=F ) {
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

ERGM.generate.samples = function(n.nodes, n.samples, theta, use.pkg=T, model='ET') {
  if (use.pkg) {
    if (model=='ET') {
      net = formula('network(n.nodes, directed=F) ~ edges + triangles')
      G.samples = simulate(net, nsim=n.samples, coef=theta, sequential=F)
    }
    if (model=='triad') {
      net = formula('network(n.nodes, directed=F) ~ edges + triangles + twopath')
      G.samples = simulate(net, nsim=n.samples, coef=theta, sequential=F)
    }
    if (model=='edges') {
      net = formula('network(n.nodes, directed=F) ~ edges')
      G.samples = simulate(net, nsim=n.samples, coef=theta, sequential=F)
    }
    G.samples = lapply(G.samples, as.matrix)
  }
  # use our own broken method
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

# General utility functions -----------------------------------------------

`%+%` <- function (x, y) {
  paste(x, y, sep = "")
}
avg.over.list = function( l ) {
  s = 0
  for( i in l) {
    s = s + i
  }
  return( as.matrix(s/length(l)))
}

simple.lr = function( i ) {
  return(1/i)
}

parameter.lr = function( i, a = 1) {
  return(a/i)
}

# Functions to calculate and plot MSE -------------------------------------

mse = function(sim, obs) {
  return(sum((sim - obs)^2) / length(sim))
}

get.pred.thetas = function(res) {
  t(do.call(cbind, lapply(res, function(l) {
    # get last parameter for each theta chain
    l[[length(l)]]
  })))
}

get.first.thetas = function(res, n.iters=NULL) {
  if (is.null(n.iters)) {
    n.iters = length(res[[1]]) - 1
  }
  p = nrow(res[[1]][[1]])
  thetas = do.call(cbind, lapply(res, function(l) {
    # get first 5 parameter iterations for each theta chain
    x = l[[2]]
    for (i in 3:(n.iters+1)) {
      x = cbind(x, l[[i]])
    }
    return(t(x))
  }))
  ans = vector('list', p)
  for (i in 1:p) {
    ans[[i]] = thetas[, seq(i, ncol(thetas), p)]
  }
  return(ans)
}

get.mse.var = function(thetas_list, theta.actual, theta.names=c('edges','triangles')) {
  n.iters = nrow(thetas_list[[1]])
  p = length(thetas_list)
  stats = data.frame()
  for (i in 1:p) {
    thetas = thetas_list[[i]]
    mses = numeric(n.iters)
    variances = numeric(n.iters)
    for (j in 1:n.iters) {
      mses[j] = mse(thetas[j, ], theta.actual[i])/mse(thetas[1, ], theta.actual[i])
      variances[j] = var(thetas[j, ])/var(thetas[1, ])
    }
    stats = rbind(stats, cbind(matrix(mses), matrix(1:n.iters), 'MSE', theta.names[i]))
    stats = rbind(stats, cbind(matrix(variances), matrix(1:n.iters), 'variance', theta.names[i]))
  }
  colnames(stats) = c('ratio', 'iteration', 'stat', 'theta')
  stats$theta = as.factor(stats$theta)
  stats$stat = as.factor(stats$stat)
  stats$iteration = as.numeric(stats$iteration)
  stats$ratio = as.numeric(stats$ratio)
  return(stats)
}

plot.var = function(res, theta.actual, model='ET', ylim=c(0,1)) {
  theta.names = c('edges', 'triangles')
  if (model=='triad') {
    theta.names = c('edges', 'triangles', 'twostar')
  }
  first.10.thetas = get.first.thetas(res, n.iters=10)
  stats.10 = get.mse.var(first.10.thetas, theta.actual, theta.names)

  first.thetas = get.first.thetas(res)
  stats = get.mse.var(first.thetas, theta.actual, theta.names)

  p = ggplot(stats.10, aes(x=iteration, y=ratio, colour=stat)) + ylim(ylim) + scale_x_discrete(1:9)
  p + geom_line() + facet_grid(. ~ theta) +
    ggtitle(sprintf('Variance of parameter estimates by iteration for %s model', model))

  p = ggplot(stats, aes(x=iteration, y=ratio, colour=stat)) + ylim(ylim)
  p + geom_line() + facet_grid(. ~ theta) +
    ggtitle(sprintf('Variance of parameter estimates by iteration for %s model', model))
}