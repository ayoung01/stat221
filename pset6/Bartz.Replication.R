source('ERGM.MCMC.R')

Bartz.experiment = function(G.samples, theta.actual, theta0, verbose=F, n.reps=25, model) {
  n.samples = len(G.samples)
  if (verbose) {
    # look at graph statistics
    edges = sapply(G.samples, ERGM.edges)
    triangles = sapply(G.samples, ERGM.triangles)
    twostars = sapply(G.samples, ERGM.twostars)

    print(table(edges))
    print(table(triangles))
    print(table(twostars))
  }

  res = list()

  for (i in 1:n.reps) {
    res[[i]] = SGD.Monte.Carlo(G.samples, G.samples[[1]], theta0, ERGM.ET.ss,
                               function(i) (1/(5+i)), n.draws=10, ERGM.ET.ss.diff, use.pkg=T,
                               debug=F, verbose=verbose, model=model)
    print('Actual theta: '%+%theta.actual)
    print('Predicted theta: '%+%res[[i]][[n.samples]])
    plot(1:length(unlist(res[[i]])), unlist(res[[i]]))
  }
}

n.nodes = 20
n.samples = 1000
n.reps = 25

# ET ----------------------------------------------------------------------

theta.ET.actual = as.matrix(rnorm(2, -1, 1))
print('Actual theta: '%+%theta.ET.actual)
theta0.ET = as.matrix(c(0,0))
G.ET.samples = ERGM.generate.samples(n.nodes, n.samples, theta.ET.actual, use.pkg=T, model='ET')
# G.ET.samples.filtered = ERGM.ET.filter.degenerates(G.ET.samples)

ET.res = Bartz.experiment(G.ET.samples, theta.ET.actual, theta0.ET, verbose=F, n.reps, model='ET')


# triad -------------------------------------------------------------------

theta.triad.actual = as.matrix(rnorm(3, -1, 1))
theta0.triad = as.matrix(c(-1,-1,-1))
G.triad.samples = ERGM.generate.samples(n.nodes, n.samples, theta.triad.actual, use.pkg=T, model='triad')

triad.res = Bartz.experiment(G.triad.samples, theta.triad.actual, theta0.triad, verbose=F, model='triad')


# Calculate MSE -----------------------------------------------------------

mse = function(sim, obs) {
  return(sum((sim - obs)^2) / length(sim))
}

get.pred.thetas = function(res) {
  t(do.call(cbind, lapply(res, function(l) {
    # get last parameter for each theta chain
    l[[length(l)]]
  })))
}

ET.thetas = get.pred.thetas(ET.res)
triad.thetas = get.pred.thetas(triad.res)


