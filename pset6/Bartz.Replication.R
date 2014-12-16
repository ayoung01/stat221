source('ERGM.MCMC.R')

n.nodes = 20
n.samples = 1000
n.reps = 25

# ET ----------------------------------------------------------------------

# theta.ET.actual = as.matrix(rnorm(2, -1, 1))
# > theta.ET.actual
# [,1]
# [1,]  0.2775402
# [2,] -1.6501253
print('Actual theta: '%+%theta.ET.actual)
theta0.ET = as.matrix(c(0,0))
G.ET.samples = ERGM.generate.samples(n.nodes, n.samples, theta.ET.actual, use.pkg=T, model='ET')
# G.ET.samples.filtered = ERGM.ET.filter.degenerates(G.ET.samples)

ET.res = Bartz.experiment(G.ET.samples, theta.ET.actual, theta0.ET, verbose=F, n.reps, model='ET')
plot.var(ET.res, theta.ET.actual)

# triad -------------------------------------------------------------------

# theta.triad.actual = as.matrix(rnorm(3, -1, 1))
theta.triad.actual = c(-.5, 0.1, -0.25)
theta0.triad = as.matrix(c(0,0,0))
G.triad.samples = ERGM.generate.samples(n.nodes, n.samples, theta.triad.actual, use.pkg=T, model='triad')

triad.res = Bartz.experiment(G.triad.samples, theta.triad.actual, theta0.triad, verbose=F, n.reps, model='triad')
plot.var(triad.res, theta.triad.actual, model='triad', ylim=c(0,10))
