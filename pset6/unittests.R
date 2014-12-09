source('ERGM.MCMC.R')
library(RUnit)

# 4 two-stars, 0 triangles
G1 = rbind(c(0,1,0,0,0),
           c(1,0,0,1,0),
           c(0,0,0,1,0),
           c(0,1,1,0,1),
           c(0,0,0,1,0))

# 7 two-stars, 1 triangle
G2 = rbind(c(0,1,0,0,0),
           c(1,0,1,1,0),
           c(0,1,0,1,0),
           c(0,1,1,0,1),
           c(0,0,0,1,0))

# 10 two-stars, 2 triangles
G3 = rbind(c(0,1,0,0,0),
           c(1,0,1,1,0),
           c(0,1,0,1,1),
           c(0,1,1,0,1),
           c(0,0,1,1,0))

# 15 two-stars, 4 triangles
G4 = rbind(c(0,1,0,0,0),
           c(1,0,1,1,1),
           c(0,1,0,1,1),
           c(0,1,1,0,1),
           c(0,1,1,1,0))

checkEquals(ERGM.twostars(G1), 4)
checkEquals(ERGM.twostars(G2), 7)
checkEquals(ERGM.twostars(G3), 10)

checkEquals(ERGM.triangles(G1), 0)
checkEquals(ERGM.triangles(G2), 1)
checkEquals(ERGM.triangles(G3), 2)

# add edges
checkEquals(ERGM.edges.diff(G1, c(2,3)), 1)

checkEquals(ERGM.twostars.diff(G1, c(2,3)), 3)
checkEquals(ERGM.twostars.diff(G2, c(5,3)), 3)
checkEquals(ERGM.twostars.diff(G3, c(2,5)), 5)

checkEquals(ERGM.triangles.diff(G1, c(2,3)), 1)
checkEquals(ERGM.triangles.diff(G2, c(5,3)), 1)
checkEquals(ERGM.triangles.diff(G3, c(2,5)), 2)

# remove edges
checkEquals(ERGM.edges.diff(G1, c(3,4)), -1)

checkEquals(ERGM.twostars.diff(G2, c(2,3)), -3)
checkEquals(ERGM.twostars.diff(G3, c(5,3)), -3)
checkEquals(ERGM.twostars.diff(G4, c(2,5)), -5)

checkEquals(ERGM.triangles.diff(G2, c(2,3)), -1)
checkEquals(ERGM.triangles.diff(G3, c(5,3)), -1)
checkEquals(ERGM.triangles.diff(G4, c(2,5)), -2)
