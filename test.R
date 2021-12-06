source("weight_matrix.R")
source("ipw-point-estimates.R")
source("integrand.R")
source("utils.R")

# test dataset
library(igraph)
g = sample_gnm(n = 20, m = 10) # 20 vertices, 40 edges
plot(g, vertex = 6)
G = rep(1:4, each = 5) # 4 groups
A = sample(c(0,1),20,replace = TRUE)
Y = (1:20)/2
H = h_neighborhood(g,Y,1)

allocations = c(0.1,0.2)
w.matrix = wght_matrix(integrand,allocations, G, A)
ipw_point_estimates(H,G,A,w.matrix)


