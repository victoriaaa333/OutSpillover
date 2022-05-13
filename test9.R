# 04/15/22
# 1. continuous variable

source("weight_matrix.R")
source("ipw_point_estimates(continuous).R")
source("integrand.R")
source("utils.R")
source("bootstrap_variance.R")
source("effects.R")
source("m_variance.R")
library(igraph)
library(lme4)

aa = c()
bb = c()
cc = c()

for (i in 1:100){

  #############
  #1. Generate a graph and dataset (treatments, covariates)
  graph = make_empty_graph(n = 0, directed = FALSE)
  repeat{
    g2 = sample_gnp(100, 0.5, directed = FALSE, loops = FALSE)
    graph = disjoint_union(graph, g2)
    if (clusters(graph)$no == 50){
      break}
  }
  
  G = group_vector(graph) 
  
  # Two-stage randomization
  numerator_alpha = 0.5
  denominator_alphas = c(0.4,0.6)
  P = c(0.5,0.5) 
  P_1 = sample(c(1,2), length(unique(G)), replace = TRUE, prob = P)
  P_1 = sapply(G, function(x) P_1[x])
  P_2 = rep(NA, length(P_1))
  A = rep(NA, length(P_1))
  for (i in 1:length(P_1)) {
    P_2[i] = denominator_alphas[P_1[i]]
    A[i] = sample(c(0,1), 1, prob = c(1-P_2[i],P_2[i]))
  }
  
  X = as.matrix(runif(length(clusters(graph)$membership)))
  
  df <- cbind.data.frame(A,G,X)
  df$treated_neigh <- h_neighsum(graph, A, 1, X = X)
  df$interaction <- df$X * df$treated_neigh
  # consider the interaction term (treated_neigh dependent x)
  # regress ?
  
  #############
  # 2. Outcome model (Y, H, conditional H)
  a = 2; b = 1; c = 5; d = 7
  Y = apply(cbind(df$A, df$X, df$treated_neigh, df$interaction), 1, 
            function(x)  rnorm(1, mean = a*x[1] + b*x[2] + c*x[3] + d*x[4], sd = 1))  
  H = h_neighborhood(graph, Y, 1) 
  df$Y = Y
  df$H = H
  
  #############
  # 3. calculate the point estimates and the variances (bootstrapped and analytical)
  allocations = list(c(0.5,denominator_alphas))
  w.matrix = wght_matrix(integrand, allocations, G, A, P)
  
  a = ipw_point_estimates_cont(H, G, A, w.matrix, X = X, x0 = 0.9)$outcomes$overall[2] -
      ipw_point_estimates_cont(H, G, A, w.matrix, X = X, x0 = 0.9)$outcomes$overall[1]
  b = ipw_point_estimates_cont(H, G, A, w.matrix, X = X, x0 = 0.5)$outcomes$overall[2] -
      ipw_point_estimates_cont(H, G, A, w.matrix, X = X, x0 = 0.5)$outcomes$overall[1]
  c = ipw_point_estimates_cont(H, G, A, w.matrix, X = X, x0 = 0.1)$outcomes$overall[2] -
      ipw_point_estimates_cont(H, G, A, w.matrix, X = X, x0 = 0.1)$outcomes$overall[1]
  aa = c(aa,a)
  bb = c(bb,b)
  cc = c(cc,c)
}

# write.csv(aa, "cont(x0 = 0.9)_(nozero).csv")
# write.csv(bb, "cont(x0 = 0.5)_(nozero).csv")
# write.csv(cc, "cont(x0 = 0.1)_(nozero).csv")

aa = read.csv("cont(x0 = 0.9)_(nozero).csv")
mean(aa$x) #8.714536

bb = read.csv("cont(x0 = 0.5)_(nozero).csv")
mean(bb$x) #8.621991

cc = read.csv("cont(x0 = 0.1)_(nozero).csv")
mean(cc$x) #8.529447


###
#TODO: 1. how to define H if spillover effect is heterogeneous with a continuous X
        # (in the outcome model of interaction between X and treated neighbors).
#      2. how to define the outcome model so that the influencer effect is heterogeneous 
        # (given that we dont have terms like treated females)
        # try average X among treated neighbors, meaning people with different X will influence the others differently.
#      3. show multiple covariates result next time


