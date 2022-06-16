##############
# 2022.05.20
# 1. tests for new neighborhood conditional function
# 2. changed the result to multiple x1 and x0, 
   # calculate the coeeficient in ipw_point_estimates explicitly

source("weight_matrix.R")
source("ipw_point_estimates(continuous).R")
source("integrand.R")
source("utils.R")
source("bootstrap_variance.R")
source("effects.R")
source("m_variance.R")
library(igraph)
library(lme4)

############
# 1. heterogeneous spillover effects model.
# dependent on x1, not x0

aa1 = c()
bb1 = c()
cc1 = c()

for (i in 1:100){
  # model for spillover effect
  
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
  
  G_mat = as.matrix(G)
  X <- as.matrix(apply(G_mat, 1, function(x) rnorm(1,mean = x/51, sd = 1))) # the avg should be 0.5
  
  df <- cbind.data.frame(A,G,X)
  df$treated_neigh <- h_neighsum(graph, A, 1) #TODO: check this
  df$interaction <- df$X * df$treated_neigh
  
  #############
  # 2. Outcome model
  a = 2; b = 1; c = 5; d = 7 
  Y = apply(cbind(df$A, df$X, df$treated_neigh, df$interaction), 1, 
            function(x)  rnorm(1, mean = a*x[1] + b*x[2] + c*x[3] + d*x[4], sd = 1))  
  H = h_neighborhood(graph, Y, 1) 
  df$Y = Y
  df$H = H
  
  ## helper functions to running regression on neighX
  neighX = h_neighcov(graph, Y, X, 1) # average of X for unit j's 1-order neighbor
  h_neigh = h_neighsum(graph, A, 1) # number of treated unit j's 1-order neighbor
  neigh2_treated = h_neighofneigh(graph, A, 1, h_neigh) # average number of treated neighbors l for unit j's 1-order neighbor i
  neigh2_treated_neighX = h_neighofneigh_withcov(graph, A, 1, X, h_neigh) # average number of treated neighbors l for unit j's 1-order neighbor i, times X_i
  
  neighinfo = cbind(neighX, neigh2_treated, neigh2_treated_neighX)
  colnames(neighinfo) <- c('neighX', 'neigh2_treated', 'neigh2_treated_neighX')
  
  #############
  # 3. calculate the point estimates and the variances (bootstrapped and analytical)
  allocations = list(c(0.5,denominator_alphas))
  w.matrix = wght_matrix(integrand, allocations, G, A, P)
  
  no_con = ipw_point_estimates_cont(H, G, A, w.matrix)$outcomes$overall
  neigh_con = ipw_point_estimates_cont(H, G, A, w.matrix, 
                                       neighinfo = neighinfo, x1 = 0.1*(1:9))$outcomes$overall
  group_con = ipw_point_estimates_cont(H, G, A, w.matrix, 
                                       X = X, x0 = 0.1*(1:9))$outcomes$overall
  
  a = no_con[2] - no_con[1]
  
  b = neigh_con[,2,] - neigh_con[,1,]
  
  c = group_con[,2,] - group_con[,1,]
  
  aa1 = rbind(aa1,a)
  bb1 = rbind(bb1,b)
  cc1 = rbind(cc1,c)
}

colMeans(aa1); colMeans(bb1); colMeans(cc1)
plot(0.1*(1:9),  colMeans(bb1))

# 8.656924
# x1 = 0.1  x1 = 0.2  x1 = 0.3  x1 = 0.4  x1 = 0.5  x1 = 0.6  x1 = 0.7  x1 = 0.8  x1 = 0.9
# 5.677167  6.376542  7.075918  7.775294  8.474669  9.174045  9.873420 10.572796 11.272172
# x0 = 0.1  x0 = 0.2 x0 = 0.3 x0 = 0.4 x0 = 0.5 x0 = 0.6 x0 = 0.7 x0 = 0.8 x0 = 0.9 
# 8.719291 8.711149 8.703008 8.694866 8.686725 8.678584 8.670442 8.662301 8.654159 

# 8.872509 5.846718 8.902284
# should be 8.5, 5.7, 8.5


# TODO: Compare the variance of 2 different methods of getting conditional H
#       maybe versus the density of network.


############
# 2. heterogeneous influencer effects model.
# dependent on x0, not x1
aa2 = c()
bb2 = c()
cc2 = c()

for (i in 1:100){
  # model for influencer effect
  
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
  
  G_mat = as.matrix(G)
  X <- as.matrix(apply(G_mat, 1, function(x) rnorm(1,mean = x/51, sd = 1))) # the avg should be 0.5
  #X <-  as.matrix(rnorm(length(G),mean = 0.5, sd = 1))
  
  df <- cbind.data.frame(A,G,X)
  df$treated_neigh <- h_neighsum(graph, A, 1)
  df$interaction <- cov_neighsum(graph, A, 1, X = X) #the sum of the covariates in the TREATED units for each node's h-order neighborhood
  
  #############
  # 2. Outcome model
  a = 2; b = 1; c = 5; d = 7 
  Y = apply(cbind(df$A, df$X, df$treated_neigh, df$interaction), 1, 
            function(x)  rnorm(1, mean = a*x[1] + b*x[2] + c*x[3] + d*x[4], sd = 1))  
  H = h_neighborhood(graph, Y, 1) 
  df$Y = Y
  df$H = H
  
  ## helper functions to running regression on neighX
  # neighX = h_neighcov(graph, Y, X, 1) # average of X for unit j's 1-order neighbor
  # h_neigh = h_neighsum(graph, A, 1) # number of treated unit j's 1-order neighbor
  # neigh2_treated = h_neighofneigh(graph, A, 1, h_neigh) # average number of treated neighbors l for unit j's 1-order neighbor i
  # neigh2_treated_neighX = h_neighofneigh_withcov(graph, A, 1, X, h_neigh) # average number of treated neighbors l for unit j's 1-order neighbor i, times X_i
  # 
  # neighinfo = cbind(neighX, neigh2_treated, neigh2_treated_neighX)
  # colnames(neighinfo) <- c('neighX', 'neigh2_treated', 'neigh2_treated_neighX')
  
  #############
  # 3. calculate the point estimates and the variances (bootstrapped and analytical)
  allocations = list(c(0.5,denominator_alphas))
  w.matrix = wght_matrix(integrand, allocations, G, A, P)
  
  a = ipw_point_estimates_cont(H, G, A, w.matrix, X = X)$outcomes$overall[2] -
    ipw_point_estimates_cont(H, G, A, w.matrix, X = X)$outcomes$overall[1]
  
  # b = ipw_point_estimates_cont(H, G, A, w.matrix, 
  #                              neighinfo = neighinfo, x1 = 0.1)$outcomes$overall[2] -
  #   ipw_point_estimates_cont(H, G, A, w.matrix, 
  #                           neighinfo = neighinfo, x1 = 0.1)$outcomes$overall[1]
  
  c = ipw_point_estimates_cont(H, G, A, w.matrix, X = X, x0 = 0.1)$outcomes$overall[2] -
    ipw_point_estimates_cont(H, G, A, w.matrix, X = X, x0 = 0.1)$outcomes$overall[1]
  
  aa2 = c(aa2,a)
  #bb2 = c(bb2,b)
  cc2 = c(cc2,c)
}

mean(aa2); mean(cc2)
# 8.777794; 5.902689
