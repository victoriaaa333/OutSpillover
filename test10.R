# 05/15/22
#######
# 1. how to define the outcome model so that the influencer effect is heterogeneous 
#      (given that we dont have terms like treated females)
# continuous variable

source("weight_matrix.R")
source("ipw_point_estimates(continuous).R")
source("integrand.R")
source("utils.R")
source("bootstrap_variance.R")
source("effects.R")
source("m_variance.R")
library(igraph)
library(lme4)

aa1 = c()
bb1 = c()
cc1 = c()


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
  
  G_mat = as.matrix(G)
  X <- as.matrix(apply(G_mat, 1, function(x) rnorm(1,mean = x/51, sd = 1))) # the avg should be 0.5
  
  df <- cbind.data.frame(A,G,X)
  df$treated_neigh <- h_neighsum(graph, A, 1, X = X)
  df$interaction <- cov_neighsum(graph, A, 1, X = X)
  
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
  aa1 = c(aa1,a)
  bb1 = c(bb1,b)
  cc1 = c(cc1,c)
}

mean(aa1); mean(bb1); mean(cc1)
# 11.44207; 8.582145; 5.72222
# 5 + (0.9/ 0.5/ 0.1)* 7 = 11.3/ 8.5/ 5.7


#######
# 2. define H if spillover effect is heterogeneous with a continuous X

aa2 = c()
bb2 = c()
cc2 = c()

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
  
  G_mat = as.matrix(G)
  X <- as.matrix(apply(G_mat, 1, function(x) rnorm(1,mean = x/51, sd = 1))) # the avg should be 0.5
  
  df <- cbind.data.frame(A,G,X)
  df$treated_neigh <- h_neighsum(graph, A, 1, X = X)
  df$interaction <- df$X * df$treated_neigh
  
  #############
  # 2. Outcome model (Y, H, conditional H)
  a = 2; b = 1; c = 5; d = 7
  Y = apply(cbind(df$A, df$X, df$treated_neigh, df$interaction), 1, 
            function(x)  rnorm(1, mean = a*x[1] + b*x[2] + c*x[3] + d*x[4], sd = 1))  
  H = h_neighborhood(graph, Y, 1) 
  H_cont = h_neighborhood_cont(graph, Y, 1, X, x1 = 0.1)
  df$Y = Y
  df$H = H
  df$H_cont = H_cont
  
  #############
  # 3. calculate the point estimates and the variances (bootstrapped and analytical)
  allocations = list(c(0.5,denominator_alphas))
  w.matrix = wght_matrix(integrand, allocations, G, A, P)
  
  a = ipw_point_estimates_cont(H, G, A, w.matrix, X = X)$outcomes$overall[2] -
    ipw_point_estimates_cont(H, G, A, w.matrix, X = X)$outcomes$overall[1]
  b = ipw_point_estimates_cont(H_cont, G, A, w.matrix, X = X)$outcomes$overall[2] -
    ipw_point_estimates_cont(H_cont, G, A, w.matrix, X = X)$outcomes$overall[1]
  c = ipw_point_estimates_cont(H, G, A, w.matrix, X = X, x0 = 0.1)$outcomes$overall[2] -
    ipw_point_estimates_cont(H, G, A, w.matrix, X = X, x0 = 0.1)$outcomes$overall[1]
  aa2 = c(aa2,a)
  bb2 = c(bb2,b)
  cc2 = c(cc2,c)
}

#write.csv(aa2, "heter_H(regular).csv")
#write.csv(bb2, "heter_H(H_cont).csv")
#write.csv(cc2, "heter_H(x0 = 0.1).csv")

mean(aa2); mean(bb2); mean(cc2)
# 8.343125; 5.786585;  8.262785
# 5 + 0.5 * 7, 5 + 0.1 *7, 5 + 0.5 * 7

