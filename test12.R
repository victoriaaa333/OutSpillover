# 2022.06.27
# 1. tests for multiple X and neighX

source("weight_matrix.R")
source("ipw_point_estimates(mixed variables).R")
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
dd1 = c()

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
  X1 <- apply(G_mat, 1, function(x) rnorm(1,mean = x/51, sd = 1)) # the avg should be 0.5
  X2 <- sample(c("M", "F"), size = length(A), replace = TRUE)
  X3 <- rnorm(length(A),mean = 0.5, sd = 1)
  X4 <- sample(c("Y", "N"), size = length(A), replace = TRUE)
  X <- cbind(X1, X2, X3)
  X_type <- c("N", "C", "N") # indicating whether the covariate is numerical or categorical
  x0 <- as.matrix(c( 0.1, "M", 0.1))
  x1 <- as.matrix(c( 0.1, "M", 0.1))
  #x1p <- as.matrix(c( 0.1, "F", 0.1))
  X_num <- apply(X[, X_type == "N"], 2, as.numeric)
  X_cat <- as.matrix(X[, X_type == "C"])
  # X <- cbind(X1, X2)
  # X_type <- c("N", "C")
  # x0 <- as.matrix(c(0.1,"F"))
  # x1 <- as.matrix(c(0.1,"F"))
  # X_num <- apply(as.matrix(X[, X_type == "N"]), 2, as.numeric)
  # 
  df <- cbind.data.frame(A,G,X)
  df$treated_neigh <- h_neighsum(graph, A, 1) 
  df$interaction1 <- X_num[,1] * df$treated_neigh
  df$interaction2 <- X_num[,2] * df$treated_neigh
  
  #############
  # 2. Outcome model
  a = 2; b = 5; c = 7; d = 9
  Y = apply(cbind(df$A, df$treated_neigh, df$interaction1, df$interaction2), 1, #X_num,
            function(x)  rnorm(1, mean = a*x[1] + b*x[2] + c*x[3] + d*x[4], sd = 1))  
  H = h_neighborhood(graph, Y, 1) 
  H_M =  h_neighborhood(graph, Y, 1, X_cat, "M") 
  df$Y = Y
  df$H = H
  df$H_M = H_M
  #TODO: h_neighsum, h_neighofneigh needed X?
  ## helper functions to running regression on neighX
  neighX = h_neighcov(graph, Y, 1, X, X_type, x1) # average of X for unit j's 1-order neighbor
  #h_neigh = h_neighsum(graph, A, 1, X, X_type, x1) # number of treated unit j's 1-order neighbor
  h_neigh = h_neighsum(graph, A, 1)
  neigh2_treated = h_neighofneigh(graph, A, 1, h_neigh, 
                                  X, X_type, x1) # average number of treated neighbors l for unit j's 1-order neighbor i
  neigh2_treated_neighX = h_neighofneigh_withcov(graph, A, 1, X, X_type, x1, h_neigh) # average number of treated neighbors l for unit j's 1-order neighbor i, times X_i
  
  neighinfo = list(neigh2_treated, neighX, neigh2_treated_neighX)
  names(neighinfo) <- c('neigh2_treated', 'neighX', 'neigh2_treated_neighX')
  
  
  # neighXp = h_neighcov(graph, Y, 1, X, X_type, x1p) # average of X for unit j's 1-order neighbor
  # h_neigh = h_neighsum(graph, A, 1)
  # neigh2_treatedp = h_neighofneigh(graph, A, 1, h_neigh, 
  #                                 X, X_type, x1p) # average number of treated neighbors l for unit j's 1-order neighbor i
  # neigh2_treated_neighXp = h_neighofneigh_withcov(graph, A, 1, X, X_type, x1p, h_neigh) # average number of treated neighbors l for unit j's 1-order neighbor i, times X_i
  # 
  # neighinfop = list(neigh2_treatedp, neighXp, neigh2_treated_neighXp)
  # names(neighinfop) <- c('neigh2_treated', 'neighX', 'neigh2_treated_neighX')
  # 
  #############
  # 3. calculate the point estimates and the variances (bootstrapped and analytical)
  allocations = list(c(0.5,denominator_alphas),c(0.4,denominator_alphas))
  w.matrix = wght_matrix(integrand, allocations, G, A, P)
  
  no_con = ipw_point_estimates_mixed(H, G, A, w.matrix)$outcomes$overall
  neigh_con = ipw_point_estimates_mixed(H_M, G, A, w.matrix, 
                                       neighinfo = neighinfo, x1 = x1, X_type = X_type)$outcomes$overall
  # neigh_conp = ipw_point_estimates_mixed(H, G, A, w.matrix, 
  #                                       neighinfo = neighinfop, x1 = x1p, X_type = X_type)$outcomes$overall
  group_con = ipw_point_estimates_mixed(H, G, A, w.matrix, 
                                       X = X, x0 = x0, X_type = X_type)$outcomes$overall
  
  a = no_con[2] - no_con[1]
  
  b = neigh_con[,2,] - neigh_con[,1,]
  
  c = group_con[,2,] - group_con[,1,]
  
 # d = neigh_conp[,2,] - neigh_conp[,1,]
  
  aa1 = rbind(aa1,a)
  bb1 = rbind(bb1,b)
  cc1 = rbind(cc1,c)
  #dd1 = rbind(dd1,d)
}

colMeans(aa1); colMeans(bb1); colMeans(cc1)
# 13.53382 6.816162 13.85307

# 5 + 7 * 0.5 + 9*0.5 = 13
# 5 + 7 * 0.1 + 9*0.1 = 6.6



# 11.86494 6.250981 11.57024

# colMeans(aa1); colMeans(bb1); colMeans(cc1)
# [1] 8.361157
# [1] 6.793112
# [1] 7.956626

# [1] 13.17149
# [1] 10.35134
# [1] 12.35708

# 0.1, M, 0.1
# [1] 12.21216
# [1] 11.15084
# [1] 12.35119
