source("weight_matrix.R")
source("ipw_point_estimate_tests.R")
source("integrand.R")
source("utils.R")
source("bootstrap_variance.R")
source("effects.R")
source("m_variance.R")
source("regression_variance.R")
source("regression_utils.R")
source("regression_utils_neigh.R")

library(igraph)
library(lme4)

aa = c()
bb = c()
#############
for (i in 1:100) {
  #1. Generate a graph and dataset (treatments, covariates)
  graph = make_empty_graph(n = 0, directed = FALSE)
  repeat{
    g2 = sample_gnp(100, 0.5, directed = FALSE, loops = FALSE)
    graph = disjoint_union(graph, g2)
    if (clusters(graph)$no == 50){
      break}
  }
  G = components(graph)$membership
  
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
  
  X <- cbind(X2, X3, X1)
  X_type <- c("C", "N", "N")
  x0 <- as.matrix(c("M", 0.1, 0.2))
  x1 <- x0
  x1_num <- c(0.1, 0.2)
  
  X_num <- apply(as.matrix(X[, X_type == "N"]), 2, as.numeric)
  X_cat <- as.matrix(X[, X_type == "C"])
  
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
  H_M =  h_neighborhood(graph, Y, 1, X_cat, c("M")) 
  df$Y = Y
  df$H = H
  df$H_M = H_M
  
  # 2.1 neighinfo
  neighX = h_neighcov(graph, 1, X, X_type, x1) # average of X for unit j's 1-order neighbor
  neighinfo = list(neighX)
  names(neighinfo) <- c('neighX')
  
  #############
  # 3. calculate the point estimates and the variances (bootstrapped and analytical)
  allocations = list(c(0.5,denominator_alphas))
  w.matrix = wght_matrix(plain_integrand, allocations, G, A, P)
  
  point_estimates_n <- ipw_point_estimates_mixed_test4(H_M, G, A, w.matrix, 
                                                       neighinfo = neighinfo, x1 = x1, 
                                                       X_type = X_type,  Con_type = "neigh")
  point_estimates_g <- ipw_point_estimates_mixed_test4(H, G, A, w.matrix, 
                                                       X = X, x0 = x0, 
                                                       X_type = X_type, Con_type = "group")
  
  a = ipw_regression_variance(H, w.matrix, point_estimates_g, effect_type ='contrast',
                          marginal = FALSE, allocation1 = allocations[1], allocation2 = allocations[1],
                          X = X, x0 = x0, X_type = X_type)
  # point_estimates_g$outcomes$overall
  # a = ipw_regression_variance(H, w.matrix, point_estimates_g, effect_type ='outcome', 
  #                             marginal = FALSE, allocation1 = allocations[1], allocation2 = allocations[1], 
  #                             X = X, x0 = x0, X_type = X_type)
  
  aa = rbind(aa, a)
  b = ipw_regression_variance_neigh(H_M, w.matrix, point_estimates_n, effect_type ='contrast', 
                                marginal = FALSE, allocation1 = allocations[1], allocation2 = allocations[1], 
                                neighinfo = neighinfo, x1_num = x1_num)
  #point_estimates_n$outcomes$overall
  bb = rbind(bb,b)
}
sd(aa$estimate)
mean(aa$std.error)

sd(bb$estimate)
mean(bb$std.error)

# > sd(aa$estimate)
# [1] 5.339242
# > mean(aa$std.error)
# [1] 3.688199
# > sd(bb$estimate)
# [1] 3.436127
# > mean(bb$std.error)
# [1] 3.220899
