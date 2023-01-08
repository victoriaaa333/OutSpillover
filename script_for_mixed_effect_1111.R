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

aa1 = c()
for (i in 1:30) {
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
  X2 <- sample(c("M", "F"), size = length(A), replace = TRUE)
  X3 <- rnorm(length(A),mean = 0.5, sd = 1)
  X1 <- apply(G_mat, 1, function(x) rnorm(1,mean = x/51, sd = 1)) # the avg should be 0.5
  
  X <- cbind(X2, X3, X1)
  X_type <- c("C", "N", "N")
  X_num <- apply(as.matrix(X[, X_type == "N"]), 2, as.numeric)
  
  df <- cbind.data.frame(A,G,X)
  df$treated_neigh <- h_neighsum(graph, A, 1) 
  df$interaction1 <- X_num[,1] * df$treated_neigh
  df$interaction2 <- cov_neighsum(graph, A, 1, X = X_num[,1])
  df$interaction3 <- df$interaction1 * df$interaction2
  
  #############
  # 2. Spillover model
  a = 2; b = 5; c = 7; d = 9; e = 11
  Y = apply(cbind(df$A, df$treated_neigh, df$interaction1, df$interaction2, df$interaction3), 1, #X_num,
            function(x)  rnorm(1, mean = a*x[1] + b*x[2] + c*x[3] + d*x[4] + e*x[5], sd = 1))  
  H = h_neighborhood(graph, Y, 1) 
  df$Y = Y
  df$H = H
  
  #############
  # 3. calculate the point estimates and the variances (bootstrapped and analytical)
  allocations = list(c(0.5,denominator_alphas))
  w.matrix = wght_matrix(plain_integrand, allocations, G, A, P)
  
  point_estimates <- ipw_point_estimates_mixed_test4(H, G, A, w.matrix, Con_type = "No-con")
  point_estimates$outcomes$overall
  
  a = ipw_m_variance(w.matrix, point_estimates, effect_type ='contrast',
                     marginal = FALSE, allocation1 = allocations[1], 
                     allocation2 = allocations[1])
  # ipw_effect_calc(w.matrix, point_estimates, effect_type ='outcome', 
  #                 marginal = TRUE, allocation1, allocation2)[2][[1]])
  
  aa1 = rbind(aa1, a)
}

sd(aa1$estimate)
mean(aa1$std.error)
mean(aa1$estimate)

# > sd(aa1$estimate)
# [1] 39.58468
# > mean(aa1$std.error)
# [1] 18.28958
# > mean(aa1$estimate)
# [1] -152.1459

5 + (7+9)/2




