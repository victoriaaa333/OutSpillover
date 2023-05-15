source("utils/weight_matrix.R")
source("point_estimates/point_estimates.R")
source("utils/integrand.R")
source("utils/utils.R")
source("variances/bootstrap_variance.R")
source("variances/m_variance.R")
source("variances/regression_variance.R")
source("variances/regression_utils.R")
source("variances/regression_utils_neigh.R")
source("utils/mixed_effects.R")

library(igraph)
library(lme4)
library(foreach)
library(doMC)

result_c1 <- foreach(i = 1:500, .combine="c") %do% {
  
  ##########
  #1. Generate a graph and dataset (treatments, covariates)
  graph = make_empty_graph(n = 0, directed = FALSE)
  repeat{
    g2 = sample_gnp(100, 0.5, directed = FALSE, loops = FALSE)
    graph = disjoint_union(graph, g2)
    if (clusters(graph)$no == 200){
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
  X1 <- sample(c("M", "F"), size = length(A), replace = TRUE)
  X2 <- sample(c("Y", "N"), size = length(A), replace = TRUE)
  X <- cbind(X1, X2)
  X_type <- c("C", "C")
  
  x0 <- as.matrix(c("M", "Y"))
  x1 <- x0
  
  X_cat <- as.matrix(X[, X_type == "C"])
  X_num_1 <- ifelse(X_cat[,1] == "M", 1, 0)
  X_num_2 <- ifelse(X_cat[,2] == "Y", 1, 0)
  X_num <- as.matrix(cbind(X_num_1, X_num_2))
  
  df <- cbind.data.frame(A,G,X)
  df$treated_neigh <- h_neighsum(graph, A, 1) 
  df$interaction1 <- X_num[,1]  * df$treated_neigh
  df$interaction2 <- X_num[,2]  * df$treated_neigh
  
  ##########
  # 2. Outcome model
  a = 0; b = 1; c = 1; d = 2
  Y = apply(cbind(df$A, df$treated_neigh, df$interaction1, df$interaction2), 1, #X_num,
            function(x)  rnorm(1, mean = a*x[1] + b*x[2] + c*x[3] + d*x[4], sd = 1))  
  H = h_neighborhood(graph, Y, 1) 
  H_M =  h_neighborhood(graph, Y, 1, X_cat, c("M", "Y")) 
  df$Y = Y
  df$H = H
  df$H_M = H_M
  
  ##########
  # 3. calculate the point estimates and the variances (bootstrapped and analytical)
  allocations = list(c(0.5,denominator_alphas))
  w.matrix = wght_matrix(plain_integrand, allocations, G, A, P)
  
  point_estimates <- ipw_point_estimates_mixed_test4(H, G, A, w.matrix, 
                                                     Con_type = "No-con")
  
  point_estimates_g <- ipw_point_estimates_mixed_test5(H, G, A, w.matrix, 
                                                       X = X, x0 = x0, 
                                                       X_type = X_type, Con_type = "group")
  
  point_estimates_n <- ipw_point_estimates_mixed_test4(H_M, G, A, w.matrix, Con_type = "No-con")
  
  point_estimates_m <- ipw_point_estimates_mixed_test5(H_M, G, A, w.matrix, 
                                                       X = X, x0 = x0,
                                                       X_type = X_type, Con_type = "group")
  
  a = ipw_m_variance(w.matrix, point_estimates, effect_type ='contrast',
                     marginal = FALSE, allocation1 = allocations[1], 
                     allocation2 = allocations[1])  
  
  b = ipw_m_variance_groups(w.matrix, point_estimates_g, effect_type ='contrast',
                            marginal = FALSE, allocation1 = allocations[1], 
                            allocation2 = allocations[1])
  
  c = ipw_m_variance(w.matrix, point_estimates_n, effect_type ='contrast',
                     marginal = FALSE, allocation1 = allocations[1], 
                     allocation2 = allocations[1])
  
  d = ipw_m_variance_groups(w.matrix, point_estimates_m, effect_type ='contrast',
                            marginal = FALSE, allocation1 = allocations[1], 
                            allocation2 = allocations[1])
  
  output = list(list(nocon = a, inf = b, sp = c, mixed = d))
}

saveRDS(result_c1, "cluster_results/sp_model_cat_var.RDS")
