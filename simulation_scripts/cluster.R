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


#### influencer effect model ####

### 4. mixed variable ###
aa4 = c()
bb4 = c()
cc4 = c()

for (i in 1:500) {
  ##########
  #1. Generate a graph and dataset (treatments, covariates)
  graph = make_empty_graph(n = 0, directed = FALSE)
  repeat{
    g2 = sample_gnp(200, 0.5, directed = FALSE, loops = FALSE)
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
  X2 <- rnorm(length(A),mean = 0.5, sd = 0.1)
  
  X <- cbind(X1, X2)
  X_type <- c("C", "N")
  x0 <- as.matrix(c("M", 0.1))
  x1 <- x0
  x1_num <- c(0.1)
  
  X_cat <- as.matrix(X[, X_type == "C"])
  X_num <- as.matrix(cbind(ifelse(X_cat[,1] == "M", 1, 0), X2))
  
  
  df <- cbind.data.frame(A,G,X)
  df$treated_neigh <- h_neighsum(graph, A, 1) 
  df$interaction1 <- cov_neighsum(graph, A, 1, X = X_num[,1]) #, X_cat, "M"
  df$interaction2 <- cov_neighsum(graph, A, 1, X = X_num[,2]) #, X_cat, "M"
  
  ##########
  # 2. Outcome model
  a = 0; b = 3; c = 1; d = 2
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
  
  ##########
  # 3. calculate the point estimates and the variances (bootstrapped and analytical)
  allocations = list(c(0.5,denominator_alphas))
  w.matrix = wght_matrix(plain_integrand, allocations, G, A, P)
  
  point_estimates <- ipw_point_estimates_mixed_test4(H, G, A, w.matrix, 
                                                     Con_type = "No-con")
  point_estimates_n <- ipw_point_estimates_mixed_test4(H_M, G, A, w.matrix, 
                                                       neighinfo = neighinfo, x1 = x1, 
                                                       X_type = X_type,  Con_type = "neigh")
  point_estimates_g <- ipw_point_estimates_mixed_test4(H, G, A, w.matrix, 
                                                       X = X, x0 = x0, 
                                                       X_type = X_type, Con_type = "group")
  
  a = ipw_m_variance(w.matrix, point_estimates, effect_type ='contrast',
                     marginal = FALSE, allocation1 = allocations[1], 
                     allocation2 = allocations[1])  
  
  b = ipw_regression_variance(H, w.matrix, point_estimates_g, A, 
                              effect_type ='contrast',marginal = FALSE, 
                              allocation1 = allocations[1], allocation2 = allocations[1],
                              X = X, x0 = x0, X_type = X_type)
  
  c = ipw_regression_variance_neigh(H_M, w.matrix, point_estimates_n, A, 
                                    effect_type ='contrast', marginal = FALSE, 
                                    allocation1 = allocations[1], allocation2 = allocations[1], 
                                    neighinfo = neighinfo, x1_num = x1_num)
  
  aa4 = rbind(aa4, a)
  bb4 = rbind(bb4, b)
  cc4 = rbind(cc4, c)
 
}

#### results from influencer outcome model w/ mixed vars ####
# inf_model_mixed_var
results_nocon = as.data.frame(matrix(unlist(lapply(results, function(l) l$nocon)), nrow = 500, byrow = TRUE))
colnames(results_nocon) <- colnames(results[[1]]$nocon)
result_stats(results_nocon, 3 + 1 * 0.5 + 2 * 0.5)

results_inf = as.data.frame(matrix(unlist(lapply(results, function(l) l$inf)), nrow = 500, byrow = TRUE))
colnames(results_inf) <- colnames(results[[1]]$inf)
result_stats(results_inf, 3 + 1 + 2 * 0.1)

results_sp = as.data.frame(matrix(unlist(lapply(results, function(l) l$sp)), nrow = 500, byrow = TRUE))
colnames(results_sp) <- colnames(results[[1]]$sp)
result_stats(results_sp,  3 + 1 * 0.5 + 2 * 0.5)

results_mixed = as.data.frame(matrix(unlist(lapply(results, function(l) l$mixed)), nrow = 500, byrow = TRUE))
colnames(results_mixed) <- colnames(results[[1]]$mixed)
result_stats(results_mixed, 1.5 + 1 * 0.75 + 2 * 0.25)

