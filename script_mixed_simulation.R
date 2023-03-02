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
source("mixed_effects.R")

library(igraph)
library(lme4)

aa4 = c()
bb4 = c()
est4g = c()
est4n = c()


aa5 = c()
for (i in 1:20) {
  ##########
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
  X1 <- sample(c("M", "F"), size = length(A), replace = TRUE)
  X2 <- rnorm(length(A),mean = 0.5, sd = 1)
  
  X <- cbind(X1, X2)
  X_type <- c("C", "N")
  x0 <- as.matrix(c("M", 0.1))
  x1 <- x0
  x1_num <- c(0.1)
  
  X_cat <- as.matrix(X[, X_type == "C"])
  X_num <- as.matrix(cbind(X2,  ifelse(X_cat[,1] == "M", 1, 0)))
  
  
  df <- cbind.data.frame(A,G,X)
  df$treated_neigh <- h_neighsum(graph, A, 1) 
  df$interaction1 <- cov_neighsum(graph, A, 1, X = X_num[,1]) #, X_cat, "M"
  df$interaction2 <- cov_neighsum(graph, A, 1, X = X_num[,2]) #, X_cat, "M"
  df$interaction3 <- X_num[,1] * df$treated_neigh
  df$interaction4 <- X_num[,2] * df$treated_neigh
  #xdf$interaction5 <- df$interaction1 * df$interaction3
  ##########
  # 2. Outcome model
  a = 2; b = 4; c = 7; d = 9; e = 10; f = 12
  Y = apply(cbind(df$A, df$treated_neigh, df$interaction1, df$interaction2, 
                  df$interaction3, df$interaction4), 1, #X_num,
            function(x)  rnorm(1, mean = a*x[1] + b*x[2] + c*x[3] + 
                                 d*x[4] + e*x[5] + f*x[6], sd = 1))  
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
  
  point_estimates <- ipw_point_estimates_mixed_test4(H_M, G, A, w.matrix, 
                                                       X = X, x0 = x0, neighinfo = neighinfo, x1= x1,
                                                       X_type = X_type, Con_type = "mixed")
  
  aa5 = c(aa5, point_estimates$outcomes$overall[1]-
                       point_estimates$outcomes$overall[2])
  # a = ipw_regression_variance_neigh(H_M, w.matrix, point_estimates_n, A, 
  #                                   effect_type ='contrast', marginal = FALSE, 
  #                                   allocation1 = allocations[1], allocation2 = allocations[1], 
  #                                   neighinfo = neighinfo, x1_num = x1_num)
  # b = ipw_regression_variance(H, w.matrix, point_estimates_g, A, 
  #                             effect_type ='contrast',marginal = FALSE, 
  #                             allocation1 = allocations[1], allocation2 = allocations[1],
  #                             X = X, x0 = x0, X_type = X_type)
  # aa4 = rbind(aa4, a)
  # bb4 = rbind(bb4, b)
  # 
  # point_estimates_n2 <- ipw_point_estimates_mixed_test5(H_M, G, A, w.matrix, 
  #                                                       neighinfo = neighinfo, x1 = x1, 
  #                                                       X_type = X_type,  Con_type = "neigh")
  # point_estimates_g2 <- ipw_point_estimates_mixed_test5(H, G, A, w.matrix, 
  #                                                       X = X, x0 = x0, 
  #                                                       X_type = X_type, Con_type = "group")
  # est4g = c(est4g,
  #           point_estimates_g2$outcomes$overall[1]-
  #             point_estimates_g2$outcomes$overall[2])
  # est4n = c(est4n,
  #           point_estimates_n2$outcomes$overall[1]-
  #             point_estimates_n2$outcomes$overall[2])
}

saveRDS(aa5, "../kaggle/working/mixed.RDS")

sd(aa4$estimate) # 5.651792
mean(aa4$std.error)# 4.805365
mean(aa4$estimate) # -13.62658
mean(est4n) # -14.18967
# - (5 + 0.5*7 + 0.5*9) = -13
saveRDS(aa4, "../kaggle/working/inf_mixed1.RDS")
saveRDS(est4n, "../kaggle/working/inf_mixed1_pest.RDS")

sd(bb4$estimate) #3.581574
mean(bb4$std.error) # 2.419118
mean(bb4$estimate) #-8.896595
mean(est4g) # -7.232731
# - (5 + 0.5*7 + 0.1*9) = -7.5
saveRDS(bb4, "../kaggle/working/inf_mixed2.RDS")
saveRDS(est4g, "../kaggle/working/inf_mixed2_pest.RDS")





