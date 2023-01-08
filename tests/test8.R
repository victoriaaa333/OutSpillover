# 04/08/22
# 1. added treated friends of friends to the outcome model


source("weight_matrix.R")
source("ipw_point_estimates.R")
source("integrand.R")
source("utils.R")
source("bootstrap_variance.R")
source("effects.R")
source("m_variance.R")
library(igraph)

######
# adding the number of treated friends of friends in the outcome model

aa = c()
bb = c()
cc = c()

for (i in 1:100) {
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
  
  X = as.matrix(sample(c("M", "F"), length(clusters(graph)$membership), 
                       replace = TRUE, prob = c(0.5, 0.5)))
  
  df <- cbind.data.frame(A,G,X)
  df$treated_male_neigh <- h_neighsum(graph, A, 1, X = X, x1 = "M")
  df$treated_female_neigh <- h_neighsum(graph, A, 1, X = X, x1 = "F")
  df$treated_neigh <- h_neighsum(graph, A, 1, X = X)
  df$treated_neighofneigh <- h_neighsum(graph, A, 2, X = X)
  df$X_numeric <- ifelse(df$X == "F", 1, 0)
  df$interaction <- df$X_numeric * df$treated_neigh
  
  
  #############
  # 2. Outcome model (Y, H, conditional H)
  a = 2; b = 1; c = 5; d = 7; e = 10
  Y = apply(cbind(df$A, df$X_numeric, df$treated_neigh, df$interaction, df$treated_neighofneigh), 1, 
            function(x)  rnorm(1, mean = a*x[1] + b*x[2] + c*x[3] + d*x[4] + e*x[5], sd = 1))  
  H = h_neighborhood(graph, Y, 1) 
  H_F = h_neighborhood(graph, Y, 1, X, x1 = "F") 
  H_M = h_neighborhood(graph, Y, 1, X, x1 = "M") 
  df$Y = Y
  df$H = H
  df$H_F = H_F
  df$H_M = H_M
  
  #############
  # 3. calculate the point estimates and the variances (bootstrapped and analytical)
  allocations = list(c(0.5,denominator_alphas))
  w.matrix = wght_matrix(integrand, allocations, G, A, P)
  
  point_estimates  <- ipw_point_estimates(H, G, A, w.matrix)
  point_estimatesF  <- ipw_point_estimates(H_F, G, A, w.matrix)
  point_estimatesM  <- ipw_point_estimates(H_M, G, A, w.matrix)
  point_estimates1  <- ipw_point_estimates(H, G, A, w.matrix, X = X, x0 = "F")
  point_estimates2  <- ipw_point_estimates(H, G, A, w.matrix, X = X, x0 = "M")
  
  aa <- rbind(aa,
              c(point_estimates$outcomes$overall[2] -  point_estimates$outcomes$overall[1],
                point_estimatesF$outcomes$overall[2] -  point_estimatesF$outcomes$overall[1],
                point_estimatesM$outcomes$overall[2] -  point_estimatesM$outcomes$overall[1],
                point_estimates1$outcomes$overall[2] -  point_estimates1$outcomes$overall[1],
                point_estimates2$outcomes$overall[2] -  point_estimates2$outcomes$overall[1])
  )
  
  alphas   <- dimnames(w.matrix)[[length(dim(w.matrix))]]
  allocation1 <- alphas[1]
  m_var <- ipw_effect_calc(w.matrix, point_estimates, effect_type ='contrast', 
                           marginal = FALSE, allocation1, allocation2 = allocation1)[2][[1]]
  m_varF <- ipw_effect_calc(w.matrix, point_estimatesF, effect_type ='contrast', 
                            marginal = FALSE, allocation1, allocation2 = allocation1)[2][[1]]
  m_varM <- ipw_effect_calc(w.matrix, point_estimatesM, effect_type ='contrast', 
                            marginal = FALSE, allocation1, allocation2 = allocation1)[2][[1]]
  m_var1 <- ipw_effect_calc(w.matrix, point_estimates1, effect_type ='contrast', 
                            marginal = FALSE, allocation1, allocation2 = allocation1)[2][[1]]
  m_var2 <- ipw_effect_calc(w.matrix, point_estimates2, effect_type ='contrast', 
                            marginal = FALSE, allocation1, allocation2 = allocation1)[2][[1]]
  
  bb = rbind(bb,
             c(m_var, m_varF, m_varM, m_var1, m_var2))
  
  bootstrap_avg <- BootVar(df, 0.5, denominator_alphas, "H")
  boot_var <- sd(apply(bootstrap_avg, 3, function(x) x[2] - x[1]))
  bootstrap_avgF <- BootVar(df, 0.5, denominator_alphas, "H_F")
  boot_varF <- sd(apply(bootstrap_avgF, 3, function(x) x[2] - x[1]))
  bootstrap_avgM <- BootVar(df, 0.5, denominator_alphas, "H_M")
  boot_varM <- sd(apply(bootstrap_avgM, 3, function(x) x[2] - x[1]))
  bootstrap_avg1 <- BootVar(df, 0.5, denominator_alphas, "H", x0 = "F")
  boot_var1 <- sd(apply(bootstrap_avg1, 3, function(x) x[2] - x[1]))
  bootstrap_avg2 <- BootVar(df, 0.5, denominator_alphas, "H", x0 = "M")
  boot_var2 <- sd(apply(bootstrap_avg2, 3, function(x) x[2] - x[1]))
  cc = rbind(cc,
             c(boot_var, boot_varF, boot_varM, boot_var1, boot_var2))
}


########
# multiple covariates


aa = c()
bb = c()
cc = c()

for (i in 1:100) {
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
  
  X = as.matrix(cbind(sample(c("M", "F"), length(clusters(graph)$membership), 
                       replace = TRUE, prob = c(0.5, 0.5)),
                      sample(c("Y", "N"), length(clusters(graph)$membership), 
                             replace = TRUE, prob = c(0.5, 0.5))))
  
  df <- cbind.data.frame(A,G,X)
  df$treated_neigh <- h_neighsum(graph, A, 1, X = X)
  df$X_numeric <- ifelse(X[,1] == "F", 1, 0)
  df$interaction <- df$X_numeric * df$treated_neigh
  
  a = 2; b = 1; c = 5; d = 7
  Y = apply(cbind(df$A, df$X_numeric, df$treated_neigh, df$interaction), 1, 
            function(x)  rnorm(1, mean = a*x[1] + b*x[2] + c*x[3] + d*x[4], sd = 1))  
  
  H = h_neighborhood(graph, Y, 1) 
  H_F = h_neighborhood(graph, Y, 1, X, x1 = c("F","Y")) 
  H_M = h_neighborhood(graph, Y, 1, X, x1 = c("M","Y")) 
  df$Y = Y
  df$H = H
  df$H_F = H_F
  df$H_M = H_M
  
  #############
  # 3. calculate the point estimates and the variances (bootstrapped and analytical)
  allocations = list(c(0.5,denominator_alphas))
  w.matrix = wght_matrix(integrand, allocations, G, A, P)
  
  point_estimates  <- ipw_point_estimates(H, G, A, w.matrix)
  point_estimatesF  <- ipw_point_estimates(H_F, G, A, w.matrix)
  point_estimatesM  <- ipw_point_estimates(H_M, G, A, w.matrix)
  point_estimates1  <- ipw_point_estimates(H, G, A, w.matrix, X = X, x0 = "F")
  point_estimates2  <- ipw_point_estimates(H, G, A, w.matrix, X = X, x0 = "M")
  
  aa <- rbind(aa,
              c(point_estimates$outcomes$overall[2] -  point_estimates$outcomes$overall[1],
                point_estimatesF$outcomes$overall[2] -  point_estimatesF$outcomes$overall[1],
                point_estimatesM$outcomes$overall[2] -  point_estimatesM$outcomes$overall[1],
                point_estimates1$outcomes$overall[2] -  point_estimates1$outcomes$overall[1],
                point_estimates2$outcomes$overall[2] -  point_estimates2$outcomes$overall[1])
  )
  
  alphas   <- dimnames(w.matrix)[[length(dim(w.matrix))]]
  allocation1 <- alphas[1]
  m_var <- ipw_effect_calc(w.matrix, point_estimates, effect_type ='contrast', 
                           marginal = FALSE, allocation1, allocation2 = allocation1)[2][[1]]
  m_varF <- ipw_effect_calc(w.matrix, point_estimatesF, effect_type ='contrast', 
                            marginal = FALSE, allocation1, allocation2 = allocation1)[2][[1]]
  m_varM <- ipw_effect_calc(w.matrix, point_estimatesM, effect_type ='contrast', 
                            marginal = FALSE, allocation1, allocation2 = allocation1)[2][[1]]
  m_var1 <- ipw_effect_calc(w.matrix, point_estimates1, effect_type ='contrast', 
                            marginal = FALSE, allocation1, allocation2 = allocation1)[2][[1]]
  m_var2 <- ipw_effect_calc(w.matrix, point_estimates2, effect_type ='contrast', 
                            marginal = FALSE, allocation1, allocation2 = allocation1)[2][[1]]
  
  bb = rbind(bb,
             c(m_var, m_varF, m_varM, m_var1, m_var2))
  
  bootstrap_avg <- BootVar(df, 0.5, denominator_alphas, "H")
  boot_var <- sd(apply(bootstrap_avg, 3, function(x) x[2] - x[1]))
  bootstrap_avgF <- BootVar(df, 0.5, denominator_alphas, "H_F")
  boot_varF <- sd(apply(bootstrap_avgF, 3, function(x) x[2] - x[1]))
  bootstrap_avgM <- BootVar(df, 0.5, denominator_alphas, "H_M")
  boot_varM <- sd(apply(bootstrap_avgM, 3, function(x) x[2] - x[1]))
  bootstrap_avg1 <- BootVar(df, 0.5, denominator_alphas, "H", x0 = "F")
  boot_var1 <- sd(apply(bootstrap_avg1, 3, function(x) x[2] - x[1]))
  bootstrap_avg2 <- BootVar(df, 0.5, denominator_alphas, "H", x0 = "M")
  boot_var2 <- sd(apply(bootstrap_avg2, 3, function(x) x[2] - x[1]))
  cc = rbind(cc,
             c(boot_var, boot_varF, boot_varM, boot_var1, boot_var2))
}

