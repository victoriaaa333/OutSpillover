# 03/06/22
# 1. updated outcome models

source("weight_matrix.R")
source("ipw_point_estimates.R")
source("integrand.R")
source("utils.R")
source("bootstrap_variance.R")
source("effects.R")
source("m_variance.R")
library(igraph)

aa = c()
bb = c()
cc = c()

for (i in 1:10) {
  #############
  #1. Generate a graph and dataset (treatments, covariates )
  graph = make_empty_graph(n = 0, directed = FALSE)
  repeat{
    g2 = sample_gnp(100, 0.5, directed = FALSE, loops = FALSE)
    graph = disjoint_union(graph, g2)
    if (clusters(graph)$no == 50){
      break}
    }
  
  G = group_vector(graph) 
  
  denominator_alphas = 0.5; P = 1
  A = sample(c(0,1), length(clusters(graph)$membership), 
             replace = TRUE, prob = c(1 - denominator_alphas, denominator_alphas))
  X = as.matrix(sample(c("M", "F"), length(clusters(graph)$membership), 
                       replace = TRUE, prob = c(1 - denominator_alphas, denominator_alphas)))
  
  df <- cbind.data.frame(A,G,X)
  df$treated_male_neigh <- h_neighsum(graph, A, 1, X = X, x1 = "M")
  df$treated_female_neigh <- h_neighsum(graph, A, 1, X = X, x1 = "F")
  df$treated_neigh <- h_neighsum(graph, A, 1, X = X)
  df$X_numeric <- ifelse(df$X == "F", 1, 0)
  df$interaction <- df$X_numeric * df$treated_neigh
  
  
  #############
  # 2. Outcome model (Y, H, conditional H)
  a = 2; b = 1; c = 5; d = 7
  Y = apply(cbind(df$A, df$X_numeric, df$treated_neigh, df$interaction), 1, 
            function(x)  rnorm(1, mean = a*x[1] + b*x[2] + c*x[3] + d*x[4], sd = 1))  
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
  # bb <- rbind(bb,
  #             c(point_estimates$marginal_outcomes$overall,
  #               point_estimatesF$marginal_outcomes$overall, 
  #               point_estimatesM$marginal_outcomes$overall, 
  #               point_estimates1$marginal_outcomes$overall,
  #               point_estimates2$marginal_outcomes$overall)
  # )
  
}

aa = read.csv("aa.csv")
bb = read.csv("bb.csv")
cc = read.csv("cc.csv")

colMeans(aa)
# c = 5; d = 7
# 9.079188 12.772647  5.336617  8.278655 10.084902
# theoretical value:
# 8.5      12         5         8.5      8.5
colMeans(bb)
# 6.000006 8.442657 3.522500 8.513539 8.657067
colMeans(cc)
# 5.940747 8.408501 3.502335 8.458320 8.649868
apply(aa, 2, sd)
# 6.401325 9.039709 3.762130 9.182009 9.054657

# colMeans(aa)
# a = 2; b = 1; c = 5; d = 7
# c + 0.5d, c + d, c, c+0.5d, c+0.5d
# 7.299035 ((c+(c+d))/2) 10.32799 (c+d) 6.66149 (c)  7.988405 ((c+(c+d))/2 = 8.5) 
# b = 0 8.835030 12.544559  8.640701  9.236120
# 9.261170 12.961077  5.437511  8.422770  9.990158

# TODO
# run multiple time simulations
# check estimands and replacement the outcome models H -> mu
# calculate the average bias (average of point-estimates minus the truth)
# calculate the empirical SD of causal effect (point-estimates)
# within each simulation, calculate the bootstrapped variance and the analytical variance, and take average across simulations
# the average should be close to the empirical SD


# TODO
# 1. two-stage randomization 
# 2. simulate 500 times
