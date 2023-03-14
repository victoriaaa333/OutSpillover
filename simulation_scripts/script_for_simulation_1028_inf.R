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

# influencer effect model
aa1 = c()
### 1. plain case

for (i in 1:20) {
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
  df$interaction1 <- cov_neighsum(graph, A, 1, X = X_num[,1])
  df$interaction2 <- cov_neighsum(graph, A, 1, X = X_num[,2])
  
  
  #############
  # 2. Spillover model
  a = 2; b = 5; c = 7; d = 9
  Y = apply(cbind(df$A, df$treated_neigh, df$interaction1, df$interaction2), 1, #X_num,
            function(x)  rnorm(1, mean = a*x[1] + b*x[2] + c*x[3] + d*x[4], sd = 1))  
  H = h_neighborhood(graph, Y, 1) 
  df$Y = Y
  df$H = H
  
  #############
  # 3. calculate the point estimates and the variances (bootstrapped and analytical)
  allocations = list(c(0.5,denominator_alphas))
  w.matrix = wght_matrix(plain_integrand, allocations, G, A, P)
  
  point_estimates1 <- ipw_point_estimates_mixed_test4(H, G, A, w.matrix, Con_type = "No-con")
  #point_estimates2 <- ipw_point_estimates_mixed_test5(H, G, A, w.matrix, Con_type = "No-con")
  
  a = ipw_m_variance(w.matrix, point_estimates1, effect_type ='contrast',
                     marginal = FALSE, allocation1 = allocations[1], 
                     allocation2 = allocations[1])
  
  aa1 = rbind(aa1, a)
}

sd(aa1$estimate)
#3.177133
mean(aa1$std.error)
#3.363055
mean(aa1$estimate)
#-12.84796
#-(5 + 0.5*(7+9)) = -13


### 2. binary variable
# in this case, no need for neighinfo, just pass in H_M
aa21 = c()
bb21 = c()
aa22 = c()
bb22 = c()

for (i in 1:20) {
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
  X4 <- sample(c("Y", "N"), size = length(A), replace = TRUE)
  X <- cbind(X2, X4)
  X_type <- c("C", "C")
  X_cat <- as.matrix(X[, X_type == "C"])
  X_num_1 <- ifelse(X_cat[,1] == "M", 1, 0)
  X_num_2 <- ifelse(X_cat[,2] == "N", 1, 0)
  X_num <- as.matrix(cbind(X_num_1, X_num_2))
  
  df <- cbind.data.frame(A,G,X)
  df$treated_neigh <- h_neighsum(graph, A, 1) 
  df$interaction1 <- cov_neighsum(graph, A, 1, X_num[,1])
  df$interaction2 <- cov_neighsum(graph, A, 1, X_num[,2])
  
  #############
  # 2. Spillover model
  a = 2; b = 5; c = 7; d = 9
  Y = apply(cbind(df$A, df$treated_neigh, df$interaction1, df$interaction2), 1, #X_num,
            function(x)  rnorm(1, mean = a*x[1] + b*x[2] + c*x[3] + d*x[4], sd = 1))  
  H = h_neighborhood(graph, Y, 1) 
  H_M =  h_neighborhood(graph, Y, 1, X_cat, c("M","N")) 
  df$Y = Y
  df$H = H
  df$H_M = H_M
  x0 <- as.matrix(c("M", "N"))
  x1 <- x0
  
  #############
  # 3. calculate the point estimates and the variances (bootstrapped and analytical)
  allocations = list(c(0.5,denominator_alphas))
  w.matrix = wght_matrix(plain_integrand, allocations, G, A, P)
  
  point_estimates_n1 <- ipw_point_estimates_mixed_test4(H_M, G, A, w.matrix, Con_type = "No-con")
  point_estimates_g1 <- ipw_point_estimates_mixed_test4(H, G, A, w.matrix, 
                                                       X = X, x0 = x0, 
                                                       X_type = X_type, Con_type = "group")
  
  point_estimates_n2 <- ipw_point_estimates_mixed_test5(H_M, G, A, w.matrix, Con_type = "No-con")
  point_estimates_g2 <- ipw_point_estimates_mixed_test5(H, G, A, w.matrix, 
                                                       X = X, x0 = x0, 
                                                       X_type = X_type, Con_type = "group")
  
  a1 = ipw_m_variance(w.matrix, point_estimates_n1, effect_type ='contrast',
                     marginal = FALSE, allocation1 = allocations[1], 
                     allocation2 = allocations[1])
  a2 = ipw_m_variance_groups(w.matrix, point_estimates_n2, effect_type ='contrast',
                     marginal = FALSE, allocation1 = allocations[1], 
                     allocation2 = allocations[1])
  point_estimates_n1$outcomes$overall
  point_estimates_n2$outcomes$overall
  
  b1 = ipw_m_variance(w.matrix, point_estimates_g1, effect_type ='contrast',
                     marginal = FALSE, allocation1 = allocations[1], 
                     allocation2 = allocations[1])
  b2 = ipw_m_variance_groups(w.matrix, point_estimates_g2, effect_type ='contrast',
                     marginal = FALSE, allocation1 = allocations[1], 
                     allocation2 = allocations[1])
  point_estimates_g1$outcomes$overall
  point_estimates_g2$outcomes$overall
  # the above is equal to  colMeans(point_estimates_g2$outcomes$groups[,,,1])
  
  # ipw_effect_calc(w.matrix, point_estimates, effect_type ='outcome', 
  #                 marginal = TRUE, allocation1, allocation2)[2][[1]])
  
  aa21 = rbind(aa21, a1)
  bb21 = rbind(bb21, b1)
  aa22 = rbind(aa22, a2)
  bb22 = rbind(bb22, b2)
}

sd(aa21$estimate) #3.704291
mean(aa21$std.error)# 3.45546
mean(aa21$estimate) # -13.92418
# - (5 + 0.5*7 + 0.5*9) = -13
sd(bb21$estimate) #2.803222
mean(bb21$std.error) # 1.425992
mean(bb21$estimate) #-20.40872
# - (5 + 7 + 9) = -21

sd(aa22$estimate) #3.704291
mean(aa22$std.error)# 3.45546
mean(aa22$estimate) # -13.92418
# - (5 + 0.5*7 + 0.5*9) = -13
sd(bb22$estimate) #2.803222
mean(bb22$std.error) # 1.425992
mean(bb22$estimate) #-20.40872
# - (5 + 7 + *9) = -21

#  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 
# 12 12 15  9  6 23 17 16 11 13  7 12 14  8 15 10  6 
# 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 
# 11  6 13 17 12  8 19 15  8  8 20 15 12 24 14 18 21 
# 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 
# 15  9 10  6 17 20 15 17 11 17 19 13 16  7 12  9 

# > sd(aa21$estimate) #3.704291
# [1] 3.892123
# > mean(aa21$std.error)# 3.45546
# [1] 3.684921
# > mean(aa21$estimate) # -13.92418
# [1] -13.854
# > # - (5 + 0.5*7 + 0.5*9) = -13
#   > sd(bb21$estimate) #2.803222
# [1] 2.870147
# > mean(bb21$std.error) # 1.425992
# [1] 5.28311
# > mean(bb21$estimate) #-20.40872
# [1] -19.36857
# > sd(aa22$estimate) #3.704291
# [1] 3.892123
# > mean(aa22$std.error)# 3.45546
# [1] 3.684921
# > mean(aa22$estimate) # -13.92418
# [1] -13.854
# > # - (5 + 0.5*7 + 0.5*9) = -13
#   > sd(bb22$estimate) #2.803222
# [1] 5.02214
# > mean(bb22$std.error) # 1.425992
# [1] 5.205775
# > mean(bb22$estimate) #-20.40872
# [1] -21.35115

### 3. continuous variable
aa31 = c()
aa32 = c()
bb31 = c()
bb32 = c()
est3g = c()
est3n = c()

for (i in 1:20) {
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
  X3 <- rnorm(length(A),mean = 0.5, sd = 1)
  X1 <- rnorm(length(A),mean = 0.5, sd = 0.5)
  #X1 <- apply(G_mat, 1, function(x) rnorm(1,mean = x/51, sd = 1))
  X <- cbind(X3, X1)
  X_type <- c("N", "N")
  X_num <- apply(as.matrix(X[, X_type == "N"]), 2, as.numeric)
  
  df <- cbind.data.frame(A,G,X)
  df$treated_neigh <- h_neighsum(graph, A, 1) 
  df$interaction1 <- cov_neighsum(graph, A, 1, X = X_num[,1])
  df$interaction2 <- cov_neighsum(graph, A, 1, X = X_num[,2])
  #df$interaction <- h_neighcov(graph, 1, X, X_type, x1)
  
  #############
  # 2. Spillover model
  a = 2; b = 5; c = 7; d = 9
  Y = apply(cbind(df$A, df$treated_neigh, df$interaction1, df$interaction2), 1, #X_num,
            function(x)  rnorm(1, mean = a*x[1] + b*x[2] + c*x[3] + d*x[4], sd = 1))  
  H = h_neighborhood(graph, Y, 1) 
  df$Y = Y
  df$H = H
  x0 <- as.matrix(c(0.1, 0.2))
  x1 <- x0
  x1_num <- c(0.1, 0.2)
  
  neighX = h_neighcov(graph, 1, X, X_type, x1) # average of X for unit j's 1-order neighbor
  neighinfo = list(neighX)
  names(neighinfo) <- c('neighX')
  
  #############
  # 3. calculate the point estimates and the variances (bootstrapped and analytical)
  allocations = list(c(0.5,denominator_alphas))
  w.matrix = wght_matrix(plain_integrand, allocations, G, A, P)
  
  point_estimates_n1 <- ipw_point_estimates_mixed_test4(H, G, A, w.matrix, 
                                                       neighinfo = neighinfo, x1 = x1, 
                                                       X_type = X_type,  Con_type = "neigh")
  point_estimates_g1 <- ipw_point_estimates_mixed_test4(H, G, A, w.matrix, 
                                                       X = X, x0 = x0, 
                                                       X_type = X_type, Con_type = "group")
  
  point_estimates_n2 <- ipw_point_estimates_mixed_test5(H, G, A, w.matrix, 
                                                       neighinfo = neighinfo, x1 = x1, 
                                                       X_type = X_type,  Con_type = "neigh")
  point_estimates_g2 <- ipw_point_estimates_mixed_test5(H, G, A, w.matrix, 
                                                       X = X, x0 = x0, 
                                                       X_type = X_type, Con_type = "group")
  
  a1 = ipw_regression_variance_neigh(H, w.matrix, point_estimates_n1, A, effect_type ='contrast', 
                                    marginal = FALSE, allocation1 = allocations[1], allocation2 = allocations[1], 
                                    neighinfo = neighinfo, x1_num = x1_num)
  a2 = ipw_regression_variance_neigh(H, w.matrix, point_estimates_n2, A, effect_type ='contrast', 
                                     marginal = FALSE, allocation1 = allocations[1], allocation2 = allocations[1], 
                                     neighinfo = neighinfo, x1_num = x1_num)
  
  point_estimates_g1$outcomes$overall
  point_estimates_g2$outcomes$overall
  point_estimates_n1$outcomes$overall
  point_estimates_n2$outcomes$overall
  
  b1 = ipw_regression_variance(H, w.matrix, point_estimates_g1, A, effect_type ='contrast',
                              marginal = FALSE, allocation1 = allocations[1], allocation2 = allocations[1],
                              X = X, x0 = x0, X_type = X_type)
  b2 = ipw_regression_variance(H, w.matrix, point_estimates_g2, A, effect_type ='contrast',
                              marginal = FALSE, allocation1 = allocations[1], allocation2 = allocations[1],
                              X = X, x0 = x0, X_type = X_type)
  
  point_estimates_g2$outcomes$overall
  point_estimates_n2$outcomes$overall
  
  est3g = c(est3g,
            point_estimates_g2$outcomes$overall[1]-
              point_estimates_g2$outcomes$overall[2])
  est3n = c(est3n,
            point_estimates_n2$outcomes$overall[1]-
              point_estimates_n2$outcomes$overall[2])
  
  # ipw_effect_calc(w.matrix, point_estimates, effect_type ='outcome', 
  #                 marginal = TRUE, allocation1, allocation2)[2][[1]])
  
  aa31 = rbind(aa31, a1)
  bb31 = rbind(bb31, b1)
  aa32 = rbind(aa32, a2)
  bb32 = rbind(bb32, b2)
}

# neigh
sd(aa31$estimate) #6.712766
mean(aa31$std.error)# 6.281982
mean(aa31$estimate) # -11.77797
#-( 5 + 0.5*(7+9)) = -13

sd(aa32$estimate) #6.712766
mean(aa32$std.error)# 6.281982
mean(aa32$estimate) # -11.77797

mean(est3n)

# group
sd(bb31$estimate) # 2.943855
mean(bb31$std.error) # 2.761648
mean(bb31$estimate) # -7.72814
# - (5 + 0.1*7 + 0.2*9) = -7.5

sd(bb32$estimate) #6.562281
mean(bb32$std.error)# 6.400047
mean(bb32$estimate) # -12.66101

mean(est3g)

### 4. mixed variable
aa41 = c()
bb41 = c()
aa42 = c()
bb42 = c()

#############
for (i in 1:20) {
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
  df$interaction1 <- cov_neighsum(graph, A, 1, X = X_num[,1]) #, X_cat, "M"
  df$interaction2 <- cov_neighsum(graph, A, 1, X = X_num[,2]) #, X_cat, "M"
  
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
  
  point_estimates_n1 <- ipw_point_estimates_mixed_test4(H_M, G, A, w.matrix, 
                                                       neighinfo = neighinfo, x1 = x1, 
                                                       X_type = X_type,  Con_type = "neigh")
  point_estimates_g1 <- ipw_point_estimates_mixed_test4(H, G, A, w.matrix, 
                                                       X = X, x0 = x0, 
                                                       X_type = X_type, Con_type = "group")
  point_estimates_n2 <- ipw_point_estimates_mixed_test5(H_M, G, A, w.matrix, 
                                                        neighinfo = neighinfo, x1 = x1, 
                                                        X_type = X_type,  Con_type = "neigh")
  point_estimates_g2 <- ipw_point_estimates_mixed_test5(H, G, A, w.matrix, 
                                                        X = X, x0 = x0, 
                                                        X_type = X_type, Con_type = "group")
  
  a1 = ipw_regression_variance_neigh(H_M, w.matrix, point_estimates_n1, A, effect_type ='contrast', 
                                    marginal = FALSE, allocation1 = allocations[1], allocation2 = allocations[1], 
                                    neighinfo = neighinfo, x1_num = x1_num)
  a2 = ipw_regression_variance_neigh(H_M, w.matrix, point_estimates_n2, A, effect_type ='contrast', 
                                     marginal = FALSE, allocation1 = allocations[1], allocation2 = allocations[1], 
                                     neighinfo = neighinfo, x1_num = x1_num)
  
  aa41 = rbind(aa41, a1)
  #point_estimates_n$outcomes$overall
  b1 = ipw_regression_variance(H, w.matrix, point_estimates_g1, A, effect_type ='contrast',
                              marginal = FALSE, allocation1 = allocations[1], allocation2 = allocations[1],
                              X = X, x0 = x0, X_type = X_type)
  b2 = ipw_regression_variance(H, w.matrix, point_estimates_g2, A, effect_type ='contrast',
                              marginal = FALSE, allocation1 = allocations[1], allocation2 = allocations[1],
                              X = X, x0 = x0, X_type = X_type)
  
  bb4 = rbind(bb4,b)
}
sd(aa4$estimate) #5.835166
mean(aa4$std.error)# 5.224558
mean(aa4$estimate) # -12.02852
# - (5 + 0.5*7 + 0.5*9) = -13
sd(bb4$estimate) #5.719304
mean(bb4$std.error) # 3.683535
mean(bb4$estimate) #-7.306903
# - (5 + 0.1*7 + 0.2*9) = -7.5


# > sd(aa4$estimate) #5.835166
# [1] 6.346423
# > mean(aa4$std.error)# 5.224558
# [1] 5.152636
# > mean(aa4$estimate) # -12.02852
# [1] -12.23907
# > # - (5 + 0.5*7 + 0.5*9) = -13
#   > sd(bb4$estimate) #5.719304
# [1] 5.023914
# > mean(bb4$std.error) # 3.683535
# [1] 3.230766
# > mean(bb4$estimate) #-7.306903
# [1] -8.722717
# > 
