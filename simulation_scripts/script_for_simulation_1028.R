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

# spillover effect model
aa1 = c()
aa1_est1 = c()
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
  df$interaction1 <- X_num[,1] * df$treated_neigh
  df$interaction2 <- X_num[,2] * df$treated_neigh
  
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
  
  point_estimates <- ipw_point_estimates_mixed_test5(H, G, A, w.matrix, Con_type = "No-con")
  a = ipw_m_variance_groups(w.matrix, point_estimates, effect_type ='contrast',
                              marginal = FALSE, allocation1 = allocations[1], 
                              allocation2 = allocations[1])
  aa1 = rbind(aa1, a)
  aa1_est1 = c(aa1_est1,
               point_estimates$outcomes$overall[1] - 
                 point_estimates$outcomes$overall[2])
}

sd(aa1$estimate)
#3.904736
mean(aa1$std.error)
#4.245152
mean(aa1$estimate)
#-12.69773
mean(aa1_est1)
#-(5 + 0.5*(7+9)) = -13


### 2. binary variable
# in this case, no need for neighinfo, just pass in H_M
aa2 = c()
bb2 = c()
cc2 = c()

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
  X2 <- sample(c("M", "F"), size = length(A), replace = TRUE)
  X4 <- sample(c("Y", "N"), size = length(A), replace = TRUE)
  X <- cbind(X2, X4)
  X_type <- c("C", "C")
  X_cat <- as.matrix(X[, X_type == "C"])
  
  df <- cbind.data.frame(A,G,X)
  df$treated_neigh <- h_neighsum(graph, A, 1) 
  df$interaction1 <- ifelse(X_cat[,1] == "M", 1, 0) * df$treated_neigh
  df$interaction2 <- ifelse(X_cat[,2] == "N", 1, 0) * df$treated_neigh
    
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
  
  point_estimates_n <- ipw_point_estimates_mixed_test4(H_M, G, A, w.matrix, Con_type = "No-con")
  point_estimates_g1 <- ipw_point_estimates_mixed_test4(H, G, A, w.matrix, 
                                                       X = X, x0 = x0, 
                                                       X_type = X_type, Con_type = "group")
  point_estimates_g2 <- ipw_point_estimates_mixed_test5(H, G, A, w.matrix, 
                                                       X = X, x0 = x0, 
                                                       X_type = X_type, Con_type = "group")
  
  a = ipw_m_variance(w.matrix, point_estimates_n, effect_type ='contrast',
                     marginal = FALSE, allocation1 = allocations[1], 
                     allocation2 = allocations[1])
  point_estimates_n$outcomes$overall
  
  b = ipw_m_variance(w.matrix, point_estimates_g1, effect_type ='contrast',
                     marginal = FALSE, allocation1 = allocations[1], 
                     allocation2 = allocations[1])
  
  c = ipw_m_variance_groups(w.matrix, point_estimates_g2, effect_type ='contrast',
                     marginal = FALSE, allocation1 = allocations[1], 
                     allocation2 = allocations[1])
  point_estimates_g1$outcomes$overall
  point_estimates_g2$outcomes$overall
  
  aa2 = rbind(aa2, a)
  bb2 = rbind(bb2, b)
  cc2 = rbind(cc2, c)
}

sd(aa2$estimate) # 6.058569
mean(aa2$std.error)# 5.103824
mean(aa2$estimate) # -20.76101
# - (5 + 7 + 9) = -21

sd(bb2$estimate) # 3.686918
mean(bb2$std.error) # 4.590228
mean(bb2$estimate) # -13.49831

sd(cc2$estimate) # 5.283076
mean(cc2$std.error) # 4.480577
mean(cc2$estimate) # -13.3053

# > sd(aa2$estimate) #4.945716
# [1] 6.06548
# > mean(aa2$std.error)# 5.058368
# [1] 5.419865
# > mean(aa2$estimate) # -20.28974
# [1] -21.22276
# > # - (5 + 7 + 9) = -21
#   > sd(bb2$estimate) #2.93329
# [1] 0.9345009
# > mean(bb2$std.error) # 1.330891
# [1] 0.909375
# > mean(bb2$estimate) #-12.77848
# [1] -6.504326
### TODO: some problem with the lmlist model??


# > sd(aa2$estimate) # 6.058569
# [1] 5.415669
# > mean(aa2$std.error)# 5.103824
# [1] 5.150655
# > mean(aa2$estimate) # -20.76101
# [1] -20.60806
# > sd(bb2$estimate) # 3.686918
# [1] 3.04871
# > mean(bb2$std.error) # 4.590228
# [1] 4.48398
# > mean(bb2$estimate) # -13.49831
# [1] -13.72745
# > sd(cc2$estimate) # 5.283076
# [1] 5.133854
# > mean(cc2$std.error) # 4.480577
# [1] 4.411645
# > mean(cc2$estimate) # -13.3053
# [1] -13.27939

### 3. continuous variable
aa3 = c()
bb3 = c()
est3g = c()
est3n = c()

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
  X3 <- rnorm(length(A),mean = 0.5, sd = 1)
  X1 <- apply(G_mat, 1, function(x) rnorm(1,mean = x/51, sd = 1))
  X <- cbind(X3, X1)
  X_type <- c("N", "N")
  X_num <- apply(as.matrix(X[, X_type == "N"]), 2, as.numeric)
  
  df <- cbind.data.frame(A,G,X)
  df$treated_neigh <- h_neighsum(graph, A, 1) 
  df$interaction1 <- X_num[,1] * df$treated_neigh
  df$interaction2 <- X_num[,2] * df$treated_neigh  
  
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
  
  a = ipw_regression_variance_neigh(H, w.matrix, point_estimates_n1, A, effect_type ='contrast', 
                                    marginal = FALSE, allocation1 = allocations[1], allocation2 = allocations[1], 
                                    neighinfo = neighinfo, x1_num = x1_num)
  
  point_estimates_n1$outcomes$overall
  point_estimates_n2$outcomes$overall
  
  b = ipw_regression_variance(H, w.matrix, point_estimates_g1, A, effect_type ='contrast',
                              marginal = FALSE, allocation1 = allocations[1], allocation2 = allocations[1],
                              X = X, x0 = x0, X_type = X_type)
 
  point_estimates_g1$outcomes$overall
  point_estimates_g2$outcomes$overall
  
  aa3 = rbind(aa3, a)
  bb3 = rbind(bb3, b)
  
  est3g = c(est3g,
            point_estimates_g2$outcomes$overall[1]-
              point_estimates_g2$outcomes$overall[2])
  est3n = c(est3n,
            point_estimates_n2$outcomes$overall[1]-
              point_estimates_n2$outcomes$overall[2])
}

sd(aa3$estimate) #4.865189
mean(aa3$std.error)# 4.278735
mean(aa3$estimate) # -7.550644
# - (5 + 0.1*7 + 0.2*9) = -7.5
mean(est3n) # -13.50154

sd(bb3$estimate) #3.487326
mean(bb3$std.error) # 3.240342
mean(bb3$estimate) #-13.09059
#-(5 + 0.5*(7+9)) = -13
mean(est3g) # -12.72239

### 4. mixed variable
aa4 = c()
bb4 = c()
est4g = c()
est4n = c()

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
  #X1 <- apply(G_mat, 1, function(x) rnorm(1,mean = x/51, sd = 1)) # the avg should be 0.5
  X1 <- rnorm(length(A),mean = 0.5, sd = 1)
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
  
  a = ipw_regression_variance_neigh(H_M, w.matrix, point_estimates_n, A, effect_type ='contrast', 
                                    marginal = FALSE, allocation1 = allocations[1], allocation2 = allocations[1], 
                                    neighinfo = neighinfo, x1_num = x1_num)
  aa4 = rbind(aa4, a)

  b = ipw_regression_variance(H, w.matrix, point_estimates_g, A, effect_type ='contrast',
                              marginal = FALSE, allocation1 = allocations[1], allocation2 = allocations[1],
                              X = X, x0 = x0, X_type = X_type)
  bb4 = rbind(bb4,b)
  
  point_estimates_n2 <- ipw_point_estimates_mixed_test5(H_M, G, A, w.matrix, 
                                                       neighinfo = neighinfo, x1 = x1, 
                                                       X_type = X_type,  Con_type = "neigh")
  point_estimates_g2 <- ipw_point_estimates_mixed_test5(H, G, A, w.matrix, 
                                                       X = X, x0 = x0, 
                                                       X_type = X_type, Con_type = "group")
  
  est4g = c(est4g,
            point_estimates_g2$outcomes$overall[1]-
              point_estimates_g2$outcomes$overall[2])
  est4n = c(est4n,
            point_estimates_n2$outcomes$overall[1]-
              point_estimates_n2$outcomes$overall[2])
}
sd(aa4$estimate) #3.907932
mean(aa4$std.error)# 3.338361
mean(aa4$estimate) # -7.562572
# - (5 + 0.1*7 + 0.2*9) = -7.5
sd(bb4$estimate) #6.40128
mean(bb4$std.error) # 4.167996
mean(bb4$estimate) #-13.31909
# - (5 + 0.5*7 + 0.5*9) = -13

# > sd(aa4$estimate) #3.907932
# [1] 3.662001
# > mean(aa4$std.error)# 3.338361
# [1] 3.251725
# > mean(aa4$estimate) # -7.562572
# [1] -7.842786
# > # - (5 + 0.1*7 + 0.2*9) = -7.5
#   > sd(bb4$estimate) #6.40128
# [1] 5.269664
# > mean(bb4$std.error) # 4.167996
# [1] 4.158461
# > mean(bb4$estimate) #-13.31909
# [1] -13.69261
mean(est4g)
mean(est4n)
