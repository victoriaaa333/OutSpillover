# test for the difference between point estimates with or without extra terms
# test for weights inside or outside the response

source("weight_matrix.R")
source("ipw_point_estimate_tests.R")
source("ipw_point_estimates_tests2.R")
#source("ipw_point_estimates_mixed.R")
source("integrand.R")
source("utils.R")
source("bootstrap_variance.R")
source("effects.R")
source("m_variance.R")
source("regression_variance.R")
source("regression_utils.R")
library(igraph)
library(lme4)

aa = c()
bb = c()
cc = c()
dd = c()
ee = c()
ff = c()

for (k in 1:100) {
  graph = make_empty_graph(n = 0, directed = FALSE)
  repeat{
    g2 = sample_gnp(100, 0.5, directed = FALSE, loops = FALSE)
    graph = disjoint_union(graph, g2)
    if (clusters(graph)$no == 50){
      break}
  }
  #G = group_vector(graph) 
  G = rep(1:50, each = 100)
  #############
  # 1. Two-stage randomization
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
  X <- cbind(X1, X2, X3)
  X_type <- c("N", "C", "N")
  x0 <- as.matrix(c( 0.1, "M", 0.1))
  x1 <- as.matrix(c( 0.1, "M", 0.1))
  # X <- cbind(X1, X3)
  # X_type <- c("N", "N") 
  # x0 <- as.matrix(c( 0.1, 0.1))
  # x1 <- as.matrix(c( 0.1, 0.1))
  
  X_num <- apply(X[, X_type == "N"], 2, as.numeric)
  X_cat <- as.matrix(X[, X_type == "C"])
  
  df <- cbind.data.frame(A,G,X)
  df$treated_neigh <- h_neighsum(graph, A, 1) 
  df$interaction1 <- X_num[,1] * df$treated_neigh
  df$interaction2 <- X_num[,2] * df$treated_neigh
  
  #############
  # 2. Outcome model
  #a = 2; b = 5; c = 7; d = 9
  Y = apply(cbind(df$A, df$treated_neigh, df$interaction1, df$interaction2), 1, #X_num,
            function(x)  rnorm(1, mean = 2*x[1] + 5*x[2] + 7*x[3] + 9*x[4], sd = 1))  
  H = h_neighborhood(graph, Y, 1) 
  H_M =  h_neighborhood(graph, Y, 1, X_cat, "M") 
  df$Y = Y
  df$H = H
  df$H_M = H_M
  neighX = h_neighcov(graph, 1, X, X_type, x1) # average of X for unit j's 1-order neighbor
  h_neigh = h_neighsum(graph, A, 1)
  neigh2_treated = h_neighofneigh(graph, A, 1, h_neigh, 
                                  X, X_type, x1) 
  # average number of treated neighbors l for unit j's 1-order neighbor i
  neigh2_treated_neighX = h_neighofneigh_withcov(graph, A, 1, X, X_type, x1, h_neigh) 
  # average number of treated neighbors l for unit j's 1-order neighbor i, times X_i
  
  neighinfo = list(neigh2_treated, neighX, neigh2_treated_neighX)
  names(neighinfo) <- c('neigh2_treated', 'neighX', 'neigh2_treated_neighX')
  
  #############
  # 3. calculate the point estimates and the variances (bootstrapped and analytical)
  allocations = list(c(0.5,denominator_alphas))
  w.matrix = wght_matrix(plain_integrand, allocations, G, A, P)
  
  no_con = ipw_point_estimates_mixed_test(H, G, A, w.matrix)$outcomes$overall
  #1. running a regression with weights in the outcomes, by group
  neigh_con = ipw_point_estimates_mixed_test(H_M, G, A, w.matrix,
                                             neighinfo = neighinfo, x1 = x1, 
                                             X_type = X_type, Con_type = "neigh")$outcomes$overall
  #2. running a regression with weights outside the outcomes, overall estimates
  neigh_con2 = ipw_point_estimates_mixed_test2(H_M, G, A, w.matrix,
                                  neighinfo = neighinfo, x1 = x1, 
                                  X_type = X_type, Con_type = "neigh")$outcomes$overall
  #3. running a regression with weights in the outcomes, by group, with no additional terms 
  neigh_con3 = ipw_point_estimates_mixed_test3(H_M, G, A, w.matrix,
                                               neighinfo = neighinfo, x1 = x1,
                                               X_type = X_type, Con_type = "neigh")$outcomes$overall
  #4. running a regression with weights outside the outcomes, with no additional terms, (both overall and group estimate work)
  neigh_con4 = ipw_point_estimates_mixed_test4(H_M, G, A, w.matrix,
                                               neighinfo = neighinfo, x1 = x1,
                                               X_type = X_type, Con_type = "neigh")$outcomes$overall
  neigh_con5 = ipw_point_estimates_mixed_test5(H_M, G, A, w.matrix,
                                               neighinfo = neighinfo, x1 = x1,
                                               X_type = X_type, Con_type = "neigh")$outcomes$overall
  
  # TODO: H_M must be align with x1? If that is the case we should calculate H within the function
  # group_con = ipw_point_estimates_mixed_test(H, G, A, w.matrix,
  #                                       X = X, x0 = x0, X_type = X_type,
  #                                       Con_type = "group")$outcomes$overall
  # 
  a = no_con[2] - no_con[1]
  b = neigh_con[,2,] - neigh_con[,1,]
  c = neigh_con2[,2,] - neigh_con2[,1,]
  d = neigh_con3[,2,] - neigh_con3[,1,]
  e = neigh_con4[,2,] - neigh_con4[,1,]
  f = neigh_con5[,2,] - neigh_con5[,1,]
  aa = rbind(aa,a)
  bb = rbind(bb,b)
  cc = rbind(cc,c)
  dd = rbind(dd,d)
  ee = rbind(ee,e)
  ff = rbind(ff,f)
}

mean(aa)
mean(bb)
mean(cc)
mean(dd)
mean(ee)
mean(ff)
# 5 + 7 * 0.5 + 9*0.5 = 13
# 5 + 7 * 0.1 + 9*0.1 = 6.6

# > mean(aa)
# [1] 13.66901
# > mean(bb)
# [1] 6.696426
# > mean(cc)
# [1] 19.4465
# > mean(dd)
# [1] 11.23212
# > mean(ee)
# [1] 6.57031

# bias
# 0.669
# 0.010
# 12.847
# 4.632
# -0.030

# When weights in the response, considering extra terms lead to correct model
# When weights pass in as a parameter in 'lm' function, not considering extra terms lead to correct model


# After simulation for 400 times:
# > mean(aa)
# [1] 13.33095 bias: 0.33
# > mean(bb)
# [1] 6.688555 bias: 0.08       
# > mean(cc)
# [1] 19.28595 bias: 12.68
# > mean(dd)
# [1] 11.91444 bias: 5.31
# > mean(ee)
# [1] 6.798934 bias: 0.14
# > mean(ff)
# [1] 13.04659 bias: 6.44

