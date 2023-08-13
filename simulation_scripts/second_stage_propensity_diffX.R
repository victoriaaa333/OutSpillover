setwd("~/OutSpillover")

source("utils/weight_matrix.R")
source("utils/integrand.R")
source("utils/utils.R")
source("utils/mixed_effects.R")

source("point_estimates/point_estimates.R")
source("variances/bootstrap_variance.R")
source("variances/regression_variance.R")
source("variances/regression_utils.R")
source("variances/regression_utils_neigh.R")

source("propensity/propensity.R")
source("propensity/point_estimate_propensity.R")
source("propensity/propensity_utils.R")
source("propensity/propensity_effect_calc.R")
source("point_estimates/point_estimates_utils.R")
source("propensity/second_stage_utils.R")

library(igraph)
library(lme4)
library(boot)

result = c()
result_G = c()
result_N = c()

result_G1 = c()
result_G2 = c()

for (i in 1:200) {
  
  # library(foreach)
  # library(doRNG)
  # library(doMC)
  # registerDoMC(min(detectCores() - 1, 6))
  
  library(parallel)
  library(MASS)
  numCores <- detectCores()
  
  #result <- foreach(i = 1:10, .combine="c") %dorng% {
  
  #trials <- seq(1, 10)
  #boot_fx <- function(trial) {
  ##########
  #1. Generate a graph and dataset (treatments, covariates)
  graph = make_empty_graph(n = 0, directed = FALSE)
  noc = 50
  ss = 100
  
  repeat{
    g2 = sample_gnp(ss, 0.5, directed = FALSE, loops = FALSE)
    graph = disjoint_union(graph, g2)
    if (clusters(graph)$no == noc){
      break}
  }
  G = rep(1:noc, each = ss)
  
  numerator_alpha = 0.5
  denominator_alphas = c(0.4,0.6)
  P = c(0.5,0.5) 
  # Assume random effect ~ N(0, 0.2^2) for each group
  ranef <- rnorm(noc, mean = 0, sd = 0.2)
  group_ef<- rep(ranef, each = ss)
  
  # assign to group 1 (0.4) and group 2 (0.6) randomly
  G1 = sample(c(1,2), length(unique(G)), replace = TRUE, prob = P)
  G1_alpha = denominator_alphas[G1]
  P1 = sapply(G, function(x) G1_alpha[x])
  
  # for each group with assigned treatment probability, the probability of treatment 
  P2 = rep(NA, length(P1))
  A = rep(NA, length(P1))
  # Assume there are 2 covariates
  X1 <- sample(c("M", "F"), size = length(A), replace = TRUE) # M is 1, F is 0
  X1_num <- ifelse(X1 == "M", 1, 0)
  X2 <- rnorm(length(A), mean = 0.5, sd = 1)
  X <- cbind(X1, X2)
  X_type <- c("C", "N")
  
  
  # different cov for conditional X
  X3 <- sample(c("Y", "N"), size = length(A), replace = TRUE)
  X4 <- rnorm(length(A), mean = 0.5, sd = 1)
  # A ～ X1 + X2 + P_2[i] + intercept
  
  
  for (i in 1:length(P1)) {
    # Base denominator assigned to each group
    if(P1[i] == 0.4){
      # mean of X1_num, X2, group_ef = .5, .5, 0
      # intercept = logit(0.4) - 1
      P2[i] = inv.logit(X1_num[i] + X2[i] + group_ef[i] - 1.4 - max(ranef[G1 == 1]))
    }else{
      P2[i] = inv.logit(X1_num[i] + X2[i] + group_ef[i] - 0.6 - max(ranef[G1 == 2]))
    }
  }
  
  df <- as.data.frame((cbind(G, P2)))
  
  for (i in 1:dim(df)[1]) {
    df$A[i] <- sample(c(1, 0), size = 1, prob = c(P2[i], 1-P2[i]))
  }
  A <- df$A
  
  x0 <- as.matrix(c("M", 0.1))
  x1 <- x0
  x1_num <- c(0.1)
  
  X_num <- apply(as.matrix(X[, X_type == "N"]), 2, as.numeric)
  X_cat <- as.matrix(X[, X_type == "C"])
  
  df$treated_neigh <- h_neighsum(graph, df$A, 1) 
  df$interaction1 <- ifelse(X_cat[,1] == "M", 1, 0)  * df$treated_neigh
  df$interaction2 <- X_num[,1] * df$treated_neigh
  
  ##########
  # 2. Outcome model
  aa = 2; bb = 5; cc = 7; dd = 9
  Y = apply(cbind(df$A, df$treated_neigh, df$interaction1, df$interaction2), 1, function(x)  rnorm(1, mean = aa*x[1] + bb*x[2] + cc*x[3] + dd*x[4], sd = 1))  
  H = h_neighborhood(graph, Y, 1) 
  H_M =  h_neighborhood(graph, Y, 1, X_cat, c("M")) 
  H_Y =  h_neighborhood(graph, Y, 1, as.matrix(X3), c("Y")) 
  df$Y = Y
  df$H = H
  df$H_M = H_M
  
  # 2.1 neighinfo
  neighX = h_neighcov(graph, 1, X, X_type, x1) # average of X for unit j's 1-order neighbor
  neighX_diff = h_neighcov(graph, 1,  cbind(X3, X4), X_type, as.matrix(c("Y", 0.1)))
  neighinfo = list(neighX)
  names(neighinfo) <- c('neighX')
  neighinfo_diff = list(neighX_diff)
  names(neighinfo_diff) <- c('neighX')
  
  df$X1_num <- X1_num
  df$X2 <- X2
  formula <- H | A  ~  X1_num + X2 + (1|G) | G
  df1 <- df[P1 == 0.4,]
  df2 <- df[P1 == 0.6,]
  
  X <- cbind(1, X1_num, X2)
  allocations = list(c(numerator_alpha, denominator_alphas))
  first_assignments = G1 - 1
  
  # error handling
  try({
    parameters1 <- unlist(propensity_parameter(formula, df1)[1])
    parameters2 <- unlist(propensity_parameter(formula, df2)[1])
    
    parameters = c(parameters1, parameters2)
    
    # no condition
    obj <- ipw_propensity_variance_second(parameters = parameters,
                                          allocations = allocations,
                                          causal_estimation_options = 
                                            list(variance_estimation = 'robust'),
                                          integrate_allocation = FALSE,
                                          H = H, propensity_X = X, P = P,
                                          A = A, G = G, first_assignments = first_assignments,
                                          effect_type = "contrast")
    result = rbind(result, obj)
    
    # conditional on group parameters (X3, X4)
    cond_X <- cbind(X3, X4)
    obj_G <- ipw_propensity_variance_second(parameters = parameters,
                                            allocations = allocations,
                                            causal_estimation_options = 
                                              list(variance_estimation = 'robust'),
                                            integrate_allocation = FALSE,
                                            propensity_X = X,
                                            H = H, P = P,
                                            A = A, G = G, 
                                            first_assignments = first_assignments,
                                            effect_type = "contrast", 
                                            X = cond_X,
                                            x0 = as.matrix(c("Y", 0.1)), # x0
                                            X_type = c("C", "N"),
                                            Con_type = "group")
    
    obj_G1 <- ipw_propensity_variance_second(parameters = parameters,
                                             allocations = allocations,
                                             causal_estimation_options = 
                                               list(variance_estimation = 'robust'),
                                             integrate_allocation = FALSE,
                                             propensity_X = X,
                                             H = H, P = P,
                                             A = A, G = G, 
                                             first_assignments = first_assignments,
                                             effect_type = "contrast", 
                                             X = cbind(X4),
                                             x0 = as.matrix(c(0.1)), 
                                             X_type = c("N"),
                                             Con_type = "group")
    
    obj_G2 <- ipw_propensity_variance_second(parameters = parameters,
                                             allocations = allocations,
                                             causal_estimation_options = 
                                               list(variance_estimation = 'robust'),
                                             integrate_allocation = FALSE,
                                             propensity_X = X,
                                             H = H, P = P,
                                             A = df$A, G = df$G, 
                                             first_assignments = first_assignments,
                                             effect_type = "contrast", 
                                             X = cbind(X3),
                                             x0 = as.matrix(c("Y")), 
                                             X_type = c("C"),
                                             Con_type = "group")
    
    result_G = rbind(result_G, obj_G)
    result_G1 = rbind(result_G1, obj_G1)
    result_G2 = rbind(result_G2, obj_G2)
    
    # conditional on neighbor parameters, remember to use H_M
    obj_N <- ipw_propensity_variance_second(parameters = parameters,
                                            allocations = allocations,
                                            causal_estimation_options = 
                                              list(variance_estimation = 'robust'),
                                            integrate_allocation = FALSE,
                                            H = H_Y, propensity_X = X, P = P, 
                                            A = A, G = G, 
                                            first_assignments = first_assignments,
                                            effect_type = "contrast", 
                                            neighinfo = neighinfo_diff,
                                            x1 = as.matrix(x1_num), # x1
                                            X_type = c("N"), # X_type
                                            Con_type = "neigh")
    result_N = rbind(result_N, obj_N)
  }, silent = TRUE)
  
  #output = list(list(nocon = obj, inf = obj_G, sp = obj_N))
}

# system.time({
#   results <- mclapply(trials, boot_fx, mc.cores = numCores)
#   # results <- lapply(trials, boot_fx)
# })

saveRDS(result, "cluster_results/second_propensity_nocon_diffX 1.RDS")
saveRDS(result_G, "cluster_results/second_propensity_group_diffX 1.RDS")
saveRDS(result_N, "cluster_results/second_propensity_neigh_diffX 1.RDS")
saveRDS(result_G1, "cluster_results/second_propensity_group1_diffX 1.RDS")
saveRDS(result_G2, "cluster_results/second_propensity_group2_diffX 1.RDS")
# (new intercept)

sd(result$estimate)
mean(result$std.error)