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

source("propensity/propensity_weighted_regression.R")

library(igraph)
library(lme4)
library(boot)

result = c()
result_G = c()
result_N = c()

result_G1 = c()
result_G2 = c()

for (i in 1:100) {
  
  ##########
  #1. Generate a graph and dataset (treatments, covariates)
  graph = make_empty_graph(n = 0, directed = FALSE)
  noc = 100
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
  group_ef <- rep(ranef, each = ss)
  
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
  
  # A ～ X1 + X2 + P_2[i] + intercept
  for (i in 1:length(P1)) {
    # Base denominator assigned to each group
    if(P1[i] == 0.4){
      P2[i] = inv.logit(X1_num[i] + X2[i] + group_ef[i] - 1.4)  
    }else{
      P2[i] = inv.logit(X1_num[i] + X2[i] + group_ef[i] - 0.6) 
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
  df$interaction1 <- cov_neighsum(graph, A, 1, X = ifelse(X_cat[,1] == "M", 1, 0))
  df$interaction2 <- cov_neighsum(graph, A, 1, X = X_num[,1])
  
  ##########
  # 2. Outcome model
  aa = 2; bb = 5; cc = 7; dd = 9
  Y = apply(cbind(df$A, df$treated_neigh, df$interaction1, df$interaction2), 1, function(x)  rnorm(1, mean = aa*x[1] + bb*x[2] + cc*x[3] + dd*x[4], sd = 1))  
  H = h_neighborhood(graph, Y, 1) 
  H_M =  h_neighborhood(graph, Y, 1, X_cat, c("M")) 
  df$Y = Y
  df$H = H
  df$H_M = H_M
  
  # 2.1 neighinfo
  neighX = h_neighcov(graph, 1, X, X_type, x1) # average of X for unit j's 1-order neighbor
  neighinfo = list(neighX)
  names(neighinfo) <- c('neighX')
  
  df$X1_num <- X1_num
  df$X2 <- X2
  df1 <- df[P1 == 0.4,]
  df2 <- df[P1 == 0.6,]
  
  X <- cbind(1, X1_num, X2)
  allocations = list(c(numerator_alpha, denominator_alphas))
  first_assignments = G1 - 1
  
  # error handling
  try({
    formula <- H | A  ~  X1_num + X2 + (1|G) | G
    parameters1 <- unlist(propensity_parameter(formula, df1)[1])
    parameters2 <- unlist(propensity_parameter(formula, df2)[1])
    
    parameters = c(parameters1, parameters2)
    
    # no condition
    # obj <- ipw_propensity_variance_second(parameters = parameters,
    #                                       allocations = allocations,
    #                                       causal_estimation_options = 
    #                                         list(variance_estimation = 'robust'),
    #                                       integrate_allocation = FALSE,
    #                                       H = H, propensity_X = X, P = P,
    #                                       A = A, G = G, first_assignments = first_assignments,
    #                                       effect_type = "outcome")
    
    #result = rbind(result, obj)
    
    # conditional on group parameters
    cond_X <- cbind(X1, X2)
    obj_G <- ipw_propensity_variance_regression(parameters = parameters,
                                                allocations = allocations,
                                                causal_estimation_options = 
                                                  list(variance_estimation = 'robust'),
                                                integrate_allocation = FALSE,
                                                propensity_X = X,
                                                H = H, P = P,
                                                A = A, G = G, 
                                                first_assignments = first_assignments,
                                                effect_type = "outcome", 
                                                X = cond_X,
                                                x0 = x0,
                                                X_type = c("C", "N"),
                                                Con_type = "group")
    result_G = rbind(result_G, obj_G)
    
    obj_G1 <- ipw_propensity_variance_regression(parameters = parameters,
                                                 allocations = allocations,
                                                 causal_estimation_options =
                                                   list(variance_estimation = 'robust'),
                                                 integrate_allocation = FALSE,
                                                 propensity_X = X,
                                                 H = H, P = P,
                                                 A = A, G = G,
                                                 first_assignments = first_assignments,
                                                 effect_type = "outcome",
                                                 X = cbind(X2),
                                                 x0 = as.matrix(c(0.1)),
                                                 X_type = c("N"),
                                                 Con_type = "group")
    result_G1 = rbind(result_G1, obj_G1)
    
    # obj_G2 <- ipw_propensity_variance_second(parameters = parameters,
    #                                       allocations = allocations,
    #                                       causal_estimation_options = 
    #                                         list(variance_estimation = 'robust'),
    #                                       integrate_allocation = FALSE,
    #                                       H = H, propensity_X = X, P = P,
    #                                       A = A, G = G, first_assignments = first_assignments,
    #                                       effect_type = "contrast",
    #                                       X = cbind(X1),
    #                                       x0 = as.matrix(c("M")),
    #                                       X_type = c("C"),
    #                                       Con_type = "group")
    # 
    
    # weight_args <- append(append(integrand_args, integrate_args),
    #                       list(integrand   = propensity_integrand,
    #                            allocations = allocations,
    #                            propensity_X = X, A = A, G = G, P = P,
    #                            parameters = parameters,
    #                            runSilent  = TRUE, #BB 2015-06-23
    #                            integrate_allocation = FALSE
    #                       ))
    # weightd <- do.call(wght_deriv_array_second, args = append(weight_args, grad_args))
    # weights <- do.call(wght_matrix_second, args = weight_args)
    # pest <- ipw_point_estimates_propensity(H = H_M, G = G, A = A, weights, X = NULL, x0 = NULL, 
    #                                neighinfo = neighinfo, x1= as.matrix(x1_num), X_type = c("N"),
    #                                Con_type = "neigh")
    
    formula <- H_M | A  ~  X1_num + X2 + (1|G) | G
    parameters1 <- unlist(propensity_parameter(formula, df1)[1])
    parameters2 <- unlist(propensity_parameter(formula, df2)[1])
    parameters = c(parameters1, parameters2)
    
    obj_N <- ipw_propensity_variance_regression(parameters = parameters,
                                                allocations = allocations,
                                                causal_estimation_options =
                                                  list(variance_estimation = 'robust'),
                                                integrate_allocation = FALSE,
                                                H = H_M, propensity_X = X, P = P,
                                                A = A, G = G,
                                                first_assignments = first_assignments,
                                                effect_type = "outcome",
                                                neighinfo = neighinfo,
                                                x1 = as.matrix(x1_num), # x1
                                                X_type = c("N"), # X_type
                                                Con_type = "neigh")
    result_N = rbind(result_N, obj_N)
  }, silent = TRUE)
}

saveRDS(result_G, "cluster_results/second_propensity_outcome_group.RDS")
#saveRDS(result, "cluster_results/second_propensity_outcome_nocon.RDS")
saveRDS(result_N, "cluster_results/second_propensity_outcome_neigh.RDS")
saveRDS(result_G1, "cluster_results/second_propensity_outcome_group1.RDS")
