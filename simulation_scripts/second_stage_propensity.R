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

library(igraph)
library(lme4)
library(boot)

result = c()
for (i in 1:50) {
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
  
  # A ～ X1 + X2 + P_2[i] + intercept
  for (i in 1:length(P1)) {
    # Base denominator assigned to each group
    if(P1[i] == 0.4){
      P2[i] = inv.logit(X1_num[i] + X2[i] + group_ef[i] -2)
    }else{
      P2[i] = inv.logit(X1_num[i] + X2[i] + group_ef[i] -1)
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
  df$Y = Y
  df$H = H
  df$H_M = H_M
  
  # 2.1 neighinfo
  neighX = h_neighcov(graph, 1, X, X_type, x1) # average of X for unit j's 1-order neighbor
  neighinfo = list(neighX)
  names(neighinfo) <- c('neighX')
  
  df$X1_num <- X1_num
  df$X2 <- X2
  formula <- H | A  ~  X1_num + X2 + (1|G) | G
  df1 <- df[P1 == 0.4,]
  df2 <- df[P1 == 0.6,]
  
  X <-cbind(1, X1_num, X2)
  allocations = list(c(numerator_alpha, denominator_alphas))

# error handling
  try({
    parameters1 <- unlist(propensity_parameter(formula, df1)[1])
    parameters2 <- unlist(propensity_parameter(formula, df2)[1])
  
    parameters = c(parameters1, parameters2)
    obj <- ipw_propensity_variance_second(parameters,
                                   allocations,
                                   causal_estimation_options = 
                                     list(variance_estimation = 'robust'),
                                   integrate_allocation = FALSE,
                                   H = H, X = X, P = P,
                                   A = A, G = G, effect_type = "contrast")
    result = rbind(result, obj)
    }, silent = TRUE)
}

sd(result$estimate)
mean(result$std.error)
