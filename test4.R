# 01/29/22
# 1. plot out the variance and distribution
# 2. run multiple times for 2-stage
# 3. compute the variance for m-estimator

source("weight_matrix.R")
source("ipw_point_estimates.R")
source("integrand.R")
source("utils.R")
source("bootstrap_variance.R")
source("effects.R")
source("m_variance.R")

library(igraph)

twostage0 = c()
twostage_ipw = c()
mvar_ipw = c()
for (i in 1:10) {
  graph = make_empty_graph(n = 0, directed = FALSE)
  repeat{
    g2 = sample_gnp(100, 0.5, directed = FALSE, loops = FALSE)
    graph = disjoint_union(graph, g2)
    if (clusters(graph)$no == 50){
      break
    }
  }
  G = group_vector(graph) 
  
  # Two-stage when h = 1
  numerator_alpha = 0.5
  denominator_alphas = c(0.4,0.6)
  P = c(0.6,0.4)
  a = 5; b = 2

  #randomly two-stage assignment
  P_1 = sample(c(1,2), length(unique(G)), replace = TRUE, prob = P)
  P_1 = sapply(G, function(x) P_1[x])
  P_2 = rep(NA, length(P_1))
  
  A = rep(NA, length(P_1))
  for (i in 1:length(P_1)) {
    P_2[i] = denominator_alphas[P_1[i]]
    A[i] = sample(c(0,1), 1, prob = c(1-P_2[i],P_2[i]))
  }
  
  df <- cbind.data.frame(A,G)
  df$neighbor_sum <- h_neighsum(graph, A, 1)
  
  Y = apply(cbind(df$A,df$neighbor_sum), 1, 
            function(x)  rnorm(1, mean = a*x[1] +  b*x[2], sd = 1))  
  # TODO: interaction between x[2] and covariates; change x[2] to number of treated neighbors with certain characteristics.
  H = h_neighborhood(graph,Y,0) 
  df$Y = Y
  df$H = H
  #head(df)
  
  allocations = list(c(numerator_alpha,denominator_alphas))
  w.matrix = wght_matrix(integrand, allocations, G, A, P)
  
  point_estimates  <- ipw_point_estimates(H, G, A, w.matrix)
  twostage_ipw <- c(twostage_ipw,
                    point_estimates$marginal_outcomes$overall)
  
  alphas   <- dimnames(w.matrix)[[length(dim(w.matrix))]]
  allocation1 <- alphas[1]
  allocation2 <- NA
  
  # twostage0 = c(twostage0,
  #               population_direct_effect(ipw_point_estimates(df$H, df$G, df$A,w.matrix)))
  mvar_ipw <- c(mvar_ipw,
                ipw_effect_calc(w.matrix, point_estimates, effect_type ='outcome', 
                                marginal = FALSE, allocation1, allocation2)[2][[1]])
                
}
mean(mvar_ipw) #12.06922
sd(twostage_ipw) #12.55968


onestage_ipw = c()
mvar_ipw = c()
for (i in 1:10) {
  graph = make_empty_graph(n = 0, directed = FALSE)
  repeat{
    g2 = sample_gnp(100, 0.5, directed = FALSE, loops = FALSE)
    graph = disjoint_union(graph, g2)
    if (clusters(graph)$no == 50){
      break
    }
  }
  G = group_vector(graph) 
  
  denominator_alphas = 0.5
  P = 1
  A = sample(c(0,1), length(clusters(graph)$membership), 
             replace = TRUE, prob = c(1 - denominator_alphas, denominator_alphas))
  
  df <- cbind.data.frame(A,G)
  df$neighbor_sum <- h_neighsum(graph, A, 1)
  
  a = 5; b = 2
  Y = apply(cbind(df$A,df$neighbor_sum), 1, 
            function(x)  rnorm(1, mean = a*x[1] +  b*x[2], sd = 1))  
  
  H = h_neighborhood(graph,Y,1) 
  df$Y = Y
  df$H = H
  head(df)
  
  allocations = list(c(0.5,denominator_alphas))
  w.matrix = wght_matrix(integrand, allocations, G, A, P)
  
  na.group = rep(NA, length(unique(df$G)))
  na.group[df[!is.na(df$H),]$G] = 1
  G_info = cbind(as.vector(table(G)),na.group)
  
  point_estimates  <- ipw_point_estimates(H, G, A, w.matrix)
  onestage_ipw <- c(onestage_ipw,
                    point_estimates$marginal_outcomes$overall)
  
  alphas   <- dimnames(w.matrix)[[length(dim(w.matrix))]]
  allocation1 <- alphas[1]
  allocation2 <- NA
  mvar_ipw <- c(mvar_ipw,
                ipw_effect_calc(w.matrix, point_estimates, effect_type ='outcome', 
                                marginal = TRUE, allocation1, allocation2)[2][[1]])
}
mean(mvar_ipw) #0.7525083
sd(onestage_ipw) #0.7682715

##################### different allocations
onestage_ipw = c()
mvar_ipw = c()
for (i in 1:10) {
  graph = make_empty_graph(n = 0, directed = FALSE)
  repeat{
    g2 = sample_gnp(100, 0.5, directed = FALSE, loops = FALSE)
    graph = disjoint_union(graph, g2)
    if (clusters(graph)$no == 50){
      break
    }
  }
  G = group_vector(graph) 
  
  denominator_alphas = 0.5
  P = 1
  A = sample(c(0,1), length(clusters(graph)$membership), 
             replace = TRUE, prob = c(1 - denominator_alphas, denominator_alphas))
  
  df <- cbind.data.frame(A,G)
  df$neighbor_sum <- h_neighsum(graph, A, 1)
  
  a = 5; b = 2
  Y = apply(cbind(df$A,df$neighbor_sum), 1, 
            function(x)  rnorm(1, mean = a*x[1] +  b*x[2], sd = 1))  
  
  H = h_neighborhood(graph,Y,1) 
  df$Y = Y
  df$H = H
  head(df)
  
  allocations = list(c(0.5,denominator_alphas), c(0.4, denominator_alphas))
  w.matrix = wght_matrix(integrand, allocations, G, A, P)
  
  na.group = rep(NA, length(unique(df$G)))
  na.group[df[!is.na(df$H),]$G] = 1
  G_info = cbind(as.vector(table(G)),na.group)
  
  point_estimates  <- ipw_point_estimates(H, G, A, w.matrix)
  onestage_ipw <- c(onestage_ipw,
                    point_estimates$marginal_outcomes$overall)
  
  alphas   <- dimnames(w.matrix)[[length(dim(w.matrix))]]
  allocation1 <- alphas[1]
  allocation2 <- alphas[2]
  mvar_ipw <- c(mvar_ipw,
                ipw_effect_calc(w.matrix, point_estimates, effect_type ='contrast', 
                                marginal = TRUE, allocation1, allocation2)[2][[1]]
  )
}
mean(mvar_ipw) #0.7525083
sd(onestage_ipw) #0.7682715


####
# 1. conditional on influencers, notice the variance needs to be changed
# 2. check for two stage result


