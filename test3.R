# 01/15/22
## 1. run multiple times
## 2. change graph and outcome models
## 3. see variance
source("weight_matrix.R")
source("ipw_point_estimates.R")
source("integrand.R")
source("utils.R")
source("bootstrap_variance.R")
source("effects.R")
source("m_variance.R")

library(igraph)

onestage1 = c()
for (i in 1:100) {
  graph = make_empty_graph(n = 0, directed = FALSE)
  repeat{
    g2 = sample_gnp(100, 0.5, directed = FALSE, loops = FALSE)
    graph = disjoint_union(graph, g2)
    if (clusters(graph)$no == 50){
      break
    }
  }
  G = group_vector(graph) 
  
  ## One-stage when h = 1
  denominator_alphas = 0.5
  P = 1
  A = sample(c(0,1), length(clusters(graph)$membership), 
             replace = TRUE, prob = c(1 - denominator_alphas, denominator_alphas))
  
  df <- cbind.data.frame(A,G)
  df$neighbor_sum <- h_neighsum(graph, A, 1)
  # df$neighbor_sum <- rep(NA, length(A))
  # for (i in 1:length(A)) {
  #   df$neighbor_sum[i] <- sum(A[G == G[i]], na.rm = TRUE) - A[i]
  # }
  
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
  
  #true_population_effect(G_info, 0.5, a, b)
  onestage1 = c(onestage1,
               population_direct_effect(ipw_point_estimates(df$H, df$G, df$A,w.matrix)))
  
  
}
mean(onestage1) # 1.845043
hist(onestage1)

onestage0 = c()
for (i in 1:100) {
  graph = make_empty_graph(n = 0, directed = FALSE)
  repeat{
    g2 = sample_gnp(100, 0.5, directed = FALSE, loops = FALSE)
    graph = disjoint_union(graph, g2)
    if (clusters(graph)$no == 50){
      break
    }
  }
  G = group_vector(graph) 
  
  ## One-stage when h = 0
  denominator_alphas = 0.5
  P = 1
  A = sample(c(0,1), length(clusters(graph)$membership), 
             replace = TRUE, prob = c(1 - denominator_alphas, denominator_alphas))
  
  df <- cbind.data.frame(A,G)
  df$neighbor_sum <- h_neighsum(graph, A, 1)
  
  a = 5; b = 2
  Y = apply(cbind(df$A,df$neighbor_sum), 1, 
            function(x)  rnorm(1, mean = a*x[1] +  b*x[2], sd = 1))  
  
  H = h_neighborhood(graph,Y,0) 
  df$Y = Y
  df$H = H
  head(df)
  
  allocations = list(c(0.5,denominator_alphas))
  w.matrix = wght_matrix(integrand, allocations, G, A, P)
  
  na.group = rep(NA, length(unique(df$G)))
  na.group[df[!is.na(df$H),]$G] = 1
  G_info = cbind(as.vector(table(G)),na.group)
  
  #true_population_effect(G_info, 0.5, a, b)
  onestage0 = c(onestage0,
               population_direct_effect(ipw_point_estimates(df$H, df$G, df$A,w.matrix)))
}
mean(onestage0) # 5.198845
hist(onestage0)

onestage2 = c()
for (i in 1:100) {
  graph = make_empty_graph(n = 0, directed = FALSE)
  repeat{
    g2 = sample_gnp(100, 0.5, directed = FALSE, loops = FALSE)
    graph = disjoint_union(graph, g2)
    if (clusters(graph)$no == 50){
      break
    }
  }
  G = group_vector(graph) 
  
  ## One-stage when h = 1
  denominator_alphas = 0.5
  P = 1
  A = sample(c(0,1), length(clusters(graph)$membership), 
             replace = TRUE, prob = c(1 - denominator_alphas, denominator_alphas))
  
  df <- cbind.data.frame(A,G)
  
  df$neighbor_sum <- h_neighsum(graph, A, 1)
  df$neighbor_sum2 <- h_neighsum(graph, A, 2)
  
  a = 5; b = 2; c = 10
  Y = apply(cbind(df$A,df$neighbor_sum, df$neighbor_sum2), 1, 
            function(x)  rnorm(1, mean = a*x[1] +  b*x[2] + c*x[3], sd = 1))  
  
  H = h_neighborhood(graph,Y,2) 
  df$Y = Y
  df$H = H
  head(df)
  
  allocations = list(c(0.5,denominator_alphas))
  w.matrix = wght_matrix(integrand, allocations, G, A, P)
  
  na.group = rep(NA, length(unique(df$G)))
  na.group[df[!is.na(df$H),]$G] = 1
  G_info = cbind(as.vector(table(G)),na.group)
  
  #true_population_effect(G_info, 0.5, a, b)
  onestage2 = c(onestage2,
                population_direct_effect(ipw_point_estimates(df$H, df$G, df$A,w.matrix)))
}#10.51043

# check for h = 2, it should be 0 checked
# number of treated neighbors of neighbors for h = 2 checked

# plot out the variance and distribution
# run multiple times for 2-stage

# compute the variance for m-estimator (check for package),
# compare with the monte-carlo estimator (empirical variance) from the distribution of point-estimate 

# conditional effects

