# 02/06/22
# 1. change the outcome model to include interaction between x[2] and covariates; 
# change x[2] to number of treated neighbors with certain characteristics.
# 2. complete randomization check
# 3. check for different combination of allocation and treatments


source("weight_matrix.R")
source("ipw_point_estimates.R")
source("integrand.R")
source("utils.R")
source("bootstrap_variance.R")
source("effects.R")
source("m_variance.R")
library(igraph)
#################
# complete randomization design

# package for complete randomization
library(randomizr)

aa = c()
bb = c()
for (i in 1:50) {
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
  #A = complete_ra(N = length(clusters(graph)$membership))
  A = sample(c(0,1), length(clusters(graph)$membership), 
             replace = TRUE, prob = c(1 - denominator_alphas, denominator_alphas))
  X = as.matrix(sample(c("M", "F"), length(clusters(graph)$membership), 
                       replace = TRUE, prob = c(1 - denominator_alphas, denominator_alphas)))
  
  df <- cbind.data.frame(A,G,X)
  df$treated_male_neigh <- h_neighsum(graph, A, 1, X = X, x1 = "M")
  df$treated_female_neigh <- h_neighsum(graph, A, 1, X = X, x1 = "F")
  #checked that the sum of above is the same as h_neighsum(graph, A, 1)
  df$X_numeric <- ifelse(df$X == "F", 1, 0)
  
  a = 2; b = 0; c = 5
  Y = apply(cbind(df$A, df$X_numeric, df$treated_female_neigh), 1, 
            function(x)  rnorm(1, mean = a*x[1] + b*x[2] + c*x[3], sd = 1))  
  # x[2] change to the number of females neighbors 
  
  H = h_neighborhood(graph, Y, 1) 
  H_F = h_neighborhood(graph, Y, 1, X, x1 = "F") 
  df$Y = Y
  df$H = H
  head(df)
  
   
  allocations = list(c(0.5,denominator_alphas))
  w.matrix = wght_matrix(integrand, allocations, G, A, P)
  
  na.group = rep(NA, length(unique(df$G)))
  na.group[df[!is.na(df$H),]$G] = 1
  G_info = cbind(as.vector(table(G)),na.group)
  
  point_estimates  <- ipw_point_estimates(H, G, A, w.matrix)
  point_estimates1  <- ipw_point_estimates(H, G, A, w.matrix, X = X, x0 = "M")
  point_estimatesF  <- ipw_point_estimates(H_F, G, A, w.matrix)
  
  aa <- rbind(aa,
              c(point_estimates$outcomes$overall[2] -  point_estimates$outcomes$overall[1],
                point_estimatesF$outcomes$overall[2] -  point_estimatesF$outcomes$overall[1],
                point_estimates1$outcomes$overall[2] -  point_estimates1$outcomes$overall[1])
              )
  bb <- rbind(bb,
              c(
                point_estimates$marginal_outcomes$overall,
                point_estimatesF$marginal_outcomes$overall, #should be smaller
                point_estimates1$marginal_outcomes$overall # should be nearly the same
              )
              )
}

# colMeans(aa)  # b = 0
#               # x0 = NULL, x1 = NULL, when h = 1, 7.167375 (c/2)
#               # x0 = NULL, x1 = "F", when h = 1, 7.181179 (c/2)
#               # x0 = "F", x1 = NULL, when h = 1, 7.234829 (c/2)

colMeans(aa)  
# [1] 2.9674848 2.9647760 0.3055283

colMeans(bb) 
# 63.90863    63.30353    63.27293 



# b = 0 
#7.167375 7.181179 7.234829
# 188.4293    186.3833    188.4651 
# b= 10
#6.872444 6.848343 7.006094
#193.2092    196.3681    193.2459 

# b = 20
#6.965107 6.944501 7.281182
#198.3724    206.4918    198.4122 

# increasing b will increase marginal effect for both "overall" population (by b/2), and 
# only consider female when calculating H (by b)