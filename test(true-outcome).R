source("weight_matrix.R")
source("ipw_point_estimates.R")
source("integrand.R")
source("utils.R")
source("bootstrap_variance.R")

group_vector <- function(g){
  number_G = length(groups(components(g)))
  G = rep(NA,number_G)
  for (i in 1:number_G) {
    G[unlist(groups(components(g))[[i]])] = i # vertex indices in group i
  }
  return(G)
}

## test dataset
library(igraph)
g = sample_gnm(n = 200, m = 150) # 200 vertices, 150 edges

plot(g, vertex = 6)
G = group_vector(g) 

## One-stage
denominator_alphas = 0.5
#P = 1
A = sample(c(0,1), 200, replace = TRUE, prob = c(0.5,0.5))

df <- cbind.data.frame(A,G)

# number of treated h-order neighbors for each node
df$neighbor_sum <- rep(NA, length(A))
for (i in 1:length(A)) {
  df$neighbor_sum[i] <- sum(A[G == G[i]], na.rm = TRUE) - A[i]
}

Y = apply(cbind(df$A,df$neighbor_sum), 1, 
          function(x)  rnorm(1, mean = x[1] + 2*x[2], sd = 1))  
H = h_neighborhood(g,Y,0) 
df$Y = Y
df$H = H
#df

allocations = list(c(0.5,denominator_alphas))
w.matrix = wght_matrix(integrand, allocations, G, A, P)

## true individual effect
# z_jk is the treatment for the unit jk itself
# n_k is the number of elements in the cluster k
# alpha is the probability of treatment being 1
# m is the magnitude of influence on neighborhood
true_ind_effect <- function(z_jk, n_k, alpha, a, b){
  sum((a * z_jk + b * 0:(n_k - 1)) *
      dbinom(0:(n_k - 1), size = (n_k - 1), prob = alpha))
} #(1-alpha)^(n_k - 1) * (alpha/(1-alpha))^(0:(n_k - 1)) * 2^(n_k - 1)

# Number if elements in each groups
# df$count_G <-  rep(NA, length(A))
# for (i in 1:length(A)) {
#   df$count_G[i] <- sum(G == G[i], na.rm = TRUE)
# }

true_cluster_effect <- function(g_info, alpha, a, b){
  # number of units in that cluster
  count_g = g_info[1] 
  # whether this cluster has no units with h-order neighbors,
  # if so, set the cluster effect to NA.
  is.h = g_info[2] 
  c(true_ind_effect(0, count_g, alpha, a, b), true_ind_effect(1, count_g, alpha, a, b)) * is.h
}

true_population_effect <- function(G_info, alpha, a, b){
  #sapply(count_G, true_cluster_effect, alpha, m)
  
  # G_info contains information about whether this cluster is NA or not
  rowMeans(apply(G_info, 1, true_cluster_effect, alpha, a, b), na.rm = TRUE)
}

na.group = rep(NA, length(unique(df$G)))
na.group[df[!is.na(df$H),]$G] = 1
G_info = cbind(as.vector(table(G)),na.group)
#apply(G_info, 1, true_cluster_effect, 0.5, 2)
true_population_effect(G_info, 0.5, 5, 2)
#  9.642857 10.642857

ipw_point_estimates(df$H, df$G, df$A,w.matrix)$outcomes$overall
# 8.985511 10.483468 

bootstrap_avg = BootVar(df, 0.5, denominator_alphas)
apply(bootstrap_avg, 1, mean) # bootstrapped population mean
# 9.050644 10.548419 

apply(bootstrap_avg, 1, sd) # bootstrapped population std
# 7.852672 7.313605

#TODO:
# 1) 5000 units 50 clusters
# 2) check for the diff-in-means between treated and untreated; check for 2-stage
# 3) check for h = 0 with that package, should be 1
# 
# > c(mean(df[df$A == 0,]$H), mean(df[df$A == 1,]$H))
# [1] 93.54618 98.14634
# > bootstrap_avg = BootVar(df, 0.5, denominator_alphas)
# > ipw_point_estimates(df$H, df$G, df$A,w.matrix)$outcomes$overall
# 0        1 
# 94.86333 96.31396 

mean(df[df$G == 1 & df$A == 0,]$H) # 88.11291
