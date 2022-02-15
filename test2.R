# 01/08/22
# 1) 5000 units 50 clusters
# 2) check for the diff-in-means between treated and untreated; check for 2-stage
# 3) check for h = 0 with that package
source("weight_matrix.R")
source("ipw_point_estimates.R")
source("integrand.R")
source("utils.R")
source("bootstrap_variance.R")

library(igraph)
graph = make_empty_graph(n = 0, directed = FALSE)

##TODO
repeat{
  #g = sample_gnm(n = 120, m = 120)
  g2 = sample_gnp(100, 0.5, directed = FALSE, loops = FALSE)
  #cls = clusters(g)
  #g2 <- delete_vertices(g, V(g)[cls$membership %in% which(cls$csize <= 80)])
  clusters(g2)$csize
  if (clusters(g2)$no == 1){
    graph = disjoint_union(graph, g2)
  }
  if (clusters(graph)$no == 50){
    break
  }
}
#sample_gnp(100, 0.5, directed = FALSE, loops = FALSE)

G = group_vector(graph) 

## One-stage when h = 0
denominator_alphas = 0.5
P = 1
A = sample(c(0,1), length(clusters(graph)$membership), 
           replace = TRUE, prob = c(1 - denominator_alphas, denominator_alphas))

df <- cbind.data.frame(A,G)

# number of treated h-order neighbors for each node
df$neighbor_sum <- rep(NA, length(A))
for (i in 1:length(A)) {
  df$neighbor_sum[i] <- sum(A[G == G[i]], na.rm = TRUE) - A[i]
}

a = 5; b = 2
Y = apply(cbind(df$A,df$neighbor_sum), 1, 
          function(x)  rnorm(1, mean = a*x[1] +  b*x[2], sd = 1))  
H = h_neighborhood(graph,Y,1) 
df$Y = Y
df$H = H
head(df)

allocations = list(c(0.4,denominator_alphas),c(0.5,denominator_alphas))
w.matrix = wght_matrix(integrand, allocations, G, A, P)

na.group = rep(NA, length(unique(df$G)))
na.group[df[!is.na(df$H),]$G] = 1
G_info = cbind(as.vector(table(G)),na.group)

true_population_effect(G_info, 0.5, a, b)
# 76.752 81.752
# 95.06 100.06

ipw_point_estimates(df$H, df$G, df$A,w.matrix)$outcomes$overall
#  h = 0            0         1
# c(0.4, 0.5) 64.08977 72.02618
# c(0.5, 0.5) 95.88820 101.59236
#  h = 1            0         1
# c(0.4, 0.5) 65.16808  70.45095
# c(0.5, 0.5) 97.39337 100.08929


# sum(ind_est) is the same as treatment mean for each group
#c(mean(df[df$A == 0,]$H), mean(df[df$A == 1,]$H))

bootstrap_avg = BootVar(df, 0.5, denominator_alphas)
apply(bootstrap_avg, 1, mean) # bootstrapped population mean
# 63.67563 71.59888  
# 94.96673 100.88503

apply(bootstrap_avg, 1, sd) # bootstrapped population std
# 25.59980 25.14421  
# 0.6863778 2.9587914 

library(inferference)

test_propensity <- function(b, X, A, parameters, 
                                 P = 1){
  h = as.vector(parameters)
  prob_alphas = unlist(lapply(h, function(alpha) prod(alpha^A * (1-alpha)^(1-A))))
  out = sum(prob_alphas * P)
  out = prod(h^A * (1-h)^(1-A)) * dnorm(b) 
  out
}

example1 <- interference(
       formula = H | A  ~ 1 | G, 
       propensity_integrand = 'test_propensity',
       data = df,
       model_method = 'oracle', 
       model_options = list(fixed.effects = 0.5, random.effects = NULL),
       allocations = c(0.4, 0.5),
       integrate_allocation = FALSE,
       causal_estimation_options = list(variance_estimation = 'naive'),
       conf.level = .9
       )

print(example1)

# direct_effect(example1)
# Direct Effects when h = 0
# alpha1 trt1 alpha2 trt2 estimate std.error conf.low conf.high
# 0.4    0    0.4    1   -7.936     1.456   -10.33   -5.5413
# 0.5    0    0.5    1   -5.704     3.024   -10.68   -0.7305

# Direct Effects when h = 1
# alpha1 trt1 alpha2 trt2 estimate std.error conf.low conf.high
# 0.4    0    0.4    1   -5.283     2.037   -8.634    -1.932
# 0.5    0    0.5    1   -2.696     3.037   -7.692     2.300

## Two-stage when h = 0
numerator_alpha = 0.4
denominator_alphas = c(0.4,0.6)
P = c(0.6,0.4)
a = 5; b = 2
G_info = cbind(as.vector(table(G)),1)
true_population_effect(G_info, numerator_alpha, a, b)
# 76.08 81.08

ipw_2stage = c()

for (k in 1:100) {
  P_1 = sample(c(1,2), length(unique(G)), replace = TRUE, prob = P)
  P_1 = sapply(G, function(x) P_1[x])
  P_2 = rep(NA, length(P_1))
  
  A = rep(NA, length(P_1))
  for (i in 1:length(P_1)) {
    P_2[i] = denominator_alphas[P_1[i]]
    A[i] = sample(c(0,1), 1, prob = c(1-P_2[i],P_2[i]))
  }
  
  df <- cbind.data.frame(A,G)
  df$neighbor_sum <- rep(NA, length(A))
  for (i in 1:length(A)) {
    df$neighbor_sum[i] <- sum(A[G == G[i]], na.rm = TRUE) - A[i]
  }
  
  Y = apply(cbind(df$A,df$neighbor_sum), 1, 
            function(x)  rnorm(1, mean = a*x[1] +  b*x[2], sd = 1))  
  H = h_neighborhood(graph,Y,0) 
  df$Y = Y
  df$H = H
  #head(df)
  
  allocations = list(c(numerator_alpha,denominator_alphas))
  w.matrix = wght_matrix(integrand, allocations, G, A, P)
  ipw_2stage = rbind(ipw_2stage,
                     as.vector(ipw_point_estimates(df$H, df$G, df$A,w.matrix)$outcomes$overall))
}
mean(ipw_2stage[,1]) #76.32368
mean(ipw_2stage[,2]) #81.9589


## 1. run multiple times
## 2. change graph and outcome models
## 3. see variance

# causal effect depends on hypothetical treatment allocation
# denominator is the propensity score, different experimental design 
# lead to different estimators, resulting in different denominators.

