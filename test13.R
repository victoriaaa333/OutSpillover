source("weight_matrix.R")
source("ipw_point_estimates(mixed variables).R")
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



# model for spillover effect

#############
#1. Generate a graph and dataset (treatments, covariates)
graph = make_empty_graph(n = 0, directed = FALSE)
repeat{
  g2 = sample_gnp(100, 0.5, directed = FALSE, loops = FALSE)
  graph = disjoint_union(graph, g2)
  if (clusters(graph)$no == 50){
    break}
}

#G = group_vector(graph) 
G = components(graph)$membership

# Two-stage randomization
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
X3 <- rnorm(length(A),mean = 0.5, sd = 1)
# X <- cbind(X1, X3)
# X_type <- c("N", "N") # indicating whether the covariate is numerical or categorical
# x0 <- as.matrix(c( 0.1, 0.1))

X2 <- sample(c("M", "F"), size = length(A), replace = TRUE)
X4 <- sample(c("Y", "N"), size = length(A), replace = TRUE)
X <- cbind(X2, X3, X4)
X_type <- c("C", "N", "C")
x0 <- as.matrix(c("M", 0.1, "Y"))
x1 <- x0
x1_num <- c(0.1)

X_num <- apply(as.matrix(X[, X_type == "N"]), 2, as.numeric)
X_cat <- as.matrix(X[, X_type == "C"])

df <- cbind.data.frame(A,G,X)
df$treated_neigh <- h_neighsum(graph, A, 1) 
#df$interaction1 <- X_num[,1] * df$treated_neigh
#df$interaction2 <- X_num[,2] * df$treated_neigh
df$interaction1 <- ifelse(X_cat[,1] == "M", 1, 0) * df$treated_neigh
df$interaction2 <- ifelse(X_cat[,2] == "N", 1, 0) * df$treated_neigh

#############
# 2. Outcome model
a = 2; b = 5; c = 7; d = 9
Y = apply(cbind(df$A, df$treated_neigh, df$interaction1, df$interaction2), 1, #X_num,
          function(x)  rnorm(1, mean = a*x[1] + b*x[2] + c*x[3] + d*x[4], sd = 1))  
H = h_neighborhood(graph, Y, 1) 
H_M =  h_neighborhood(graph, Y, 1, X_cat, c("M", "Y")) 
df$Y = Y
df$H = H
df$H_M = H_M

# 2.1 neighinfo
neighX = h_neighcov(graph, Y, 1, X, X_type, x1) # average of X for unit j's 1-order neighbor
h_neigh = h_neighsum(graph, A, 1)
neigh2_treated = h_neighofneigh(graph, A, 1, h_neigh, 
                                X, X_type, x1) # average number of treated neighbors l for unit j's 1-order neighbor i
neigh2_treated_neighX = h_neighofneigh_withcov(graph, A, 1, X, X_type, x1, h_neigh) # average number of treated neighbors l for unit j's 1-order neighbor i, times X_i

neighinfo = list(neigh2_treated, neighX, neigh2_treated_neighX)
names(neighinfo) <- c('neigh2_treated', 'neighX', 'neigh2_treated_neighX')



#############
# 3. calculate the point estimates and the variances (bootstrapped and analytical)
allocations = list(c(0.5,denominator_alphas))
w.matrix = wght_matrix(integrand, allocations, G, A, P)


# 3.0 no-condition
X_non = cbind(X2, X4)
x0_non = as.matrix(c("M", "Y"))
X_type_non = c("C","C")
point_estimates = ipw_point_estimates_mixed(H, G, A, w.matrix, 
                                            X = X_non, x0 = x0_non,
                                            X_type = X_type_non)
point_estimates$outcomes$overall
ipw_m_variance(w.matrix, point_estimates, effect_type ='contrast', 
              marginal = FALSE, allocation1 = allocations[1], allocation2 = allocations[1])
               

# 3.1. group avg conditional
point_estimates_g = ipw_point_estimates_mixed(H, G, A, w.matrix, 
                                            X = X, x0 = x0, X_type = X_type)
point_estimates_g$outcomes$overall
# 0        1
# c(0.5, 0.4, 0.6) 313.073 320.3357

# ipw_m_variance(w.matrix, point_estimates, effect_type ='outcome', 
#                marginal = FALSE, allocation1 = allocations[1])

ipw_regression_variance(H, w.matrix, point_estimates_g, effect_type ='outcome', 
                        marginal = FALSE, allocation1 = allocations[1], 
                        X = X, x0 = x0, X_type = X_type)

# estimate std.error conf.low conf.high
#  313.073  5.848343 301.6105  324.5356

# ipw_m_variance(w.matrix, point_estimates, effect_type ='contrast', 
#                marginal = FALSE, allocation1 = allocations[1], allocation2 = allocations[1])

ipw_regression_variance(H, w.matrix, point_estimates_g, effect_type ='contrast', 
                        marginal = FALSE, allocation1 = allocations[1], allocation2 = allocations[1], 
                        X = X, x0 = x0, X_type = X_type)

# estimate std.error conf.low conf.high
# -7.262691  0.444655 -8.134199 -6.391184

# 3.2. H conditional
point_estimates_neigh =  ipw_point_estimates_mixed(H_M, G, A, w.matrix, 
                                             neighinfo = neighinfo, x1 = x1, X_type = X_type)

point_estimates_neigh$outcomes$overall
#                         0        1
# c(0.5, 0.4, 0.6) 301.7175 339.3596

# ipw_m_variance(w.matrix, point_estimates_neigh, effect_type ='outcome', 
#               marginal = FALSE, allocation1 = allocations[1])

ipw_regression_variance_neigh(H_M, w.matrix, point_estimates_neigh, effect_type ='outcome', 
                        marginal = FALSE, allocation1 = allocations[1], 
                        neighinfo = neighinfo, x1_num = x1_num)

# estimate std.error conf.low conf.high
# 301.7175 0.1200248 301.4822  301.9527
# 339.3596 0.1037242 339.1563  339.5629


# ipw_m_variance(w.matrix, point_estimates_neigh, effect_type ='contrast', 
#                marginal = FALSE, allocation1 = allocations[1], allocation2 = allocations[1])

ipw_regression_variance_neigh(H_M, w.matrix, point_estimates_neigh, effect_type ='contrast', 
                              marginal = FALSE, allocation1 = allocations[1], allocation2 = allocations[1], 
                              neighinfo = neighinfo, x1_num = x1_num)
# estimate std.error conf.low conf.high
# -37.64217 0.007102515 -37.65609 -37.62825
