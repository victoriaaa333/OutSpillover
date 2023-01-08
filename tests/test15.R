#### 2022.08
# 1. test for estimated propensity score

source("propensity.R")
source("integrand.R")
source("weight_matrix.R")
source("utils.R")
source("bootstrap_variance.R")
source("effects.R")
source("m_variance.R")
source("regression_variance.R")
source("regression_utils.R")
library(igraph)
library(lme4)

# 0.1 first generating the model for spillover effect

graph = make_empty_graph(n = 0, directed = FALSE)
repeat{
  g2 = sample_gnp(100, 0.5, directed = FALSE, loops = FALSE)
  graph = disjoint_union(graph, g2)
  if (clusters(graph)$no == 50){
    break}
}
G = group_vector(graph) 

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
X <- cbind(X1, X3)
X_type <- c("N", "N") # indicating whether the covariate is numerical or categorical
x0 <- as.matrix(c( 0.1, 0.1))
x1 <- as.matrix(c( 0.1, 0.1))
X_num <- apply(X[, X_type == "N"], 2, as.numeric)
df <- cbind.data.frame(A,G,X)
df$treated_neigh <- h_neighsum(graph, A, 1) 
df$interaction1 <- X_num[,1] * df$treated_neigh
df$interaction2 <- X_num[,2] * df$treated_neigh

#############
# 2. Outcome model
aa = 2; bb = 5; cc = 7; dd = 9
Y = apply(cbind(df$A, df$treated_neigh, df$interaction1, df$interaction2), 1, #X_num,
          function(x)  rnorm(1, mean = aa*x[1] + bb*x[2] + cc*x[3] + dd*x[4], sd = 1))  
H = h_neighborhood(graph, Y, 1) 
df$Y = Y
df$H = H

############
#3. neighbor information 
neighX = h_neighcov(graph, Y, 1, X, X_type, x1) # average of X for unit j's 1-order neighbor
h_neigh = h_neighsum(graph, A, 1)
neigh2_treated = h_neighofneigh(graph, A, 1, h_neigh, 
                                X, X_type, x1) # average number of treated neighbors l for unit j's 1-order neighbor i
neigh2_treated_neighX = h_neighofneigh_withcov(graph, A, 1, X, X_type, x1, h_neigh) # average number of treated neighbors l for unit j's 1-order neighbor i, times X_i

neighinfo = list(neigh2_treated, neighX, neigh2_treated_neighX)
names(neighinfo) <- c('neigh2_treated', 'neighX', 'neigh2_treated_neighX')




# 0.2 Necessary bits for propensity score function
formula <- H | A  ~ X1 + X3 + (1|G) | G
#b <- 0.1

# 1. test integrand function
parameters <- unlist(propensity_parameter(formula, df)[1])
X_wint <- cbind(1, df$X1, df$X3)
# logit_integrand(b, X[1:100,], A[1:100], parameters, 
#                 allocation = A[1:100], randomization = 1)
plain_integrand(A[1:100], denominator_alphas = denominator_alphas, P = P)

# 2. test weights function
#wght_calc(integrand, numerator_alpha, 
#          denominator_alphas, parameters = parameters, ...)


## NOTICE:
# If logit_integrand, X needs to incorporate intercept
# TODO: modify this 
allocations = list(c(0.5,denominator_alphas))
w.matrix = wght_matrix(integrand =logit_integrand, 
                        allocations, 
                        G, A, P,
                        X = X_wint, 
                        parameters = parameters,
                        randomizations = c(0.4, 0.6),
                        runSilent = TRUE)

# 3. how to incorporate the second-stage denominator (0.4, 0.6) in propensity score

