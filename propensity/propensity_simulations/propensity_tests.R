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

library(igraph)
library(lme4)

##########
#1. Generate a graph and dataset (treatments, covariates)
graph = make_empty_graph(n = 0, directed = FALSE)
repeat{
  g2 = sample_gnp(50, 0.5, directed = FALSE, loops = FALSE)
  graph = disjoint_union(graph, g2)
  if (clusters(graph)$no == 20){
    break}
}
G = group_vector(graph) 

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
X1 <- sample(c("M", "F"), size = length(A), replace = TRUE)
X2 <- rnorm(length(A),mean = 0.5, sd = 1)

X <- cbind(X1, X2)
X_type <- c("C", "N")
x0 <- as.matrix(c("M", 0.1))
x1 <- x0
x1_num <- c(0.1)

X_num <- apply(as.matrix(X[, X_type == "N"]), 2, as.numeric)
X_cat <- as.matrix(X[, X_type == "C"])

df <- cbind.data.frame(A,G,X)
df$treated_neigh <- h_neighsum(graph, A, 1) 
df$interaction1 <-  ifelse(X_cat[,1] == "M", 1, 0)  * df$treated_neigh
df$interaction2 <- X_num[,1] * df$treated_neigh

##########
# 2. Outcome model
aa = 2; bb = 5; cc = 7; dd = 9
Y = apply(cbind(df$A, df$treated_neigh, df$interaction1, df$interaction2), 1, #X_num,
          function(x)  rnorm(1, mean = aa*x[1] + bb*x[2] + cc*x[3] + dd*x[4], sd = 1))  
H = h_neighborhood(graph, Y, 1) 
H_M =  h_neighborhood(graph, Y, 1, X_cat, c("M")) 
df$Y = Y
df$H = H
df$H_M = H_M

# 2.1 neighinfo
neighX = h_neighcov(graph, 1, X, X_type, x1) # average of X for unit j's 1-order neighbor
neighinfo = list(neighX)
names(neighinfo) <- c('neighX')

##########
# 3. calculate the point estimates and the variances (bootstrapped and analytical)
allocations = list(c(0.5,denominator_alphas))
w.matrix = wght_matrix(plain_integrand, allocations, G, A, P)

#######
# propensity score
formula <- H | A  ~  X2  + (1|G) | G
parameters <- unlist(propensity_parameter(formula, df)[1])
wght_matrix(parameters = parameters, integrand = logit_integrand, 
            G = G, A = A, P = P, X = as.matrix(X2), allocations = allocations)
# weight is ppp/PrA (Pi(A_i,\alpha)/Pr(A_i | X+i, \phi_hat))
wght_deriv_array(parameters, logit_integrand, G, A, P, X = as.matrix(X2), allocations = allocations)

ipw_propensity_variance(parameters,
                        allocations,
                        causal_estimation_options = 
                          list(variance_estimation = 'robust'),
                        integrate_allocation = FALSE,
                        H, as.matrix(X2), A, G,
                        effect_type = "outcome")
