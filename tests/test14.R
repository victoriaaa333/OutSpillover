##
# 22.07.27
# 1. outcome model with both influencer and spillover effect

# neighinfo has to be non-NULL, or refer to the function group_coefs_oncont2 and group_means_oncont3

mixed_coef <- function(weights_df, H, G, X, X_type, x0, neighinfo, A, a){
  if (is.null(neighinfo)) stop("neighinfo cannot be NULL, 
                               refer to the function group_coef")
  
  trt_cond <- which(A == a)
  #### information needed for conditional H ####
  neigh2_treated = neighinfo$neigh2_treated
  neighX = neighinfo$neighX
  neigh2_treated_neighX = neighinfo$neigh2_treated_neighX
  colnames(neighX) <- paste("neighX", colnames(neighX), sep = "_")
  colnames(neigh2_treated_neighX) <- paste("neigh2_treated_neighX", 
                                           colnames(neigh2_treated_neighX), sep = "_")
  neigh_names <- c("neigh2_treated", colnames(neighX), colnames(neigh2_treated_neighX))
  
  #### select out categorical values for conditional group avg ####
  # X_cat <- as.matrix(X[, X_type == "C"])
  # x0_cat <- x0[X_type == "C"]
  # X_num <- apply(as.matrix(X[, X_type == "N"]), 2, as.numeric)
  # num_names <- colnames(X)[X_type == "N"]
  # x0_num <- as.numeric(x0[X_type == "N"])
  
  if (sum(X_type == "C") > 0){
    X_cat <- as.matrix(X[, X_type == "C"])
    x0_cat <- x0[X_type == "C"]
    ind_cond <- which(apply(X_cat, 1, function(x) prod(x == x0_cat)) == 1)
  }else{
    ind_cond <- 1:length(neigh2_treated)}
  
  
  if (sum(X_type == "N") > 0){
    X_num <- apply(as.matrix(X[, X_type == "N"]), 2, as.numeric)
    num_names <- colnames(X)[X_type == "N"]
    #x0_num <- as.numeric(x0[X_type == "N"])
    group_df <- as.data.frame(cbind(weights_df, H, G, X_num, neigh2_treated, neighX, neigh2_treated_neighX))
    colnames(group_df) <- c('w', 'H', 'G', num_names, neigh_names)
    fits_df <- group_df[intersect(ind_cond, trt_cond), ]
    fits <- lmList(as.formula(paste("H ~ ", paste(num_names, collapse= "+ "), "+", 
                                    paste(colnames(fits_df)[-1:-(3+length(num_names))], collapse= "+ "), "| G")),
                   weights = w, data = fits_df)
    fits2 <- lmList(neigh2_treated ~ 1 | G, data=fits_df)
    cond_coef <- cbind(coef(fits), coef(fits2))

    overall_fits <- lm(as.formula(paste("H ~ ", paste(num_names, collapse= "+ "), "+", 
                                        paste(colnames(fits_df)[-1:-(3+length(num_names))], collapse= "+ "))),
                       weights = w, data=fits_df)
    overall_fits2 <- lm(neigh2_treated ~ 1, data=fits_df)
    overall_coef <- c(coef(overall_fits), coef(overall_fits2))

    coef <- list(cond_coef, overall_coef)
  }else{
    group_df <- as.data.frame(cbind(weights_df, H, G, neigh2_treated, neighX, neigh2_treated_neighX))
    colnames(group_df) <- c('w', 'H', 'G', neigh_names)
    fits_df <- group_df[intersect(ind_cond, trt_cond), ]
    fits <- lmList(as.formula(paste("H ~ ", paste(colnames(fits_df)[-1:-3], 
                                                       collapse= "+ "), "| G")), weights = w, data = fits_df)
    fits2 <- lmList(neigh2_treated ~ 1 | G, data=fits_df)
    cond_coef <- cbind(coef(fits), coef(fits2))
    
    overall_fits <- lm(as.formula(paste("H ~ ",  paste(colnames(fits_df)[-1:-3], 
                                                       collapse= "+ "))), weights = w, data=fits_df)
    overall_fits2 <- lm(neigh2_treated ~ 1, data=fits_df)
    overall_coef <- c(coef(overall_fits), coef(overall_fits2))
  
    coef <- list(cond_coef, overall_coef)
    }
  coef
}


mixed_means <- function(overall_coef, X_type, x0, x1){
  x0_num <- as.numeric(x0[X_type == "N"])
  x1_num <- as.numeric(x1[X_type == "N"])
  lenn <- length(x1_num)
  leng <- length(x0_num)
  cond_coefs <- matrix(rep(overall_coef, N), nrow = N, byrow = TRUE)
  
  if (leng > 0){
    cond_group_means <- cond_coefs[,1] + cond_coefs[,2:(1+leng)] %*% as.matrix(x0_num) +
      cond_coefs[,(2+leng)] * cond_coefs[,dim(cond_coefs)[2]] +
      as.matrix(cond_coefs[,(3+leng):(2+leng+lenn)]) %*% as.matrix(x1_num) + 
      as.matrix(cond_coefs[,(3+leng+lenn):(2+leng+2*lenn)]) %*% as.matrix(x1_num) * cond_coefs[,dim(cond_coefs)[2]]
  }else{
    cond_group_means <- cond_coefs[,1] + cond_coefs[,2] * cond_coefs[,dim(cond_coefs)[2]] +
      as.matrix(cond_coefs[,3:(2+lenn)]) %*% as.matrix(x1_num) + 
      as.matrix(cond_coefs[,(3+lenn):(2+2*lenn)]) %*% as.matrix(x1_num) * cond_coefs[,dim(cond_coefs)[2]]
  }
  cond_group_means
}


source("ipw_point_estimates_mixed.R")


non = c()
grp = c()
nei = c()
mix = c()

for (p in 1:100) {
  
#############
#1. Generate a graph and dataset (treatments, covariates)
graph = make_empty_graph(n = 0, directed = FALSE)
repeat{
  g2 = sample_gnp(100, 0.5, directed = FALSE, loops = FALSE)
  graph = disjoint_union(graph, g2)
  if (clusters(graph)$no == 50){
    break}
}

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
X2 <- sample(c("M", "F"), size = length(A), replace = TRUE)
X4 <- sample(c("Y", "N"), size = length(A), replace = TRUE)
X <- cbind(X1,X2,X3)
X_type <- c("N","C","N")
x0 <- as.matrix(c(0.1, "M", 0.2))
x1 <- x0
#x1_num <- c(0.1)

X_num <- apply(as.matrix(X[, X_type == "N"]), 2, as.numeric)
X_cat <- as.matrix(X[, X_type == "C"])

df <- cbind.data.frame(A,G,X)
df$treated_neigh <- h_neighsum(graph, A, 1) 
df$interaction1 <- ifelse(X_cat[,1] == "M", 1, 0) * df$treated_neigh
df$interaction2 <- X_num[,1] * df$treated_neigh
df$interaction3 <- cov_neighsum(graph, A, 1, X = X_num[,2]) 

# spillover effect with categorical condition
# spillover effect with numerical regression
# influencer effect with numerical regression

#############
# 2. Outcome model
a = 2; b = 5; c = 7; d = 9; e = 11
Y = apply(cbind(df$A, df$treated_neigh, df$interaction1, df$interaction2, df$interaction3), 1, 
          function(x)  rnorm(1, mean = a*x[1] + b*x[2] + c*x[3] + d*x[4] + e*x[5], sd = 1))  
H = h_neighborhood(graph, Y, 1) 
H_M =  h_neighborhood(graph, Y, 1, X_cat, c("M")) 
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

# original estimates without any condition
point_estimates = ipw_point_estimates_mixed2(H, G, A, w.matrix)
point_estimates$outcomes$overall

# estimates with only group condition
point_estimates2 = ipw_point_estimates_mixed2(H, G, A, w.matrix, 
                                              X = X, x0 = x0, 
                                              X_type = X_type, Con_type = "group")
point_estimates2$outcomes$overall
point_estimates2$outcomes$overall_coefG

# estimates with only neigh condition
point_estimates3 = ipw_point_estimates_mixed2(H_M, G, A, w.matrix, 
                                              neighinfo = neighinfo, x1 = x1, 
                                              X_type = X_type, Con_type = "neigh")
point_estimates3$outcomes$overall
point_estimates3$outcomes$overall_coefH

# estimates with mixed condition
point_estimates4 = ipw_point_estimates_mixed2(H_M, G, A, w.matrix, 
                                              X = X, x0 = x0, 
                                              neighinfo = neighinfo, x1 = x1,
                                              X_type = X_type, Con_type = "mixed")
point_estimates4$outcomes$overall
point_estimates4$outcomes$overall_coefM


non = c(non, point_estimates$outcomes$overall[2] - point_estimates$outcomes$overall[1])
grp = c(grp, point_estimates2$outcomes$overall[2] - point_estimates2$outcomes$overall[1] )
nei = c(nei, point_estimates3$outcomes$overall[2] - point_estimates3$outcomes$overall[1])
mix = c(mix, point_estimates4$outcomes$overall[2] - point_estimates4$outcomes$overall[1])
}

mean(non); mean(grp); mean(nei); mean(mix)

#TODO: write down the expected values with conditional $X_jk$ and $X_ik$

# > 5 + 7 * 0.5 + 9 *0.5 + 11 * 0.5
# [1] 13
# > 5 + 7 * 1 + 9 * 0.1 + 11 * 0.5
# [1] 12.9
# > 5 + 7 * 0.5 + 9 *0.5 + 11 * 0.2
# [1] 7.2
# > 5 + 7 * 1 + 9 * 0.1 + 11 * 0.2
# [1] 15.1

# point_estimates2$outcomes$overall_coefH 
# #[1] "(Intercept)" "X3" "neigh2_treated" "neighX_" "neigh2_treated_neighX_" "(Intercept)"           
# point_estimates2$outcomes$overall_coefG
# point_estimatesg = ipw_point_estimates_mixed2(H, G, A, w.matrix, 
#                                              X = X, x0 = x0, X_type = X_type)
# point_estimatesg$outcomes$overall_coefH 
# point_estimatesg$outcomes$overall_coefG
