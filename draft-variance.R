# a = 1
a = 1
weights_trt <- array(dim= c(length(A), p, k))

for(pp in 1:p){
  for (kk in 1:k) {
    weights_trt[ , pp, kk] <- apply(as.array(G), 1, function(x) weights[,pp,kk][x])
    weights_trt[ , pp, kk] <- weights_trt[ , pp, kk] * (A == a)/
      (numerator_alphas[kk]^a * (1-numerator_alphas[kk])^(1-a)) # add some comments
  }
}

# Compute estimates
ind_est <- apply(weights_trt, 2:3, function(x) x * H) 
ind_est_df <- ind_est[ , , 1]

X_num <- apply(as.matrix(X[, X_type == "N"]), 2, as.numeric)
group_df <- as.data.frame(cbind(ind_est_df, G, X_num, A))
colnames(group_df) <- c('ind_est', 'G', 'X1', 'X3', 'A')

group_dfa <- as.data.frame(group_df[group_df$A == a,]) # group_df <- as.data.frame(group_df[group_df$ind_est != 0,]) 
fits <- lmList(ind_est ~ X1 + X3 | G, data=group_dfa) # 
head(coef(fits))

#fits2 <- lmList(ind_est ~ X1 + X3 + A + X1*A + X3*A | G, data = group_df)
#head(coef(fits2))


ipw_point_estimates_mixed(H, G, A, w.matrix, 
                          X = X, x0 = x0, X_type = X_type)$outcomes$overall_coefG

ipw_point_estimates_mixed(H, G, A, w.matrix, 
                          X = X, x0 = x0, X_type = X_type)$outcomes$grp_coefG

#r_ij = Y_ij - grp_coefG * X_ij

# for each group, and a = 1

a = 1
numerator_alphas = 0.5; kk = 1

X_cat <- as.matrix(X[, X_type == "C"])
X_num <- apply(as.matrix(X[, X_type == "N"]), 2, as.numeric)
num_names <- colnames(X)[X_type == "N"]
x0_cat <- x0[X_type == "C"]
x0_num <- as.numeric(x0[X_type == "N"])

grp_coef <- ipw_point_estimates_mixed(H, G, A, w.matrix, 
                          X = X, x0 = x0, X_type = X_type)$outcomes$grp_coefG

trt_ind <- (A == a) # A = a
if (sum(X_type == "C") > 0){
  cat_ind <- apply(X_cat, 1, function(x) prod(x == x0_cat)) == 1
}else{
  cat_ind <- rep(1,dim(X_cat)[1])} #X = X0 categorical

X_reg <- cbind(G, 1, X_num)
X_coef <- grp_coef[, , , (a+1)] # pick out the coefficients for treatment a
X_fit <- apply(X_reg, 1, function(x) sum(x[-1] * X_coef[x[1],]))
X_fit_trt <- X_fit * trt_ind
res_ij <- ind_est_df - X_fit_trt
res_data <- as.data.frame(cbind(res_ij, G))
res_group <- aggregate(res_data$res_ij, list(res_data$G), FUN=mean) # mean of residual for each group i
var <-  mean((res_group$x - mean(res_group$x))^2) 


# w.matrix * (A == a)/
#   (numerator_alphas[kk]^a * (1-numerator_alphas[kk])^(1-a))
