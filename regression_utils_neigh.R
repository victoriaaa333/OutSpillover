# we only need to consider numerical variance in this case, if categorical,
# we select those out when calculating H.
point_estimates <- ipw_point_estimates_mixed(H_M, G, A, w.matrix, 
                                             neighinfo = neighinfo, x1 = x1, X_type = X_type)
reg_coef <- point_estimates[[fff]]$overall_coefH

var_para <- function(reg_coef, a1, t1, neighinfo, H, weights_ind, A, G){
  neigh2_treated = neighinfo$neigh2_treated
  neighX = neighinfo$neighX
  neigh2_treated_neighX = neighinfo$neigh2_treated_neighX
  colnames(neighX) <- paste("neighX", colnames(neighX), sep = "_")
  colnames(neigh2_treated_neighX) <- paste("neigh2_treated_neighX", colnames(neigh2_treated_neighX), sep = "_")
  len_h <- dim(neighinfo$neighX)[2]
  
  N <- length(unique(G))
  oval_coef_a1 <- reg_coef[, a1, t1][1:(2 + 2* len_h)]
  X_reg <- cbind(1, neigh2_treated, neighX, neigh2_treated_neighX)
  X_fit <- apply(X_reg, 1, function(x) sum(x * oval_coef_a1))
  res_ij <- (H - X_fit) # \mu - \betaX
  trt_ind <- ifelse((A == t1), 1, NA)
  weight_t1 <- weights_ind[, , a1, t1] * trt_ind
  psi_mat <- as.matrix(replicate(dim(X_reg)[2], (weight_t1 * res_ij))) * X_reg # extend to the dimension of X_reg
  psi_data <- as.data.frame(cbind(G, psi_mat))
  psi_group <- aggregate(psi_data[,-1], list(psi_data$G), FUN = mean, na.rm = TRUE)
  colnames(psi_group) <- c("G", "intcp", "neigh2_treated", 
                          colnames(neighX), colnames(neigh2_treated_neighX))
  
  psid_mat <- as.matrix(replicate(dim(X_reg)[2], (weights_ind[, , a1, t1]))) * X_reg
  psid_data <- as.data.frame(cbind(G, A, psid_mat, X_reg))
  psid_mat_grp <- array(dim = c((2 + 2* len_h), (2 + 2* len_h), N)) 
  
  for (w in 1:N){
    psid_data_grp <- psid_data[G == w,]
    psid_mat_grp[, , w] <- t(as.matrix(psid_data_grp[,3:(4+2*len_h)])) %*% 
      as.matrix(psid_data_grp[,(5+2*len_h):(6+4*len_h)])/sum(psid_data_grp$A == t1)
   }
  
  outcomes = list(psi_group, psid_mat_grp)
  outcomes
}

sig_outcome <- function(psi_group, psid_mat_grp){
  N <- length(psi_group[,1])
  V_mat <- crossprod(as.matrix(psi_group[,-1]))/N #rowMeans(apply(as.matrix(psi_group[,-1]), 1, function(x) x %*% t(x)) ) 
  U_mat <- apply(psid_mat_grp, 1:2, mean)
  U_mat_inv <- solve(U_mat)
  var_mat <- U_mat_inv %*% V_mat %*% t(U_mat_inv)/N
  var_mat
}

sig_effect  <- function(psi_group1, psid_mat_grp1, psi_group2, psid_mat_grp2, len_h){
  N <- length(psi_group1[,1])
  psi_group <- cbind(psi_group1[,-1], psi_group2[,-1])
  V_mat <- crossprod(as.matrix(psi_group))/N
  U_mat <- rbind(cbind(apply(psid_mat_grp1, 1:2, mean), matrix(0,2+2*len_h,2+2*len_h)),
                 cbind(matrix(0,2+2*len_h,2+2*len_h), apply(psid_mat_grp2, 1:2, mean)))
  U_mat_inv <- solve(U_mat)
  var_mat <- U_mat_inv %*% V_mat %*% t(U_mat_inv)/N
  var_mat
} 

