
var_para_mixed <- function(reg_coef, a, t, mixed_info, H, weights_ind, A, G, cat_ind){
  N <- length(unique(G))
  oval_coef_a <- reg_coef[, a, t]
  len_m <- length(oval_coef_a)
  
  X_fit <- apply(mixed_info, 1, function(x) sum(x * oval_coef_a))
  res_ij <- (H - X_fit) # \mu - \betaX
  trt_ind <- ifelse((A == t), 1, NA)
  cross_ind <- trt_ind * cat_ind
  weight_t <- weights_ind[, , a, t] * cross_ind # select A = t and x0 = some categorical value
  psi_mat <- as.matrix(replicate(dim(mixed_info)[2], (weight_t * res_ij))) * mixed_info # extend to the dimension of X_reg
  psi_data <- as.data.frame(cbind(G, psi_mat))
  psi_group <- aggregate(psi_data[,-1], list(psi_data$G), FUN = mean, na.rm = TRUE)
  
  psid_mat <- as.matrix(replicate(dim(mixed_info)[2], (weights_ind[, , a, t]))) * mixed_info
  psid_data <- as.data.frame(cbind(G, A, psid_mat, mixed_info, cat_ind))
  psid_mat_grp <- array(dim = c(len_m, len_m, N)) 
  
  for (w in 1:N){
    psid_data_grp <- psid_data[psid_data$G == w & !is.na(psid_data$cat_ind),]
    psid_mat_grp[, , w] <- t(as.matrix(psid_data_grp[,3:(2+len_m)])) %*% 
      as.matrix(psid_data_grp[,(3+len_m):(2+2*len_m)])/sum(psid_data_grp$A == t)
    }
  
  outcomes = list(psi_group, psid_mat_grp, X_fit * cross_ind)
  outcomes
}

sig_outcome_mixed <- function(psi_group, psid_mat_grp){
  N <- length(psi_group[,1])
  V_mat <- crossprod(as.matrix(psi_group[,-1]))/N 
  U_mat <- apply(psid_mat_grp, 1:2, mean, na.rm = TRUE)
  U_mat_inv <- solve(U_mat)
  var_mat <- U_mat_inv %*% V_mat %*% t(U_mat_inv)/N
  var_mat
}

sig_effect_mixed  <- function(psi_group1, psid_mat_grp1, psi_group2, psid_mat_grp2, len_m){
  N <- length(psi_group1[,1])
  psi_group <- cbind(psi_group1[,-1], psi_group2[,-1])
  V_mat <- crossprod(as.matrix(psi_group))/N
  U_mat <- rbind(cbind(apply(psid_mat_grp1, 1:2, mean, na.rm = TRUE), matrix(0,len_m,len_m)),
                 cbind(matrix(0,len_m,len_m), apply(psid_mat_grp2, 1:2, mean, na.rm = TRUE)))
  U_mat_inv <- solve(U_mat)
  var_mat <- U_mat_inv %*% V_mat %*% t(U_mat_inv)/N
  var_mat
}  

var_outcome_mixed <- function(G, A,  var_mat, cat_ind, mixed_info, X_mean = NULL){ 
  if(is.null(X_mean)){
    X_data <- as.data.frame(cbind(G, A, mixed_info))
    X_data_t1 <- X_data[cat_ind == 1,]
    X_data_group <- aggregate(X_data_t1[,-1:-2], list(X_data_t1$G),
                              FUN = mean, na.rm = TRUE)
    X_mean <- t(colMeans(X_data_group)[-1])
  }
  X_mean <- t(X_mean)
  ave <- X_mean %*% var_mat %*% t(X_mean)
  ave
}

var_effect_mixed <- function(G, A, var_mat, cat_ind, mixed_info, X_mean = NULL){ #X_reg, t1, t2, cat_ind
  if (is.null(X_mean)){
    X_data <- as.data.frame(cbind(G, A, mixed_info))
    X_data_t1 <- X_data[cat_ind == 1,]#[A == t1, ]
    X_data_group_t1 <- aggregate(X_data_t1[,-1:-2], list(X_data_t1$G),
                                 FUN = mean, na.rm = TRUE)
    X_mean_t1 <- colMeans(X_data_group_t1)[-1]
    X_mean_diff <- t(c(X_mean_t1, -X_mean_t1))
    }
  X_mean_diff <- t(c(X_mean, -X_mean))
  ave <- X_mean_diff %*% var_mat %*% t(X_mean_diff)#/length(unique(G))
  ave
}


