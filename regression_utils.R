var_para <- function(reg_coef, a1, t1, X_reg, H, weights_ind, A, G, cat_ind){
  len_n <- dim(X_reg)[2] - 1
  N <- length(unique(G))
  oval_coef_a1 <- reg_coef[, a1, t1]
  #X_reg <- cbind(1, X_num)
  X_fit <- apply(X_reg, 1, function(x) sum(x * oval_coef_a1))
  res_ij <- (H - X_fit) # \mu - \betaX
  trt_ind <- ifelse((A == t1), 1, NA)
  cross_ind <- trt_ind * cat_ind
  weight_t1 <- weights_ind[, , a1, t1] * cross_ind # select A = t1 and x0 = some categorical value
  psi_mat <- as.matrix(replicate(dim(X_reg)[2], (weight_t1 * res_ij))) * X_reg # extend to the dimension of X_reg
  psi_data <- as.data.frame(cbind(G, psi_mat))
  psi_group <- aggregate(psi_data[,-1], list(psi_data$G), FUN = mean, na.rm = TRUE)
  
  #weight_set0 <- weights_ind[, , a1, t1]
  #weight_set0[is.na(cat_ind)] = 0
  #psid_mat <- as.matrix(replicate(dim(X_reg)[2], weight_set0)) * X_reg
  psid_mat <- as.matrix(replicate(dim(X_reg)[2], (weights_ind[, , a1, t1]))) * X_reg
  psid_data <- as.data.frame(cbind(G, A, psid_mat, X_reg))
  psid_mat_grp <- array(dim = c(len_n+1, len_n+1, N)) 
  
  for (w in 1:N){
    psid_data_grp <- psid_data[G == w & !is.na(cat_ind),]
    psid_mat_grp[, , w] <- t(as.matrix(psid_data_grp[,3:(3+len_n)])) %*% 
      as.matrix(psid_data_grp[,(4+len_n):(4+2*len_n)])/sum(psid_data_grp$A == t1)
    # psid_mat_grp[, , w] <- t(as.matrix(psid_data_grp[,3:(3+len_n)])) %*% 
    #    as.matrix(psid_data_grp[,(4+len_n):(4+2*len_n)])/dim(psid_data_grp)[1]
  }
  
  outcomes = list(psi_group, psid_mat_grp, X_fit * cross_ind)
  outcomes
}

# psi_group1 = var_para(reg_coef, a1, t1, X_num, H, weights_ind, A, G)[[1]]
# psid_mat_grp1 = var_para(reg_coef, a1, t1, X_num, H, weights_ind, A, G)[[2]]

sig_outcome <- function(psi_group, psid_mat_grp){
  N <- length(psi_group[,1])
  V_mat <- crossprod(as.matrix(psi_group[,-1]))/N #rowMeans(apply(as.matrix(psi_group[,-1]), 1, function(x) x %*% t(x)) ) 
  U_mat <- apply(psid_mat_grp, 1:2, mean)
  U_mat_inv <- solve(U_mat)
  var_mat <- U_mat_inv %*% V_mat %*% t(U_mat_inv)/N
  var_mat
}

sig_effect  <- function(psi_group1, psid_mat_grp1, psi_group2, psid_mat_grp2, len_n){
  N <- length(psi_group1[,1])
  psi_group <- cbind(psi_group1[,-1], psi_group2[,-1])
  V_mat <- crossprod(as.matrix(psi_group))/N
  U_mat <- rbind(cbind(apply(psid_mat_grp1, 1:2, mean), matrix(0,len_n+1,len_n+1)),
                 cbind(matrix(0,len_n+1,len_n+1), apply(psid_mat_grp2, 1:2, mean)))
  U_mat_inv <- solve(U_mat)
  var_mat <- U_mat_inv %*% V_mat %*% t(U_mat_inv)/N
  var_mat
}  

var_outcome <- function(G, A,  var_mat, X_mean = NULL){ # X_reg, t1, cat_ind, 
  if(is.null(X_mean)){
    #X_reg <- cbind(1, X_num)
    X_data <- as.data.frame(cbind(G, A, X_reg))
    X_data_t1 <- X_data[cat_ind == 1,]#[A == t1, ]
    X_data_group <- aggregate(X_data_t1[,-1:-2], list(X_data_t1$G),
                              FUN = mean, na.rm = TRUE)
    X_mean <- t(colMeans(X_data_group)[-1])
  }
  X_mean <- t(X_mean)
  ave <- X_mean %*% var_mat %*% t(X_mean)
  ave
}

var_effect <- function(G, A, var_mat,  X_mean = NULL){ #X_reg, t1, t2, cat_ind
  if (is.null(X_mean)){
    #X_reg <- cbind(1, X_num)
    X_data <- as.data.frame(cbind(G, A, X_reg))
    X_data_t1 <- X_data[cat_ind == 1,]#[A == t1, ]
    X_data_group_t1 <- aggregate(X_data_t1[,-1:-2], list(X_data_t1$G),
                                  FUN = mean, na.rm = TRUE)
    X_mean_t1 <- colMeans(X_data_group_t1)[-1]
    X_mean_diff <- t(c(X_mean_t1, -X_mean_t1))
    # X_data_t1 <- X_data[cat_ind == 1,]#[A == t1, ]
    # X_data_group_t1 <- aggregate(X_data_t1[,-1:-2], list(X_data_t1$G),
    #                              FUN = mean, na.rm = TRUE)
    # X_mean_t1 <- colMeans(X_data_group_t1)[-1]
    # 
    # X_data_t2 <- X_data[cat_ind == 1,]#[A == t2, ]
    # X_data_group_t2 <- aggregate(X_data_t2[,-1:-2], list(X_data_t2$G),
    #                              FUN = mean, na.rm = TRUE)
    # X_mean_t2 <- colMeans(X_data_group_t2)[-1]
    
  }
  X_mean_diff <- t(c(X_mean, -X_mean))
  ave <- X_mean_diff %*% var_mat %*% t(X_mean_diff)/length(unique(G))
  ave
}