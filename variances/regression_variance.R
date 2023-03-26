source("variances/regression_utils.R")
source("variances/regression_utils_neigh.R")
source("variances/regression_utils_mixed.R")

# no implementation of marginal yet
# regression variance for conditional within neighbors
# for this function, input is H_M, x1_num is the condition then
ipw_regression_variance_neigh <- function(H,
                                    weights, 
                                    point_estimates, 
                                    A,
                                    effect_type, 
                                    marginal,
                                    allocation1, 
                                    allocation2 = NA,
                                    t1 = 0, # default of outcome variance is control group
                                    t2 = 1, # default of contrast is t1(0) - t2(1)
                                    neighinfo = NULL,
                                    x1_num = NULL,
                                    rescale.factor = 1,
                                    conf.level = 0.95,
                                    print = FALSE){  
  
  ## Necessary bits ##
  N  <- dim(weights)[1] 
  p  <- dim(weights)[2] 
  k  <- length(allocations)
  l  <- dim(point_estimates$outcomes$overall)[2]
  a1 <- as.character(allocation1)
  a2 <- as.character(allocation2)
  t1 <- as.character(t1)
  t2 <- as.character(t2)
  #t1 <- as.character(trt.lvl1)
  #t2 <- as.character(trt.lvl2)
  
  len_h <- dim(neighinfo$neighX)[2]
  
  fff <- ifelse(marginal == TRUE, 'marginal_outcomes', 'outcomes')
  #ind_est_df <- point_estimates[[fff]]$weighted_ind
  reg_coef <- point_estimates[[fff]]$overall_coefH
  #reg_coef <- point_estimates[[fff]]$grp_coefH
  weights_ind <- point_estimates[[fff]]$weights_ind  
  X_mean <- c(1, as.numeric(x1_num))
  
  if(effect_type == 'contrast'){ 
    trt_ind_t1 <- ifelse((A == t1), 1, NA)
    trt_ind_t2 <- ifelse((A == t2), 1, NA)
    if(marginal == TRUE){
      # pe          <- oal[a1] - oal[a2]
      # pe_grp_diff <- (grp[ , a1] - oal[a1]) - (grp[, a2] - oal[a2])
      # if(variance_estimation == 'robust'){
      #   U_pe_grp    <- Ugrp[ , , a1] - Ugrp[ , , a2]
      # }
    } else {
      #reg_coef, a1, t1, neighinfo, H, weights_ind, A, G
      psi_group1 = var_para_neigh(reg_coef, a1, t1, neighinfo, H, weights_ind, A, G)[[1]]
      psid_mat_grp1 = var_para_neigh(reg_coef, a1, t1, neighinfo, H, weights_ind, A, G)[[2]]
      
      psi_group2 = var_para_neigh(reg_coef, a2, t2, neighinfo, H, weights_ind, A, G)[[1]]
      psid_mat_grp2 = var_para_neigh(reg_coef, a2, t2, neighinfo, H, weights_ind, A, G)[[2]]
      
      var_mat = sig_effect_neigh(psi_group1, psid_mat_grp1, psi_group2, psid_mat_grp2, len_h)
      ave <- var_effect_neigh(G, A, var_mat, neighinfo, X_mean)
      
      # X_fit1 = var_para(reg_coef, a1, t1, X_reg, H, weights_ind, A, G, cat_ind)[[3]] #TODO: all categorical, reg_coef is uncorrect?
      # X_fit2 = var_para(reg_coef, a2, t2, X_reg, H, weights_ind, A, G, cat_ind)[[3]]
      # X_data = as.data.frame(cbind(G, X_fit1, X_fit2))
      # X_group = aggregate(X_data[,-1], list(X_data$G), FUN=mean, na.rm = TRUE) # mean of residual for each group i
      pe <- sum(X_mean * reg_coef[,a1,t1]) - sum(X_mean * reg_coef[,a2,t2])
      #pe <- sum(X_mean * colMeans(reg_coef)[,a1,t1]) - sum(X_mean * colMeans(reg_coef)[,a2,t2])
    }
  } else {
    # default is control-outcome
    if(marginal == TRUE){
      # pe          <- oal[a1] 
      # pe_grp_diff <- (grp[ , a1] - oal[a1])
      # if(variance_estimation == 'robust'){
      #   U_pe_grp    <- Ugrp[ , , a1]
      # }
    } else {
      psi_group = var_para_neigh(reg_coef, a1, t1, neighinfo, H, weights_ind, A, G)[[1]]
      psid_mat_grp = var_para_neigh(reg_coef, a1, t1, neighinfo, H, weights_ind, A, G)[[2]]
      var_mat = sig_outcome_neigh(psi_group, psid_mat_grp)
      ave <- var_outcome_neigh(G, A, var_mat, neighinfo, X_mean)
      
      # X_fit = var_para(reg_coef, a1, t1, X_reg, H, weights_ind, A, G, cat_ind)[[3]]
      # X_data = as.data.frame(cbind(G, X_fit))
      # X_group = aggregate(X_data$X_fit, list(X_data$G), FUN=mean, na.rm = TRUE) # mean of residual for each group i

      pe <- sum(X_mean * reg_coef[,a1,t1])
      #pe <- sum(X_mean * colMeans(reg_coef)[,a1,t1])
    }
  }
  
  #ave <- (1/(N^2)) * (sum((pe_grp_diff)^2, na.rm = T)) * rescale.factor^2
  #ave <- X_mean %*% var_mat %*% t(X_mean) * rescale.factor^2 # variance of outcome
  
  ## Confidence Intervals ##
  qq <- qnorm(conf.level + (1 - conf.level)/2)
  me <- qq * sqrt(ave)
  
  ## Prepare Output ##
  pe <- pe * rescale.factor
  
  if(print == TRUE){
    toprint <- paste0('Estimate: ', round(pe, 2), ' ',
                      conf.level*100, '% CI: (', 
                      round(pe - me, 2), ', ', round(pe + me, 2), ')' )
    print(toprint)
  }
  
  out <- data.frame(estimate = pe,
                    std.error = sqrt(ave), 
                    conf.low = pe - me, conf.high = pe + me)
  return(out)
}


ipw_regression_variance <- function(H,
                                    weights, 
                                    point_estimates, 
                                    A,
                                    effect_type, 
                                    marginal,
                                    allocation1, 
                                    allocation2 = NA,
                                    #trt.lvl1 = 0, 
                                    #trt.lvl2 = 1, 
                                    t1 = 0, # default of outcome variance is control group
                                    t2 = 1, # default of contrast is t1(0) - t2(1)
                                    X = NULL,
                                    X_type = NULL, 
                                    x0 = NULL,
                                    rescale.factor = 1,
                                    conf.level = 0.95,
                                    print = FALSE){  
  
  ## Necessary bits ##
  N  <- dim(weights)[1] 
  p  <- dim(weights)[2] 
  k  <- length(allocations)
  l  <- dim(point_estimates$outcomes$overall)[2]
  a1 <- as.character(allocation1)
  a2 <- as.character(allocation2)
  t1 <- as.character(t1)
  t2 <- as.character(t2)
  #t1 <- as.character(trt.lvl1)
  #t2 <- as.character(trt.lvl2)
  
  ## group coefficients conditional on x0_cat ##
  X_cat <- as.matrix(X[, X_type == "C"])
  x0_cat <- x0[X_type == "C"]
  X_num <- apply(as.matrix(X[, X_type == "N"]), 2, as.numeric)
  len_n <- sum(X_type == "N")
  X_mean <- c(1, as.numeric(x0[X_type == "N"]))
  ## categorical indicator ##
  if (sum(X_type == "C") > 0){
    cat_ind <- apply(X_cat, 1, function(x) ifelse(prod(x == x0_cat), 1, NA)) 
  }else{
    cat_ind <- rep(1,dim(X_cat)[1])} #X = X0 categorical
  
  ## group coefficients for marginal / non-marginal (a = 0 / 1 separately) ##
  #grp_coef <- ipw_point_estimates_mixed(H, G, A, w.matrix, X = X, x0 = x0, X_type = X_type)$outcomes$grp_coefG
  
  fff <- ifelse(marginal == TRUE, 'marginal_outcomes', 'outcomes')
  ind_est_df <- point_estimates[[fff]]$weighted_ind
  reg_coef <- point_estimates[[fff]]$overall_coefG
  weights_ind <- point_estimates[[fff]]$weights_ind  
  if(effect_type == 'contrast'){ 
    X_reg <- as.matrix(cbind(rep(1,dim(X)[1]), X_num))
    trt_ind_t1 <- ifelse((A == t1), 1, NA)
    trt_ind_t2 <- ifelse((A == t2), 1, NA)
    if(marginal == TRUE){
      # pe          <- oal[a1] - oal[a2]
      # pe_grp_diff <- (grp[ , a1] - oal[a1]) - (grp[, a2] - oal[a2])
      # if(variance_estimation == 'robust'){
      #   U_pe_grp    <- Ugrp[ , , a1] - Ugrp[ , , a2]
      # }
    } else {
      #X_fit1 = apply(X_reg, 1, function(x) sum(x * reg_coef[,a1,t1])) *trt_ind_t1 #TODO: all categorical, reg_coef is uncorrect?
      #X_fit2 = apply(X_reg, 1, function(x) sum(x * reg_coef[,a2,t2])) *trt_ind_t2
      #X_data = as.data.frame(cbind(G, X_fit1, X_fit2))
      #X_group = aggregate(X_data[,-1], list(X_data$G), FUN=mean, na.rm = TRUE) # mean of residual for each group i
      #pe <- colMeans(X_group)[2] -  colMeans(X_group)[3]
      psi_group1 = var_para(reg_coef, a1, t1, X_reg, H, weights_ind, A, G, cat_ind)[[1]]
      psid_mat_grp1 = var_para(reg_coef, a1, t1, X_reg, H, weights_ind, A, G, cat_ind)[[2]]
      
      psi_group2 = var_para(reg_coef, a2, t2, X_reg, H, weights_ind, A, G, cat_ind)[[1]]
      psid_mat_grp2 = var_para(reg_coef, a2, t2, X_reg, H, weights_ind, A, G, cat_ind)[[2]]
      
      var_mat = sig_effect(psi_group1, psid_mat_grp1, psi_group2, psid_mat_grp2, len_n)
      #ave <- var_effect(G, A, X_reg, t1, t2, var_mat, cat_ind)
      ave <- var_effect(G, A, var_mat, cat_ind, X_reg, X_mean = X_mean)
      
      X_fit1 = var_para(reg_coef, a1, t1, X_reg, H, weights_ind, A, G, cat_ind)[[3]] #TODO: all categorical, reg_coef is uncorrect?
      X_fit2 = var_para(reg_coef, a2, t2, X_reg, H, weights_ind, A, G, cat_ind)[[3]]
      X_data = as.data.frame(cbind(G, X_fit1, X_fit2))
      X_group = aggregate(X_data[,-1], list(X_data$G), FUN=mean, na.rm = TRUE) # mean of residual for each group i
      #pe <- colMeans(X_group)[2] -  colMeans(X_group)[3]
      pe <- sum(X_mean * reg_coef[,a1,t1]) - sum(X_mean * reg_coef[,a2,t2])
    }
  } else {
    # default is control-outcome
    X_reg <- as.matrix(cbind(rep(1,dim(X)[1]), X_num))
    #trt_ind <- ifelse((A == t1), 1, NA) 
    if(marginal == TRUE){
      # pe          <- oal[a1] 
      # pe_grp_diff <- (grp[ , a1] - oal[a1])
      # if(variance_estimation == 'robust'){
      #   U_pe_grp    <- Ugrp[ , , a1]
      # }
    } else {
      # X_fit = apply(X_reg, 1, function(x) sum(x * reg_coef[,a1,t1]))
      # X_data = as.data.frame(cbind(G, X_fit))
      # X_group = aggregate(X_data$X_fit, list(X_data$G), FUN=mean, na.rm = TRUE) # mean of residual for each group i
      # pe <- mean(X_group$x)
      
      psi_group = var_para(reg_coef, a1, t1, X_reg, H, weights_ind, A, G, cat_ind)[[1]]
      psid_mat_grp = var_para(reg_coef, a1, t1, X_reg, H, weights_ind, A, G, cat_ind)[[2]]
      var_mat = sig_outcome(psi_group, psid_mat_grp)
      #ave <- var_outcome(G, A, X_reg, t1, var_mat, cat_ind)
      ave <- var_outcome(G, A, var_mat, cat_ind, X_reg, X_mean = X_mean)
      
      X_fit = var_para(reg_coef, a1, t1, X_reg, H, weights_ind, A, G, cat_ind)[[3]]
      X_data = as.data.frame(cbind(G, X_fit))
      X_group = aggregate(X_data$X_fit, list(X_data$G), FUN=mean, na.rm = TRUE) # mean of residual for each group i
      #pe <- mean(X_group$x)
      pe <- sum(X_mean * reg_coef[,a1,t1])
    }
  }
  
  #ave <- (1/(N^2)) * (sum((pe_grp_diff)^2, na.rm = T)) * rescale.factor^2
  #ave <- X_mean %*% var_mat %*% t(X_mean) * rescale.factor^2 # variance of outcome
  
  ## Confidence Intervals ##
  qq <- qnorm(conf.level + (1 - conf.level)/2)
  me <- qq * sqrt(ave)
  
  ## Prepare Output ##
  pe <- pe * rescale.factor
  
  if(print == TRUE){
    toprint <- paste0('Estimate: ', round(pe, 2), ' ',
                      conf.level*100, '% CI: (', 
                      round(pe - me, 2), ', ', round(pe + me, 2), ')' )
    print(toprint)
  }
  
  out <- data.frame(estimate = pe,
                    std.error = sqrt(ave), 
                    conf.low = pe - me, conf.high = pe + me)
  return(out)
}


# ipw_regression_variance <- function(weights, 
#                             point_estimates, 
#                             effect_type, 
#                             marginal,
#                             allocation1, 
#                             allocation2 = NA,
#                             trt.lvl1 = 0, 
#                             trt.lvl2 = 1, 
#                             t1 = 0, # default of outcome variance is control group
#                             t2 = 1, # default of contrast is t1(0) - t2(1)
#                             X = NULL,
#                             X_type = NULL, 
#                             x0 = NULL,
#                             rescale.factor = 1,
#                             conf.level = 0.95,
#                             print = FALSE){  
#   
#   ## Necessary bits ##
#   N  <- dim(weights)[1] 
#   p  <- dim(weights)[2] 
#   k  <- length(allocations)
#   l  <- dim(point_estimates$outcomes$overall)[2]
#   a1 <- as.character(allocation1)
#   a2 <- as.character(allocation2)
#   t1 <- as.character(trt.lvl1)
#   t2 <- as.character(trt.lvl2)
#   
#   ## group coefficients conditional on x0_cat ##
#   X_cat <- as.matrix(X[, X_type == "C"])
#   x0_cat <- x0[X_type == "C"]
#   X_num <- apply(as.matrix(X[, X_type == "N"]), 2, as.numeric)
#   len_n <- sum(X_type == "N")
#   
#   ## categorical indicator ##
#   if (sum(X_type == "C") > 0){
#     cat_ind <- apply(X_cat, 1, function(x) ifelse(prod(x == x0_cat), 1, NA)) 
#   }else{
#     cat_ind <- rep(1,dim(X_cat)[1])} #X = X0 categorical
#   
#   ## group coefficients for marginal / non-marginal (a = 0 / 1 separately) ##
#   #grp_coef <- ipw_point_estimates_mixed(H, G, A, w.matrix, X = X, x0 = x0, X_type = X_type)$outcomes$grp_coefG
# 
#   fff <- ifelse(marginal == TRUE, 'marginal_outcomes', 'outcomes')
#   #oal  <- point_estimates[[fff]]$overall_coefG
#   #grp  <- point_estimates[[fff]]$grp_coefG
#   ind_est_df <- point_estimates[[fff]]$weighted_ind
#   reg_coef <- point_estimates[[fff]]$overall_coefG
#   weights_ind <- point_estimates[[fff]]$weights_ind  
#   if(effect_type == 'contrast'){ 
#     #X_reg <- cbind(G, 1, X_num)
#     X_reg <- cbind(1, X_num)
#     if(marginal == TRUE){
#       # pe          <- oal[a1] - oal[a2]
#       # pe_grp_diff <- (grp[ , a1] - oal[a1]) - (grp[, a2] - oal[a2])
#       # if(variance_estimation == 'robust'){
#       #   U_pe_grp    <- Ugrp[ , , a1] - Ugrp[ , , a2]
#       # }
#     } else {
#       #grp_coef_a1 <- grp[, , a1, t1]
#       #X_fit_a1 <- apply(X_reg, 1, function(x) sum(x[-1] * grp_coef_a1[x[1],]))
#       oval_coef_a1 <- reg_coef[, a1, t1]
#       X_fit_a1 <- apply(X_reg, 1, function(x) sum(x * oval_coef_a1))
#       trt_ind_a1 <- (A == t1) 
#       X_fit_trt_a1 <- X_fit_a1 * trt_ind_a1
#       res_ij_a1 <- (ind_est_df[, , a1, t1] - X_fit_trt_a1) 
#       res_data_a1 <- as.data.frame(cbind(res_ij_a1, G))
#       res_group_a1 <- aggregate(res_data_a1$res_ij_a1, list(res_data_a1$G), FUN=mean, na.rm = TRUE)
#       
#       grp_coef_a2 <- grp[, , a2, t2]
#       X_fit_a2 <- apply(X_reg, 1, function(x) sum(x[-1] * grp_coef_a2[x[1],]))
#       trt_ind_a2 <- (A == t2) 
#       X_fit_trt_a2 <- X_fit_a2 * trt_ind_a2
#       res_ij_a2 <- (ind_est_df[, , a2, t2] - X_fit_trt_a2)
#       res_data_a2 <- as.data.frame(cbind(res_ij_a2, G))
#       res_group_a2 <- aggregate(res_data_a2$res_ij_a2, list(res_data_a2$G), FUN=mean, na.rm = TRUE)
#       
#       pe <- mean(res_group_a1$x) - mean(res_group_a2$x)
#       pe_grp_diff <- (res_group_a1$x - mean(res_group_a1$x)) - (res_group_a2$x - mean(res_group_a2$x))
#       
#       # pe          <- oal[a1, t1] - oal[a2, t2] 
#       # pe_grp_diff <- (grp[ , a1, t1] - oal[a1, t1]) - (grp[ , a2, t2] - oal[a2, t2])
#       # if(variance_estimation == 'robust'){
#       #   U_pe_grp    <- Ugrp[ , , a1, t1] - Ugrp[ , , a2, t2]
#       # }
#     }
#   } else {
#     # default is control-outcome
#     X_reg <- cbind(1, X_num)
#     trt_ind <- ifelse((A == t1), 1, NA) 
#     if(marginal == TRUE){
#       # pe          <- oal[a1] 
#       # pe_grp_diff <- (grp[ , a1] - oal[a1])
#       # if(variance_estimation == 'robust'){
#       #   U_pe_grp    <- Ugrp[ , , a1]
#       # }
#     } else {
#       oval_coef_a1 <- reg_coef[, a1, t1]
#       X_fit <- apply(X_reg, 1, function(x) sum(x * oval_coef_a1))
#       res_ij <- (H - X_fit) # \mu - \betaX
#       weight_t1 <- weights_ind[, , a1, t1] * trt_ind
#       psi_mat <- as.matrix(replicate(dim(X_reg)[2], (weight_t1 * res_ij))) * X_reg # extend to the dimension of X_reg
#       psi_data <- as.data.frame(cbind(G, psi_mat))
#       psi_group <- aggregate(psi_data[,-1], list(psi_data$G), FUN = mean, na.rm = TRUE)
#       V_mat <- crossprod(as.matrix(psi_group[,-1]))/N
#       
#       psid_mat <- as.matrix(replicate(dim(X_reg)[2], (weights_ind[, , a1, t1]))) * X_reg
#       psid_data <- as.data.frame(cbind(G, A, psid_mat, X_reg))
#       psid_mat_grp <- array(dim = c(len_n+1, len_n+1, N)) 
#       
#       for (w in 1:N){
#         psid_data_grp <- psid_data[G == w,]
#         psid_mat_grp[, , w] <- t(as.matrix(psid_data_grp[,3:(3+len_n)])) %*% 
#            as.matrix(psid_data_grp[,(4+len_n):(4+2*len_n)])/sum(psid_data_grp$A == t1)
#         # psid_mat_grp[, , w] <- t(as.matrix(psid_data_grp[,3:(3+len_n)])) %*% 
#         #    as.matrix(psid_data_grp[,(4+len_n):(4+2*len_n)])/dim(psid_data_grp)[1]
#         
#       }
#       
#       U_mat <- apply(psid_mat_grp, 1:2, mean)
#       U_mat_inv <- solve(U_mat)
#       var_mat <- U_mat_inv %*% V_mat %*% t(U_mat_inv)/N
#       
#       
#       X_data <- as.data.frame(cbind(G, A, X_reg))
#       X_data_t1 <- X_data[A == t1, ]
#       X_data_group <- aggregate(X_data_t1[,-1:-2], list(X_data_t1$G), 
#                                 FUN = mean, na.rm = TRUE)
#       X_mean <- t(colMeans(X_data_group)[-1])
#       #var_fitted <- X_mean %*% var_mat %*% t(X_mean)
# 
#       # res_data <- as.data.frame(cbind(res_ij, G))
#       # res_group <- aggregate(res_data$res_ij, list(res_data$G), FUN=mean, na.rm = TRUE) # mean of residual for each group i
#       # pe <- mean(res_group$x)
#       # pe_grp_diff <- res_group$x - pe
#       #var <-  mean((res_group$x - mean(res_group$x))^2) 
#       # pe          <- ifelse(k == 1, oal[t1], oal[a1, t1]) 
#       # pe_grp_diff <- ifelse(k == 1, (grp[ , t1] - oal[t1]),
#       #                       (grp[ , a1, t1] - oal[a1, t1]))
#       # if(variance_estimation == 'robust'){
#       #   U_pe_grp    <- Ugrp[ , , a1, t1]
#       # }
#       
#     }
#   }
#   
#   #ave <- (1/(N^2)) * (sum((pe_grp_diff)^2, na.rm = T)) * rescale.factor^2
#   ave <- X_mean %*% var_mat %*% t(X_mean) * rescale.factor^2 # variance of outcome
#   
#   ## Confidence Intervals ##
#   qq <- qnorm(conf.level + (1 - conf.level)/2)
#   me <- qq * sqrt(ave)
#   
#   ## Prepare Output ##
#   pe <- pe * rescale.factor
#   
#   if(print == TRUE){
#     toprint <- paste0('Estimate: ', round(pe, 2), ' ',
#                       conf.level*100, '% CI: (', 
#                       round(pe - me, 2), ', ', round(pe + me, 2), ')' )
#     print(toprint)
#   }
#   
#   out <- data.frame(estimate = pe,
#                     std.error = sqrt(ave), 
#                     conf.low = pe - me, conf.high = pe + me)
#   return(out)
# }

# X = NULL,
# X_type = NULL, 
# x0 = NULL,
#qnames <- paste0("x0 = ", toString(x0))

ipw_m_variance <- function(weights, 
                           point_estimates, 
                           effect_type, 
                           marginal,
                           allocation1, 
                           allocation2 = NA,
                           trt.lvl1 = 0, 
                           trt.lvl2 = 1, 
                           t1 = 0, # default of outcome variance is control group
                           t2 = 1, # default of contrast is t1(0) - t2(1)
                           rescale.factor = 1,
                           conf.level = 0.95,
                           print = FALSE){  
  
  ## Necessary bits ##
  N  <- dim(weights)[1] 
  p  <- dim(weights)[2] 
  k  <- length(allocations)
  l  <- dim(point_estimates$outcomes$overall)[2]
  a1 <- as.character(allocation1)
  a2 <- as.character(allocation2)
  t1 <- as.character(trt.lvl1)
  t2 <- as.character(trt.lvl2)
  
  fff <- ifelse(marginal == TRUE, 'marginal_outcomes', 'outcomes')
  
  oal  <- point_estimates[[fff]]$overall
  grp  <- point_estimates[[fff]]$groups
  #oal  <- apply(grp, 2:4, mean, na.rm = TRUE)
    
  if(effect_type == 'contrast'){ # TODO: change the name and check for a1 a2 to be the same allocation.
    if(marginal == TRUE){
      pe          <- oal[a1,] - oal[a2,]
      pe_grp_diff <- (grp[ , a1,] - oal[a1,]) - (grp[, a2,] - oal[a2,])
      # if(variance_estimation == 'robust'){
      #   U_pe_grp    <- Ugrp[ , , a1] - Ugrp[ , , a2]
      # }
    } else {
        if (k ==1){
        pe <-  oal[a1, t1,] - oal[a2, t2,]
        pe_grp_diff <- (grp[ , a1, t1,] - oal[a1, t1,]) - (grp[ , a2, t2,] - oal[a2, t2,])
      }else{
        pe <-  oal[a1, t1] - oal[a2, t2]
        pe_grp_diff <- (grp[ , a1, t1] - oal[a1, t1]) - (grp[ , a2, t2] - oal[a2, t2])
      }
      # pe          <- ifelse(k == 1, oal[a1, t1,] - oal[a2, t2,],
      #                       oal[a1, t1] - oal[a2, t2])
      # pe_grp_diff <- ifelse(k == 1, (grp[ , a1, t1,] - oal[a1, t1,]) - (grp[ , a2, t2,] - oal[a2, t2,]),
      #                       (grp[ , a1, t1] - oal[a1, t1]) - (grp[ , a2, t2] - oal[a2, t2]))
      
      # if(variance_estimation == 'robust'){
      #   U_pe_grp    <- Ugrp[ , , a1, t1] - Ugrp[ , , a2, t2]
      # }
    }
  } else {
    # default is control-outcome
    if(marginal == TRUE){
      pe          <- oal[a1] 
      pe_grp_diff <- (grp[ , a1] - oal[a1])
      # if(variance_estimation == 'robust'){
      #   U_pe_grp    <- Ugrp[ , , a1]
      # }
    } else {
      if (k == 1){
        pe          <- oal[a1, t1, ]
        pe_grp_diff <- grp[ , a1, t1, ] - oal[a1, t1, ]
      }else{
        pe          <- oal[a1, t1]
        pe_grp_diff <- grp[ , a1, t1] - oal[a1, t1]
      }
      # if(variance_estimation == 'robust'){
      #   U_pe_grp    <- Ugrp[ , , a1, t1]
      # }
    }
  }
  
  ave <- (1/(N^2)) * (sum((pe_grp_diff)^2, na.rm = T)) * rescale.factor^2
  
  ## Confidence Intervals ##
  qq <- qnorm(conf.level + (1 - conf.level)/2)
  me <- qq * sqrt(ave)
  
  ## Prepare Output ##
  pe <- pe * rescale.factor
  
  if(print == TRUE){
    toprint <- paste0('Estimate: ', round(pe, 2), ' ',
                      conf.level*100, '% CI: (', 
                      round(pe - me, 2), ', ', round(pe + me, 2), ')' )
    print(toprint)
  }
  
  out <- data.frame(estimate = pe,
                    std.error = sqrt(ave), 
                    conf.low = pe - me, conf.high = pe + me)
  return(out)
}

# ipw_m_variance <- Vectorize(ipw_m_variance, 
#                              vectorize.args = c("allocation1", 'allocation2', 'trt.lvl1', 'trt.lvl2', 
#                                                 'marginal', 'effect_type'),
#                              SIMPLIFY = TRUE)

ipw_m_variance_groups <- function(weights, 
                           point_estimates, 
                           effect_type, 
                           marginal,
                           allocation1, 
                           allocation2 = NA,
                           trt.lvl1 = 0, 
                           trt.lvl2 = 1, 
                           t1 = 0, # default of outcome variance is control group
                           t2 = 1, # default of contrast is t1(0) - t2(1)
                           rescale.factor = 1,
                           conf.level = 0.95,
                           print = FALSE){  
  
  ## Necessary bits ##
  N  <- dim(weights)[1] 
  p  <- dim(weights)[2] 
  k  <- length(allocations)
  l  <- dim(point_estimates$outcomes$overall)[2]
  a1 <- as.character(allocation1)
  a2 <- as.character(allocation2)
  t1 <- as.character(trt.lvl1)
  t2 <- as.character(trt.lvl2)
  
  fff <- ifelse(marginal == TRUE, 'marginal_outcomes', 'outcomes')
  
  #oal  <- point_estimates[[fff]]$overall
  grp  <- point_estimates[[fff]]$groups
  oal  <- apply(grp, 2:4, mean, na.rm = TRUE)
  
  if(effect_type == 'contrast'){ # TODO: change the name and check for a1 a2 to be the same allocation.
    if(marginal == TRUE){
      pe          <- oal[a1,] - oal[a2,]
      pe_grp_diff <- (grp[ , a1,] - oal[a1,]) - (grp[, a2,] - oal[a2,])
      # if(variance_estimation == 'robust'){
      #   U_pe_grp    <- Ugrp[ , , a1] - Ugrp[ , , a2]
      # }
    } else {
      if (k ==1){
        pe <-  oal[a1, t1,] - oal[a2, t2,]
        pe_grp_diff <- (grp[ , a1, t1,] - oal[a1, t1,]) - (grp[ , a2, t2,] - oal[a2, t2,])
      }else{
        pe <-  oal[a1, t1] - oal[a2, t2]
        pe_grp_diff <- (grp[ , a1, t1] - oal[a1, t1]) - (grp[ , a2, t2] - oal[a2, t2])
      }
      # pe          <- ifelse(k == 1, oal[a1, t1,] - oal[a2, t2,],
      #                       oal[a1, t1] - oal[a2, t2])
      # pe_grp_diff <- ifelse(k == 1, (grp[ , a1, t1,] - oal[a1, t1,]) - (grp[ , a2, t2,] - oal[a2, t2,]),
      #                       (grp[ , a1, t1] - oal[a1, t1]) - (grp[ , a2, t2] - oal[a2, t2]))
      
      # if(variance_estimation == 'robust'){
      #   U_pe_grp    <- Ugrp[ , , a1, t1] - Ugrp[ , , a2, t2]
      # }
    }
  } else {
    # default is control-outcome
    if(marginal == TRUE){
      pe          <- oal[a1] 
      pe_grp_diff <- (grp[ , a1] - oal[a1])
      # if(variance_estimation == 'robust'){
      #   U_pe_grp    <- Ugrp[ , , a1]
      # }
    } else {
      if (k == 1){
        pe          <- oal[a1, t1, ]
        pe_grp_diff <- grp[ , a1, t1, ] - oal[a1, t1, ]
      }else{
        pe          <- oal[a1, t1]
        pe_grp_diff <- grp[ , a1, t1] - oal[a1, t1]
      }
      # if(variance_estimation == 'robust'){
      #   U_pe_grp    <- Ugrp[ , , a1, t1]
      # }
    }
  }
  
  ave <- (1/(N^2)) * (sum((pe_grp_diff)^2, na.rm = T)) * rescale.factor^2
  
  ## Confidence Intervals ##
  qq <- qnorm(conf.level + (1 - conf.level)/2)
  me <- qq * sqrt(ave)
  
  ## Prepare Output ##
  pe <- pe * rescale.factor
  
  if(print == TRUE){
    toprint <- paste0('Estimate: ', round(pe, 2), ' ',
                      conf.level*100, '% CI: (', 
                      round(pe - me, 2), ', ', round(pe + me, 2), ')' )
    print(toprint)
  }
  
  out <- data.frame(estimate = pe,
                    std.error = sqrt(ave), 
                    conf.low = pe - me, conf.high = pe + me)
  return(out)
}




ipw_regression_variance_mixed <- function(H,
                                    weights, 
                                    point_estimates, 
                                    A,
                                    effect_type, 
                                    marginal,
                                    allocation1, 
                                    allocation2 = NA,
                                    t1 = 0, 
                                    t2 = 1, 
                                    X = NULL,
                                    X_type = NULL, 
                                    x0 = NULL,
                                    neighinfo = NULL,
                                    x1_num = NULL,
                                    rescale.factor = 1,
                                    conf.level = 0.95,
                                    print = FALSE){  
  
  ## Necessary bits ##
  N  <- dim(weights)[1] 
  p  <- dim(weights)[2] 
  k  <- length(allocations)
  l  <- dim(point_estimates$outcomes$overall)[2]
  a1 <- as.character(allocation1)
  a2 <- as.character(allocation2)
  t1 <- as.character(t1)
  t2 <- as.character(t2)
  
  ## group coefficients conditional on x0_cat ##
  X_cat <- as.matrix(X[, X_type == "C"])
  x0_cat <- x0[X_type == "C"]
  
  X_num <- apply(as.matrix(X[, X_type == "N"]), 2, as.numeric)
  X_inter <- column_multiply(X_num, neighinfo$neighX)
  X_neighX <- X_inter[[1]]
  mixed_info <- cbind(1, X_num, neighinfo$neighX, X_neighX)
  len_m <- dim(mixed_info)[2]
  
  x0_num <- as.numeric(x0[X_type == "N"])
  
  X_inter_num <- c()
  for (j in 1:length(x0_num)) {
    X_inter_num <- c(X_inter_num,
                     apply(as.array(x1_num), 1, function(x) x*x0_num[j]))}
  X_mean <- c(1, x0_num, x1_num, X_inter_num)
  
  ## categorical indicator ##
  if (sum(X_type == "C") > 0){
    cat_ind <- apply(X_cat, 1, function(x) ifelse(prod(x == x0_cat), 1, NA)) 
  }else{
    cat_ind <- rep(1,dim(X_cat)[1])} 
  
  ## group coefficients for marginal / non-marginal (a = 0 / 1 separately) ##
  fff <- ifelse(marginal == TRUE, 'marginal_outcomes', 'outcomes')
  ind_est_df <- point_estimates[[fff]]$weighted_ind
  reg_coef <- point_estimates[[fff]]$overall_coefM
  weights_ind <- point_estimates[[fff]]$weights_ind  
  
  if(effect_type == 'contrast'){ 
    X_reg <- as.matrix(cbind(rep(1,dim(X)[1]), X_num))
    trt_ind_t1 <- ifelse((A == t1), 1, NA)
    trt_ind_t2 <- ifelse((A == t2), 1, NA)
    if(marginal == TRUE){
      #
    } else {
      psi_group1 = var_para_mixed(reg_coef, a1, t1, mixed_info, H, weights_ind, A, G, cat_ind)[[1]]
      psid_mat_grp1 = var_para_mixed(reg_coef, a1, t1, mixed_info, H, weights_ind, A, G, cat_ind)[[2]]
      
      psi_group2 = var_para_mixed(reg_coef, a2, t2, mixed_info, H, weights_ind, A, G, cat_ind)[[1]]
      psid_mat_grp2 = var_para_mixed(reg_coef, a2, t2, mixed_info, H, weights_ind, A, G, cat_ind)[[2]]
      
      var_mat = sig_effect_mixed(psi_group1, psid_mat_grp1, psi_group2, psid_mat_grp2, len_m)
      ave <- var_effect_mixed(G, A, var_mat, cat_ind, mixed_info, X_mean = X_mean)
      
      X_fit1 = var_para_mixed(reg_coef, a1, t1, mixed_info, H, weights_ind, A, G, cat_ind)[[3]] 
      X_fit2 = var_para_mixed(reg_coef, a2, t2, mixed_info, H, weights_ind, A, G, cat_ind)[[3]]
      X_data = as.data.frame(cbind(G, X_fit1, X_fit2))
      X_group = aggregate(X_data[,-1], list(X_data$G), FUN=mean, na.rm = TRUE) # mean of residual for each group i
      pe <- sum(X_mean * reg_coef[,a1,t1]) - sum(X_mean * reg_coef[,a2,t2])
    }
  } else {
    # default is control-outcome
    if(marginal == TRUE){
      #
    } else {
      psi_group = var_para_mixed(reg_coef, a1, t1, mixed_info, H, weights_ind, A, G, cat_ind)[[1]]
      psid_mat_grp = var_para_mixed(reg_coef, a1, t1, mixed_info, H, weights_ind, A, G, cat_ind)[[2]]
      var_mat = sig_outcome_mixed(psi_group, psid_mat_grp)
      ave <- var_outcome_mixed(G, A, var_mat, cat_ind, mixed_info, X_mean = X_mean)
      
      X_fit = var_para_mixed(reg_coef, a1, t1, mixed_info, H, weights_ind, A, G, cat_ind)[[3]]
      X_data = as.data.frame(cbind(G, X_fit))
      X_group = aggregate(X_data$X_fit, list(X_data$G), FUN=mean, na.rm = TRUE) # mean of residual for each group i
      pe <- sum(X_mean * reg_coef[,a1,t1])
    }
  }
  
  ## Confidence Intervals ##
  qq <- qnorm(conf.level + (1 - conf.level)/2)
  me <- qq * sqrt(ave)
  
  ## Prepare Output ##
  pe <- pe * rescale.factor
  
  if(print == TRUE){
    toprint <- paste0('Estimate: ', round(pe, 2), ' ',
                      conf.level*100, '% CI: (', 
                      round(pe - me, 2), ', ', round(pe + me, 2), ')' )
    print(toprint)
  }
  
  out <- data.frame(estimate = pe,
                    std.error = sqrt(ave), 
                    conf.low = pe - me, conf.high = pe + me)
  return(out)
}

