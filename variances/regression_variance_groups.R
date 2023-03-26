ipw_regression_variance_neigh_groups <- function(H,
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


var_para_neigh_groups <- function(reg_coef, a1, t1, neighinfo, H, weights_ind, A, G){
  neighX = neighinfo$neighX
  colnames(neighX) <- paste("neighX", colnames(neighX), sep = "_")
  len_h <- dim(neighinfo$neighX)[2]
  
  N <- length(unique(G))
  oval_coef_a1 <- reg_coef[, a1, t1][1:(1+len_h)]
  X_reg <- cbind(1, neighX)
  X_fit <- apply(X_reg, 1, function(x) sum(x * oval_coef_a1))
  res_ij <- (H - X_fit) # \mu - \betaX
  trt_ind <- ifelse((A == t1), 1, NA)
  weight_t1 <- weights_ind[, , a1, t1] * trt_ind
  psi_mat <- as.matrix(replicate(dim(X_reg)[2], (weight_t1 * res_ij))) * X_reg # extend to the dimension of X_reg
  psi_data <- as.data.frame(cbind(G, psi_mat))
  psi_group <- aggregate(psi_data[,-1], list(psi_data$G), FUN = mean, na.rm = TRUE)
  colnames(psi_group) <- c("G", "intcp", colnames(neighX))
  
  psid_mat <- as.matrix(replicate(dim(X_reg)[2], (weights_ind[, , a1, t1]))) * X_reg
  psid_data <- as.data.frame(cbind(G, A, psid_mat, X_reg))
  psid_mat_grp <- array(dim = c((1 + len_h), (1 + len_h), N)) 
  
  for (w in 1:N){
    psid_data_grp <- psid_data[G == w,]
    psid_mat_grp[, , w] <- t(as.matrix(psid_data_grp[,3:(3+len_h)])) %*% 
      as.matrix(psid_data_grp[,(4+len_h):(4+2*len_h)])/sum(psid_data_grp$A == t1)
  }
  
  outcomes = list(psi_group, psid_mat_grp, X_fit* trt_ind)
  outcomes
}