point_estimates = ipw_point_estimates_mixed(H, G, A, w.matrix, X = X, x0 = x0, X_type = X_type)
marginal = FALSE
# no implementation of marginal or constrast yet
ipw_regression_variance <- function(weights, 
                            point_estimates, 
                            effect_type, 
                            marginal,
                            allocation1, 
                            allocation2 = NA,
                            trt.lvl1 = 0, 
                            trt.lvl2 = 1, 
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
  t1 <- as.character(trt.lvl1)
  t2 <- as.character(trt.lvl2)
  
  ## group coefficients conditional on x0_cat ##
  X_cat <- as.matrix(X[, X_type == "C"])
  x0_cat <- x0[X_type == "C"]
  X_num <- apply(as.matrix(X[, X_type == "N"]), 2, as.numeric)
  
  ## categorical indicator ##
  if (sum(X_type == "C") > 0){
    cat_ind <- apply(X_cat, 1, function(x) ifelse(prod(x == x0_cat), 1, NA)) 
  }else{
    cat_ind <- rep(1,dim(X_cat)[1])} #X = X0 categorical
  
  ## group coefficients for marginal / non-marginal (a = 0 / 1 separately) ##
  #grp_coef <- ipw_point_estimates_mixed(H, G, A, w.matrix, X = X, x0 = x0, X_type = X_type)$outcomes$grp_coefG

  fff <- ifelse(marginal == TRUE, 'marginal_outcomes', 'outcomes')
  #oal  <- point_estimates[[fff]]$overall_coefG
  #grp  <- point_estimates[[fff]]$grp_coefG
  ind_est_df <- point_estimates[[fff]]$weighted_ind
  reg_coef <- point_estimates[[fff]]$overall_coefG
    
  if(effect_type == 'contrast'){ 
    #X_reg <- cbind(G, 1, X_num)
    X_reg <- cbind(1, X_num)
    if(marginal == TRUE){
      # pe          <- oal[a1] - oal[a2]
      # pe_grp_diff <- (grp[ , a1] - oal[a1]) - (grp[, a2] - oal[a2])
      # if(variance_estimation == 'robust'){
      #   U_pe_grp    <- Ugrp[ , , a1] - Ugrp[ , , a2]
      # }
    } else {
      #grp_coef_a1 <- grp[, , a1, t1]
      #X_fit_a1 <- apply(X_reg, 1, function(x) sum(x[-1] * grp_coef_a1[x[1],]))
      oval_coef_a1 <- reg_coef[, a1, t1]
      X_fit_a1 <- apply(X_reg, 1, function(x) sum(x * oval_coef_a1))
      trt_ind_a1 <- (A == t1) 
      X_fit_trt_a1 <- X_fit_a1 * trt_ind_a1
      res_ij_a1 <- (ind_est_df[, , a1, t1] - X_fit_trt_a1) 
      res_data_a1 <- as.data.frame(cbind(res_ij_a1, G))
      res_group_a1 <- aggregate(res_data_a1$res_ij_a1, list(res_data_a1$G), FUN=mean, na.rm = TRUE)
      
      grp_coef_a2 <- grp[, , a2, t2]
      X_fit_a2 <- apply(X_reg, 1, function(x) sum(x[-1] * grp_coef_a2[x[1],]))
      trt_ind_a2 <- (A == t2) 
      X_fit_trt_a2 <- X_fit_a2 * trt_ind_a2
      res_ij_a2 <- (ind_est_df[, , a2, t2] - X_fit_trt_a2)
      res_data_a2 <- as.data.frame(cbind(res_ij_a2, G))
      res_group_a2 <- aggregate(res_data_a2$res_ij_a2, list(res_data_a2$G), FUN=mean, na.rm = TRUE)
      
      pe <- mean(res_group_a1$x) - mean(res_group_a2$x)
      pe_grp_diff <- (res_group_a1$x - mean(res_group_a1$x)) - (res_group_a2$x - mean(res_group_a2$x))
      
      # pe          <- oal[a1, t1] - oal[a2, t2] 
      # pe_grp_diff <- (grp[ , a1, t1] - oal[a1, t1]) - (grp[ , a2, t2] - oal[a2, t2])
      # if(variance_estimation == 'robust'){
      #   U_pe_grp    <- Ugrp[ , , a1, t1] - Ugrp[ , , a2, t2]
      # }
    }
  } else {
    # default is control-outcome
    X_reg <- cbind(G, 1, X_num)
    trt_ind <- ifelse((A == t1), 1, 0) 
    if(marginal == TRUE){
      # pe          <- oal[a1] 
      # pe_grp_diff <- (grp[ , a1] - oal[a1])
      # if(variance_estimation == 'robust'){
      #   U_pe_grp    <- Ugrp[ , , a1]
      # }
    } else {
      grp_coef <- grp[, , a1, t1]
      #oval_coef_a1 <- reg_coef[, a1, t1]
      X_fit <- apply(X_reg, 1, function(x) sum(x[-1] * grp_coef[x[1],]))
      X_fit_trt <- X_fit * trt_ind
      res_ij <- (ind_est_df[, , a1, t1] - X_fit_trt) 
      res_data <- as.data.frame(cbind(res_ij, G))
      res_group <- aggregate(res_data$res_ij, list(res_data$G), FUN=mean, na.rm = TRUE) # mean of residual for each group i
      pe <- mean(res_group$x)
      pe_grp_diff <- res_group$x - pe
      #var <-  mean((res_group$x - mean(res_group$x))^2) 
      # pe          <- ifelse(k == 1, oal[t1], oal[a1, t1]) 
      # pe_grp_diff <- ifelse(k == 1, (grp[ , t1] - oal[t1]),
      #                       (grp[ , a1, t1] - oal[a1, t1]))
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

# X = NULL,
# X_type = NULL, 
# x0 = NULL,
qnames <- paste0("x0 = ", toString(x0))

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
  
  if(effect_type == 'contrast'){ # TODO: change the name and check for a1 a2 to be the same allocation.
    if(marginal == TRUE){
      pe          <- oal[a1,] - oal[a2,]
      pe_grp_diff <- (grp[ , a1,] - oal[a1,]) - (grp[, a2,] - oal[a2,])
      # if(variance_estimation == 'robust'){
      #   U_pe_grp    <- Ugrp[ , , a1] - Ugrp[ , , a2]
      # }
    } else {
      pe          <- ifelse(k == 1, oal[a1, t1,] - oal[a2, t2,],
                            oal[a1, t1] - oal[a2, t2])
      pe_grp_diff <- ifelse(k == 1, (grp[ , a1, t1,] - oal[a1, t1,]) - (grp[ , a2, t2,] - oal[a2, t2,]),
                            (grp[ , a1, t1] - oal[a1, t1]) - (grp[ , a2, t2] - oal[a2, t2]))
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
