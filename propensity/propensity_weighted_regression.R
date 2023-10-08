# coefs <- out$point_estimates$outcomes$overall_coefG from second_stages_utils.R line 322
# X <- cond_X

# U21
ipw_point_estimates_propensity_regression <- function(H, G, A, weights, objs, X = NULL, x0 = NULL, 
                                           neighinfo = NULL, x1= NULL, X_type = NULL,
                                           Con_type = "No-con"){

  grps     <- dimnames(weights)[[1]]
  alphas   <- dimnames(weights)[[length(dim(weights))]]
  numerator_alphas <- as.numeric(lapply(alphas, function(l) substr(l[1],3,5)))
  trt_lvls <- sort(unique(A)) # binary if A = 1/0
  
  N <- length(grps)
  k <- length(alphas)
  l <- length(trt_lvls)
  q <- ifelse(!is.null(x1), dim(x1)[2],
              ifelse(!is.null(x0), dim(x0)[2], 1))
  
  # TODO: what about condition on both?
  if (!is.null(neighinfo)){
    qnames <- paste0("x1 = ", toString(x1))
  }else if (!is.null(X)){
    qnames <- paste0("x0 = ", toString(x0))
  }else{
    qnames <- "No-con"
  }
  
  # Use to handle the array of weight derivatives when provided ::
  p <- ifelse(is.matrix(weights), 1, dim(weights)[2]) 
  
  ## Generalize function to work on arrays. Add a dimension to matrices. ##
  if(is.matrix(weights)){
    weights <- array(c(weights, 1), dim=c(N, 1, k))
    predictors <- NULL
  } else {
    #predictors <- dimnames(weights)[[2]]
    predictors <- NULL
  }
  
  out <- list()
  
  # for neighinfo, the dimension is also sum(X_type == "N")
  len_g <- ifelse(is.null(X), 1, sum(X_type == "N"))
  len_h <- ifelse(is.null(neighinfo), 1, dim(neighinfo$neighX)[2])
  
  ## CALCULATE OUTCOME ESTIMATES PER TREATMENT LEVEL####
  
  hold_grp <- array(dim = c(N, p, k, l, q), dimnames = list(grps, predictors, 
                                                            alphas, trt_lvls, qnames))
  hold_oal <- array(dim = c(p, k, l, q),
                    dimnames = list(predictors, alphas, trt_lvls, qnames))
  
  weights_ind <- array(dim = c(length(A), p, k, l),
                       dimnames = list(NULL, NULL, alphas, trt_lvls))
   
  reg_est_grps <- array(dim = c(N, p, len_g+1, k, l, q),
                           dimnames = list(grps, NULL, NULL,
                                           alphas, trt_lvls, qnames)) # TODO: rename dimnames
  
  reg_est_overall <- array(dim = c(p, len_g+1, k, l, q),
                           dimnames = list(NULL, NULL,
                                           alphas, trt_lvls, qnames)) # TODO: rename dimnames
  
  for(ll in 1:l){    
    a <- trt_lvls[ll]
    
    # Modify weights per treatment level
    weights_trt <- array(dim= c(length(A), p, k))
    
    for(pp in 1:p){
      for (kk in 1:k) {
        weights_trt[ , pp, kk] <- apply(as.array(G), 1, function(x) weights[,pp,kk][x])
        weights_trt[ , pp, kk] <- weights_trt[ , pp, kk] * (A == a)/
          (numerator_alphas[kk]^a * (1-numerator_alphas[kk])^(1-a)) # indicator, \pi(i-j)
      }
    }
    
    
    # Compute estimates
    ind_est <- apply(weights_trt, 2:3, function(x) x * H) 
    
    grp_est <- array(dim= c(N, p, k, q))
    grp_est_overall <- array(dim= c(N, p, k, q))
    
    # 081423: For weighted regression, we replace H with X*(H - X \beta) or neighX *(H - neighX \beta)
    if (Con_type == "group"){
      # We have indicator A = a in weights_trt, so we can use beta under a for all units, the result is still 0
      if (sum(X_type == "C") > 0){
        X_cat <-  as.matrix(X[, X_type == "C"])
        cond_ind <- ifelse(X_cat == x0[X_type == "C"], 1, 0) 
      } else{
        cond_ind <- rep(1, length(A))
      }
      X_num <- apply(as.matrix(X[, X_type == "N"]), 2, as.numeric)
      X_num_intercept <- cbind(1, X_num)
      
      coefs <- objs$outcomes$overall_coefG
      coef_a <- as.matrix(coefs[, , as.character(a)])
      
      # reg_a is X'*(H - X \beta) instead of H
      reg_a <- (H - (X_num_intercept %*% coef_a))
      
      reg_est <- array(dim = c(length(A), p, len_g+1, k))
      for (i in 1:length(A)) {
        reg_est[i, , ,] <- reg_a[i] * weights_trt[i,,] %*% t(X_num_intercept[i,])
      }
      
      # calculate the group average
      reg_est_grp <- array(dim = c(length(grps), p, len_g+1, k))
      for (g in 1:length(grps)) {
        ind <- intersect(which(G == g), which(cond_ind == 1))
        #ind <- intersect(intersect(which(G == g), which(cond_ind == 1)), which(A == a))
        reg_est_grp[g, , ,] <- apply(reg_est[ind, , ,], 2:3, mean)
        reg_est_grps[g, , , , ll, ] <- apply(reg_est[ind, , ,], 2:3, mean)
      }
      
      reg_est_overall[, , k, ll,] <- apply(reg_est_grp, 2:3, mean) 
      #dim( out$Upart$outcomes$overall) 8 1 2 1 #  NULL "c(0.5, 0.4, 0.6)" [1] "0" "1" [1] "No-con"
    }  
    else if (Con_type == "neigh"){
      # We have indicator A = a in weights_trt, so we can use beta under a for all units, the result is still 0
      X_num <- neighinfo$neighX
      X_num_intercept <- cbind(1, X_num)
      
      coefs <- objs$outcomes$overall_coefH
      coef_a <- as.matrix(coefs[, , as.character(a)])
      
      # reg_a is X'*(H - X \beta) instead of H
      reg_a <- (H - (X_num_intercept %*% coef_a))
      
      reg_est <- array(dim = c(length(A), p, len_g+1, k))
      for (i in 1:length(A)) {
        reg_est[i, , ,] <- reg_a[i] * weights_trt[i,,] %*% t(X_num_intercept[i,])
        }
      
      # calculate the group average
      reg_est_grp <- array(dim = c(length(grps), p, len_g+1, k))
      for (g in 1:length(grps)) {
        ind <- which(G == g)
        #ind <- intersect(which(G == g), which(A == a))
        reg_est_grp[g, , ,] <- apply(reg_est[ind, , ,], 2:3, mean)
        reg_est_grps[g, , , , ll, ] <- apply(reg_est[ind, , ,], 2:3, mean)
        }
      
      reg_est_overall[, , k, ll,] <- apply(reg_est_grp, 2:3, mean) 
      #dim( out$Upart$outcomes$overall) 8 1 2 1 #  NULL "c(0.5, 0.4, 0.6)" [1] "0" "1" [1] "No-con"
      }
    }
  
  ## DONE ####
  return(list(U_grps = reg_est_grps, U_oal = reg_est_overall))
}

# U22      
# the derivative to \beta, replace X'*(H - X \beta) with X'* - X, still use weights instead of weightd
ipw_point_estimates_propensity_beta <- function(H, G, A, weights, objs, X = NULL, x0 = NULL, 
                                                      neighinfo = NULL, x1= NULL, X_type = NULL,
                                                      Con_type = "No-con"){
  
  grps     <- dimnames(weights)[[1]]
  alphas   <- dimnames(weights)[[length(dim(weights))]]
  numerator_alphas <- as.numeric(lapply(alphas, function(l) substr(l[1],3,5)))
  trt_lvls <- sort(unique(A)) # binary if A = 1/0
  
  N <- length(grps)
  k <- length(alphas)
  l <- length(trt_lvls)
  q <- ifelse(!is.null(x1), dim(x1)[2],
              ifelse(!is.null(x0), dim(x0)[2], 1))
  
  # TODO: what about condition on both?
  if (!is.null(neighinfo)){
    qnames <- paste0("x1 = ", toString(x1))
  }else if (!is.null(X)){
    qnames <- paste0("x0 = ", toString(x0))
  }else{
    qnames <- "No-con"
  }
  
  # Use to handle the array of weight derivatives when provided ::
  p <- ifelse(is.matrix(weights), 1, dim(weights)[2]) 
  
  ## Generalize function to work on arrays. Add a dimension to matrices. ##
  if(is.matrix(weights)){
    weights <- array(c(weights, 1), dim=c(N, 1, k))
    predictors <- NULL
  } else {
    #predictors <- dimnames(weights)[[2]]
    predictors <- NULL
  }
  
  out <- list()
  
  len_g <- ifelse(is.null(X), 1, sum(X_type == "N"))
  len_h <- ifelse(is.null(neighinfo), 1, dim(neighinfo$neighX)[2])
  
  ## CALCULATE OUTCOME ESTIMATES PER TREATMENT LEVEL####
  
  hold_grp <- array(dim = c(N, p, k, l, q), dimnames = list(grps, predictors, 
                                                            alphas, trt_lvls, qnames))
  hold_oal <- array(dim = c(p, k, l, q),
                    dimnames = list(predictors, alphas, trt_lvls, qnames))
  
  weights_ind <- array(dim = c(length(A), p, k, l),
                       dimnames = list(NULL, NULL, alphas, trt_lvls))
  
  reg_est_grps <- array(dim = c(N, len_g+1, len_g+1, k, l, q),
                        dimnames = list(grps, NULL, NULL,
                                        alphas, trt_lvls, qnames)) # TODO: rename dimnames
  
  reg_est_overall <- array(dim = c(len_g+1, len_g+1, k, l, q), 
                           dimnames = list(NULL, NULL, 
                                           alphas, trt_lvls, qnames)) # TODO: rename dimnames
  

  for(ll in 1:l){    
    a <- trt_lvls[ll]
    
    # Modify weights per treatment level
    weights_trt <- array(dim= c(length(A), p, k))
    
    for(pp in 1:p){
      for (kk in 1:k) {
        weights_trt[ , pp, kk] <- apply(as.array(G), 1, function(x) weights[,pp,kk][x])
        weights_trt[ , pp, kk] <- weights_trt[ , pp, kk] * (A == a)/
          (numerator_alphas[kk]^a * (1-numerator_alphas[kk])^(1-a)) # indicator, \pi(i-j)
        }
      }

    # Compute estimates
    ind_est <- apply(weights_trt, 2:3, function(x) x * H) 
    
    grp_est <- array(dim= c(N, p, k, q))
    grp_est_overall <- array(dim= c(N, p, k, q))
    
    #  X'*(H - X \beta) with X'* - X
    if (Con_type == "group"){
      # We have indicator A = a in weights_trt, so we can use beta under a for all units, the result is still 0
      if (sum(X_type == "C") > 0){
        X_cat <-  as.matrix(X[, X_type == "C"])
        cond_ind <- ifelse(X_cat == x0[X_type == "C"], 1, 0) 
      } else{
        cond_ind <- rep(1, length(A))
      }
      
      X_num <- apply(as.matrix(X[, X_type == "N"]), 2, as.numeric)
      X_num_intercept <- cbind(1, X_num)
      
      coefs <- objs$outcomes$overall_coefG
      coef_a <- as.matrix(coefs[, , as.character(a)])
      
      # reg_a is - X instead of H
      reg_a <- -1*X_num_intercept
      
      reg_est <- array(dim = c(length(A), len_g+1, len_g+1, k))
      for (i in 1:length(A)) {
        reg_est[i, , ,] <- weights_trt[i,,] * reg_a[i,]  %*% t(X_num_intercept[i,])
      }
      
      reg_est_grp <- array(dim = c(length(grps), len_g+1, len_g+1, k))
      #cond_ind <- ifelse(X_cat == x0[X_type == "C"], 1, 0) 
      
      # calculate the group average, select only the unit that categorical variable = x0
      for (g in 1:length(grps)) {
        ind <- intersect(which(G == g), which(cond_ind == 1))
        #ind <- intersect(intersect(which(G == g), which(cond_ind == 1)), which(A == a))
        reg_est_grp[g, , ,] <- apply(reg_est[ind, , ,], 2:3, mean)
        
        reg_est_grps[g, , , , ll, ] <- apply(reg_est[ind, , ,], 2:3, mean)
      }
      
      reg_est_overall[, , k, ll,] <- apply(reg_est_grp, 2:3, mean) 
    }
    else if (Con_type == "neigh"){
      # We have indicator A = a in weights_trt, so we can use beta under a for all units, the result is still 0
      X_num <- neighinfo$neighX
      X_num_intercept <- cbind(1, X_num)
      
      coefs <- objs$outcomes$overall_coefH
      coef_a <- as.matrix(coefs[, , as.character(a)])
      
      # reg_a is - X instead of H
      reg_a <- -1*X_num_intercept
      
      reg_est <- array(dim = c(length(A), len_g+1, len_g+1, k))
      for (i in 1:length(A)) {
        reg_est[i, , ,] <- weights_trt[i,,] * reg_a[i,]  %*% t(X_num_intercept[i,])
      }
      
      reg_est_grp <- array(dim = c(length(grps), len_g+1, len_g+1, k))

      for (g in 1:length(grps)) {
        ind <- which(G == g)
        #ind <- intersect(which(G == g), which(A == a))
        reg_est_grp[g, , ,] <- apply(reg_est[ind, , ,], 2:3, mean)
        
        reg_est_grps[g, , , , ll, ] <- apply(reg_est[ind, , ,], 2:3, mean)
      }
      
      reg_est_overall[, , k, ll,] <- apply(reg_est_grp, 2:3, mean) 
      }
  }
  ## DONE ####
  return(reg_est_grps)
}

# V matrix, use weights here
ipw_point_estimates_propensity_Vmatrix <- function(H, G, A, weights, objs, X = NULL, x0 = NULL, 
                                                neighinfo = NULL, x1= NULL, X_type = NULL,
                                                Con_type = "No-con"){
  
  grps     <- dimnames(weights)[[1]]
  alphas   <- dimnames(weights)[[length(dim(weights))]]
  numerator_alphas <- as.numeric(lapply(alphas, function(l) substr(l[1],3,5)))
  trt_lvls <- sort(unique(A)) # binary if A = 1/0
  
  N <- length(grps)
  k <- length(alphas)
  l <- length(trt_lvls)
  q <- ifelse(!is.null(x1), dim(x1)[2],
              ifelse(!is.null(x0), dim(x0)[2], 1))
  
  if (!is.null(neighinfo)){
    qnames <- paste0("x1 = ", toString(x1))
  }else if (!is.null(X)){
    qnames <- paste0("x0 = ", toString(x0))
  }else{
    qnames <- "No-con"
  }
  
  # Use to handle the array of weight derivatives when provided ::
  p <- ifelse(is.matrix(weights), 1, dim(weights)[2]) 
  
  ## Generalize function to work on arrays. Add a dimension to matrices. ##
  if(is.matrix(weights)){
    weights <- array(c(weights, 1), dim=c(N, 1, k))
    predictors <- NULL
  } else {
    predictors <- NULL
  }
  
  out <- list()
  
  len_g <- ifelse(is.null(X), 1, sum(X_type == "N"))
  len_h <- ifelse(is.null(neighinfo), 1, dim(neighinfo$neighX)[2])
  
  ## CALCULATE OUTCOME ESTIMATES PER TREATMENT LEVEL####
  
  hold_grp <- array(dim = c(N, p, k, l, q), dimnames = list(grps, predictors, 
                                                            alphas, trt_lvls, qnames))
  hold_oal <- array(dim = c(p, k, l, q),
                    dimnames = list(predictors, alphas, trt_lvls, qnames))
  
  weights_ind <- array(dim = c(length(A), p, k, l),
                       dimnames = list(NULL, NULL, alphas, trt_lvls))
  
  reg_est_grps <- array(dim = c(length(grps), len_g+1, k, l, q),
                        dimnames = list(NULL, NULL, 
                                        alphas, trt_lvls, qnames)) # TODO: rename dimnames
  
  for(ll in 1:l){    
    a <- trt_lvls[ll]
    
    # Modify weights per treatment level
    weights_trt <- array(dim= c(length(A), p, k))
    
    for(pp in 1:p){
      for (kk in 1:k) {
        weights_trt[ , pp, kk] <- apply(as.array(G), 1, function(x) weights[,pp,kk][x])
        weights_trt[ , pp, kk] <- weights_trt[ , pp, kk] * (A == a)/
          (numerator_alphas[kk]^a * (1-numerator_alphas[kk])^(1-a)) # indicator, \pi(i-j)
      }
    }
    
    
    # Compute estimates
    ind_est <- apply(weights_trt, 2:3, function(x) x * H) 
    
    grp_est <- array(dim= c(N, p, k, q))
    grp_est_overall <- array(dim= c(N, p, k, q))
    
    if (Con_type == "group"){
      # We have indicator A = a in weights_trt, so we can use beta under a for all units, the result is still 0
      if (sum(X_type == "C") > 0){
        X_cat <-  as.matrix(X[, X_type == "C"])
        cond_ind <- ifelse(X_cat == x0[X_type == "C"], 1, 0) 
      } else{
        cond_ind <- rep(1, length(A))
      }
      
      X_num <- apply(as.matrix(X[, X_type == "N"]), 2, as.numeric)
      X_num_intercept <- cbind(1, X_num)
      
      coefs <- objs$outcomes$overall_coefG
      coef_a <- as.matrix(coefs[, , as.character(a)])
      
      # reg_a is X' (H - X \beta)
      reg_a <-  (H - (X_num_intercept %*% coef_a))
      
      reg_est <- array(dim = c(length(A), len_g+1, k))
      for (i in 1:length(A)) {
        reg_est[i, ,] <- weights_trt[i,,] * reg_a[i,]  %*% t(X_num_intercept[i,])
      }
      #cond_ind <- ifelse(X_cat == x0[X_type == "C"], 1, 0) 
      
      #reg_est_grp <- array(dim = c(length(grps), len_g+1, k))
      # calculate the group average
      for (g in 1:length(grps)) {
        ind <- intersect(which(G == g), which(cond_ind == 1))
        #ind <- intersect(intersect(which(G == g), which(cond_ind == 1)), which(A == a))
        #reg_est_grps[g, , , as.character(a), ] <- apply(reg_est[ind, ,], 2, mean)
        reg_est_grps[g, , , ll, ] <- apply(reg_est[ind, ,], 2, mean)
        }    
    }
    else if (Con_type == "neigh"){
      # We have indicator A = a in weights_trt, so we can use beta under a for all units, the result is still 0
      X_num <- neighinfo$neighX
      X_num_intercept <- cbind(1, X_num)
      
      coefs <- objs$outcomes$overall_coefH
      coef_a <- as.matrix(coefs[, , as.character(a)])
      
      # reg_a is X' (H - X \beta)
      reg_a <-  (H - (X_num_intercept %*% coef_a))
      
      reg_est <- array(dim = c(length(A), len_g+1, k))
      for (i in 1:length(A)) {
        reg_est[i, ,] <- weights_trt[i,,] * reg_a[i,]  %*% t(X_num_intercept[i,])
      }
      
      # calculate the group average
      for (g in 1:length(grps)) {
        ind <- which(G == g)
        #ind <- intersect(which(G == g), which(A == a))
        reg_est_grps[g, , , as.character(a), ] <- apply(reg_est[ind, ,], 2, mean)
      }    
    }
  }
  
  ## DONE ####
  return(reg_est_grps)
}

#beta_est_eqs <- reg_est_grps

V_matrix_regression <- function(scores, 
                                beta_est_eqs, 
                                allocation1, 
                                trt.lvl1, 
                                allocation2 = NA, 
                                trt.lvl2    = NA, 
                                effect_type, 
                                marginal){
  ## Necessary bits ##
  N  <- dim(scores)[1]
  p  <- dim(scores)[2]
  a1 <- allocation1
  a2 <- allocation2
  t1 <- trt.lvl1
  t2 <- trt.lvl2
  
  ## Grab the last element of the psi(O, theta) vector: psi_a, alpha ##
  fff <- ifelse(marginal == TRUE, 'marginal_outcomes', 'outcomes')
  
  if(effect_type == 'contrast'){   
    if(marginal == TRUE){
      #xx <- (hold_grp[ , a1,] - hold_oal[, a1]) - (hold_grp[, a2,] - hold_oal[, a2])
    } else {
      # xx <- (hold_grp[, a1, t1,] - hold_oal[a1, t1,]) -
      #   (hold_grp[ , a2, t2,] - hold_oal[a2, t2,])
      xx <- cbind(beta_est_eqs[, , a1, t1,], beta_est_eqs[, , a2, t2,])
    }
  } 
  else if(effect_type == 'outcome'){
    if(marginal == TRUE){
      #xx <- hold_grp[ , a1,] - hold_oal[a1, ]
    } else {
      #xx <- hold_grp[ , a1, t1,] - hold_oal[a1, t1,]
      xx <- beta_est_eqs[, , a1, t1,]
    }
  }
  
  ee <- cbind(scores, xx)
  V <- crossprod(ee)/N
  V
}



# causal_estimation_options = list(variance_estimation = 'robust')
# # dots <- list(X = cond_X,
# #              x0 = x0,
# #              X_type = c("C", "N"),
# #              Con_type = "group")
# 
# dots <- list(neighinfo = neighinfo,
#              x1 = as.matrix(x1_num),
#              X_type = c("N"),
#              Con_type = "neigh")

ipw_interference_regression <- function(propensity_integrand,
                                    loglihood_integrand = propensity_integrand,
                                    allocations,
                                    H, propensity_X, A, B = A, G, P,
                                    parameters,
                                    first_assignments, 
                                    variance_estimation,
                                    runSilent   = TRUE, 
                                    integrate_allocation,
                                    ...)
{
  dots <- list(...)
  
  ## Warnings ##
  
  #### Arguments Necessary for Causal Estimation Functions ####
  integrand_args <- get_args(FUN = propensity_integrand, args_list = dots)
  point_est_args <- get_args(FUN = ipw_point_estimates_mixed_test4, args_list = dots)
  loglihood_args <- get_args(FUN = loglihood_integrand, args_list = dots)
  grad_args      <- get_args(FUN = numDeriv::grad, args_list = dots)
  integrate_args <- get_args(FUN = stats::integrate, args_list = dots)
  
  weight_args <- append(append(integrand_args, integrate_args),
                        list(integrand   = propensity_integrand, 
                             allocations = allocations, 
                             propensity_X = propensity_X, A = A, G = G, P = P,
                             parameters = parameters,
                             runSilent  = runSilent, #BB 2015-06-23
                             integrate_allocation = integrate_allocation,
                             first_assignments = first_assignments #100123
                        ))
  #score_args <- append(weight_args, list(first_assignments = first_assignments))
  
  #### Prepare output ####
  out <- list()  
  
  ## Compute Weights ##
  # 100123: change wght_matrix_second to wght_matrix_second_prop
  weights <- do.call(wght_matrix_second_prop, args = weight_args)
  
  if(variance_estimation == 'robust'){
    weightd <- do.call(wght_deriv_array_second_prop, args = append(weight_args, grad_args))
    out$weightd <- weightd
  }
  #   U11 <- do.call(score_matrix_deriv, args = append(score_args, grad_args))
  #   out$U11 <- U11
  #   }
  
  
  #### COMPUTE ESTIMATES AND OUTPUT ####
  estimate_args <- append(point_est_args, list(H = H, G = G, A = A))#, list(Y = Y, G = G, A = A)
  point_args    <- append(estimate_args, list(weights = weights))
  
  
  #### Calculate output ####
  #out$point_estimates <- do.call(ipw_point_estimates_mixed_test4, args = point_args)
  out$point_estimates <- do.call(ipw_point_estimates_propensity, args = point_args)
  
  if(variance_estimation == 'robust'){
    U_args     <- append(estimate_args, list(weights = weightd, objs =  out$point_estimates))
    # 081723: for U22, we need weights instead of weightd
    beta_args     <- append(point_args, list(objs =  out$point_estimates))
    
    sargs      <- append(append(loglihood_args, grad_args), integrate_args)
    score_args <- append(sargs, list(integrand = loglihood_integrand,
                                     propensity_X = propensity_X, G = G, P = P,
                                     A = B, # Use B for treatment in scores
                                     parameters = parameters,
                                     first_assignments = first_assignments,
                                     runSilent  = runSilent #BB 2015-06-23
    ))
    
    U11 <- do.call(score_matrix_deriv, args = append(score_args, grad_args))
    out$U11 <- U11
    
    # set randomization scheme to 1 for scores for logit_integrand
    score_args$randomization <- 1
    
    # 081823: add U22 and the last two columns for V matrix
    out$Upart <- do.call(ipw_point_estimates_propensity_regression, args = U_args)
    #out$Upart           <- do.call(ipw_point_estimates_propensity, args = U_args)
    out$U22_grps <- do.call(ipw_point_estimates_propensity_beta, args = beta_args)
    
    out$V_grps <- do.call(ipw_point_estimates_propensity_Vmatrix, args = beta_args)
    
    out$scores          <- do.call(score_matrix_second, args = score_args)
  } 
  
  out$weights <- weights

  return(out)
}

#allocations = list(c(numerator_alpha, denominator_alphas))

ipw_effect_calc_regression <- function(obj, 
                                   weights, #added parameter, w.matrix
                                   variance_estimation,
                                   P,
                                   propensity_X, 
                                   alpha1, 
                                   trt.lvl1, 
                                   alpha2 = NA, 
                                   trt.lvl2 = NA,
                                   effect_type,
                                   marginal,
                                   rescale.factor = 1,
                                   conf.level = 0.95,
                                   print = FALSE,
                                   X, x0,
                                   neighinfo, x1,
                                   X_type, Con_type)
{
  allocations <- dimnames(obj$weights)[[2]] 
  
  ## Warnings ##
  # Print error if either estimates with alpha1 have been computed 
  # or a constrast is being estimated when estimates for alpha2
  # have not been computed
  if (!(alpha1 %in% allocations) | (effect_type == 'contrast' & !(alpha2 %in% allocations))){
    stop(paste('At least one of the chosen coverage levels has not been estimated.\n',
               'Select from the following: \n', 
               paste(allocations, collapse = ' ')))
  }
  
  ## Necessary bits ##
  N  <- dim(obj$weights)[1] 
  p  <- dim(obj$scores)[2] 
  k  <- length(allocations)
  l  <- dim(obj$point_estimates$outcomes$overall)[2]
  a1 <- as.character(alpha1)
  a2 <- as.character(alpha2)
  t1 <- as.character(trt.lvl1)
  t2 <- as.character(trt.lvl2)
  
  fff <- ifelse(marginal == TRUE, 'marginal_outcomes', 'outcomes')
  
  oal  <- obj$point_estimates[[fff]]$overall
  grp  <- obj$point_estimates[[fff]]$groups
  
  if(variance_estimation == 'robust'){
    #Uoal <- obj$Upart[[fff]]$overall 
    #Ugrp <- obj$Upart[[fff]]$groups
    Uoal <- obj$Upart$U_oal 
    Ugrp <- obj$Upart$U_grps
    U22_grps <- obj$U22_grps
    V_grps <- obj$V_grps
      
    # TODO: change this part
    # Cludgy workaround for case of 1 fixed effect: add dimension to Ugrp array #
    if(p == 1){
      names <- dimnames(Ugrp)
      if(marginal == TRUE){
        Ugrp <-  array(c(Ugrp[1:N, ], 1, Ugrp[, 1:k]),
                       dim=c(N, 1, k),
                       dimnames = list(names[[1]], 'Intercept', names[[2]]))
      } else {
        Ugrp <-  array(c(Ugrp[1:N, , ], 1, Ugrp[, 1:k , 1:l]),
                       dim=c(N, 1, k, l),
                       dimnames = list(names[[1]], 'Intercept', names[[2]], names[[3]]))
      }
    }
  }
  
  if(effect_type == 'contrast'){
    if(marginal == TRUE){
      pe          <- oal[a1] - oal[a2]
      pe_grp_diff <- (grp[ , a1] - oal[a1]) - (grp[, a2] - oal[a2])
      if(variance_estimation == 'robust'){
        U_pe_grp    <- Ugrp[ , , a1] - Ugrp[ , , a2]
      }
    } else {
      # pe          <- oal[a1, t1, ] - oal[a2, t2, ]
      # pe_grp_diff <- (grp[ , a1, t1, ] - oal[a1, t1, ]) - (grp[ , a2, t2, ] - oal[a2, t2, ])
      pe          <- oal[,a1, t1, ] - oal[,a2, t2, ]
      pe_grp_diff <- (grp[, , a1, t1, ] - oal[,a1, t1, ]) - (grp[, , a2, t2, ] - oal[,a2, t2, ])
      
      # if(variance_estimation == 'robust'){
      #   #U_pe_grp    <- as.matrix(Ugrp[, , a1, t1, ] - Ugrp[, , a2, t2, ])
      #   U_pe_grp <- Ugrp[ , , , a1, t1,] - Ugrp[ , , , a2, t2,]
      #   U_beta_grp <- U22_grps[ , , , a1, t1,] - U22_grps[ , , , a2, t2,]
      # }
    }
  } else {
    if(marginal == TRUE){
      pe          <- oal[a1] 
      pe_grp_diff <- (grp[ , a1] - oal[a1])
      if(variance_estimation == 'robust'){
        U_pe_grp    <- Ugrp[ , , a1]
      }
    } else {
      pe          <- oal[, a1, t1, ] 
      pe_grp_diff <- (grp[ , , a1, t1, ] - oal[, a1, t1, ])
      if(variance_estimation == 'robust'){
        U_pe_grp <- Ugrp[ , , , a1, t1,]
        U_beta_grp <- U22_grps[ , , , a1, t1,]
      }
    }
  }
  
  #### VARIANCE ESTIMATION ####
  if(variance_estimation == 'robust'){
    # partial U matrix
    if (effect_type == 'outcome'){
    if(p == 1){
      U21 <- sum(-U_pe_grp)/N
      U22 <- sum(-U_beta_grp)/N
    } else {
      U21 <- (t(as.matrix(apply(-U_pe_grp, 2:3, sum, na.rm = T))))/N # times -1
      U22 <- (t(as.matrix(apply(-U_beta_grp, 2:3, sum, na.rm = T))))/N# times -1
    }
    
    # U matrix
    U11 <- obj$U11 #( it is already -1 * deriv)
    U <- cbind(rbind(U11, U21), rbind(matrix(0, nrow = dim(U11)[1], ncol = dim(U22)[1]), U22))
    
    # V matrix
    V <- V_matrix_regression(scores = obj$scores, 
                             beta_est_eqs = obj$V_grps, 
                         allocation1 = a1, allocation2 = a2, 
                         trt.lvl1 = t1, trt.lvl2 = t2, 
                         effect_type = effect_type, marginal = marginal)
    
      
    ## Sandwich Variance Estimate ##
    inv_U = solve(U) 
    sigma = inv_U %*% V %*% t(inv_U)
    vdim <- dim(V)[1]
    
    # variance for the bottom right part (estimate)
    #beta_ave <- sigma[(p+1):vdim, (p+1):vdim]/N
    V11 <- V[1:dim(U11)[1], 1:dim(U11)[1]]
    V21 <- V[(dim(U11)[1]+1):vdim, 1:dim(U11)[1]]
    V22 <- V[(dim(U11)[1]+1):vdim, (dim(U11)[1]+1):vdim]
    beta_ave <- ((solve(U22)%*%U21 - 2*solve(U22)%*%V21) %*% t(solve(U22)%*% U21 %*% solve(U11)) + solve(U22) %*% V22 %*% t(solve(U22)))/N
    }
    else if (effect_type == 'contrast'){
      U_pe_grp1 <- Ugrp[ , , , a1, t1,]
      U_pe_grp2 <- Ugrp[ , , , a2, t2,]
      U_beta_grp1 <- U22_grps[ , , , a1, t1,] 
      U_beta_grp2 <- U22_grps[ , , , a2, t2,]
      
      if(p == 1){
        U21_1 <- sum(-U_pe_grp1)/N
        U21_2 <- sum(-U_pe_grp2)/N
        U22_1 <- sum(-U_beta_grp1)/N
        U22_2 <- sum(-U_beta_grp2)/N
      } else {
        U21_1 <- (t(as.matrix(apply(-U_pe_grp1, 2:3, sum, na.rm = T))))/N # times -1
        U22_1 <- (t(as.matrix(apply(-U_beta_grp1, 2:3, sum, na.rm = T))))/N# times -1
        U21_2 <- (t(as.matrix(apply(-U_pe_grp2, 2:3, sum, na.rm = T))))/N # times -1
        U22_2 <- (t(as.matrix(apply(-U_beta_grp2, 2:3, sum, na.rm = T))))/N# times -1
      }
      U11 <- obj$U11
      U <- cbind(rbind(U11, U21_1, U21_2), 
                 rbind(matrix(0, nrow = dim(U11)[1], ncol = dim(U22_1)[1]),
                       U22_1, matrix(0, nrow = dim(U22_1)[1], ncol = dim(U22_1)[1])),
                 rbind(matrix(0, nrow = dim(U11)[1], ncol = dim(U22_2)[1]),
                       matrix(0, nrow = dim(U22_2)[1], ncol = dim(U22_2)[1]), U22_2))
      
      # reselect U21 and U22 from the big U now
      U21 <- rbind(U21_1, U21_2)
      U22 <- U[(dim(U11)[1]+1):dim(U)[1],(dim(U11)[1]+1):dim(U)[1]]
      
      V <- V_matrix_regression(scores = obj$scores, 
                               beta_est_eqs = obj$V_grps, 
                               allocation1 = a1, allocation2 = a2, 
                               trt.lvl1 = t1, trt.lvl2 = t2, 
                               effect_type = effect_type, marginal = marginal)
      
      ## Sandwich Variance Estimate ##
      inv_U = solve(U) 
      sigma = inv_U %*% V %*% t(inv_U)
      vdim <- dim(V)[1]
      
      V11 <- V[1:dim(U11)[1], 1:dim(U11)[1]]
      V21 <- V[(dim(U11)[1]+1):vdim, 1:dim(U11)[1]]
      V22 <- V[(dim(U11)[1]+1):vdim, (dim(U11)[1]+1):vdim]
      beta_ave <- ((solve(U22)%*%U21 - 2*solve(U22)%*%V21) %*% t(solve(U22)%*% U21 %*% solve(U11)) + solve(U22) %*% V22 %*% t(solve(U22)))/N
    }
  } else if(variance_estimation == 'naive'){
    ave <- (1/(N^2)) * (sum((pe_grp_diff)^2, na.rm = T)) * rescale.factor^2
  }
  
  if(Con_type == "group"){
    X_num <- X[, X_type == "N"]
    X0_num <- x0[X_type == "N"]
    X_num_intercept <- as.matrix(apply(as.matrix(c(1, X0_num)), 1, as.numeric))
    if (effect_type == 'outcome'){
    ave <- t(X_num_intercept) %*% beta_ave %*% X_num_intercept
    } else if (effect_type == 'contrast'){
      ave <- t(c(X_num_intercept, -1*X_num_intercept)) %*% beta_ave %*% c(X_num_intercept, -1*X_num_intercept)
    }
  }
  else if (Con_type == "neigh"){
    X_num <- neighinfo$neighX
    X1_num <- x1[X_type == "N"]
    X_num_intercept <- as.matrix(apply(as.matrix(c(1, X1_num)), 1, as.numeric))
    if (effect_type == 'outcome'){
      ave <- t(X_num_intercept) %*% beta_ave %*% X_num_intercept
    } else if (effect_type == 'contrast'){
      ave <- t(c(X_num_intercept, -1*X_num_intercept)) %*% beta_ave %*% c(X_num_intercept, -1*X_num_intercept)
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
                    conf.low = pe - me, 
                    conf.high = pe + me)
  return(out)
}

ipw_propensity_variance_regression <- function(parameters,
                                               allocations,
                                               H, propensity_X, A, G, P, 
                                               first_assignments, 
                                               causal_estimation_options = 
                                                 list(variance_estimation = 'robust'),
                                               integrate_allocation = FALSE,
                                               effect_type = effect_type, #or contrast
                                               propensity_integrand = logit_integrand_second,
                                               ...){
  
  out <- list()
  dots <- list(...)
  ipw_args <- append(append(dots, causal_estimation_options),
                     list(propensity_integrand = logit_integrand_second, 
                          loglihood_integrand  = propensity_integrand,
                          allocations          = allocations,
                          parameters           = parameters,
                          runSilent            = TRUE, 
                          integrate_allocation = integrate_allocation,
                          H = H, propensity_X = propensity_X, 
                          A = A, B = A, G = G, P = P, first_assignments = first_assignments))
  ipw <- do.call(ipw_interference_regression, args = ipw_args)
  out <- append(out, ipw)
  estimate_args <- list(obj = ipw,
                        variance_estimation = causal_estimation_options$variance_estimation,
                        causal_estimation_options$variance_estimation,
                        propensity_X = propensity_X, 
                        P           = P, 
                        alpha1      = allocations[1],
                        trt.lvl1    = 1,
                        alpha2      = allocations[1],
                        trt.lvl2    = 0,
                        marginal    = FALSE,
                        effect_type = effect_type,
                        rescale.factor = 1,
                        conf.level = 0.95,
                        print = FALSE)
  regression_args <- append(estimate_args, dots)
  est <- do.call(ipw_effect_calc_regression, args = regression_args)
  return(est)
}


# this function is to run regression for different cluster size
ipw_point_estimates_propensity2 <- function(H, G, A, weights, X = NULL, x0 = NULL, 
                                           neighinfo = NULL, x1= NULL, X_type = NULL,
                                           Con_type = "No-con"){
  
  
  grps     <- dimnames(weights)[[1]]
  alphas   <- dimnames(weights)[[length(dim(weights))]]
  numerator_alphas <- as.numeric(lapply(alphas, function(l) substr(l[1],3,5)))
  trt_lvls <- sort(unique(A)) # binary if A = 1/0
  
  N <- length(grps)
  k <- length(alphas)
  l <- length(trt_lvls)
  q <- ifelse(!is.null(x1), dim(x1)[2],
              ifelse(!is.null(x0), dim(x0)[2], 1))
  
  # TODO: what about condition on both?
  if (!is.null(neighinfo)){
    qnames <- paste0("x1 = ", toString(x1))
  }else if (!is.null(X)){
    qnames <- paste0("x0 = ", toString(x0))
  }else{
    qnames <- "No-con"
  }
  
  # Use to handle the array of weight derivatives when provided ::
  p <- ifelse(is.matrix(weights), 1, dim(weights)[2]) 
  
  ## Generalize function to work on arrays. Add a dimension to matrices. ##
  if(is.matrix(weights)){
    weights <- array(c(weights, 1), dim=c(N, 1, k))
    predictors <- NULL
  } else {
    #predictors <- dimnames(weights)[[2]]
    predictors <- NULL
  }
  
  out <- list()
  
  
  len_g <- ifelse(is.null(X), 1, sum(X_type == "N"))
  len_h <- ifelse(is.null(neighinfo), 1, dim(neighinfo$neighX)[2])
  
  ## CALCULATE OUTCOME ESTIMATES PER TREATMENT LEVEL####
  
  hold_grp <- array(dim = c(N, p, k, l, q), dimnames = list(grps, predictors, 
                                                            alphas, trt_lvls, qnames))
  hold_oal <- array(dim = c(p, k, l, q),
                    dimnames = list(predictors, alphas, trt_lvls, qnames))
  
  hold_grp_coefM <- array(dim = c(N, (1+len_g+len_h+len_g*len_h)*p, k, l), dimnames = list(grps, predictors,
                                                                                           alphas, trt_lvls))
  hold_grp_coefH <- array(dim = c(N, (1+len_h)*p, k, l), dimnames = list(grps, predictors,
                                                                         alphas, trt_lvls))
  hold_grp_coefG <- array(dim = c(N, (1+len_g)*p, k, l), dimnames = list(grps, predictors,
                                                                         alphas, trt_lvls))
  hold_oal_coefM <- array(dim = c((1+len_g+len_h+len_g*len_h)*p, k, l),
                          dimnames = list(predictors, alphas, trt_lvls))
  hold_oal_coefH <- array(dim = c((1+len_h)*p, k, l),
                          dimnames = list(predictors, alphas, trt_lvls))
  hold_oal_coefG <- array(dim = c((1+len_g)*p, k, l),
                          dimnames = list(predictors, alphas, trt_lvls))
  
  weights_ind <- array(dim = c(length(A), p, k, l),
                       dimnames = list(NULL, NULL, alphas, trt_lvls))
  
  
  for(ll in 1:l){    
    a <- trt_lvls[ll]
    
    # Compute means per group
    #Hbar_trt <- group_means(H = H, A = A, G = G, a = a)
    
    # Modify weights per treatment level
    weights_trt <- array(dim= c(length(A), p, k))
    
    for(pp in 1:p){
      for (kk in 1:k) {
        weights_trt[ , pp, kk] <- apply(as.array(G), 1, function(x) weights[,pp,kk][x])
        weights_trt[ , pp, kk] <- weights_trt[ , pp, kk] * (A == a)/
          (numerator_alphas[kk]^a * (1-numerator_alphas[kk])^(1-a)) # add some comments
      }
    }
    
    #weights_ind[, pp, kk, ll] <- weights_trt[ , pp, kk]
    
    # Compute estimates
    ind_est <- apply(weights_trt, 2:3, function(x) x * H) 
    coef_est_M <- array(dim= c(N, (1+len_g+len_h+len_g*len_h)*p, k))
    coef_est_H <- array(dim= c(N, (1+len_h)*p, k))
    coef_est_G <- array(dim= c(N, (1+len_g)*p, k))
    
    ova_coef_est_M <- array(dim= c(1, (1+len_g+len_h+len_g*len_h)*p, k))
    ova_coef_est_H <- array(dim= c(1, (1+len_h)*p, k))
    ova_coef_est_G <- array(dim= c(1, (1+len_g)*p, k))
    grp_est <- array(dim= c(N, p, k, q))
    grp_est_overall <- array(dim= c(N, p, k, q))
    
    ##############
    # Calculate the coefficients for conditional H and group means
    for (j in 1:k){
      ind_est_df <- ind_est[ , , j]
      weights_df <- weights_trt[ , , j]
      if (Con_type == "mixed"){
        M_list <- mixed_coef(weights_df, H, G, X, X_type, x0, neighinfo, A, a)
        coef_est_M[ , , j] <- as.matrix(M_list[[1]])
        ova_coef_est_M[ , , j] <- as.matrix(M_list[[2]])
      }else if (Con_type == "neigh"){
        H_list <- neigh_coefs_oncont5(weights_df, H, G, neighinfo, A, a)
        coef_est_H[ , , j] <- as.matrix(H_list[[1]])
        ova_coef_est_H[ , , j] <- as.matrix(H_list[[2]])
      }else if (Con_type == "group"){
        G_list <- group_coefs_oncont1(weights_df, H, G, X, X_type, x0, A, a)
        coef_est_G[ , , j] <- as.matrix(G_list[[1]])
        ova_coef_est_G[ , , j] <- as.matrix(G_list[[2]])
      }
    }
    
    for (j in 1:k) {
      ind_est_df <- as.matrix(ind_est[ , , j])
      weights_est_df <- weights_trt[ , , j]
      
      if (Con_type == "mixed"){
        grp_est_overall[ , p, j, ] <- apply(as.matrix(x1), 2, mixed_means_overall,
                                            overall_coef = ova_coef_est_M[ , , j],
                                            X_type = X_type, x0 = x0, N = N)
        grp_est[ , p, j, ] <- apply(as.matrix(x1), 2, mixed_means_group,
                                    cond_coefs = coef_est_M[ , , j],
                                    X_type = X_type, x0 = x0)
      }else if (Con_type == "neigh"){
        grp_est_overall[ , p, j, ] <- apply(as.matrix(x1), 2, neigh_means_oncont4,
                                            overall_coef = ova_coef_est_H[ , , j],
                                            X_type = X_type, N = N)
        grp_est[ , p, j, ] <- apply(as.matrix(x1), 2, neigh_means_oncont3,
                                    cond_coefs = coef_est_H[ , , j],
                                    X_type = X_type)
      }else if (Con_type == "group"){
        grp_est_overall[ , p, j, ] <- apply(as.matrix(x0), 2, group_means_oncont1,
                                            overall_coef = ova_coef_est_G[ , , j],
                                            X_type = X_type, N = N)
        grp_est[ , p, j, ] <- apply(as.matrix(x0), 2, group_means_oncont2,
                                    cond_coefs = coef_est_G[ , , j],
                                    X_type = X_type)
      }else{
        grp_est_overall[, , j, ] <- apply(ind_est_df, 2, group_means_propensity, A, G, a)
        grp_est[, , j, ] <- apply(ind_est_df, 2, group_means_propensity, A, G, a)
        #grp_est_overall[, p, j, ] <- group_means_null(ind_est_df, G, A, a)
        #grp_est[, p, j, ] <- group_means_null(ind_est_df, G, A, a)
      }
    }
    
    oal_est <- apply(grp_est_overall, 2:4, mean, na.rm = TRUE)
    
    hold_grp[ , , , ll, ] <- grp_est
    hold_oal[ , , ll, ]   <- oal_est
    
    hold_grp_coefM[ , , , ll] <- coef_est_M
    hold_grp_coefH[ , , , ll] <- coef_est_H
    hold_grp_coefG[ , , , ll] <- coef_est_G
    
    hold_oal_coefM[ , , ll] <- ova_coef_est_M
    hold_oal_coefH[ , , ll] <- ova_coef_est_H
    hold_oal_coefG[ , , ll] <- ova_coef_est_G
  }
  
  if (k == 1){
    out$outcomes$groups <- array(hold_grp, dim = c(N, p, k, l, q), 
                                 dimnames = list(grps, predictors, alphas, trt_lvls, qnames))
    #hold_grp: dim = c(N, p, k, l, q), dimnames = list(grps, predictors, alphas, trt_lvls, qnames)
    out$outcomes$overall <- array(hold_oal, dim = c(p, k, l, q),
                                  dimnames = list(predictors, alphas, trt_lvls, qnames))
    
    out$outcomes$grp_coefM <- array(hold_grp_coefM, dim = c(N, (1+len_h+len_g+len_h*len_g)*p, k, l),
                                    dimnames = list(grps, predictors, alphas, trt_lvls))
    out$outcomes$grp_coefH <- array(hold_grp_coefH, dim = c(N, (1+len_h)*p, k, l),
                                    dimnames = list(grps, predictors, alphas, trt_lvls))
    out$outcomes$grp_coefG <- array(hold_grp_coefG, dim = c(N, (1+len_g)*p, k, l),
                                    dimnames = list(grps, predictors, alphas, trt_lvls))
    
    out$outcomes$overall_coefM <- array(hold_oal_coefM, dim = c((1+len_h+len_g+len_h*len_g)*p, k, l),
                                        dimnames = list(predictors, alphas, trt_lvls))
    out$outcomes$overall_coefH <- array(hold_oal_coefH, dim = c((1+len_h)*p, k, l),
                                        dimnames = list(predictors, alphas, trt_lvls))
    out$outcomes$overall_coefG <- array(hold_oal_coefG, dim = c((1+len_g)*p, k, l),
                                        dimnames = list(predictors, alphas, trt_lvls))
    
    out$outcomes$weights_ind <- array(weights_ind, dim = c(length(A), p, k, l),
                                      dimnames = list(NULL, NULL, alphas, trt_lvls))
  }else{
    out$outcomes <- list(groups = drop(hold_grp), 
                         overall = drop(hold_oal),
                         grp_coefM = drop(hold_grp_coefM),
                         grp_coefH = drop(hold_grp_coefH),
                         grp_coefG = drop(hold_grp_coefG),
                         overall_coefM = drop(hold_oal_coefM),
                         overall_coefH = drop(hold_oal_coefH),
                         overall_coefG = drop(hold_oal_coefG),
                         weights_ind = weights_ind
    )
  }
  
  ## DONE ####
  return(out)
}

