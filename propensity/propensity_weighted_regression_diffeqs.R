# change the last different estimating equations
# V matrix, use weights here
ipw_point_estimates_propensity_Vmatrix2 <- function(H, G, A, weights, objs, X = NULL, x0 = NULL, 
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

      # calculate the group average
      for (g in 1:length(grps)) {
        ind <- intersect(which(G == g), which(cond_ind == 1))
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
        reg_est_grps[g, , , as.character(a), ] <- apply(reg_est[ind, ,], 2, mean)
      }    
    }
  }
  
  ## DONE ####
  reg_est_grps <- (reg_est_grps[, , , 1, ] + reg_est_grps[, , , 2, ])
  return(reg_est_grps)
}

V_matrix_regression2 <- function(scores, 
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
      xx <- beta_est_eqs
    }
  } 
  else if(effect_type == 'outcome'){
    if(marginal == TRUE){
      #xx <- hold_grp[ , a1,] - hold_oal[a1, ]
    } else {
      xx <- beta_est_eqs[, , a1, t1,]
    }
  }
  
  ee <- cbind(scores, xx)
  V <- crossprod(ee)/N
  V
}

ipw_interference_regression2 <- function(propensity_integrand,
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
                             integrate_allocation = integrate_allocation
                        ))
  #score_args <- append(weight_args, list(first_assignments = first_assignments))
  
  #### Prepare output ####
  out <- list()  
  
  ## Compute Weights ##
  weights <- do.call(wght_matrix_second, args = weight_args)
  
  if(variance_estimation == 'robust'){
    weightd <- do.call(wght_deriv_array_second, args = append(weight_args, grad_args))
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
    out$U22_grps <- do.call(ipw_point_estimates_propensity_beta, args = beta_args)
    
    out$V_grps <- do.call(ipw_point_estimates_propensity_Vmatrix2, args = beta_args)
    
    out$scores          <- do.call(score_matrix_second, args = score_args)
  } 
  
  out$weights <- weights
  
  return(out)
}

ipw_effect_calc_regression2 <- function(obj, 
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
      pe          <- oal[,a1, t1, ] - oal[,a2, t2, ]
      pe_grp_diff <- (grp[, , a1, t1, ] - oal[,a1, t1, ]) - (grp[, , a2, t2, ] - oal[,a2, t2, ])
    }
  } else {
    if(marginal == TRUE){
      pe          <- oal[a1] 
      pe_grp_diff <- (grp[ , a1] - oal[a1])
    } else {
      pe          <- oal[, a1, t1, ] 
      pe_grp_diff <- (grp[ , , a1, t1, ] - oal[, a1, t1, ])
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
      V <- V_matrix_regression2(scores = obj$scores, 
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
      U_pe_grp <- Ugrp[ , , , a1, t1,] +  Ugrp[ , , , a2, t2,]
      U_beta_grp <- U22_grps[ , , , a1, t1,] + U22_grps[ , , , a2, t2,]
      
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

ipw_propensity_variance_regression2 <- function(parameters,
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
  ipw <- do.call(ipw_interference_regression2, args = ipw_args)
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
  est <- do.call(ipw_effect_calc_regression2, args = regression_args)
  return(est)
}


