# coefs <- out$point_estimates$outcomes$overall_coefG from second_stages_utils.R line 322
# X <- cond_X

# U21
ipw_point_estimates_propensity_regression <- function(H, G, A, weights, coefs, X = NULL, x0 = NULL, 
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
      X_cat <- as.matrix(X[, X_type == "C"])
      X_num <- apply(as.matrix(X[, X_type == "N"]), 2, as.numeric)
      X_num_intercept <- cbind(1, X_num)
      coef_a <- as.matrix(coefs[, , as.character(a)])
      
      # reg_a is X'*(H - X \beta) instead of H
      reg_a <- (H - (X_num_intercept %*% coef_a))
      
      reg_est <- array(dim = c(length(A), p, len_g+1, k))
      for (i in 1:length(A)) {
        reg_est[i, , ,] <- reg_a[i] * weights_trt[i,,] %*% t(X_num_intercept[i,])
      }
      cond_ind <- ifelse(X_cat == x0[X_type == "C"], 1, 0) 
      
      reg_est_grp <- array(dim = c(length(grps), p, len_g+1, k))
      # calculate the group average
      for (g in 1:length(grps)) {
        ind <- intersect(which(G == g), which(cond_ind == 1))
        reg_est_grp[g, , ,] <- apply(reg_est[ind, , ,], 2:3, mean)
      }
      
      reg_est_overall[, , k, ll,] <- apply(reg_est_grp, 2:3, mean) 
      #dim( out$Upart$outcomes$overall) 8 1 2 1 #  NULL "c(0.5, 0.4, 0.6)" [1] "0" "1" [1] "No-con"
      }
    }
  
  ## DONE ####
  return(reg_est_overall)
}

# U22      
# the derivative to \beta, replace X'*(H - X \beta) with X'* - X, still use weights instead of weightd
ipw_point_estimates_propensity_beta <- function(H, G, A, weights, coefs, X = NULL, x0 = NULL, 
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
      X_cat <- as.matrix(X[, X_type == "C"])
      X_num <- apply(as.matrix(X[, X_type == "N"]), 2, as.numeric)
      X_num_intercept <- cbind(1, X_num)
      coef_a <- as.matrix(coefs[, , as.character(a)])
      
      # reg_a is - X instead of H
      reg_a <- -1*X_num_intercept
      
      reg_est <- array(dim = c(length(A), len_g+1, len_g+1, k))
      for (i in 1:length(A)) {
        reg_est[i, , ,] <- weights_trt[i,,] * reg_a[i,]  %*% t(X_num_intercept[i,])
      }
      cond_ind <- ifelse(X_cat == x0[X_type == "C"], 1, 0) 
      
      reg_est_grp <- array(dim = c(length(grps), len_g+1, len_g+1, k))
      # calculate the group average
      for (g in 1:length(grps)) {
        ind <- intersect(which(G == g), which(cond_ind == 1))
        reg_est_grp[g, , ,] <- apply(reg_est[ind, , ,], 2:3, mean)
      }
      
      reg_est_overall[, , k, ll,] <- apply(reg_est_grp, 2:3, mean) 
      #dim( out$Upart$outcomes$overall) 8 1 2 1 #  NULL "c(0.5, 0.4, 0.6)" [1] "0" "1" [1] "No-con"
    }
  }
  
  ## DONE ####
  return(reg_est_overall)
}

# V matrix, use weights here
ipw_point_estimates_propensity_Vmatrix <- function(H, G, A, weights, coefs, X = NULL, x0 = NULL, 
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
      X_cat <- as.matrix(X[, X_type == "C"])
      X_num <- apply(as.matrix(X[, X_type == "N"]), 2, as.numeric)
      X_num_intercept <- cbind(1, X_num)
      coef_a <- as.matrix(coefs[, , as.character(a)])
      
      # reg_a is X' (H - X \beta)
      reg_a <-  (H - (X_num_intercept %*% coef_a))
      
      reg_est <- array(dim = c(length(A), len_g+1, k))
      for (i in 1:length(A)) {
        reg_est[i, ,] <- weights_trt[i,,] * reg_a[i,]  %*% t(X_num_intercept[i,])
      }
      cond_ind <- ifelse(X_cat == x0[X_type == "C"], 1, 0) 
      
      #reg_est_grp <- array(dim = c(length(grps), len_g+1, k))
      # calculate the group average
      for (g in 1:length(grps)) {
        ind <- intersect(which(G == g), which(cond_ind == 1))
        reg_est_grps[g, , , as.character(a), ] <- apply(reg_est[ind, ,], 2, mean)
        }    
      }
  }
  
  ## DONE ####
  return(reg_est_grps)
}

beta_est_eqs <- reg_est_grps

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
      xx <- beta_est_eqs[, , a1, t1,] -  beta_est_eqs[, , a2, t2,]
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
                                   print = FALSE)
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
    Uoal <- obj$Upart[[fff]]$overall 
    Ugrp <- obj$Upart[[fff]]$groups
    
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
      pe          <- oal[a1, t1, ] - oal[a2, t2, ]
      pe_grp_diff <- (grp[ , a1, t1, ] - oal[a1, t1, ]) - (grp[ , a2, t2, ] - oal[a2, t2, ])
      
      # pe          <- oal[, a1, t1, ] - oal[, a2, t2,]
      # pe_grp_diff <- (grp[, , a1, t1, ] - oal[,a1, t1,]) - (grp[ , , a2, t2,] - oal[, a2, t2,])
      
      if(variance_estimation == 'robust'){
        U_pe_grp    <- as.matrix(Ugrp[, , a1, t1, ] - Ugrp[, , a2, t2, ])
        #as.matrix(Ugrp[ , a1, t1, ] - Ugrp[ , a2, t2, ])
      }
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
        U_pe_grp    <- Ugrp[ , , a1, t1, ]
      }
    }
  }
  
  #### VARIANCE ESTIMATION ####
  if(variance_estimation == 'robust'){
    # partial U matrix
    if(p == 1){
      U21 <- sum(-U_pe_grp)/N
    } else {
      U21 <- (t(as.matrix(apply(-U_pe_grp, 2, sum, na.rm = T))))/N
    }
    
    # U11
    U11 <- obj$U11
    U <- cbind(rbind(U11, U21), c(rep(0, dim(U11)[1]), -1))
    # V matrix
    V <- V_matrix_second(scores = obj$scores, 
                         point_estimates = obj$point_estimates, 
                         allocation1 = a1, allocation2 = a2, 
                         trt.lvl1 = t1, trt.lvl2 = t2, 
                         effect_type = effect_type, marginal = marginal)
    
    vdim <- dim(V)[1]
    
    V21 <- V[vdim, 1:(vdim - 1)] # Last row, up to last column
    V11 <- V[1:(vdim - 1), 1:(vdim - 1)] # up to last row, up to last column
    V22 <- V[vdim, vdim] # bottom right element
    
    ## Sandwich Variance Estimate ##
    inv_U = solve(U) 
    sigma = inv_U %*% V %*% t(inv_U)
    
    #ave <- sigma[dim(sigma)[1], dim(sigma)[2]]/N * rescale.factor^2
    ave <- ((U21 - 2*V21) %*% solve(V11) %*% t(U21) + V22)/N * rescale.factor^2
    
  } else if(variance_estimation == 'naive'){
    ave <- (1/(N^2)) * (sum((pe_grp_diff)^2, na.rm = T)) * rescale.factor^2
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