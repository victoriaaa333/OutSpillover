# tests for group_estimates and direct overall

ipw_point_estimates_mixed_test4 <- function(H, G, A, weights, X = NULL, x0 = NULL, 
                                            neighinfo = NULL, x1= NULL, X_type = NULL,
                                            Con_type = "No-con"){
  ## Necessary Bits ##
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
    predictors <- dimnames(weights)[[2]]
  }
  
  out <- list()
  
  ## CALCULATE MARGINAL ESTIMATES ####
  # Hbar <- group_means(H = H, A = A, G = G,
  #                     X = X, x0 = x0, a = NA) # average for each groups, marginalized over all treatments.
  # 
  # grp_est <- apply(weights, 2:3, function(x) x * Hbar) 
  # dimnames(grp_est) <- list(grps, predictors, alphas)
  # 
  # oa_est <- apply(grp_est, 2:3, mean, na.rm = TRUE)
  # 
  # if (k == 1){
  #   grp_est1 <- as.matrix(drop(grp_est))
  #   colnames(grp_est1) <- alphas
  #   out$marginal_outcomes$groups <- grp_est1
  # }else{
  #   out$marginal_outcomes$groups <- drop(grp_est)
  # }
  # out$marginal_outcomes$overall <- drop(oa_est)
  
  len_g <- ifelse(is.null(X), 1, sum(X_type == "N"))
  len_h <- ifelse(is.null(neighinfo), 1, dim(neighinfo$neighX)[2])
  
  ## CALCULATE OUTCOME ESTIMATES PER TREATMENT LEVEL####
  
  hold_grp <- array(dim = c(N, p, k, l, q), dimnames = list(grps, predictors, 
                                                            alphas, trt_lvls, qnames))
  hold_oal <- array(dim = c(p, k, l, q),
                    dimnames = list(predictors, alphas, trt_lvls, qnames))
  
  hold_grp_coefM <- array(dim = c(N, (1+len_h+len_g)*p, k, l), dimnames = list(grps, predictors, 
                                                                               alphas, trt_lvls))
  hold_grp_coefH <- array(dim = c(N, (1+len_h)*p, k, l), dimnames = list(grps, predictors, 
                                                                         alphas, trt_lvls))
  hold_grp_coefG <- array(dim = c(N, (1+len_g)*p, k, l), dimnames = list(grps, predictors, 
                                                                         alphas, trt_lvls))
  hold_oal_coefM <- array(dim = c((1+len_h+len_g)*p, k, l),
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
          (numerator_alphas[kk]^a * (1-numerator_alphas[kk])^(1-a)) 
      }
    }
    
    weights_ind[, pp, kk, ll] <- weights_trt[ , pp, kk]
    
    # Compute estimates
    ind_est <- apply(weights_trt, 2:3, function(x) x * H) 
    coef_est_M <- array(dim= c(N, (1+len_h+len_g)*p, k))
    coef_est_H <- array(dim= c(N, (1+len_h)*p, k))
    coef_est_G <- array(dim= c(N, (1+len_g)*p, k))
    
    ova_coef_est_M <- array(dim= c(1, (1+len_h+len_g)*p, k))
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
        # TODO: check this
        H_list <- neigh_coefs_oncont4(weights_df, H, G, neighinfo, A, a)
        #H_list <- neigh_coefs_oncont_old(ind_est_df, H, G, neighinfo, A, a)
        coef_est_H[ , , j] <- as.matrix(H_list[[1]])
        ova_coef_est_H[ , , j] <- as.matrix(H_list[[2]])
      }else if (Con_type == "group"){
        G_list <- group_coefs_oncont1(weights_df, H, G, X, X_type, x0, A, a)
        #G_list <- group_coefs_oncont2(weights_df, H, G, X, X_type, x0, A, a)
        coef_est_G[ , , j] <- as.matrix(G_list[[1]])
        ova_coef_est_G[ , , j] <- as.matrix(G_list[[2]])
      }
    }
    
    for (j in 1:k) {
      ind_est_df <- ind_est[ , , j]
      if (Con_type == "mixed"){
        grp_est[ , p, j, ] <- apply(as.matrix(x1), 2, mixed_means, 
                                    overall_coef = ova_coef_est_M[ , , j], 
                                    X_type = X_type, x0 = x0) 
        #TODO: fill out the function for mixed means 
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
        grp_est_overall[, p, j, ] <- group_means_null(ind_est_df, G, A, a)
        grp_est[, p, j, ] <- group_means_null(ind_est_df, G, A, a)
      }
    }
    
    
    oal_est <- apply(grp_est, 2:4, mean, na.rm = TRUE)
    #oal_est <- apply(grp_est_overall, 2:4, mean, na.rm = TRUE)
    
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
    out$outcomes$groups <- array(hold_grp, dim = c(N, k, l, q), 
                                 dimnames = list(grps, alphas, trt_lvls, qnames))
    
    out$outcomes$overall <- array(hold_oal, dim = c(k, l, q),
                                  dimnames = list(alphas, trt_lvls, qnames))
    
    out$outcomes$grp_coefM <- array(hold_grp_coefM, dim = c(N, (1+len_h+len_g)*p, k, l), 
                                    dimnames = list(grps, predictors, alphas, trt_lvls))
    out$outcomes$grp_coefH <- array(hold_grp_coefH, dim = c(N, (1+len_h)*p, k, l), 
                                    dimnames = list(grps, predictors, alphas, trt_lvls))
    out$outcomes$grp_coefG <- array(hold_grp_coefG, dim = c(N, (1+len_g)*p, k, l),
                                    dimnames = list(grps, predictors, alphas, trt_lvls))
    
    out$outcomes$overall_coefM <- array(hold_oal_coefM, dim = c((1+len_h+len_g)*p, k, l),
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
                         weights_ind = weights_ind)
  }
  
  ## DONE ####
  return(out)
}

# weights outside without extra terms
#function(weights_df, H, G, neighinfo, A, a)
neigh_coefs_oncont4 <- function(weights_df, H, G, neighinfo, A, a){
  neighX = neighinfo$neighX
  colnames(neighX) <- paste("neighX", colnames(neighX), sep = "_")
  
  if (!is.null(neighinfo) && !is.null(x1)){
    group_df <- as.data.frame(cbind(weights_df, H, G, A, neighX))
    colnames(group_df)[1] <- 'w'
    
    fits_df <- as.data.frame(group_df[group_df$A == a,])
    fits <- lmList(as.formula(paste("H ~ ", paste(colnames(fits_df)[-1:-4], 
                                                  collapse= "+ "), "| G")), 
                   data = fits_df, weights = w)
    cond_coef <- as.matrix(coef(fits))
    
    overall_fits <- lm(as.formula(paste("H ~ ", paste(colnames(fits_df)[-1:-4],
                                                      collapse= "+ "))), 
                       data = fits_df, weights = w)
    overall_coef <- as.matrix(coef(overall_fits))
    coef <- list(cond_coef, overall_coef, colnames(coef(fits)))
    
  }else{
    group_df <- as.data.frame(cbind(weights_df, H, G, A))
    colnames(group_df) <- c('w', 'H', 'G',' A')
    fits_df <- as.data.frame(group_df[group_df$A == a,])
    
    fits <- lmList(H ~ 1 | G, data=fits_df, weights = w)
    cond_coef <- coef(fits)[,1]
    
    overall_fits <- lm(H ~ 1, data=fits_df, weights = w)
    overall_coef <- coef(overall_fits)[[1]]
    
    coef <- list(cond_coef, overall_coef, colnames(coef(fits)))
  }
  coef
}

neigh_means_oncont4 <-function(overall_coef, X_type, x1, N){
  x1_num <- as.numeric(x1[X_type == "N"])
  lenn <- length(x1_num)
  cond_coefs <- matrix(rep(overall_coef, N), nrow = N, byrow = TRUE)
  if (lenn > 0){
    cond_group_means <- cond_coefs[,1] + as.matrix(cond_coefs[,-1])%*%as.matrix(x1_num)
  }else{
    cond_group_means <- cond_coefs
  }
  cond_group_means
}