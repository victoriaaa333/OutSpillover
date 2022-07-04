#-----------------------------------------------------------------------------#
# Calculate IPW point estimates
#  
# @param H vector of average of h-neighborhood
# @param G vector of group assignments
# @param A vector of treatment assignments
# @param weights weight matrix/array to use from \code{\link{wght_matrix}}
# @param X matrix of covariates
# @param x0 vector of specific covariates conditioned on 
# @return list containing point estimates for marginal outcomes and estimates
# per treatment level
# @export
#-----------------------------------------------------------------------------#

# At the moment, we assume weights is a matrix, and type I randomization

ipw_point_estimates_mixed <- function(H, G, A, weights, X = NULL, x0 = NULL, 
                                     neighinfo = NULL, x1= NULL, X_type = NULL){
  ## Necessary Bits ##
  grps     <- dimnames(weights)[[1]]
  alphas   <- dimnames(weights)[[length(dim(weights))]]
  numerator_alphas <- as.numeric(lapply(alphas, function(l) substr(l[1],3,5)))
  trt_lvls <- sort(unique(A)) # binary if A = 1/0
  #x1_lvls <- sort(unique(x1))
  #x0_lvls <- sort(unique(x0))
  
  N <- length(grps)
  k <- length(alphas)
  l <- length(trt_lvls)
  q <- ifelse(!is.null(x1), dim(x1)[2], #TODO: change this to dimension
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
  # 
  
  len_g <- ifelse(is.null(X), 1, sum(X_type == "N"))
  len_h <- ifelse(is.null(neighinfo), 1, dim(neighinfo$neighX)[2])
  ## CALCULATE OUTCOME ESTIMATES PER TREATMENT LEVEL####
  
  hold_grp <- array(dim = c(N, p, k, l, q), dimnames = list(grps, predictors, 
                                                            alphas, trt_lvls, qnames))
  hold_oal <- array(dim = c(p, k, l, q),
                    dimnames = list(predictors, alphas, trt_lvls, qnames))
  
  hold_grp_coefH <- array(dim = c(N, (3+2*len_h)*p, k, l), dimnames = list(grps, predictors, 
                                                                 alphas, trt_lvls))
  hold_grp_coefG <- array(dim = c(N, (1+len_g)*p, k, l), dimnames = list(grps, predictors, 
                                                                 alphas, trt_lvls))
  hold_oal_coefH <- array(dim = c((3+2*len_h)*p, k, l),
                          dimnames = list(predictors, alphas, trt_lvls))
  hold_oal_coefG <- array(dim = c((1+len_g)*p, k, l),
                          dimnames = list(predictors, alphas, trt_lvls))

  weighted_ind <- array(dim = c(length(A), 1, k, l), dimnames = list(NULL, NULL, alphas, trt_lvls))
  
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
    
    # Compute estimates
    ind_est <- apply(weights_trt, 2:3, function(x) x * H) 
    coef_est_H <- array(dim= c(N, (3+2*len_h)*p, k))
    coef_est_G <- array(dim= c(N, (1+len_g)*p, k))
    ova_coef_est_H <- array(dim= c(1, (3+2*len_h)*p, k))
    ova_coef_est_G <- array(dim= c(1, (1+len_g)*p, k))
    grp_est <- array(dim= c(N, p, k, q))
    
    ##############
    # Calculate the coefficients for conditional H and group means
    for (j in 1:k){
      ind_est_df <- ind_est[ , , j]
      weighted_ind[, , j, ll] <- ind_est_df
      if (!is.null(neighinfo)){
        H_list <- neigh_coefs_oncont2(ind_est_df, G, neighinfo, A, a)
        coef_est_H[ , , j] <- as.matrix(H_list[[1]])
        ova_coef_est_H[ , , j] <- as.matrix(H_list[[2]])
        predictors <- names(H_list[[2]])
      }
      if (!is.null(X)){
        G_list <- group_coefs_oncont2(ind_est_df, G, X, x0, A, a)
        coef_est_G[ , , j] <- as.matrix(G_list[[1]])
        ova_coef_est_G[ , , j] <- as.matrix(G_list[[2]])
        predictors <- names(G_list[[2]])
      }
    }
    
    for (j in 1:k) {
      ind_est_df <- ind_est[ , , j]
      if (!is.null(x1)){
        #grp_est[, p, j] <- neigh_means_oncont2(coef_est_H[ , , j], x1)
        grp_est[ , p, j, ] <- apply(as.matrix(x1), 2, neigh_means_oncont2, 
                                    cond_coef = coef_est_H[ , , j])
      }else if (!is.null(x0)){
        #grp_est[, p, j] <- group_means_oncont2(coef_est_G[ , , j], x0)
        grp_est[ , p, j, ] <- apply(as.matrix(x0), 2, group_means_oncont2, 
                                    cond_coef = coef_est_G[ , , j])
      }else{
        grp_est[, p, j, ] <- group_means_null(ind_est_df, G, A, a)
      }
    }
    
    
    oal_est <- apply(grp_est, 2:4, mean, na.rm = TRUE)
    
    hold_grp[ , , , ll, ] <- grp_est
    hold_oal[ , , ll, ]   <- oal_est
    
    hold_grp_coefH[ , , , ll] <- coef_est_H
    hold_grp_coefG[ , , , ll] <- coef_est_G
    hold_oal_coefH[ , , ll] <- ova_coef_est_H
    hold_oal_coefG[ , , ll] <- ova_coef_est_G
  }
  
  if (k == 1){
    out$outcomes$groups <- array(hold_grp, dim = c(N, k, l, q), 
                                 dimnames = list(grps, alphas, trt_lvls, qnames))
    
    out$outcomes$overall <- array(hold_oal, dim = c(k, l, q),
                                  dimnames = list(alphas, trt_lvls, qnames))
    
    # out$outcomes$grp_coefH <- hold_grp_coefH
    # 
    # out$outcomes$grp_coefG <- hold_grp_coefG
    # 
    # out$outcomes$overall_coefH <- hold_oal_coefH
    # 
    # out$outcomes$overall_coefG <- hold_oal_coefG
    
    out$outcomes$grp_coefH <- array(hold_grp_coefH, dim = c(N, (3+2*len_h)*p, k, l), 
                                    dimnames = list(grps, predictors, alphas, trt_lvls))
    out$outcomes$grp_coefG <- array(hold_grp_coefG, dim = c(N, (1+len_g)*p, k, l),
                                    dimnames = list(grps, predictors, alphas, trt_lvls))
    out$outcomes$overall_coefH <- array(hold_oal_coefH, dim = c((3+2*len_h)*p, k, l),
                            dimnames = list(predictors, alphas, trt_lvls))
    out$outcomes$overall_coefG <- array(hold_oal_coefG, dim = c((1+len_g)*p, k, l),
                            dimnames = list(predictors, alphas, trt_lvls))
    out$outcomes$weighted_ind <- array(weighted_ind, dim = c(length(A), 1, k, l), 
                                  dimnames = list(NULL, NULL, alphas, trt_lvls))
    # out$outcomes$grp_coefH <- array(hold_grp_coefH, dim = c(N, k, l), 
    #                                 dimnames = list(grps, alphas, trt_lvls))
    # 
    # out$outcomes$grp_coefG <- array(hold_grp_coefG, dim = c(N, k, l), 
    #                                 dimnames = list(grps, alphas, trt_lvls))
    # 
    # out$outcomes$overall_coefH <- array(hold_oal_coefH, dim = c(k, l),
    #                                     dimnames = list(alphas, trt_lvls))
    # 
    # out$outcomes$overall_coefG <- array(hold_oal_coefG, dim = c(k, l),
    #                                     dimnames = list(alphas, trt_lvls))
    
  }else{
    out$outcomes <- list(groups = drop(hold_grp), 
                         overall = drop(hold_oal),
                         grp_coefH = drop(hold_grp_coefH),
                         grp_coefG = drop(hold_grp_coefG),
                         overall_coefH = drop(hold_oal_coefH),
                         overall_coefG = drop(hold_oal_coefG),
                         weighted_ind = weighted_ind)
  }
  
  ## DONE ####
  return(out)
}

group_means_null <- function(ind_est_df, G, A, a){
    group_df <- as.data.frame(cbind(ind_est_df, G, A))
    colnames(group_df) <- c('ind_est', 'G','A')
    group_df <- as.data.frame(group_df[group_df$A == a,])
    fits <- lmList(ind_est ~ 1 | G, data=group_df) 
    cond_group_means <- coef(fits)[,1]
    cond_group_means <- as.matrix(cond_group_means)
    cond_group_means
}

# TODO: Output the coefficient so that we dont need to run regression everytime we change value of x1
neigh_means_oncont <- function(ind_est_df, G, neighinfo, x1, A, a){
  neigh2_treated = neighinfo$neigh2_treated
  neighX = neighinfo$neighX
  neigh2_treated_neighX = neighinfo$neigh2_treated_neighX
  colnames(neighX) <- paste("neighX", colnames(neighX), sep = "_")
  colnames(neigh2_treated_neighX) <- paste("neigh2_treated_neighX", colnames(neigh2_treated_neighX), sep = "_")
  
  if (!is.null(neighinfo) && !is.null(x1)){
    group_df <- as.data.frame(cbind(ind_est_df, G, neighX, 
                                    A, neigh2_treated, neigh2_treated_neighX))
    #colnames(group_df) <- c('ind_est', 'G', 'neighX', 
    #                        'A', 'neigh2_treated', 'neigh2_treated_neighX')
    colnames(group_df)[1] <- 'ind_est'
    group_df <- as.data.frame(group_df[group_df$A == a,])
    
    fits <- lmList(ind_est ~ neighX + neigh2_treated + neigh2_treated_neighX | G, data=group_df) 
    fits2 <- lmList(neigh2_treated ~ 1 | G, data=group_df)
    
    cond_group_means <- coef(fits)[,1] + coef(fits)[,2] * x1 + 
      coef(fits)[,3] * coef(fits2)[,1] + coef(fits)[,4] * coef(fits2)[,1] * x1
  }else{
    group_df <- as.data.frame(cbind(ind_est_df, G, A))
    colnames(group_df) <- c('ind_est', 'G', 'A')
    group_df <- as.data.frame(group_df[group_df$A == a,])
    fits <- lmList(ind_est ~ 1 | G, data=group_df)
    cond_group_means <- coef(fits)[,1]
  }
  
  cond_group_means <- as.matrix(cond_group_means)
}



### Output the coef
neigh_coefs_oncont2 <- function(ind_est_df, G, neighinfo, A, a){
  neigh2_treated = neighinfo$neigh2_treated
  neighX = neighinfo$neighX
  neigh2_treated_neighX = neighinfo$neigh2_treated_neighX
  colnames(neighX) <- paste("neighX", colnames(neighX), sep = "_")
  colnames(neigh2_treated_neighX) <- paste("neigh2_treated_neighX", colnames(neigh2_treated_neighX), sep = "_")
  
  if (!is.null(neighinfo) && !is.null(x1)){
    group_df <- as.data.frame(cbind(ind_est_df, G, A, neigh2_treated, neighX, 
                                    neigh2_treated_neighX))
    #colnames(group_df) <- c('ind_est', 'G', 'neighX', 
    #                        'A', 'neigh2_treated', 'neigh2_treated_neighX')
    colnames(group_df)[1] <- 'ind_est'
    group_df <- as.data.frame(group_df[group_df$A == a,])
    fits <- lmList(as.formula(paste("ind_est ~ ", paste(colnames(group_df)[-1:-3], collapse= "+ "), "| G")),
                   data = group_df)
    #fits <- lmList(ind_est ~ neighX + neigh2_treated + neigh2_treated_neighX | G, data=group_df) 
    fits2 <- lmList(neigh2_treated ~ 1 | G, data=group_df)
    cond_coef <- cbind(coef(fits), coef(fits2))
    
    #overall_fits <- lm(ind_est ~ neighX + neigh2_treated + neigh2_treated_neighX, data=group_df)
    overall_fits <- lm(as.formula(paste("ind_est ~ ", paste(colnames(group_df)[-1:-3], collapse= "+ "))),
                   data = group_df)
    overall_fits2 <- lm(neigh2_treated ~ 1, data=group_df)
    #overall_coef <- c(coef(overall_fits)[[1]], coef(overall_fits)[[2]], coef(overall_fits)[[3]],
    #                  coef(overall_fits)[[4]], coef(overall_fits2)[[1]])
    overall_coef <- c(coef(overall_fits), coef(overall_fits2))
    coef <- list(cond_coef, overall_coef, colnames(coef(fits)))
    
  }else{
    group_df <- as.data.frame(cbind(ind_est_df, G, A))
    colnames(group_df) <- c('ind_est', 'G','A')
    group_df <- as.data.frame(group_df[group_df$A == a,])
    
    fits <- lmList(ind_est ~ 1 | G, data=group_df)
    cond_coef <- coef(fits)[,1]
    
    overall_fits <- lm(ind_est ~ 1, data=group_df)
    overall_coef <- coef(overall_fits)[[1]]
    
    coef <- list(cond_coef, overall_coef,colnames(coef(fits)))
  }
  coef
}

neigh_means_oncont2 <- function(cond_coef, x1){
  x1_num <- as.numeric(x1[X_type == "N"])
  lenn <- length(x1_num)
  
  if (lenn > 0){
    #      apply(cond_coef[,3:(2+lenn)],2, function(x) x* x1_num)
    cond_group_means <- cond_coef[,1] + cond_coef[,2] * cond_coef[,dim(cond_coef)[2]] +
      as.matrix(cond_coef[,3:(2+lenn)]) %*% as.matrix(x1_num) + 
      as.matrix(cond_coef[,(3+lenn):(2+2*lenn)]) %*% as.matrix(x1_num) * cond_coef[,dim(cond_coef)[2]]
  }else{
    cond_group_means <- cond_coef[,1]
  }
  cond_group_means
}


group_coefs_oncont2 <- function(ind_est_df, G, X, x0, A, a){
  
  X_cat <- as.matrix(X[, X_type == "C"])
  X_num <- apply(as.matrix(X[, X_type == "N"]), 2, as.numeric)
  num_names <- colnames(X)[X_type == "N"]
  x0_cat <- x0[X_type == "C"]
  x0_num <- as.numeric(x0[X_type == "N"])
  
  trt_cond <- which(A == a)
  
  if (sum(X_type == "C") > 0){
    ind_cond <- which(apply(X_cat, 1, function(x) prod(x == x0_cat)) == 1)
  }else{
    ind_cond <- 1:dim(X_cat)[1]}
  
  if (sum(X_type == "N") > 0){
    group_df <- as.data.frame(cbind(ind_est_df, G, X_num))
    colnames(group_df) <- c('ind_est', 'G', num_names)
    fits_df <- group_df[intersect(ind_cond,trt_cond), ]
    fits <- lmList(as.formula(paste("ind_est ~ ", paste(num_names, collapse= "+ "), "| G")),
                   data = fits_df)
    
    #cond_group_means <- coef(fits)[,1] + as.matrix(coef(fits)[,-1])%*%as.matrix(x0_num)
    cond_coef <- as.matrix(coef(fits))
    overall_fits <- lm(as.formula(paste("ind_est ~ ", paste(num_names, collapse= "+ "))),
                       data=group_df)
    overall_coef <- as.vector(coef(overall_fits))
    coef <- list(cond_coef, overall_coef,colnames(coef(fits)))
  }else{
    group_df <- as.data.frame(cbind(ind_est_df, G))
    colnames(group_df)[1] <- 'ind_est' 
    fits_df <- group_df[intersect(ind_cond,trt_cond), ]
    fits <- lmList(ind_est ~ 1 | G, data=fits_df) 
    
    #cond_group_means <- coef(fits)[,1]
    cond_coef <- coef(fits)[,1]
    overall_fits <- lm(ind_est ~ 1, data=group_df)
    overall_coef <- coef(overall_fits)[[1]]
    coef <- list(cond_coef, overall_coef,colnames(coef(fits)))
  }
  
  coef
  
}

group_means_oncont2 <- function(cond_coef, x0){
  #x0_cat <- x0[X_type == "C"]
  x0_num <- as.numeric(x0[X_type == "N"])
  
  if (length(x0_num) > 0){
    cond_group_means <- cond_coef[,1] + as.matrix(cond_coef[,-1])%*%as.matrix(x0_num)
  }else{
    cond_group_means <- cond_coef[,1]
  }
  cond_group_means
}
