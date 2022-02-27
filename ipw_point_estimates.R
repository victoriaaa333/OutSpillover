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

ipw_point_estimates <- function(H, G, A, weights, X = NULL, x0 = NULL){
  ## Necessary Bits ##
  grps     <- dimnames(weights)[[1]]
  alphas   <- dimnames(weights)[[length(dim(weights))]]
  numerator_alphas <- as.numeric(lapply(alphas, function(l) substr(l[1],3,5)))
  trt_lvls <- sort(unique(A)) # binary if A = 1/0
  grp_counts = as.vector(table(G))
  N <- length(grps)
  k <- length(alphas)
  l <- length(trt_lvls)
  
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
  Hbar <- group_means(H = H, A = A, G = G,
                      X = X, x0 = x0, a = NA) # average for each groups, marginalized over all treatments.
  
  grp_est <- apply(weights, 2:3, function(x) x * Hbar) 
  dimnames(grp_est) <- list(grps, predictors, alphas)
  
  oa_est <- apply(grp_est, 2:3, mean, na.rm = TRUE)
  
  if (k == 1){
    grp_est1 <- as.matrix(drop(grp_est))
    colnames(grp_est1) <- alphas
    out$marginal_outcomes$groups <- grp_est1
  }else{
    out$marginal_outcomes$groups <- drop(grp_est)
  }
  out$marginal_outcomes$overall <- drop(oa_est)
  
  ## CALCULATE OUTCOME ESTIMATES PER TREATMENT LEVEL####
  
  hold_grp <- array(dim = c(N, p, k, l), dimnames = list(grps, predictors, 
                                                         alphas, trt_lvls))
  hold_oal <- array(dim = c(p, k, l),
                    dimnames = list(predictors, alphas, trt_lvls))
  
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
    grp_est <- array(dim= c(N, p, k))
    grp_est[, p, ] <- apply(ind_est, 3, group_means, A, G, X, x0, a) 
    oal_est <- apply(grp_est, 2:3, mean, na.rm = TRUE)
    
    hold_grp[ , , , ll] <- grp_est
    hold_oal[ , , ll]   <- oal_est
  }
  if (k == 1){
    out$outcomes$groups <- array(hold_grp, dim = c(N, k, l), 
                                 dimnames = list(grps, alphas, trt_lvls))
    
    out$outcomes$overall <- array(hold_oal, dim = c(k, l),
                                    dimnames = list(alphas, trt_lvls))
  }else{
    out$outcomes <- list(groups = drop(hold_grp), 
                         overall = drop(hold_oal))
  }
  
  ## DONE ####
  return(out)
}
