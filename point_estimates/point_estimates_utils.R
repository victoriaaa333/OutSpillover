######
# functions used in test4 and test5 (group avg vs one reg)
#####
# 1. non-condition
group_means_null <- function(ind_est_df, G, A, a){
  group_df <- as.data.frame(cbind(ind_est_df, G, A))
  colnames(group_df) <- c('ind_est', 'G','A')
  group_df <- as.data.frame(group_df[group_df$A == a,])
  fits <- lmList(ind_est ~ 1 | G, data=group_df) 
  cond_group_means <- coef(fits)[,1]
  cond_group_means <- as.matrix(cond_group_means)
  cond_group_means
}

#####
# 2. group avg
## 1. conditional within groups
group_coefs_oncont1 <- function(weights_df, H, G, X, X_type, x0, A, a){
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
    group_df <- as.data.frame(cbind(weights_df, H, G, X_num))
    colnames(group_df) <- c('w', 'H', 'G', num_names)
    group_df$weighted_H <- group_df$w * group_df$H
    fits_df <- group_df[intersect(ind_cond,trt_cond), ]
    # fits <- lmList(as.formula(paste("H ~ ", paste(num_names, collapse= "+ "), "| G")),
    #                weights = w, data = fits_df)
    fits <- lmList(as.formula(paste("weighted_H ~ ", paste(num_names, 
                                collapse= "+ "), "| G")), data = fits_df)
    cond_coef <- as.matrix(coef(fits))
    overall_fits <- lm(as.formula(paste("H ~ ", paste(num_names, collapse= "+ "))),
                       weights = w, data=fits_df)
    overall_coef <- as.vector(coef(overall_fits))
    
    # TODO: fix this when cond_coef with dim > 2 
    if (dim(cond_coef)[1] < length(unique(G))){
      cond_coef <- c(cond_coef, rep(NA, (length(unique(G))-
                                           dim(cond_coef)[1])))}
    pvals <- as.vector(summary(overall_fits)$coefficients[,4] )
    coef <- list(cond_coef, overall_coef,pvals)
  }else{
    group_df <- as.data.frame(cbind(weights_df, H, G))
    colnames(group_df)[1] <- 'w' 
    group_df$weighted_H <- group_df$w * group_df$H
    fits_df <- group_df[intersect(ind_cond,trt_cond), ]
    fits <- lmList(weighted_H ~ 1 | G, data=fits_df) 
    #fits <- lmList(H ~ 1 | G, weights = w, data=fits_df) 
    #cond_group_means <- coef(fits)[,1]
    cond_coef <- coef(fits)[,1]
    
    if (length(cond_coef) < length(unique(G))){
      cond_coef <- c(cond_coef, rep(NA, (length(unique(G))-
                                           length(cond_coef))))}
    
    overall_fits <- lm(H ~ 1, weights = w, data=fits_df)
    overall_coef <- coef(overall_fits)[[1]]
    pvals <- as.vector(summary(overall_fits)$coefficients[,4])
    coef <- list(cond_coef, overall_coef, pvals)
  }
  coef
}

# coef same for different groups
group_means_oncont1 <- function(overall_coef, X_type, x0, N){
  x0_num <- as.numeric(x0[X_type == "N"])
  cond_coefs <- matrix(rep(overall_coef, N), nrow = N, byrow = TRUE)
  if (length(x0_num) > 0){
    cond_group_means <- cond_coefs[,1] + as.matrix(cond_coefs[,-1])%*%as.matrix(x0_num)
  }else{
    cond_group_means <- cond_coefs
  }
  cond_group_means
}

# coef different for different groups
group_means_oncont2 <- function(cond_coefs, X_type, x0){
  x0_num <- as.numeric(x0[X_type == "N"])
  if (length(x0_num) > 0){
    cond_group_means <- cond_coefs[,1] + as.matrix(cond_coefs[,-1])%*%as.matrix(x0_num)
  }else{
    cond_group_means <- cond_coefs
  }
  cond_group_means
}

#####
# 3. neigh avg
# weights inside, because for neigh avg solely, the weight for each group
# is the same, so "weights = w" is ignored in lmlist
neigh_coefs_oncont5 <- function(weights_df, H, G, neighinfo, A, a){
  neighX = neighinfo$neighX
  colnames(neighX) <- paste("neighX", colnames(neighX), sep = "_")
  
  if (!is.null(neighinfo) && !is.null(x1)){
    group_df <- as.data.frame(cbind(weights_df, H, G, A, neighX))
    colnames(group_df)[1] <- 'w'
    
    fits_df <- as.data.frame(group_df[group_df$A == a,])
    fits_df$weighted_H <- fits_df$w * fits_df$H
    fits <- lmList(as.formula(paste("weighted_H ~ ", 
                                    paste(colnames(fits_df)[5:(4+length(colnames(neighX)))], 
                                          collapse= "+ "), "| G")), data = fits_df)
    cond_coef <- as.matrix(coef(fits))
    
    overall_fits <- lm(as.formula(paste("H ~ ", 
                                        paste(colnames(fits_df)[5:(4+length(colnames(neighX)))],
                                              collapse= "+ "))), data = fits_df, weights = w)
    overall_coef <- as.matrix(coef(overall_fits))
    coef <- list(cond_coef, overall_coef, colnames(coef(fits)))
    
  }else{
    group_df <- as.data.frame(cbind(weights_df, H, G, A))
    colnames(group_df) <- c('w', 'H', 'G',' A')
    fits_df <- as.data.frame(group_df[group_df$A == a,])
    fits_df$weighted_H <- fits_df$w * fits_df$H
    
    fits <- lmList(weighted_H ~ 1 | G, data=fits_df)
    # fits <- lmList(H ~ 1 | G, data=fits_df, weights = w)
    cond_coef <- coef(fits)[,1]
    # overall_fits <- lm(weighted_H ~ 1, data=fits_df)
    overall_fits <- lm(H ~ 1, data=fits_df, weights = w)
    overall_coef <- coef(overall_fits)[[1]]
    
    coef <- list(cond_coef, overall_coef, colnames(coef(fits)))
  }
  coef
}

# coef same for different neigh groups
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

# coef different for different neigh groups
neigh_means_oncont3 <- function(cond_coefs, X_type, x1){
  x1_num <- as.numeric(x1[X_type == "N"])
  lenn <- length(x1_num)
  if (lenn > 0){
    cond_group_means <- cond_coefs[,1] + as.matrix(cond_coefs[,-1])%*%as.matrix(x1_num)
  }else{
    cond_group_means <- cond_coefs
  }
  cond_group_means
}

###### old functions not used in point estimates now ######
# old functions not used in point estimates now
group_means_null0 <- function(ind_est_df, G, A, a){
  group_df <- as.data.frame(cbind(ind_est_df, G, A))
  colnames(group_df) <- c('ind_est', 'G','A')
  group_df <- as.data.frame(group_df[group_df$A == a,])
  fits <- lm(ind_est ~ 1, data=group_df) 
  cond_group_means <- rep(as.numeric(coef(fits)[1]), length(unique(G)))
  cond_group_means <- as.matrix(cond_group_means)
  cond_group_means
}

group_means_null_1 <- function(weights_df, H, G, A, a){
  group_df <- as.data.frame(cbind(weights_df, H, G, A))
  colnames(group_df) <- c('w', 'H', 'G', 'A')
  group_df <- as.data.frame(group_df[group_df$A == a,])
  fits <- lmList(H ~ 1 | G, data=group_df, weights = w) 
  cond_group_means <- coef(fits)[,1]
  cond_group_means <- as.matrix(cond_group_means)
  cond_group_means
}

group_means_null_0 <- function(weights_df, H, G, A, a){
  group_df <- as.data.frame(cbind(weights_df, H, G, A))
  colnames(group_df) <- c('w', 'H', 'G', 'A')
  group_df <- as.data.frame(group_df[group_df$A == a,])
  fits <- lm(H ~ 1, data=group_df, weights = w) 
  cond_group_means <- rep(as.numeric(coef(fits)[1]), length(unique(G)))
  cond_group_means <- as.matrix(cond_group_means)
  cond_group_means
}

## 1. conditional within groups
group_coefs_oncont1 <- function(weights_df, H, G, X, X_type, x0, A, a){
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
    group_df <- as.data.frame(cbind(weights_df, H, G, X_num))
    colnames(group_df) <- c('w', 'H', 'G', num_names)
    group_df$weighted_H <- group_df$w * group_df$H
    fits_df <- group_df[intersect(ind_cond,trt_cond), ]
    # fits <- lmList(as.formula(paste("H ~ ", paste(num_names, collapse= "+ "), "| G")),
    #                weights = w, data = fits_df)
    fits <- lmList(as.formula(paste("weighted_H ~ ", paste(num_names, 
                                                           collapse= "+ "), "| G")), data = fits_df)
    cond_coef <- as.matrix(coef(fits))
    overall_fits <- lm(as.formula(paste("H ~ ", paste(num_names, collapse= "+ "))),
                       weights = w, data=fits_df)
    overall_coef <- as.vector(coef(overall_fits))
    
    # TODO: fix this when cond_coef with dim > 2 
    if (dim(cond_coef)[1] < length(unique(G))){
      cond_coef <- c(cond_coef, rep(NA, (length(unique(G))-
                                           dim(cond_coef)[1])))}
    pvals <- as.vector(summary(overall_fits)$coefficients[,4] )
    coef <- list(cond_coef, overall_coef,pvals)
  }else{
    group_df <- as.data.frame(cbind(weights_df, H, G))
    colnames(group_df)[1] <- 'w' 
    group_df$weighted_H <- group_df$w * group_df$H
    fits_df <- group_df[intersect(ind_cond,trt_cond), ]
    fits <- lmList(weighted_H ~ 1 | G, data=fits_df) 
    #fits <- lmList(H ~ 1 | G, weights = w, data=fits_df) 
    #cond_group_means <- coef(fits)[,1]
    cond_coef <- coef(fits)[,1]
    
    if (length(cond_coef) < length(unique(G))){
      cond_coef <- c(cond_coef, rep(NA, (length(unique(G))-
                                           length(cond_coef))))}
    
    overall_fits <- lm(H ~ 1, weights = w, data=fits_df)
    overall_coef <- coef(overall_fits)[[1]]
    pvals <- as.vector(summary(overall_fits)$coefficients[,4])
    coef <- list(cond_coef, overall_coef, pvals)
  }
  coef
}

# weights inside
group_coefs_oncont2 <- function(ind_est_df, H, G, X, X_type, x0, A, a){
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
    group_df <- as.data.frame(cbind(ind_est_df, H, G, X_num))
    colnames(group_df) <- c('ind_est', 'H', 'G', num_names)
    fits_df <- group_df[intersect(ind_cond,trt_cond), ]
    fits <- lmList(as.formula(paste("ind_est ~ ", paste(num_names, collapse= "+ "), "| G")), data = fits_df)
    cond_coef <- as.matrix(coef(fits))
    overall_fits <- lm(as.formula(paste("ind_est ~ ", paste(num_names, collapse= "+ "))), data=fits_df)
    overall_coef <- as.vector(coef(overall_fits))
    
    # TODO: fix this when cond_coef with dim > 2 
    if (dim(cond_coef)[1] < length(unique(G))){
      cond_coef <- c(cond_coef, rep(NA, (length(unique(G))-
                                           dim(cond_coef)[1])))}
    pvals <- as.vector(summary(overall_fits)$coefficients[,4] )
    coef <- list(cond_coef, overall_coef,pvals)
    #colnames(coef(fits)))
  }else{
    group_df <- as.data.frame(cbind(weights_df, H, G))
    colnames(group_df)[1] <- 'w' 
    group_df$weighted_H <- group_df$w * group_df$H
    fits_df <- group_df[intersect(ind_cond,trt_cond), ]
    fits <- lmList(weighted_H ~ 1 | G, data=fits_df) 
    #fits <- lmList(H ~ 1 | G, weights = w, data=fits_df) 
    #cond_group_means <- coef(fits)[,1]
    cond_coef <- coef(fits)[,1]
    
    if (length(cond_coef) < length(unique(G))){
      cond_coef <- c(cond_coef, rep(NA, (length(unique(G))-
                                           length(cond_coef))))}
    
    overall_fits <- lm(H ~ 1, weights = w, data=fits_df)
    overall_coef <- coef(overall_fits)[[1]]
    pvals <- as.vector(summary(overall_fits)$coefficients[,4])
    coef <- list(cond_coef, overall_coef, pvals)
  }
  coef
}


## 2. conditional within neighbors
## neigh coefficient (weights outside) ##
# with extra terms

neigh_coefs_oncont1 <- function(weights_df, H, G, neighinfo, A, a){
  neigh2_treated = neighinfo$neigh2_treated
  neighX = neighinfo$neighX
  neigh2_treated_neighX = neighinfo$neigh2_treated_neighX
  colnames(neighX) <- paste("neighX", colnames(neighX), sep = "_")
  colnames(neigh2_treated_neighX) <- paste("neigh2_treated_neighX", colnames(neigh2_treated_neighX), sep = "_")
  
  if (!is.null(neighinfo) && !is.null(x1)){
    group_df <- as.data.frame(cbind(weights_df, H, G, A, neigh2_treated, neighX, 
                                    neigh2_treated_neighX))
    colnames(group_df)[1] <- 'w'
    
    fits_df <- as.data.frame(group_df[group_df$A == a,])
    fits <- lmList(as.formula(paste("H ~ ", paste(colnames(fits_df)[-1:-4], 
                                                  collapse= "+ "), "| G")),
                   weights = w, data = fits_df)
    fits2 <- lmList(neigh2_treated ~ 1 | G, data=fits_df)
    cond_coef <- cbind(coef(fits), coef(fits2))
    
    overall_fits <- lm(as.formula(paste("H ~ ", paste(colnames(fits_df)[-1:-4],
                                                      collapse= "+ "))),
                       weights = w, data = fits_df)
    overall_fits2 <- lm(neigh2_treated ~ 1, data=fits_df)
    
    overall_coef <- c(coef(overall_fits), coef(overall_fits2))
    coef <- list(cond_coef, overall_coef, colnames(coef(fits)))
    
  }else{
    group_df <- as.data.frame(cbind(weights_df, H, G, A))
    colnames(group_df) <- c('w', 'H', 'G',' A')
    fits_df <- as.data.frame(group_df[group_df$A == a,])
    
    fits <- lmList(H ~ 1 | G, weights = w, data=fits_df)
    cond_coef <- coef(fits)[,1]
    
    overall_fits <- lm(H ~ 1, weights = w, data=fits_df)
    overall_coef <- coef(overall_fits)[[1]]
    
    coef <- list(cond_coef, overall_coef, colnames(coef(fits)))
  }
  coef
}

neigh_coefs_oncont_old <- function(ind_est_df, H, G, neighinfo, A, a){
  neigh2_treated = neighinfo$neigh2_treated
  neighX = neighinfo$neighX
  neigh2_treated_neighX = neighinfo$neigh2_treated_neighX
  colnames(neighX) <- paste("neighX", colnames(neighX), sep = "_")
  colnames(neigh2_treated_neighX) <- paste("neigh2_treated_neighX", colnames(neigh2_treated_neighX), sep = "_")
  
  if (!is.null(neighinfo) && !is.null(x1)){
    group_df <- as.data.frame(cbind(ind_est_df, H, G, A, neigh2_treated, neighX, 
                                    neigh2_treated_neighX))
    colnames(group_df)[1] <- 'ind_est'
    
    fits_df <- as.data.frame(group_df[group_df$A == a,])
    fits <- lmList(as.formula(paste("ind_est ~ ", paste(colnames(fits_df)[-1:-4], 
                                                        collapse= "+ "), "| G")),
                   data = fits_df)
    fits2 <- lmList(neigh2_treated ~ 1 | G, data=fits_df)
    cond_coef <- cbind(coef(fits), coef(fits2))
    
    overall_fits <- lm(as.formula(paste("ind_est ~ ", paste(colnames(fits_df)[-1:-4],
                                                            collapse= "+ "))),
                       data = fits_df)
    overall_fits2 <- lm(neigh2_treated ~ 1, data=fits_df)
    
    overall_coef <- c(coef(overall_fits), coef(overall_fits2))
    coef <- list(cond_coef, overall_coef, colnames(coef(fits)))
    
  }else{
    group_df <- as.data.frame(cbind(ind_est_df, H, G, A))
    colnames(group_df) <- c('ind_est', 'H', 'G',' A')
    fits_df <- as.data.frame(group_df[group_df$A == a,])
    
    fits <- lmList(ind_est ~ 1 | G, data=fits_df)
    cond_coef <- coef(fits)[,1]
    
    overall_fits <- lm(ind_est ~ 1, data=fits_df)
    overall_coef <- coef(overall_fits)[[1]]
    
    coef <- list(cond_coef, overall_coef, colnames(coef(fits)))
  }
  coef
}

######
neigh_means_oncont1_eachg <- function(cond_coefs, X_type, x1){
  x1_num <- as.numeric(x1[X_type == "N"])
  lenn <- length(x1_num)
  #cond_coefs <- matrix(rep(overall_coef, N), nrow = N, byrow = TRUE)
  
  if (lenn > 0){
    cond_group_means <- cond_coefs[,1] + cond_coefs[,2] * cond_coefs[,dim(cond_coefs)[2]] +
      as.matrix(cond_coefs[,3:(2+lenn)]) %*% as.matrix(x1_num) + 
      as.matrix(cond_coefs[,(3+lenn):(2+2*lenn)]) %*% as.matrix(x1_num) * cond_coefs[,dim(cond_coefs)[2]]
  }else{
    cond_group_means <- cond_coefs
  }
  cond_group_means
}

neigh_means_oncont1 <- function(overall_coef, X_type, x1, N){
  x1_num <- as.numeric(x1[X_type == "N"])
  lenn <- length(x1_num)
  cond_coefs <- matrix(rep(overall_coef, N), nrow = N, byrow = TRUE)
  
  if (lenn > 0){
    cond_group_means <- cond_coefs[,1] + cond_coefs[,2] * cond_coefs[,dim(cond_coefs)[2]] +
      as.matrix(cond_coefs[,3:(2+lenn)]) %*% as.matrix(x1_num) + 
      as.matrix(cond_coefs[,(3+lenn):(2+2*lenn)]) %*% as.matrix(x1_num) * cond_coefs[,dim(cond_coefs)[2]]
  }else{
    cond_group_means <- cond_coefs
  }
  cond_group_means
}

######
# w/o extra terms

# weights outside
neigh_coefs_oncont2 <- function(weights_df, H, G, neighinfo, A, a){
  neighX = neighinfo$neighX
  colnames(neighX) <- paste("neighX", colnames(neighX), sep = "_")
  
  if (!is.null(neighinfo) && !is.null(x1)){
    group_df <- as.data.frame(cbind(weights_df, H, G, A, neighX))
    colnames(group_df)[1] <- 'weights'
    
    fits_df <- as.data.frame(group_df[group_df$A == a,])
    fits <- lmList(as.formula(paste("H ~ ", paste(colnames(fits_df)[-1:-4], 
                                                  collapse= "+ "), "| G")),
                   weights = weights, data = fits_df)
    cond_coef <- as.matrix(coef(fits))
    
    overall_fits <- lm(as.formula(paste("H ~ ", paste(colnames(fits_df)[-1:-4],
                                                      collapse= "+ "))),
                       weights = weights, data = fits_df)
    overall_coef <- as.matrix(coef(overall_fits))
    coef <- list(cond_coef, overall_coef, colnames(coef(fits)))
    
  }else{
    group_df <- as.data.frame(cbind(weights_df, H, G, A))
    colnames(group_df) <- c('weights', 'H', 'G',' A')
    fits_df <- as.data.frame(group_df[group_df$A == a,])
    
    fits <- lmList(H ~ 1 | G, weights = weights, data=fits_df)
    cond_coef <- coef(fits)[,1]
    
    overall_fits <- lm(H ~ 1, weights = weights, data=fits_df)
    overall_coef <- coef(overall_fits)[[1]]
    
    coef <- list(cond_coef, overall_coef, colnames(coef(fits)))
  }
  coef
}

# weights inside 
neigh_coefs_oncont3 <- function(ind_est_df, H, G, neighinfo, A, a){
  neighX = neighinfo$neighX
  colnames(neighX) <- paste("neighX", colnames(neighX), sep = "_")
  
  if (!is.null(neighinfo) && !is.null(x1)){
    group_df <- as.data.frame(cbind(ind_est_df, H, G, A, neighX))
    colnames(group_df)[1] <- 'ind_est'
    
    fits_df <- as.data.frame(group_df[group_df$A == a,])
    fits <- lmList(as.formula(paste("ind_est ~ ", paste(colnames(fits_df)[-1:-4], 
                                                        collapse= "+ "), "| G")), data = fits_df)
    cond_coef <- as.matrix(coef(fits))
    
    overall_fits <- lm(as.formula(paste("ind_est ~ ", paste(colnames(fits_df)[-1:-4],
                                                            collapse= "+ "))), data = fits_df)
    overall_coef <- as.matrix(coef(overall_fits))
    coef <- list(cond_coef, overall_coef, colnames(coef(fits)))
    
  }else{
    group_df <- as.data.frame(cbind(ind_est_df, H, G, A))
    colnames(group_df) <- c('ind_est', 'H', 'G',' A')
    fits_df <- as.data.frame(group_df[group_df$A == a,])
    
    fits <- lmList(ind_est ~ 1 | G, data=fits_df)
    cond_coef <- coef(fits)[,1]
    
    overall_fits <- lm(ind_est ~ 1, data=fits_df)
    overall_coef <- coef(overall_fits)[[1]]
    
    coef <- list(cond_coef, overall_coef, colnames(coef(fits)))
  }
  coef
}

neigh_means_oncont2 <- function(overall_coef, X_type, x1, N){
  x1_num <- as.numeric(x1[X_type == "N"])
  lenn <- length(x1_num)
  cond_coefs <- matrix(rep(overall_coef, N), nrow = N, byrow = TRUE)
  # TODO: maybe N (number of groups) should be an input parameter
  
  if (lenn > 0){
    cond_group_means <- cond_coefs[,1] + as.matrix(cond_coefs[,-1])%*%as.matrix(x1_num)
  }else{
    cond_group_means <- cond_coefs
  }
  cond_group_means
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

