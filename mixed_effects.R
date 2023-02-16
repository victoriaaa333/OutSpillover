column_multiply <- function(mat1, mat2){
  if(prod(dim(mat1) == dim(mat2)) == 0) stop("dimension are not the same")
  mat = c()
  mat_colnames = c()
  for (j in 1:dim(mat1)[2]){
    mat = cbind(mat, apply(mat2, 2, function(x) x*mat1[,j]))
    mat_colnames<- c(mat_colnames, paste(colnames(mat1)[j],
                                         colnames(mat2), sep = "_w_"))
  }
  list(mat,mat_colnames)
}

mixed_coef <- function(weights_df, H, G, X, X_type, x0, neighinfo, A, a){
  if (is.null(neighinfo)) stop("neighinfo cannot be NULL, 
                               refer to the function group_coef")

  #### information needed for conditional H ####
  trt_cond <- which(A == a)
  neighX = neighinfo$neighX
  colnames(neighX) <- paste("neighX", colnames(neighX), sep = "_")
  
  #### select out categorical values for conditional group avg ####
  if (sum(X_type == "C") > 0){
    X_cat <- as.matrix(X[, X_type == "C"])
    x0_cat <- x0[X_type == "C"]
    ind_cond <- which(apply(X_cat, 1, function(x) prod(x == x0_cat)) == 1)
  }else{
    ind_cond <- 1:length(neigh2_treated)}
  
  #### calculate the interactions ####
  # if there's numerical in regression for group avg
  if (sum(X_type == "N") > 0){
    X_num <- apply(as.matrix(X[, X_type == "N"]), 2, as.numeric)
    X_inter <- column_multiply(X_num, neighX)
    X_neighX <- X_inter[[1]]
    group_df <- as.data.frame(cbind(weights_df, H, G, X_num, neighX, X_neighX))
    
    num_names <- colnames(X)[X_type == "N"]
    neigh_names <- colnames(neighX)
    inter_names <- X_inter[[2]]
    colnames(group_df) <- c('w', 'H', 'G', num_names, neigh_names, inter_names)
  
  #### fit the regression ####  
    fits_df <- group_df[intersect(ind_cond, trt_cond), ]
    fits <- lmList(as.formula(paste("H ~ ", paste(colnames(fits_df)[-1:-3], collapse= "+ "), "| G")),
                   weights = w, data = fits_df)
    cond_coef <- as.matrix(coef(fits))
    overall_fits <- lm(as.formula(paste("H ~ ", paste(colnames(fits_df)[-1:-3], collapse= "+ "))),
                       weights = w, data = fits_df)
    # lm(as.formula(paste("H ~ ", paste(num_names, collapse= "+ "), "+", 
    #             paste(colnames(fits_df)[-1:-(3+length(num_names))], collapse= "+ "))), weights = w, data=fits_df)
    overall_coef <- as.matrix(coef(overall_fits))
    
    coef <- list(cond_coef, overall_coef)
  }else{
    # if only categorical value in the group avg, then there's no interaction, 
    # just regular regression on neigh factors
    group_df <- as.data.frame(cbind(weights_df, H, G, neighX))
    colnames(group_df)[1] <- 'w'
    
    #### fit the regression ####  
    fits_df <- group_df[intersect(ind_cond, trt_cond), ]
    fits_df$weighted_H <- fits_df$w * fits_df$H
    fits <- lmList(as.formula(paste("weighted_H ~ ", paste(colnames(fits_df)[c(-1:-3,
                      -dim(fits_df)[2])], collapse= "+ "), "| G")), data = fits_df)
    cond_coef <- as.matrix(coef(fits))
    overall_fits <- lm(as.formula(paste("weighted_H ~ ", paste(colnames(fits_df)[c(-1:-3,
                      -dim(fits_df)[2])], collapse= "+ "))), data = fits_df)
    overall_coef <- as.matrix(coef(overall_fits))
    coef <- list(cond_coef, overall_coef)
  }
  coef
}


mixed_means <- function(overall_coef, X_type, x0, x1){
  x0_num <- as.numeric(x0[X_type == "N"])
  x1_num <- as.numeric(x1[X_type == "N"])
  lenn <- length(x1_num)
  leng <- length(x0_num)
  cond_coefs <- matrix(rep(overall_coef, N), nrow = N, byrow = TRUE)
  
  if (leng > 0){
    cond_group_means <- cond_coefs[,1] + cond_coefs[,2:(1+leng)] %*% as.matrix(x0_num) +
      as.matrix(cond_coefs[,(2+leng):(1+leng+lenn)]) %*% as.matrix(x1_num) + 
      as.matrix(cond_coefs[,(2+leng+lenn):(3+leng+2*lenn)]) %*% c(outer(x1_num, x0_num))
  }else{
    cond_group_means <- cond_coefs[,1] + cond_coefs[,2:(1+lenn)] %*% as.matrix(x1_num)
  }
  cond_group_means
}

## TODO: integrate it into point estimates
