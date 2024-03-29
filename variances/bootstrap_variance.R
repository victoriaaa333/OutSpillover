# Assume we already have a dataframe with treatments (A), Y
# We would like to resample the clusters, and calculate its variance
#' bootstrap estimates of the population average potential outcomes for all
#' values of alpha will be returned
#' 

# TODO: the function used for point estimates might be changed

# GetBootSample <- function(dta) {
#   
#   num_clus <- max(dta$neigh)
#   boot_clusters <- sample(1 : num_clus, num_clus, replace = TRUE)
#   
#   # Binding data without accidentally merging repeated clusters.
#   boot_dta <- NULL
#   for (nn in 1 : num_clus) {
#     D <- subset(dta, neigh == boot_clusters[nn])
#     D$neigh <- nn
#     boot_dta <- rbind(boot_dta, D)
#   }
#   
#   return(list(boot_dta = boot_dta, chosen_clusters = boot_clusters))
# }

GetBootSample <- function(dta) {
  
  num_clus <- max(dta$G)
  boot_clusters <- sample(1 : num_clus, num_clus, replace = TRUE)
  
  # Binding data without accidentally merging repeated clusters.
  boot_dta <- NULL
  for (nn in 1 : num_clus) {
    D <- subset(dta, G == boot_clusters[nn])
    D$neigh <- nn
    boot_dta <- rbind(boot_dta, D)
  }
  
  return(list(boot_dta = boot_dta, chosen_clusters = boot_clusters))
}

#boot_df = GetBootSample(df)$boot_dta

# 1. boot var for influencer effect
BootVar <- function(dta, numerator_alpha, denominator_alphas, P,
                    boot_variable = "H", 
                    X_variable = "X", x0 = NULL,
                    B = 100, verbose = FALSE, return_everything = FALSE) {
  
  n_neigh <- max(dta$G)
  
  chosen_clusters <- array(NA, dim = c(n_neigh, B))
  dimnames(chosen_clusters) <- list(neigh = 1 : n_neigh, sample = 1 : B)
  
  ygroup <- array(NA, dim = c(n_neigh, 2, 1, B)) # length(alpha) = 1
  dimnames(ygroup) <- list(neigh = 1 : n_neigh, po = c('y0', 'y1'),
                           alpha = list(c(numerator_alpha, denominator_alphas)), sample = 1 : B)
  
  boots <- array(NA, dim = c(2, 1, B)) # length(alpha) = 1
  dimnames(boots) <- dimnames(ygroup)[- 1]
  
  
  for (bb in 1 : B) {
    
    if (verbose) {
      if (bb %% 100 == 0) {
        print(paste0('bootstrap sample ', bb))
      }
    }
    
    boot_dta <- GetBootSample(dta)
    chosen_clusters[, bb] <- boot_dta$chosen_clusters
    
    boot_dta <- boot_dta$boot_dta
    # neigh_ind <- lapply(1 : max(boot_dta$neigh),
    #                     function(nn) which(boot_dta$neigh == nn))
    # 
    allocations = list(c(numerator_alpha, denominator_alphas))
    w.matrix = wght_matrix(plain_integrand, allocations, boot_dta$neigh, boot_dta$A, P = P)
    names(boot_dta)[names(boot_dta) == boot_variable] <- 'boot_variable'
    
    # ygroup_boot <- ipw_point_estimates(boot_dta$boot_variable, boot_dta$neigh, 
    #                                      boot_dta$A, w.matrix,boot_dta[,X_variable], x0)$outcomes$groups
    ygroup_boot <- ipw_point_estimates_mixed_test4(boot_dta$boot_variable, boot_dta$neigh, 
                                                   boot_dta$A, w.matrix, boot_dta[,X_variable], x0)$outcomes$groups
    if (length(dim(ygroup_boot)) == 2){
        dim(ygroup_boot) <- c(dim(ygroup_boot)[1],dim(ygroup_boot)[2], 1)
    }
    ygroup[, , , bb] <- ygroup_boot
    boots[, , bb] <- apply(ygroup_boot, c(2, 3), mean, na.rm = TRUE) # take averages across all clusters
  }
    

if (return_everything) {
  return(list(boots = boots, ygroup = ygroup,
              chosen_clusters = chosen_clusters))
}

return(boots)
}

