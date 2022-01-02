# Assume we already have a dataframe with treatments (A), Y
# We would like to resample the clusters, and calculate its variance
#' bootstrap estimates of the population average potential outcomes for all
#' values of alpha will be returned
#' 

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

BootVar <- function(dta, numerator_alpha, alpha, B = 1000, 
                    verbose = FALSE, return_everything = FALSE) {
  
  n_neigh <- max(dta$G)
  
  chosen_clusters <- array(NA, dim = c(n_neigh, B))
  dimnames(chosen_clusters) <- list(neigh = 1 : n_neigh, sample = 1 : B)
  
  ygroup <- array(NA, dim = c(n_neigh, 2, length(alpha), B))
  dimnames(ygroup) <- list(neigh = 1 : n_neigh, po = c('y0', 'y1'),
                           alpha = alpha, sample = 1 : B)
  
  boots <- array(NA, dim = c(2, length(alpha), B))
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
    neigh_ind <- lapply(1 : max(boot_dta$neigh),
                        function(nn) which(boot_dta$neigh == nn))
    
    allocations = list(c(numerator_alpha,alpha))
    w.matrix = wght_matrix(integrand, allocations, boot_dta$neigh, boot_dta$A, P)
    
    ygroup_boot <- ipw_point_estimates(boot_dta$H, boot_dta$neigh, 
                                         boot_dta$A, w.matrix)$outcomes$groups
    if (length(dim(ygroup_boot)) == 2){
        dim(ygroup_boot) <- c(dim(ygroup_boot)[1],dim(ygroup_boot)[2],length(alpha))
    }
    ygroup[, , , bb] <- ygroup_boot
    boots[, , bb] <- apply(ygroup_boot, c(2, 3), mean)
    #array(NA, dim = c(n_neigh, 2, length(alpha)))
  }
    

if (return_everything) {
  return(list(boots = boots, ygroup = ygroup,
              chosen_clusters = chosen_clusters))
}

return(boots)
}

