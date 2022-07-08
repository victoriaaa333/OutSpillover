#-----------------------------------------------------------------------------#
# Calculate outcome mean per group per treatment level
# 
# @param Y vector of outcomes
# @param G vector of group assignments
# @param A vector of treatment assignments
# @param a value of treatment level, defaults to NA. NA is used for marginal
# estimates.
# @param X matrix of covariates
# @param x0 vector of specific covariates conditioned on 
# @return matrix of group means
#-----------------------------------------------------------------------------#

group_means <- function(H, A, G, X = NULL, x0 = NULL, a = NA){
  #if (dim(X)[2] != length(x0)) stop('X and x0 should be the same dimension')
  N <- length(unique(G))
  HAX <- cbind(H, A, X)
  
  vals <- by(HAX, INDICES = G, function(x){
    #n <- length(x[ , 1])
    if(is.null(X)){
      if(is.na(a)){
        mean(x[ , 1], na.rm = TRUE)
      } else {
        #trt_cond = apply(x[, 2], 1, function(x) ifelse(x == a, 1, NA))
        mean(x[ , 1] * ifelse(x[,2] == a, 1, NA) * 1, na.rm = TRUE)
      }
    }else{
      if(is.null(x0)){
        influencer_cond = rep(1, dim(x)[1])
      }else{
        influencer_cond = apply(as.matrix(x[,3:dim(HAX)[2]]), 1, 
                                function(x) ifelse(prod(x == x0), 1, NA))}
      # indicator of whether this influencer is conditional on x0, 1 or NA (remove this influencer when taking ave)
      if(is.na(a)){
        mean(as.numeric(x[ , 1]) * influencer_cond, na.rm = TRUE)
      } else {
        mean(as.numeric(x[ , 1]) * ifelse(x[, 2] == a, 1, NA) * influencer_cond, na.rm = TRUE)
      }
    }
  })
  
  out <- matrix(unlist(vals), nrow = N)
  
  return(out)
}

# only condition on neighborhoods
# group_means <- function(H, A, G, X = NULL, x0 = NULL, a = NA){
#   #if (dim(X)[2] != length(x0)) stop('X and x0 should be the same dimension')
#   N <- length(unique(G))
#   HAX <- cbind(H, A, X)
#   
#   vals <- by(HAX, INDICES = G, function(x){
#     #n <- length(x[ , 1])
#     if(is.null(X)){
#       if(is.na(a)){
#         mean(x[ , 1],na.rm = TRUE)
#       } else {
#         mean(x[ , 1] * (x[ , 2] == a) * 1,na.rm = TRUE)
#       }
#     }else{
#       if(is.na(a)){
#         mean(x[ , 1] * (x[ , 3:dim(HAX)[2]] == x0), na.rm = TRUE)
#       } else {
#         mean(x[ , 1] * (x[ , 3:dim(HAX)[2]] == x0) * (x[ , 2] == a) * 1,na.rm = TRUE)
#       }
#     }
# # TODO: change (x[ , 3:dim(HAX)[2]] == x0) to indicator
#   })
#   
#   out <- matrix(unlist(vals), nrow = N)
#   
#   return(out)
# }

# no condition
# group_means <- function(H, A, G, a = NA){
#   
#   N <- length(unique(G))
#   inds <- 1:length(A)
#   HA <- cbind(H, A)
#   
#   vals <- rep(NA, length(inds))
#   for (i in 1:length(inds)) {
#     if(is.na(a)){
#       vals[i] = mean(HA[, 1][G == G[i]], na.rm = TRUE)
#     }
#     else {
#       vals[i] = mean(HA[ , 1][G == G[i]] * 
#                        (HA[ , 2] == a)[G == G[i]] * 1,na.rm = TRUE)
#     }
#   }
#   
#   out <- matrix(unlist(vals), nrow = N)
#   
#   return(out)
# }

group_vector <- function(g){
  number_G = length(groups(components(g)))
  G = rep(NA,number_G)
  for (i in 1:number_G) {
    G[unlist(groups(components(g))[[i]])] = i # vertex indices in group i
  }
  return(G)
}

# Assume we have a dataset with the following columns:
# G vector of group assignments (clusters)
# A vector of assignments
# Y vector of outcomes
# H vector of average of h-neighborhood observed outcomes for the data point this row represents.

# We order this dataset by group.

#-----------------------------------------------------------------------------#
# Calculate h-neighborhood averages (x1 be categorical)

# The vertices in graph needs to be same order as Y
# X and x1 is the covariate matrix and the condition for neighbors
h_neighborhood <- function(graph, Y, h, X = NULL, x1 = NULL){
  h_vector = c()
  vg = V(graph)
  num_vertices = length(vg)
  if( h < 0 ) stop('h has to be non-negative')
  
  for (j in 1:num_vertices) {
    if (h == 0){
      h_neighbor = ego(graph,h,vg[j])[[1]]
    }else{
      h_neighbor = setdiff(ego(graph,h,vg[j])[[1]],
                           ego(graph,h-1,vg[j])[[1]]) 
    }
    if (!is.null(X) && !is.null(x1)){
      neigh_cond = apply(as.matrix(X[h_neighbor, ]), 1, function(x) ifelse(prod(x == x1), 1, NA))
      h_out = ifelse(length(h_neighbor) > 0, mean(Y[h_neighbor] * neigh_cond, na.rm = TRUE), NA)
    }else{
      h_out = ifelse(length(h_neighbor) > 0, mean(Y[h_neighbor]), NA)
    }
    h_vector = c(h_vector,h_out)
  }
  return(h_vector)
  # 0 if no h-neighborhood
}

# no-condition version
# h_neighborhood <- function(graph, Y, h){
#   h_vector = c()
#   vg = V(graph)
#   num_vertices = length(vg)
#   if( h < 0 ) stop('h has to be non-negative')
#   
#   for (j in 1:num_vertices) {
#     if (h == 0){
#       h_neighbor = ego(graph,h,vg[j])[[1]]
#     }else{
#       h_neighbor = setdiff(ego(graph,h,vg[j])[[1]],
#                            ego(graph,h-1,vg[j])[[1]]) 
#     }
#     h_out = ifelse(length(h_neighbor) > 0, mean(Y[h_neighbor]),NA)
#     h_vector = c(h_vector,h_out)
#   }
#   return(h_vector)
#   # 0 if no h-neighborhood
# }

# The vertices in graph needs to be same order as Y
h_counts <- function(graph, h){
  h_vector = c()
  vg = V(graph)
  num_vertices = length(vg)
  if( h < 0 ) stop('h has to be non-negative')
  
  for (j in 1:num_vertices) {
    if (h == 0){
      h_neighbor = ego(graph,h,vg[j])[[1]]
    }else{
      h_neighbor = setdiff(ego(graph,h,vg[j])[[1]],
                           ego(graph,h-1,vg[j])[[1]]) 
    }
    h_out = length(h_neighbor)
    h_vector = c(h_vector,h_out)
  }
  return(h_vector)
}

# TODO: check this function
# the number of treated units for each node's h-order neighborhood
h_neighsum <- function(graph, A, h, X = NULL, X_type = NULL, x1 = NULL){
  h_vector = c()
  vg = V(graph)
  num_vertices = length(vg)
  if( h < 0 ) stop('h has to be non-negative')
  
  if (sum(X_type == "C") > 0){
    X_cat <- as.matrix(X[, X_type == "C"])
    x1_cat <- x1[X_type == "C"]
    ind_cond <- apply(X_cat, 1, function(x) ifelse(prod(x == x1_cat), 1, NA))
  }else{
    ind_cond <- rep(1,length(A))}
  
  for (j in 1:num_vertices) {
    if (h == 0){
      h_neighbor = ego(graph,h,vg[j])[[1]]
    }else{
      h_neighbor = setdiff(ego(graph,h,vg[j])[[1]],
                           ego(graph,h-1,vg[j])[[1]]) 
    }
    
    if (!is.null(X) & !is.null(x1)){
      h_out = ifelse(length(h_neighbor) > 0, sum(A[h_neighbor] * ind_cond[h_neighbor], na.rm = TRUE), NA)
      # neigh_cond = apply(as.matrix(X[h_neighbor, ]), 1, 
      #                    function(x) ifelse(prod(x == x1), 1, NA)) # check if X is the same as x1
      # h_out = ifelse(length(h_neighbor) > 0, sum(A[h_neighbor] * neigh_cond, na.rm = TRUE), NA)
    }else{
      h_out = ifelse(length(h_neighbor) > 0, sum(A[h_neighbor]), NA)}
      #h_out = sum(A[h_neighbor])
      h_vector = c(h_vector,h_out)
    }
  return(h_vector)
}

# the mean of covariates for each node's h-order neighborhood
h_neighcov <- function(graph, A, h, X, X_type, x1){
  h_vector = c()
  vg = V(graph)
  num_vertices = length(vg)
  if( h < 0 ) stop('h has to be non-negative')
  
  X_cat <- as.matrix(X[, X_type == "C"])
  X_num <- apply(as.matrix(X[, X_type == "N"]), 2, as.numeric)
  x1_cat <- x1[X_type == "C"]
  x1_num <- as.numeric(x1[X_type == "N"])
  
  if (sum(X_type == "C") > 0){
    ind_cond <- apply(X_cat, 1, function(x) ifelse(prod(x == x1_cat), 1, NA))
  }else{
    ind_cond <- 1:dim(X)[1]}
  
  for (j in 1:num_vertices) {
    if (h == 0){
      h_neighbor = ego(graph,h,vg[j])[[1]]
    }else{
      h_neighbor = setdiff(ego(graph,h,vg[j])[[1]],
                           ego(graph,h-1,vg[j])[[1]])
      }
    if (length(h_neighbor) > 0){
      if (length(x1_num) == 1){
       h_out = mean(X_num[h_neighbor, ] * ind_cond[h_neighbor], na.rm = TRUE)
      }else{
       h_out = colMeans(X_num[h_neighbor, ] * ind_cond[h_neighbor], na.rm = TRUE)      
      }
    }else{
      h_out = NA
    }
    h_vector = rbind(h_vector,h_out)
  }
  return(h_vector)
}

true_ind_effect <- function(z_jk, n_k, alphas, a, b, P = 1){
  if( length(alphas) != length(P) ) stop('P is not the same length as alphas')
  prob_peralpha = as.vector(sapply(alphas, function(alpha)
          sum((a * z_jk + b * 0:(n_k - 1)) *
            dbinom(0:(n_k - 1), size = (n_k - 1), prob = alpha)))) 
  sum(prob_peralpha * P)
}

true_cluster_effect <- function(g_info, alphas, a, b, P = 1){
  # number of units in that cluster
  count_g = g_info[1] 
  # whether this cluster has no units with h-order neighbors,
  # if so, set the cluster effect to NA.
  is.h = g_info[2] 
  c(true_ind_effect(0, count_g, alphas, a, b, P), true_ind_effect(1, count_g, alphas, a, b, P)) * is.h
}

true_population_effect <- function(G_info, alphas, a, b, P = 1){
  # G_info contains information about whether this cluster is NA or not
  rowMeans(apply(G_info, 1, true_cluster_effect, alphas, a, b, P), na.rm = TRUE)
}

get_args <- function(FUN, args_list = NULL, ...){
  dots <- append(args_list, list(...))
  arg_names <- names(formals(match.fun(FUN)))
  
  args <- dots[arg_names]
  args[sapply(args, is.null)] <- NULL
  
  return(args)
}



####### helper function for test10
# 1. the sum of the covariates in the TREATED units for each node's h-order neighborhood
# x1 can be set to categorical value
cov_neighsum <- function(graph, A, h, X, x1 = NULL){
  h_vector = c()
  vg = V(graph)
  num_vertices = length(vg)
  if( h < 0 ) stop('h has to be non-negative')
  if(is.null(X)) stop('X should not be NULL')
  for (j in 1:num_vertices) {
    if (h == 0){
      h_neighbor = ego(graph,h,vg[j])[[1]]
    }else{
      h_neighbor = setdiff(ego(graph,h,vg[j])[[1]],
                           ego(graph,h-1,vg[j])[[1]]) 
    }
    
    if (!is.null(x1)){
      neigh_cond = apply(as.matrix(X[h_neighbor, ]), 1, 
                         function(x) ifelse(prod(x == x1), 1, NA)) # check if X is the same as x1
      h_out = ifelse(length(h_neighbor) > 0, sum(A[h_neighbor] *X[h_neighbor] * neigh_cond, na.rm = TRUE), NA)
    }else{
      h_out = ifelse(length(h_neighbor) > 0, sum(A[h_neighbor] *X[h_neighbor]), NA)}
    #h_out = sum(A[h_neighbor])
    h_vector = c(h_vector,h_out)
  }
  return(h_vector)
}

# 2. H with ONE continuous variable (run a regression)
# the number of treated units for each node's h-order neighborhood
# unit
h_neighborhood_cont <- function(graph, Y, h, X = NULL, x1 = NULL){
  h_vector = c()
  vg = V(graph)
  num_vertices = length(vg)
  if( h < 0 ) stop('h has to be non-negative')
  
  for (j in 1:num_vertices) {
    if (h == 0){
      h_neighbor = ego(graph,h,vg[j])[[1]]
    }else{
      h_neighbor = setdiff(ego(graph,h,vg[j])[[1]],
                           ego(graph,h-1,vg[j])[[1]]) 
    }
    if (!is.null(X) && !is.null(x1)){
      #neigh_cond = apply(as.matrix(X[h_neighbor, ]), 1, function(x) ifelse(prod(x == x1), 1, NA))
      #h_out = ifelse(length(h_neighbor) > 0, mean(Y[h_neighbor] * neigh_cond, na.rm = TRUE), NA)
      neigh_df <- as.data.frame(cbind(Y[h_neighbor], X[h_neighbor, ]))
      colnames(neigh_df) <- c('Y', 'X')
      fits <- lm(Y ~ X, data=neigh_df)
      h_out <- as.numeric(coef(fits)[1] + coef(fits)[2]*x1)
    }else{
      h_out = ifelse(length(h_neighbor) > 0, mean(Y[h_neighbor]), NA)
    }
    h_vector = c(h_vector,h_out)
  }
  return(h_vector)
  # 0 if no h-neighborhood
}

# 3. new way to calculate H with ONE continuous variable (run a regression of H on avg X)
# grouped
h_neighborhood_cont2 <- function(graph, Y, G, h, X, A, a, x1 = NULL){
  h_vector = c()
  vg = V(graph)
  num_vertices = length(vg)
  neigh_df = data.frame()

  if( h < 0 ) stop('h has to be non-negative')

  for (j in 1:num_vertices) {
    if (h == 0){
      h_neighbor = ego(graph,h,vg[j])[[1]]
    }else{
      h_neighbor = setdiff(ego(graph,h,vg[j])[[1]],
                           ego(graph,h-1,vg[j])[[1]])
    }
    neigh_df <- rbind(neigh_df,
                      cbind(mean(Y[h_neighbor]), mean(X[h_neighbor, ]), G[j], A[j]))
  }
  colnames(neigh_df) <- c('H', 'neigh_X', 'G', 'A')
  neigh_df <- as.data.frame(neigh_df[neigh_df$A == a,])
  fits <- lmList(H ~ neigh_X | G, data = neigh_df)
  h_out <- as.numeric(coef(fits)[,1] + coef(fits)[,2]*x1)
  h_vector <- rep(h_out, each = 100) # TODO: hardcode here, change this
  return(h_vector)
  # 0 if no h-neighborhood
}

# 3. the sum of TREATED neighbors for each node's h-order neighborhood
# x1 can be set to categorical value
h_neighofneigh_withcov <- function(graph, A, h, X, X_type, x1, h_neigh){
  h_vector = c()
  vg = V(graph)
  num_vertices = length(vg)
  if( h < 0 ) stop('h has to be non-negative')
  if(is.null(X)) stop('X should not be NULL')
  
  X_cat <- as.matrix(X[, X_type == "C"])
  X_num <- apply(as.matrix(X[, X_type == "N"]), 2, as.numeric)
  x1_cat <- x1[X_type == "C"]
  x1_num <- as.numeric(x1[X_type == "N"])
  
  if (sum(X_type == "C") > 0){
    ind_cond <- apply(X_cat, 1, function(x) ifelse(prod(x == x1_cat), 1, NA))
  }else{
    ind_cond <- rep(1, dim(X)[1])}
  
  for (j in 1:num_vertices) {
    if (h == 0){
      h_neighbor = ego(graph,h,vg[j])[[1]]
    }else{
      h_neighbor = setdiff(ego(graph,h,vg[j])[[1]],
                           ego(graph,h-1,vg[j])[[1]]) 
    }
    
    cov_neighbor = h_neighbor * ind_cond[h_neighbor]
    cov_neighbor = cov_neighbor[!is.na(cov_neighbor)]
    h_out = t((h_neigh[cov_neighbor] - A[j]))  %*% 
      (X_num[cov_neighbor,])/length(cov_neighbor)
    
    #0.1076
    # for (k in h_neighbor) {
    #   h_neigh2 = setdiff(ego(graph,h,k)[[1]],
    #                      ego(graph,h-1,k)[[1]]) 
    #   neigh2_A = c(neigh2_A, A[h_neigh2])
    # }
    
    #h_out = sum(X_num[h_neighbor,]*(h_neigh[h_neighbor] - A[j]))/length(h_neighbor)
    # h_out = t((h_neigh[h_neighbor] - A[j]) * ind_cond[h_neighbor]) %*% 
    #   (X_num[h_neighbor,])/sum(ind_cond[h_neighbor])
    h_vector = rbind(h_vector,h_out)
  }
  return(h_vector)
}

h_neighind <- function(graph, A, h, X = NULL, X_type = NULL, x1 = NULL){
  h_vector = c()
  vg = V(graph)
  num_vertices = length(vg)
  if( h < 0 ) stop('h has to be non-negative')
  
  if (sum(X_type == "C") > 0){
    X_cat <- as.matrix(X[, X_type == "C"])
    x1_cat <- x1[X_type == "C"]
    ind_cond <- apply(X_cat, 1, function(x) ifelse(prod(x == x1_cat), 1, NA))
  }else{
    ind_cond <- 1:length(A)}
  
  for (j in 1:num_vertices) {
    if (h == 0){
      h_neighbor = ego(graph,h,vg[j])[[1]]
    }else{
      h_neighbor = setdiff(ego(graph,h,vg[j])[[1]],
                           ego(graph,h-1,vg[j])[[1]]) 
    }
    
    if (!is.null(X) & !is.null(x1)){
      h_out = h_neighbor[ind_cond[h_neighbor] == 1]
      h_out = h_out[!is.na(h_out)]
      # h_out = ifelse(length(h_neighbor) > 0, sum(A[h_neighbor] * ind_cond[h_neighbor], na.rm = TRUE), NA)
      # neigh_cond = apply(as.matrix(X[h_neighbor, ]), 1, 
      #                    function(x) ifelse(prod(x == x1), 1, NA)) # check if X is the same as x1
      # h_out = ifelse(length(h_neighbor) > 0, sum(A[h_neighbor] * neigh_cond, na.rm = TRUE), NA)
    }else{
      h_out = h_neighbor
      # h_out = ifelse(length(h_neighbor) > 0, sum(A[h_neighbor]), NA)}
      # h_out = sum(A[h_neighbor])
    }
    h_vector = rbind(h_vector,h_out)
  }
  return(h_vector)
}

h_neighofneigh <- function(graph, A, h, h_neigh, 
                           X = NULL, X_type = NULL, x1 = NULL){
  h_vector = c()
  vg = V(graph)
  num_vertices = length(vg)
  if( h < 0 ) stop('h has to be non-negative')
  
  if (!is.null(X_type) && sum(X_type == "C") > 0){
    X_cat <- as.matrix(X[, X_type == "C"])
    x1_cat <- x1[X_type == "C"]
    ind_cond <- apply(X_cat, 1, function(x) ifelse(prod(x == x1_cat), 1, NA))
  }else{
    ind_cond <- rep(1, length(A))}
  
  h_vector = c()
  for (j in 1:num_vertices) {
    if (h == 0){
      h_neighbor = ego(graph,h,vg[j])[[1]]
    }else{
      h_neighbor = setdiff(ego(graph,h,vg[j])[[1]],
                           ego(graph,h-1,vg[j])[[1]]) 
    }
    cov_neighbor = h_neighbor*ind_cond[h_neighbor]
    cov_neighbor = cov_neighbor[!is.na(cov_neighbor)]
    # for (k in h_neighbor) {
    #   h_neigh2 = setdiff(ego(graph,h,k)[[1]],
    #                      ego(graph,h-1,k)[[1]]) 
    #   neigh2_A = c(neigh2_A, A[h_neigh2])
    # }
    h_out = sum(h_neigh[cov_neighbor] - A[j])/length(cov_neighbor)
    h_vector = c(h_vector,h_out)
  }
  return(h_vector)
}

# # 4. new way to calculate H with ONE continuous variable (run a regression of H on avg X)
# # grouped
# h_neighborhood_cont2 <- function(graph, Y, h, X, G, x1 = NULL){
#   h_vector = c()
#   vg = V(graph)
#   num_vertices = length(vg)
#   neigh_df = data.frame()
#   
#   if( h < 0 ) stop('h has to be non-negative')
#   
#   for (j in 1:num_vertices) {
#     if (h == 0){
#       h_neighbor = ego(graph,h,vg[j])[[1]]
#     }else{
#       h_neighbor = setdiff(ego(graph,h,vg[j])[[1]],
#                            ego(graph,h-1,vg[j])[[1]])
#     }
#     neigh_df <- rbind(neigh_df,
#                       cbind(mean(Y[h_neighbor]), mean(X[h_neighbor, ]), G[j]))
#   }
#   colnames(neigh_df) <- c('H', 'neigh_X', 'G')
#   #neigh_df <- as.data.frame(neigh_df[neigh_df$A == a,]) 
#   fits <- lmList(H ~ neigh_X | G, data = neigh_df)
#   h_out <- as.numeric(coef(fits)[,1] + coef(fits)[,2]*x1)
#   h_vector <- rep(h_out, each = 100) # TODO: hardcode here, change this
#   return(h_vector)
#   # 0 if no h-neighborhood
# }

