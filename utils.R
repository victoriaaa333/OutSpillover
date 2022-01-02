#-----------------------------------------------------------------------------#
# Calculate outcome mean per group per treatment level
# 
# @param Y vector of outcomes
# @param G vector of group assignments
# @param A vector of treatment assignments
# @param a value of treatment level, defaults to NA. NA is used for marginal
# estimates.
# @return matrix of group means
#-----------------------------------------------------------------------------#
group_means <- function(H, A, G, a = NA){
  
  N <- length(unique(G))
  HA <- cbind(H, A)
  
  vals <- by(HA, INDICES = G, function(x){
    n <- length(x[ , 1])
    
    if(is.na(a)){
      mean(x[ , 1],na.rm = TRUE)
    } else {
      mean(x[ , 1] * (x[ , 2] == a) * 1,na.rm = TRUE)
    }
  })
  
  out <- matrix(unlist(vals), nrow = N)
  
  return(out)
}

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

# Assume we have a dataset with the following columns:
# G vector of group assignments (clusters)
# A vector of assignments
# Y vector of outcomes
# H vector of average of h-neighborhood observed outcomes for the data point this row represents.

# We order this dataset by group.

#-----------------------------------------------------------------------------#
# Calculate h-neighborhood averages

# The vertices in graph needs to be same order as Y
h_neighborhood <- function(graph, Y, h){
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
    h_out = ifelse(length(h_neighbor) > 0, mean(Y[h_neighbor]),NA)
    h_vector = c(h_vector,h_out)
  }
  return(h_vector)
  # 0 if no h-neighborhood
}

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
