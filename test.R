#############################################
## 12/15/21 Yihan Bao:
## 1. generate function "group_vector" to align with actual group assignments
## 2. treatment vector A is sample from denominator alpha
## 3. generate Y that is dependent on the number of treatments in its group
## 4. 2-stage (p/alpha) in denominator

source("weight_matrix.R")
source("ipw_point_estimates.R")
source("integrand.R")
source("utils.R")

group_vector <- function(g){
  number_G = length(groups(components(g)))
  G = rep(NA,number_G)
  for (i in 1:number_G) {
    G[unlist(groups(components(g))[[i]])] = i # vertex indices in group i
  }
  return(G)
}

# test dataset
library(igraph)
g = sample_gnm(n = 20, m = 15) # 20 vertices, 15 edges
plot(g, vertex = 6)
G = group_vector(g) 
denominator_alphas = c(0.4,0.7)
P = c(0.7,0.3)

A = sample(c(0,1),20,replace = TRUE, prob = c(0.5,0.5))

df <- cbind.data.frame(A,G)
df.trt_mean <- aggregate(A, list(G), FUN=mean)
colnames(df.trt_mean) <- c("G", "group_mean")

df <- merge(x = df, y = df.trt_mean, by = "G", all.x = TRUE)
Y = apply(cbind(df$A,df$group_mean),1, function(x) rnorm(1, mean = x[1] + 2*x[2]))
H = h_neighborhood(g,Y,1) 
df$Y = Y
df$H = H
df

allocations = list(c(0.4,denominator_alphas[1]),c(0.4,denominator_alphas[2]))
w.matrix = wght_matrix(integrand,allocations, G, A, P)
ipw_point_estimates(H,G,A,w.matrix)


