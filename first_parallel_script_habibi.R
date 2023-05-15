library(foreach)
library(doParallel)
registerDoParallel(7)
#library(doRedis)
#setProgress(TRUE)
library(matlib)

invert_matrix <- function(matrix)
{
  mat_inv <- solve(matrix)
  mat_inv
}

n_matrices = 100
mat_dim = 10

print("starting non-parallel test")
non_parallel_time <- system.time(
  {
    for(i in 1:n_matrices){
      random_matrix <- matrix(rnorm(mat_dim^2), nrow=mat_dim)
      inv(random_matrix)
    }
  }
)

print(c("non-parallel time", non_parallel_time))

print("starting parallel test")
parallel_time <- system.time(
  {
    foreach(i=1:n_matrices) %dopar% {
      print(i)
      random_matrix <- matrix(rnorm(mat_dim^2), nrow=mat_dim)
      inv(random_matrix)
      }
  }
)

print(c("parallel time", parallel_time))
