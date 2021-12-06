#-----------------------------------------------------------------------------#
#' Integrand for the group-level propensity score (Type I randomization)
#' 
#' Computes the following function:
#' \deqn{\prod_{j=1}^{n} (\alpha)^{A_j}  (1 - \alpha)^{1 - A_j} 
#' @param A vector of binary treatments 
#' @param allocation The allocation strategy. Defaults to A so that is essentially 
#' ignored if allocation is not set to a value within (0, 1).
#' @return value of the integrand
#' @export
#' 
#-----------------------------------------------------------------------------#

integrand <- function(A, allocation)
{
  out <- prod(allocation^A * (1 - allocation)^(1-A))
  return(out)
}
