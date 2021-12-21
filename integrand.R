#-----------------------------------------------------------------------------#
#' Integrand for the group-level propensity score (Type I randomization)
#' 
#' Computes the following function:
#' \deqn{\sum_{i = 1}^{m} P_i \prod_{j=1}^{n} (\alpha_i)^{A_j}  (1 - \alpha_i)^{1 - A_j} 
#' @param A vector of binary treatments 
#' @param denominator_alphas The allocation strategies for denominator.
#' @param P vector of first-stage randomization, default is 1 (one-stage)
#' @return value of the integrand
#' @export
#' 
#-----------------------------------------------------------------------------#

integrand <- function(A, denominator_alphas, P = 1)
{ 
  if (length(P) != length(denominator_alphas) | sum(P) != 1) stop('P is wrong')
  prob_alphas = unlist(lapply(denominator_alphas, 
                       function(alpha) prod(alpha^A * (1-alpha)^(1-A))))
  out = sum(prob_alphas * P)
  return(out)
}
