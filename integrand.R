#-----------------------------------------------------------------------------#
#' Default integrand for the group-level propensity score
#' 
#' Computes the following function:
#' \deqn{\prod_{j=1}^{n} (r h_{j}(b))^{A_j}  (1 - r h_{j}(b))^{1 - A_j} 
#' f_b(b; \theta_b)}{ prod(r * plogis(X * fixef + b)^A * 
#' (1 - r * plogis(X * fixef+ b))^(1 - A)) * 
#' dnorm(sd = sqrt(ranef))} 
#' where \eqn{r} is the randomization scheme. \eqn{X} is the covariate(s) vectors. 
#' \eqn{fixef} is the vector of fixed effects. \eqn{b} is the random (group-level) effect.
#' \eqn{ranef} is the random effect variance. 
#' 
#' @param b vector argument of values necessary for \code{\link{integrate}}.
#' @param X n by length(fixed effects) matrix of covariates.
#' @param parameters vector of fixed effect (and random effect if applicable). 
#' Random effect should be last element in vector.
#' @param A vector of binary treatments 
#' @param allocation The allocation strategy. Defaults to A so that is essentially 
#' ignored if allocation is not set to a value within (0, 1).
#' @param randomization Randomization probability. Defaults to 1.
#' 
#' @return value of the integrand
#' @export
#' 
#-----------------------------------------------------------------------------#

integrand <- function(A, allocation)
{
  out <- prod(allocation^A * (1 - allocation)^(1-A))
  return(out)
}
