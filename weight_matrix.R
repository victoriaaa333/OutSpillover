#-----------------------------------------------------------------------------#
#' Compute IPW weight
#' 
#' Calculates the IPW for a single group. Used by \code{\link{wght_matrix}} to
#' create a matrix of weights for each group and allocation scheme.
#' 
#' If \code{allocation} is an argument in the integrand function and
#' \code{integrate_allocation == TRUE}, then the weight is calcuated as:
#' 
#' \deqn{\frac{1}{Pr(A|X)}}{1 / integrate(integrand)}
#' 
#' Otherwise, the weight is computed by:
#' \deqn{\frac{\prod_{j=1}^n \alpha^A_j (1 - \alpha)^(1-
#' A_j)}{Pr(A|X)}}{prod(allocation^A(1 - allocation)^A) / integrate(integrand)}
#' 
#' @param parameters vector of parameter values
#' @param integrand function to pass to the argument 'f' of \code{\link{integrate}}.
#' @param allocation the allocation ratio for which to compute the weight
#' @param integrate_allocation Indicator of whether the integrand function uses 
#' the allocation parameter. Defaults to TRUE.
#' @param ... other arguments passed to integrand.
#' @return scalar result of the integral
#' @export
#' @importFrom methods is
#'
#-----------------------------------------------------------------------------#

wght_calc <- function(integrand, 
                      allocation,
                      ...)
{  
  ## Necessary pieces ##
  integrand         <- match.fun(integrand)
  integrand.formals <- names(formals(integrand))
  dots              <- list(...)
  dot.names         <- names(dots)
  A                 <- dots[[match.arg('A', dot.names)]]
  
  PrA    <- integrand(A, allocation)
  ppp    <- prod(allocation^A * (1-allocation)^(1-A))
  weight <- ppp/PrA
  weight
}

#-----------------------------------------------------------------------------#
# Create a matrix of group IP weights 
#' 
#' Creates a number of groups by number of allocation schemes matrix of group weights.
#' Allocation schemes are selected by the user.
#' 
#' Groups should be numbered 1, ..., N
#' 
#' @param allocations coverage levels in [0, 1]. Can be vector.
#' @param X covariate matrix
#' @param A vector of treatment assignments
#' @param G vector of group assignments
#' @param parameters vector of parameters passed to \code{integrand}
#' @param runSilent if FALSE, errors are printed to console. Defaults to TRUE.
#' @inheritParams wght_calc
#' @return a length(unique(group)) X length(alphas) matrix of group weights
#' @export
#
#-----------------------------------------------------------------------------#

wght_matrix <- function(integrand, 
                        allocations, 
                        G, A, 
                        X = NULL, 
                        runSilent = TRUE, 
                        ...)
{
  ## Gather necessary bits ##
  p  <- ifelse(is.null(ncol(X)),0,ncol(X)) 
  aa <- sort(allocations) # Make sure alphas are sorted
  gg <- sort(unique(G))
  
  ## Compute weight for each group and allocation level ##
  if(!runSilent) print('Calculating matrix of IP weights...') 
  
  w.list <- lapply(aa, function(allocation){
    w <- by(cbind(X, A), INDICES = G, simplify = FALSE, 
            FUN = function(x) {
              wght_calc(
                        integrand  = integrand, 
                        allocation = allocation, 
                        A = x[, p+1], X = x[, 1:p])})
    as.numeric(w)
  }) 
  
  ## Reshape list into matrix ##
  w.matrix <- matrix(unlist(w.list, use.names = FALSE), 
                     ncol = length(allocations), 
                     byrow = FALSE,
                     dimnames = list(gg, aa))
  
  return(w.matrix)
}


