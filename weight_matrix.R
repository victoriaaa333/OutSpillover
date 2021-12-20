#-----------------------------------------------------------------------------#
#' Compute IPW weight
#' 
#' Calculates the IPW for a single group. Used by \code{\link{wght_matrix}} to
#' create a matrix of weights for each group and allocation scheme.
#' 
#' @param integrand function to pass to the denominator of IPW.
#' @param allocation the allocation ratio for which to compute the weight
#' @param ... other arguments passed to integrand.
#' @return scalar result of the integral
#' @export
#' @importFrom methods is
#'
#-----------------------------------------------------------------------------#

wght_calc <- function(integrand, 
                      numerator_alpha,
                      denominator_alphas,
                      ...)
{  
  ## Necessary pieces ##
  integrand         <- match.fun(integrand)
  integrand.formals <- names(formals(integrand))
  dots              <- list(...)
  dot.names         <- names(dots)
  A                 <- dots[[match.arg('A', dot.names)]]
  P                 <- dots[[match.arg('P', dot.names)]]
  
  PrA    <- integrand(A, denominator_alphas, P)
  ppp    <- prod(numerator_alpha^A * (1-numerator_alpha)^(1-A))
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
#' @param allocations a list of numerator alphas and denominator alphas.
#' @param X covariate matrix
#' @param A vector of treatment assignments
#' @param G vector of group assignments
#' @param runSilent if FALSE, errors are printed to console. Defaults to TRUE.
#' @inheritParams wght_calc
#' @return a length(unique(group)) X length(alphas) matrix of group weights
#' @export
#
#-----------------------------------------------------------------------------#

wght_matrix <- function(integrand, 
                        allocations, 
                        G, A, P,
                        X = NULL, 
                        runSilent = TRUE, 
                        ...)
{
  ## Gather necessary bits ##
  p  <- ifelse(is.null(ncol(X)),0,ncol(X)) 
  gg <- sort(unique(G))
  denominator_alphas = unlist(lapply(allocations, function (x) x[2]))
  
  ## Compute weight for each group and allocation level ##
  if(!runSilent) print('Calculating matrix of IP weights...') 
  
  w.list <- lapply(allocations, function(allocation){
    w <- by(cbind(X, A), INDICES = G, simplify = FALSE, 
            FUN = function(x) {
              print(allocation[1])
              wght_calc(
                        integrand  = integrand, 
                        numerator_alpha = allocation[1],
                        denominator_alphas = denominator_alphas, 
                        P = P,
                        A = x[, p+1], X = ifelse(p == 0, x[,p], x[, 1:p]))})
    as.numeric(w)
  }) 
  
  ## Reshape list into matrix ##
  w.matrix <- matrix(unlist(w.list, use.names = FALSE), 
                     ncol = length(allocations), 
                     byrow = FALSE,
                     dimnames = list(gg, allocations))
  
  return(w.matrix)
}


