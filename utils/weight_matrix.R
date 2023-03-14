#-----------------------------------------------------------------------------#
#' Compute IPW weight
#' 
#' Calculates the IPW for a single group. Used by \code{\link{wght_matrix}} to
#' create a matrix of weights for each group and allocation scheme.
#' 
#' @param integrand function to pass to the denominator of IPW.
#' @param numerator_alpha the allocation ratio for which to compute the weight
#' @param denominator_alphas the combinations of potential alphas
#' @param ... other arguments passed to integrand.
#' @return scalar result of the integral
#' @export
#' @importFrom methods is
#'
#-----------------------------------------------------------------------------#

wght_calc <- function(integrand, 
                      numerator_alpha,
                      denominator_alphas,
                      parameters = NULL,
                      randomizations = NULL,
                      ...)
{  
  ## Necessary pieces ##
  integrand         <- match.fun(integrand)
  integrand.formals <- names(formals(integrand))
  dots              <- list(...)
  dot.names         <- names(dots)
  A                 <- dots[[match.arg('A', dot.names)]]
  P                 <- dots[[match.arg('P', dot.names)]]
  X                 <- dots[[match.arg('X', dot.names)]]
  
  if (is.null(parameters)){
    PrA    <- integrand(A, denominator_alphas, P)
  }else{
    if(length(randomizations) != length(P)) stop("b should be equal length of P")
    
    PrA_eachb = c()
    
    # Integrate() arguments ##
    if(!'lower' %in% dot.names){
      dots$lower <- -Inf
    }
    if(!'upper' %in% dot.names){
      dots$upper <- Inf
    }
    
    for (r in randomizations) {
      int.args <- append(get_args(stats::integrate, dots),
                         list(f = integrand, parameters = parameters, 
                              randomization = r))
      args <- append(get_args(integrand, dots), int.args)
      
      ## Compute the integral ##
      # if any of the products within the integrand return Inf, then return NA
      # else return the result of integration
      
      f <- try(do.call(stats::integrate, args = args), silent = TRUE)
      PrA_eachb <- c(PrA_eachb, if(is(f, 'try-error')) NA else f$value)
    }
    
    PrA <- sum(PrA_eachb*P)
    
  }
  
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
                        parameters = NULL,
                        randomizations = NULL,
                        runSilent = TRUE, 
                        ...)
{
  ## Gather necessary bits ##
  X_col  <- ifelse(is.null(ncol(X)), 0, ncol(X)) 
  gg <- sort(unique(G))
  #denominator_alphas = unique(unlist(lapply(allocations, function (x) x[2])))
  denominator_alphas = unique(unlist(lapply(allocations, function (x) x[-1])))
  
  ## Compute weight for each group and allocation level ##
  if(!runSilent) print('Calculating matrix of IP weights...') 
  
  if(X_col == 0){
    w.list <- lapply(allocations, function(allocation){
      w <- by(cbind(X, A), INDICES = G, simplify = FALSE, 
              FUN = function(x) {
                wght_calc(
                  integrand  = integrand, 
                  numerator_alpha = allocation[1],
                  denominator_alphas = denominator_alphas, 
                  P = P,
                  A = x[, X_col+1], 
                  X = NULL, #x[,X_col]
                  parameters = parameters,
                  randomizations = randomizations)})
      as.numeric(w)})
  }else{
  w.list <- lapply(allocations, function(allocation){
    w <- by(cbind(X, A), INDICES = G, simplify = FALSE, 
            FUN = function(x) {
              wght_calc(
                        integrand  = integrand, 
                        numerator_alpha = allocation[1],
                        denominator_alphas = denominator_alphas, 
                        P = P,
                        A = x[, X_col+1], 
                        X = x[, 1:X_col],
                        parameters = parameters,
                        randomizations = randomizations)})
    as.numeric(w)})
    }
  
  ## Reshape list into matrix ##
  w.matrix <- matrix(unlist(w.list, use.names = FALSE), 
                     ncol = length(allocations), 
                     byrow = FALSE,
                     dimnames = list(gg, allocations))
  
  return(w.matrix)
}

wght_deriv_calc <- function(integrand,
                            parameters,
                            allocation,
                            integrate_allocation = TRUE,
                            ...)
{  
  ## Necessary pieces ##
  integrand <- match.fun(integrand)
  dots <- list(...)
  
  ## Integrand and  arguments ##
  int.args <- append(get_args(integrand, dots),
                     get_args(stats::integrate, dots))
  
  args <- append(append(int.args, get_args(numDeriv::grad, dots)),
                 list(func       = wght_calc, 
                      integrand  = integrand, 
                      allocation = allocation,
                      x          = parameters))
  
  dervs <- do.call(numDeriv::grad, args = args)
  
  dervs
}

wght_deriv_array <- function(parameters, 
                             integrand, 
                             allocations, 
                             G, A, P, 
                             X = NULL,
                             runSilent = TRUE, 
                             integrate_allocation = TRUE,
                             ...)
{
  ## Gather necessary bits ##
  integrand <- match.fun(integrand)
  XX <- cbind(X, A)
  #p  <- ncol(X) # number of predictors
  p <- ifelse(is.null(ncol(X)),0,ncol(X)) 
  pp <- length(parameters)
  aa <- sort(allocations) # Make sure alphas are sorted
  gg <- sort(unique(G))
  k  <- length(allocations) 
  N  <- length(unique(G))
  dots <- list(...)
  
  ## Warnings ##
  
  ## Compute weight (derivative) for each group, parameter, and alpha level ##
  if(!runSilent) print('Calculating array of IP weight derivatives...')
  
  w.list <- lapply(aa, function(allocation){
    w <- by(XX, INDICES = G, simplify = TRUE, 
            FUN = function(x) {
              wght_deriv_calc(parameters = parameters,
                              integrand  = integrand, 
                              allocation = allocation, 
                              integrate_allocation = FALSE,
                              A = x[, p+1], 
                              X = x[, 1:p], 
                              ...)})
    w2 <- matrix(unlist(w, use.names = FALSE), ncol = pp, byrow = TRUE)
    return(w2)}) 
  
  ## Reshape list into array ##
  out <- array(unlist(w.list, use.names = FALSE), 
               dim = c(N, pp, k),
               dimnames = list(gg, names(parameters), aa))
  
  return(out)
}

