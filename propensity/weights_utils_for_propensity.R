wght_calc_second_prop <- function(parameters,
                             integrand, 
                             allocation,
                             P, # first stage probability
                             propensity_X,
                             first_assignment,
                             randomization = 1,
                             integrate_allocation = FALSE,
                             ...)
{  
  ## Necessary pieces ##
  integrand         <- match.fun(integrand)
  integrand.formals <- names(formals(integrand))
  dots              <- list(...)
  dot.names         <- names(dots)
  #propensity_X      <- dots[[match.arg('propensity_X', dot.names)]]
  A                 <- dots[[match.arg('A', dot.names)]]
  #P                 <- dots[[match.arg('P', dot.names)]]
  
  numerator_alpha <- allocation[1]
  denominator_alphas <- allocation[-1]
  pp <- dim(propensity_X)[2] + 1 # length of each parameter if with random effect
  
  ## Warnings ##
  if(!'A' %in% dot.names){
    stop("The argument 'A' (treatment assignment) must be specified")
  }
  
  if (is.null(parameters)){
    PrA <- integrand(A, denominator_alphas, P)
  }else{
    
    ## Warnings ##
    if (length(parameters) != pp*length(P)) {
      stop(paste(pp, length(P), length(parameters), "length of first stage probability is not the same dim for parameter"))}
    
    PrA <- 0
    # Integrate() arguments ##
    if(!'lower' %in% dot.names){
      dots$lower <- -Inf
    }
    if(!'upper' %in% dot.names){
      dots$upper <- Inf
    }
    

    int.args <- append(get_args(stats::integrate, dots),
                       list(f = integrand, parameters = 
                              parameters[(1+pp*(first_assignment)):(pp*(first_assignment+1))],
                            propensity_X = propensity_X,
                            randomization = randomization))
    args <- append(get_args(integrand, dots), int.args)
    
    ## Compute the integral ##
    # if any of the products within the integrand return Inf, then return NA
    # else return the result of integration
    f <- try(do.call(stats::integrate, args = args), silent = TRUE)
    PrA <- PrA + if(is(f, 'try-error')) NA else f$value
  }
  
  if(integrate_allocation == TRUE){
    weight <- 1/PrA
  } else {
    ppp    <- prod(numerator_alpha^A * (1-numerator_alpha)^(1-A))
    weight <- ppp/PrA
  }
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

wght_matrix_second_prop <- function(integrand, 
                               allocations, 
                               G, A, P,
                               propensity_X = NULL, 
                               first_assignments,
                               parameters = NULL,
                               randomization = 1,
                               integrate_allocation = FALSE,
                               runSilent = TRUE,
                               ...)
{
  ## Gather necessary bits ##
  X_col  <- ifelse(is.null(ncol(propensity_X)), 0, ncol(propensity_X)) 
  gg <- sort(unique(G))
  
  ## Compute weight for each group and allocation level ##
  if(!runSilent) print('Calculating matrix of IP weights...') 
  
  if(X_col == 0){
    w.list <- lapply(allocations, function(allocation){
      w <- by(cbind(propensity_X, A, G), INDICES = G, simplify = FALSE, 
              FUN = function(x) {
                wght_calc_second_prop(
                  parameters = parameters,
                  integrand  = integrand, 
                  allocation = allocation,
                  P = P,
                  A = as.numeric(x[, X_col+1]), 
                  propensity_X = propensity_X, 
                  first_assignment = first_assignments[x[,X_col+2][1]],
                  randomization = randomization,
                  integrate_allocation = integrate_allocation)})
      as.numeric(w)})
  }else{
    w.list <- lapply(allocations, function(allocation){
      w <- by(cbind(propensity_X, A, G), INDICES = G, simplify = FALSE, 
              FUN = function(x) {
                wght_calc_second_prop(
                  parameters = parameters,
                  integrand  = integrand, 
                  allocation = allocation,
                  P = P,
                  A = as.numeric(x[, X_col+1]), 
                  propensity_X = x[, 1:X_col],
                  first_assignment =  first_assignments[x[,X_col+2][1]],
                  randomization = randomization,
                  integrate_allocation = integrate_allocation)})
      as.numeric(w)})
  }
  
  ## Reshape list into matrix ##
  w.matrix <- matrix(unlist(w.list, use.names = FALSE), 
                     ncol = length(allocations), 
                     byrow = FALSE,
                     dimnames = list(gg, allocations))
  
  return(w.matrix)
}



wght_deriv_calc_second_prop <- function(parameters,
                                   integrand,
                                   allocation,
                                   propensity_X,
                                   first_assignment,
                                   integrate_allocation = FALSE,
                                   P = P,
                                   ...)
{  
  ## Necessary pieces ##
  integrand <- match.fun(integrand)
  dots <- list(...)
  
  ## Integrand and arguments ##
  int.args <- append(get_args(integrand, dots),
                     get_args(stats::integrate, dots))
  
  args <- append(append(int.args, get_args(numDeriv::grad, dots)),
                 list(func       = wght_calc_second_prop, 
                      integrand  = integrand, 
                      allocation = allocation,
                      x          = parameters,
                      P          = P,
                      propensity_X = propensity_X,
                      first_assignment = first_assignment))
  
  dervs <- do.call(numDeriv::grad, args = args)
  dervs
}

wght_deriv_array_second_prop <- function(parameters, 
                                    integrand, 
                                    allocations, 
                                    G, A, P,
                                    propensity_X,
                                    first_assignments,
                                    integrate_allocation = FALSE,
                                    runSilent = TRUE, 
                                    ...)
{
  ## Gather necessary bits ##
  integrand <- match.fun(integrand)
  XX <- cbind(propensity_X, A, G)
  p <- ifelse(is.null(ncol(propensity_X)), 0, ncol(propensity_X)) # number of predictors
  pp <- length(parameters)
  
  ## reformatting allocation ##
  # numerator_alphas = unique(unlist(lapply(allocations, function (x) x[1])))
  # aa <- sort(allocations) # Make sure alphas are sorted
  aa <- allocations
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
              wght_deriv_calc_second_prop(parameters = parameters,
                                     integrand  = integrand, 
                                     allocation = allocation, 
                                     integrate_allocation = integrate_allocation,
                                     A = x[, p+1], 
                                     propensity_X = x[, 1:p], 
                                     P = P,
                                     first_assignment =  first_assignments[x[, p+2][1]],
                                     ...)})
    w2 <- matrix(unlist(w, use.names = FALSE), ncol = pp, byrow = TRUE)
    return(w2)}) 
  
  ## Reshape list into array ##
  out <- array(unlist(w.list, use.names = FALSE), 
               dim = c(N, pp, k),
               dimnames = list(gg, names(parameters), aa))
  
  return(out)
}
