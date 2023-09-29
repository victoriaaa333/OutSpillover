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

wght_calc_second <- function(parameters,
                      integrand, 
                      allocation,
                      P, # first stage probability
                      propensity_X,
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
    
    for (i in 1:length(P)){
      int.args <- append(get_args(stats::integrate, dots),
                         list(f = integrand, parameters = parameters[(1+pp*(i-1)):(pp*i)],
                              propensity_X = propensity_X, 
                              randomization = randomization))
      args <- append(get_args(integrand, dots), int.args)
      
      ## Compute the integral ##
      # if any of the products within the integrand return Inf, then return NA
      # else return the result of integration
      f <- try(do.call(stats::integrate, args = args), silent = TRUE)
      PrA <- PrA + if(is(f, 'try-error')) NA else P[i] * f$value 
      # add P[i] * integral under paramater[i]
    }

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

wght_matrix_second <- function(integrand, 
                        allocations, 
                        G, A, P,
                        propensity_X = NULL, 
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
      w <- by(cbind(propensity_X, A), INDICES = G, simplify = FALSE, 
              FUN = function(x) {
                wght_calc_second(
                  parameters = parameters,
                  integrand  = integrand, 
                  allocation = allocation,
                  P = P,
                  A = as.numeric(x[, X_col+1]), 
                  propensity_X = propensity_X, 
                  randomization = randomization,
                  integrate_allocation = integrate_allocation)})
      as.numeric(w)})
  }else{
    w.list <- lapply(allocations, function(allocation){
      w <- by(cbind(propensity_X, A), INDICES = G, simplify = FALSE, 
              FUN = function(x) {
                wght_calc_second(
                  parameters = parameters,
                  integrand  = integrand, 
                  allocation = allocation,
                  P = P,
                  A = as.numeric(x[, X_col+1]), 
                  propensity_X = x[, 1:X_col],
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



wght_deriv_calc_second <- function(parameters,
                            integrand,
                            allocation,
                            propensity_X,
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
                 list(func       = wght_calc_second, 
                      integrand  = integrand, 
                      allocation = allocation,
                      x          = parameters,
                      P          = P,
                      propensity_X = propensity_X))
  
  dervs <- do.call(numDeriv::grad, args = args)
  dervs
}

wght_deriv_array_second <- function(parameters, 
                             integrand, 
                             allocations, 
                             G, A, P,
                             propensity_X,
                             integrate_allocation = FALSE,
                             runSilent = TRUE, 
                             ...)
{
  ## Gather necessary bits ##
  integrand <- match.fun(integrand)
  XX <- cbind(propensity_X, A)
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
              wght_deriv_calc_second(parameters = parameters,
                              integrand  = integrand, 
                              allocation = allocation, 
                              integrate_allocation = integrate_allocation,
                              A = x[, p+1], 
                              propensity_X = x[, 1:p], 
                              P = P,
                              ...)})
    w2 <- matrix(unlist(w, use.names = FALSE), ncol = pp, byrow = TRUE)
    return(w2)}) 
  
  ## Reshape list into array ##
  out <- array(unlist(w.list, use.names = FALSE), 
               dim = c(N, pp, k),
               dimnames = list(gg, names(parameters), aa))
  
  return(out)
}

#-----------------------------------------------------------------------------#
# IPW Interference estimation 
#
# Prepares the object necessary to compute IPW effect estimates with 
# \code{\link{ipw_effect_calc}}.
# 
# @inheritParams interference
# @param Y outcome vector
# @param X covariate matrix
# @param A treatmeent vector
# @param B 'participation' vector. Defaults to A in the case there is no
# participation variable.
# @param G group assignment vector
# @param parameters a list of fixed_effects and random_effects
# @param variance_estimation currently supports 'robust' or 'naive'
# @param ... additional arguments passed to other functions such as 
# \code{\link{glmer}}, \code{\link{grad}}, and \code{integrand} or \code{likelihood}.
# @return Returns a list of overall and group-level IPW point estimates 
# (the output of \code{\link{ipw_point_estimates}}), overall and group-level IPW 
# point estimates (using the weight derivatives), scores (the output of 
# \code{\link{score_matrix}}), the computed weight matrix, and the computed 
# weight derivative array.
# @export
#-----------------------------------------------------------------------------#

ipw_interference_second <- function(propensity_integrand,
                             loglihood_integrand = propensity_integrand,
                             allocations,
                             H, propensity_X, A, B = A, G, P,
                             parameters,
                             first_assignments, 
                             variance_estimation,
                             runSilent   = TRUE, 
                             integrate_allocation,
                             ...)
{
  dots <- list(...)
  
  ## Warnings ##
  
  #### Arguments Necessary for Causal Estimation Functions ####
  integrand_args <- get_args(FUN = propensity_integrand, args_list = dots)
  #point_est_args <- get_args(FUN = ipw_point_estimates_mixed_test4, args_list = dots)
  point_est_args <- get_args(FUN = ipw_point_estimates_propensity, args_list = dots)
  loglihood_args <- get_args(FUN = loglihood_integrand, args_list = dots)
  grad_args      <- get_args(FUN = numDeriv::grad, args_list = dots)
  integrate_args <- get_args(FUN = stats::integrate, args_list = dots)
  
  weight_args <- append(append(integrand_args, integrate_args),
                        list(integrand   = propensity_integrand, 
                             allocations = allocations, 
                             propensity_X = propensity_X, A = A, G = G, P = P,
                             parameters = parameters,
                             runSilent  = runSilent, #BB 2015-06-23
                             integrate_allocation = integrate_allocation
                        ))
  #score_args <- append(weight_args, list(first_assignments = first_assignments))
  
  #### Prepare output ####
  out <- list()  
  
  ## Compute Weights ##
  weights <- do.call(wght_matrix_second, args = weight_args)
  
  if(variance_estimation == 'robust'){
    weightd <- do.call(wght_deriv_array_second, args = append(weight_args, grad_args))
    out$weightd <- weightd
    }
  #   U11 <- do.call(score_matrix_deriv, args = append(score_args, grad_args))
  #   out$U11 <- U11
  #   }
  
  
  #### COMPUTE ESTIMATES AND OUTPUT ####
  estimate_args <- append(point_est_args, list(H = H, G = G, A = A))#, list(Y = Y, G = G, A = A)
  point_args    <- append(estimate_args, list(weights = weights))
  
  
  #### Calculate output ####
  #out$point_estimates <- do.call(ipw_point_estimates_mixed_test4, args = point_args)
  out$point_estimates <- do.call(ipw_point_estimates_propensity, args = point_args)
  
  if(variance_estimation == 'robust'){
    U_args     <- append(estimate_args, list(weights = weightd))
    sargs      <- append(append(loglihood_args, grad_args), integrate_args)
    score_args <- append(sargs, list(integrand = loglihood_integrand,
                                     propensity_X = propensity_X, G = G, P = P,
                                     A = B, # Use B for treatment in scores
                                     parameters = parameters,
                                     first_assignments = first_assignments,
                                     runSilent  = runSilent #BB 2015-06-23
    ))
    
    U11 <- do.call(score_matrix_deriv, args = append(score_args, grad_args))
    out$U11 <- U11
    
    # set randomization scheme to 1 for scores for logit_integrand
    score_args$randomization <- 1
    
    out$Upart           <- do.call(ipw_point_estimates_propensity, args = U_args)
    out$scores          <- do.call(score_matrix_second, args = score_args)
  } 
  
  out$weights <- weights
  # out$variance_estimation <- variance_estimation #for use in ipw_effect_calc()
  
  return(out)
}

#-----------------------------------------------------------------------------#
#' Log Likelihood
#' 
#' Used by \code{\link{score_matrix}} to compute the log likelihood.
#' 
#' @param parameters vector of parameters passed to \code{integrand}
#' @param integrand Defaults to logit_integrand
#' @param ... additional arguments passed to \code{integrand} function.
#' @return value of log likelihood
#' @export
#' @importFrom methods is
#-----------------------------------------------------------------------------#

log_likelihood_second <- function(parameters,
                                  integrand,
                                  P, 
                                  first_assignment,
                                  # default is alpha_0 (first stage assignment)
                                  ...)
{
  ## Necessary pieces ##
  integrand <- match.fun(integrand)
  dots      <- list(...)
  dot.names <- names(dots)
  
  ## Integrate() arguments ##
  if(!'lower' %in% dot.names){
    dots$lower <- -Inf
  }
  
  if(!'upper' %in% dot.names){
    dots$upper <- Inf
  }
  
  int.args <- append(get_args(stats::integrate, dots),
                     get_args(integrand, dots))
  
  ## Calculation ##
  if (length(parameters)%%length(P) != 0) {stop("in log_likelihood_second, 
                                                length of parameters and first stage probability is not compatible")}
  if (length(P) == 1 && first_assignment > 0) {stop("in log_likelihood_second,
                                                         first_assignments should be 0 if P = 1")}
  pp <- length(parameters)/length(P)

  args <- append(int.args, list(f = integrand, parameters = 
                                  parameters[(1+pp*(first_assignment)):(pp*(first_assignment+1))]))
  attempt <- try(do.call(stats::integrate, args = args))
  val <- if(is(attempt, 'try-error')) NA else attempt$value
  likelihood = val
  
  return(log(likelihood))
}


#-----------------------------------------------------------------------------#
#' Compute scores for a single group
#' 
#' Used by \code{\link{score_matrix}} to log likelihood derivatives for
#' a single group.
#' 
#' @param parameters vector of parameters passed to \code{integrand}
#' @param integrand function to used for the integrand.
#' Defaults to \code{\link{logit_integrand}}.
#' @param hide.errors Hide errors printed from \code{\link{grad}}.
#' Defaults to true.
#' @param ... additional arguments pass to the integrand function.
#' @return length(theta) vector of scores
#' @export
#-----------------------------------------------------------------------------#

score_calc_second <- function(parameters,
                       integrand,
                       P, 
                       first_assignment, 
                       propensity_X, 
                       hide.errors = TRUE,
                       ...)
{
  ## Necessary bits ##
  integrand <- match.fun(integrand)
  dots <- list(...)
  
  ## Function arguments ##
  int.args <- append(get_args(integrand, dots),
                     get_args(stats::integrate, dots))
  fargs    <- append(int.args, get_args(numDeriv::grad, dots))
  
  args     <- append(fargs,
                     list(func = log_likelihood_second,
                          x    = parameters,
                          integrand = integrand,
                          first_assignment = first_assignment,
                          propensity_X = propensity_X, 
                          P = P))
  
  ## Compute the derivative of the log likelihood for each parameter ##
  do.call(numDeriv::grad, args = args)
}  

#-----------------------------------------------------------------------------#
#' Calculate matrix of log Likelihood derivatives
#' 
#' @param integrand function passed to \code{\link{log_likelihood}}. Defaults to
#' \code{\link{logit_integrand}}
#' @param X covariate matrix
#' @param A vector of treatment assignments
#' @param G vector of group assignments
#' @param parameters vector of parameters passed to \code{integrand}
#' @param runSilent If FALSE, prints errors to console. Defaults to TRUE.
#' @param ... additional arguments passed to \code{integrand} or \code{\link{grad}}.
#' For example, one can change the \code{method} argument in \code{grad}.
#' @return N X length(params) matrix of scores
#' @export
#-----------------------------------------------------------------------------#

score_matrix_second <- function(integrand,
                         propensity_X, A, G, P,
                         parameters,
                         first_assignments,
                         runSilent = TRUE, 
                         ...)
{
  ## Necessary bits ##
  integrand <- match.fun(integrand)
  dots <- list(...)
  XX <- cbind(propensity_X, A, first_assignments[G])
  pp <- ncol(propensity_X)
  gg <- sort(unique(G))
  
  ## Compute score for each group and parameter ##
  int.args <- append(get_args(integrand, dots),
                     get_args(stats::integrate, dots))
  fargs <- append(int.args, get_args(numDeriv::grad, dots))
  
  if(!runSilent) print("Calculating matrix of scores...")
  
  s.list <- by(XX, INDICES = G, simplify = TRUE, 
               FUN = function(xx) {
                 args <- append(fargs, 
                                list(integrand  = integrand, 
                                     parameters = parameters,
                                     A = xx[ , (pp + 1)],
                                     propensity_X = xx[ , 1:pp],
                                     P = P,
                                     first_assignment = xx[, (pp + 2)][1])) 
                 # first assignment for each group is the same, so [1]
                 do.call(score_calc_second, args = args)
               })
  
  ## Reshape list into matrix ##
  out <- matrix(unlist(s.list, use.names = FALSE), 
                ncol = length(parameters), 
                byrow = TRUE,
                dimnames = list(gg, names(parameters)))
  
  out 
}

######
# second derivative of log likelihood

score_calc_deriv <- function(parameters,
                              integrand,
                              P, 
                              first_assignment, 
                              propensity_X,
                              hide.errors = TRUE,
                              ...)
{
  ## Necessary bits ##
  integrand <- match.fun(integrand)
  dots <- list(...)
  
  ## Function arguments ##
  int.args <- append(get_args(integrand, dots),
                     get_args(stats::integrate, dots))
  fargs    <- append(int.args, get_args(numDeriv::grad, dots))
  
  pp <- length(parameters)/length(P)
  parameters_integrand <-
    parameters[(1+pp*(first_assignment)):(pp*(first_assignment+1))]
  
  args     <- append(fargs,
                     list(func = log_likelihood_second,
                          x    = parameters,
                          integrand = integrand,
                          P = P,
                          first_assignment = first_assignment,
                          #parameters_integrand = parameters_integrand,
                          #assigned_group = first_assignment,
                          propensity_X = propensity_X))
  #print(args)
  ## Compute the derivative of the log likelihood for each parameter ##
  do.call(numDeriv::hessian, args = args)
}  


score_matrix_deriv <- function(integrand,
                               propensity_X, A, G, P,
                                parameters,
                                first_assignments,
                                runSilent = TRUE, 
                                ...)
{
  ## Necessary bits ##
  integrand <- match.fun(integrand)
  dots <- list(...)
  XX <- cbind(propensity_X, A, first_assignments[G])
  pp <- ncol(propensity_X)
  gg <- sort(unique(G))
  
  ## Compute score for each group and parameter ##
  int.args <- append(get_args(integrand, dots),
                     get_args(stats::integrate, dots))
  fargs <- append(int.args, get_args(numDeriv::grad, dots))
  
  if(!runSilent) print("Calculating matrix of scores...")
  
  s.list <- by(XX, INDICES = G, simplify = TRUE, 
               FUN = function(xx) {
                 args <- append(fargs, 
                                list(parameters = parameters,
                                     integrand  = integrand, 
                                     P = P,
                                     first_assignment = as.numeric(xx[, (pp + 2)][1]),
                                     propensity_X = xx[ , 1:pp],
                                     A = xx[ , (pp + 1)])) 
                 # first assignment for each group is the same, so [1]
                 do.call(score_calc_deriv, args = args)
               })
  
  ## Reshape list into matrix ##
  U11 = matrix(0, nrow = length(parameters), ncol = length(parameters))
  
  for (i in 1:length(gg)) {
    # this is where -1 * U11 comes 
    U11 =  U11 - s.list[[i]]
    #U11 = U11 - s.list[[i]]/sum(first_assignments == first_assignments[i])
  }
  U11 = U11/length(gg)
  U11
}



############
# wrap-up function
##############
ipw_effect_calc_second <- function(obj, 
                                   weights, #added parameter, w.matrix
                                   variance_estimation,
                                   P,
                                   propensity_X, 
                                   alpha1, 
                                   trt.lvl1, 
                                   alpha2 = NA, 
                                   trt.lvl2 = NA,
                                   effect_type,
                                   marginal,
                                   rescale.factor = 1,
                                   conf.level = 0.95,
                                   print = FALSE)
{
  
  allocations <- dimnames(obj$weights)[[2]] 
  
  ## Warnings ##
  # Print error if either estimates with alpha1 have been computed 
  # or a constrast is being estimated when estimates for alpha2
  # have not been computed
  if (!(alpha1 %in% allocations) | (effect_type == 'contrast' & !(alpha2 %in% allocations))){
    stop(paste('At least one of the chosen coverage levels has not been estimated.\n',
               'Select from the following: \n', 
               paste(allocations, collapse = ' ')))
  }
  
  ## Necessary bits ##
  N  <- dim(obj$weights)[1] 
  p  <- dim(obj$scores)[2] 
  k  <- length(allocations)
  l  <- dim(obj$point_estimates$outcomes$overall)[2]
  a1 <- as.character(alpha1)
  a2 <- as.character(alpha2)
  t1 <- as.character(trt.lvl1)
  t2 <- as.character(trt.lvl2)
  
  fff <- ifelse(marginal == TRUE, 'marginal_outcomes', 'outcomes')
  
  oal  <- obj$point_estimates[[fff]]$overall
  grp  <- obj$point_estimates[[fff]]$groups
  
  if(variance_estimation == 'robust'){
    Uoal <- obj$Upart[[fff]]$overall 
    Ugrp <- obj$Upart[[fff]]$groups
    
    # Cludgy workaround for case of 1 fixed effect: add dimension to Ugrp array #
    if(p == 1){
      names <- dimnames(Ugrp)
      if(marginal == TRUE){
        Ugrp <-  array(c(Ugrp[1:N, ], 1, Ugrp[, 1:k]),
                       dim=c(N, 1, k),
                       dimnames = list(names[[1]], 'Intercept', names[[2]]))
      } else {
        Ugrp <-  array(c(Ugrp[1:N, , ], 1, Ugrp[, 1:k , 1:l]),
                       dim=c(N, 1, k, l),
                       dimnames = list(names[[1]], 'Intercept', names[[2]], names[[3]]))
      }
    }
  }
  
  if(effect_type == 'contrast'){
    if(marginal == TRUE){
      pe          <- oal[a1] - oal[a2]
      pe_grp_diff <- (grp[ , a1] - oal[a1]) - (grp[, a2] - oal[a2])
      if(variance_estimation == 'robust'){
        U_pe_grp    <- Ugrp[ , , a1] - Ugrp[ , , a2]
      }
    } else {
      # pe          <- oal[a1, t1, ] - oal[a2, t2, ]
      # pe_grp_diff <- (grp[ , a1, t1, ] - oal[a1, t1, ]) - (grp[ , a2, t2, ] - oal[a2, t2, ])
      
      pe          <- oal[, a1, t1, ] - oal[, a2, t2,]
      pe_grp_diff <- (grp[, , a1, t1, ] - oal[,a1, t1,]) - (grp[ , , a2, t2,] - oal[, a2, t2,])
      
      if(variance_estimation == 'robust'){
        U_pe_grp    <- as.matrix(Ugrp[, , a1, t1, ] - Ugrp[, , a2, t2, ])
        #as.matrix(Ugrp[ , a1, t1, ] - Ugrp[ , a2, t2, ])
      }
    }
  } else {
    if(marginal == TRUE){
      pe          <- oal[a1] 
      pe_grp_diff <- (grp[ , a1] - oal[a1])
      if(variance_estimation == 'robust'){
        U_pe_grp    <- Ugrp[ , , a1]
      }
    } else {
      pe          <- oal[, a1, t1, ] 
      pe_grp_diff <- (grp[ , , a1, t1, ] - oal[, a1, t1, ])
      if(variance_estimation == 'robust'){
        U_pe_grp    <- Ugrp[ , , a1, t1, ]
      }
    }
  }
  
  #### VARIANCE ESTIMATION ####
  if(variance_estimation == 'robust'){
    # partial U matrix
    if(p == 1){
      U21 <- sum(-U_pe_grp)/N
    } else {
      U21 <- (t(as.matrix(apply(-U_pe_grp, 2, sum, na.rm = T))))/N
    }
    
    # U11
    U11 <- obj$U11
    U <- cbind(rbind(U11, U21), c(rep(0, dim(U11)[1]), -1))
    # V matrix
    V <- V_matrix_second(scores = obj$scores, 
                  point_estimates = obj$point_estimates, 
                  allocation1 = a1, allocation2 = a2, 
                  trt.lvl1 = t1, trt.lvl2 = t2, 
                  effect_type = effect_type, marginal = marginal)
    
    vdim <- dim(V)[1]

    V21 <- V[vdim, 1:(vdim - 1)] # Last row, up to last column
    V11 <- V[1:(vdim - 1), 1:(vdim - 1)] # up to last row, up to last column
    V22 <- V[vdim, vdim] # bottom right element

    ## Sandwich Variance Estimate ##
    inv_U = solve(U) 
    sigma = inv_U %*% V %*% t(inv_U)
    
    #ave <- sigma[dim(sigma)[1], dim(sigma)[2]]/N * rescale.factor^2
    ave <- ((U21 - 2*V21) %*% solve(V11) %*% t(U21) + V22)/N * rescale.factor^2
    
  } else if(variance_estimation == 'naive'){
    ave <- (1/(N^2)) * (sum((pe_grp_diff)^2, na.rm = T)) * rescale.factor^2
  }
  
  ## Confidence Intervals ##
  qq <- qnorm(conf.level + (1 - conf.level)/2)
  me <- qq * sqrt(ave)
  
  ## Prepare Output ##
  pe <- pe * rescale.factor
  
  if(print == TRUE){
    toprint <- paste0('Estimate: ', round(pe, 2), ' ',
                      conf.level*100, '% CI: (', 
                      round(pe - me, 2), ', ', round(pe + me, 2), ')' )
    print(toprint)
  }
  
  out <- data.frame(estimate = pe,
                    std.error = sqrt(ave), 
                    conf.low = pe - me, 
                    conf.high = pe + me)
  return(out)
}

effect_grid <- function(allocations, treatments = c(0,1))
{
  marginal    <- c('TRUE', 'FALSE')
  
  # Outcomes
  g1.1 <- expand.grid(alpha1 = allocations, trt1 = treatments, 
                      alpha2 = NA, trt2 = NA,
                      marginal = FALSE, 
                      effect_type = 'outcome', effect = 'outcome',
                      stringsAsFactors = FALSE)
  g1.2 <- expand.grid(alpha1 = allocations, trt1 = NA, 
                      alpha2 = NA, trt2 = NA,
                      marginal = TRUE, 
                      effect_type = 'outcome', effect = 'outcome',
                      stringsAsFactors = FALSE)
  
  # Direct Effects
  g2 <- expand.grid(alpha1 = allocations, trt1 = treatments, 
                    alpha2 = NA, trt2 = treatments,
                    marginal = FALSE, 
                    effect_type = 'contrast', effect = 'direct',
                    stringsAsFactors = FALSE)
  g2$alpha2 <- g2$alpha1
  g2 <- g2[g2$trt1 != g2$trt2, ]
  
  # Indirect Effects
  g3 <- expand.grid(alpha1 = allocations, trt1 = treatments, 
                    alpha2 = allocations, trt2 = NA,
                    marginal = FALSE, 
                    effect_type = 'contrast', effect = 'indirect',
                    stringsAsFactors = FALSE)
  g3$trt2 <- g3$trt1
  
  # Total Effects
  g4 <- expand.grid(alpha1 = allocations, trt1 = treatments, 
                    alpha2 = allocations, trt2 = treatments,
                    marginal = FALSE, 
                    effect_type = 'contrast', effect = 'total',
                    stringsAsFactors = FALSE)
  g4 <- g4[g4$trt1 != g4$trt2, ]
  
  # Overall Effects
  g5 <- expand.grid(alpha1 = allocations, trt1 = NA, 
                    alpha2 = allocations, trt2 = NA,
                    marginal = TRUE, 
                    effect_type = 'contrast', effect = 'overall',
                    stringsAsFactors = FALSE)
  
  out <- rbind(g1.1, g1.2, g2, g3, g4, g5)
  rownames(out) <- NULL # Rownames aren't useful
  return(out) 
}

ipw_propensity_variance_second <- function(parameters,
                                    allocations,
                                    H, propensity_X, A, G, P, 
                                    first_assignments, 
                                    causal_estimation_options = 
                                      list(variance_estimation = 'robust'),
                                    integrate_allocation = FALSE,
                                    effect_type = effect_type, #or contrast
                                    propensity_integrand = logit_integrand_second,
                                    ...){
  
  out <- list()
  dots <- list(...)
  ipw_args <- append(append(dots, causal_estimation_options),
                     list(propensity_integrand = logit_integrand_second, 
                          loglihood_integrand  = propensity_integrand,
                          allocations          = allocations,
                          parameters           = parameters,
                          runSilent            = TRUE, 
                          integrate_allocation = integrate_allocation,
                          H = H, propensity_X = propensity_X, 
                          A = A, B = A, G = G, P = P, first_assignments = first_assignments))
  ipw <- do.call(ipw_interference_second, args = ipw_args)
  out <- append(out, ipw)
  estimate_args <- list(obj = ipw,
                        variance_estimation = causal_estimation_options$variance_estimation,
                        causal_estimation_options$variance_estimation,
                        propensity_X = propensity_X, 
                        P           = P, 
                        alpha1      = allocations[1],
                        trt.lvl1    = 1,
                        alpha2      = allocations[1],
                        trt.lvl2    = 0,
                        marginal    = FALSE,
                        effect_type = effect_type,
                        rescale.factor = 1,
                        conf.level = 0.95,
                        print = FALSE)
  
  est <- do.call(ipw_effect_calc_second, args = estimate_args)
  return(est)
}


V_matrix_second <- function(scores, 
                     point_estimates, 
                     allocation1, 
                     trt.lvl1, 
                     allocation2 = NA, 
                     trt.lvl2    = NA, 
                     effect_type, 
                     marginal){
  ## Necessary bits ##
  N  <- dim(scores)[1]
  p  <- dim(scores)[2]
  a1 <- allocation1
  a2 <- allocation2
  t1 <- trt.lvl1
  t2 <- trt.lvl2
  
  ## Grab the last element of the psi(O, theta) vector: psi_a, alpha ##
  fff <- ifelse(marginal == TRUE, 'marginal_outcomes', 'outcomes')
  hold_oal <- point_estimates[[fff]]$overall 
  hold_grp <- point_estimates[[fff]]$groups
  
  if(effect_type == 'contrast'){   
    if(marginal == TRUE){
      #xx <- (hold_grp[ , a1] - hold_oal[a1]) - (hold_grp[, a2] - hold_oal[a2])
      xx <- (hold_grp[ , a1,] - hold_oal[, a1]) - (hold_grp[, a2,] - hold_oal[, a2])
    } else {
      # xx <- (hold_grp[, a1, t1,] - hold_oal[a1, t1,]) -
      #   (hold_grp[ , a2, t2,] - hold_oal[a2, t2,])
      xx <- (hold_grp[, , a1, t1,] - hold_oal[, a1, t1,]) -
        (hold_grp[, , a2, t2,] - hold_oal[, a2, t2,])
    }
  } 
  else if(effect_type == 'outcome'){
    if(marginal == TRUE){
      xx <- hold_grp[ , a1,] - hold_oal[a1, ]
    } else {
      xx <- hold_grp[ , a1, t1,] - hold_oal[a1, t1,]
    }
  }
  
  ee <- cbind(scores, xx)
  V <- crossprod(ee)/N
  #ee <- ee[complete.cases(ee),]
  #V <- ee/dim(ee)[1]
  V
}

logit_integrand_second <- function(b, propensity_X, A, 
                            parameters,
                            allocation = A, 
                            randomization = 1)
{
  ## In the case of an intercept-only model, X needs to be converted to matrix
  # for the warning to work
  if(!is.matrix(propensity_X)){
    propensity_X <- as.matrix(propensity_X)
  }
  
  theta <- parameters 
  p <- ncol(propensity_X)
  #print(theta[1])
  
  ## Warnings ##
  # if(p != ncol(X)){
  #   stop('The number of fixed effect parameters is not equal to the number \n
  #        of columns in the covariate matrix')
  # }
  
  if(length(A) != nrow(propensity_X)){
    stop('Length of treatment vector is not equal to number of observations in
         propensity_X matrix')
  }
  
  # Check whether to ignore random effect
  ignore_re <- (length(theta) == p || theta[p + 1] <= 0)
  
  ## Calculations ## 
  if(ignore_re){
    pr.b <- randomization * (stats::plogis(propensity_X %*% theta[1:p]))
  } else {
    if (nrow(X)==1) {
      linpred_vec <- drop(outer(propensity_X %*% theta[1:p], b, '+'))
      ##drop() will return a vector if X has one row - BGB 2017-02-12
      ##solution: create a matrix of one row from that vector - BGB 2017-02-12
      linpred_mat <- matrix(linpred_vec, byrow=TRUE,
                            nrow=1, ncol = length(linpred_vec))
      pr.b <- randomization * (stats::plogis(linpred_mat))
      ##pr.b should not throw errors in the apply() fun below. - BGB 2017-02-12
    } else {
      pr.b <- randomization * (stats::plogis(drop(outer(propensity_X %*% theta[1:p], b, '+'))))
    }
  }
  
  hh <- as.matrix((pr.b/allocation)^A * ((1-pr.b)/(1 - allocation))^(1-A))
  
  if(ignore_re){
    # in this way dnorm integrates to one when integrating from -Inf to Inf
    out <- exp(sum(log(hh))) * stats::dnorm(b, mean=0, sd = 1) 
  } else {
    hha <- apply(hh, 2, function(x) exp(sum(log(x))))
    out <- hha * stats::dnorm(b, mean=0, sd = theta[p + 1])
  }
  return(out)
}

