logit_integrand_second_new <- function(b, propensity_X, A, 
                                       parameters,
                                       assigned_group, 
                                       allocation = A, 
                                       randomization = 1)
{
  ## In the case of an intercept-only model, X needs to be converted to matrix
  # for the warning to work
  if(!is.matrix(propensity_X)){
    propensity_X <- as.matrix(propensity_X)
  }
  
  # pass in first_assignment again
  first_assignment <- assigned_group
  
  p <- ncol(propensity_X)
  pp <- length(parameters)/(p+1)
  theta <- parameters[((p+1)*(first_assignment)+1):((p+1)*(first_assignment+1))] 

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

log_likelihood_second_new <- function(parameters,
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
  
  args <- append(int.args, list(f = integrand, parameters = parameters))
                                 # parameters[(1+pp*(first_assignment)):(pp*(first_assignment+1))]
  attempt <- try(do.call(stats::integrate, args = args))
  val <- if(is(attempt, 'try-error')) NA else attempt$value
  likelihood = val
  
  return(log(likelihood))
}


score_calc_deriv_new <- function(parameters,
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
                     list(func = log_likelihood_second_new,
                          x    = parameters,
                          integrand = integrand,
                          P = P,
                          first_assignment = first_assignment,
                          assigned_group = first_assignment,
                          propensity_X = propensity_X))
  
  ## Compute the derivative of the log likelihood for each parameter ##
  do.call(numDeriv::hessian, args = args)
}  


score_matrix_deriv_new <- function(integrand,
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
                 do.call(score_calc_deriv_new, args = args)
               })
  
  ## Reshape list into matrix ##
  U11 = matrix(0, nrow = length(parameters), ncol = length(parameters))
  
  for (i in 1:length(gg)) {
    # this is where -1 * U11 comes 
    U11 =  U11 - s.list[[i]]
  }
  U11 = U11/length(gg)
  U11
}
