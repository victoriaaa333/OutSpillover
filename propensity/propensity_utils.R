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

log_likelihood <- function(parameters,
                           integrand,
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
  args <- append(int.args, list(f = integrand, parameters = parameters))
  
  ## Calculation ##
  attempt <- try(do.call(stats::integrate, args = args))
  val <- if(is(attempt, 'try-error')) NA else attempt$value
  
  return(log(val))
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

score_calc <- function(parameters,
                       integrand,
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
                     list(func = log_likelihood,
                          x    = parameters,
                          integrand = integrand))
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

score_matrix <- function(integrand,
                         X, A, G, 
                         parameters,
                         runSilent = TRUE, 
                         ...)
{
  ## Warnings ##
  # if(length(fixed.effects) != ncol(X)){
  #   stop("The length of params is not equal to the number of predictors")
  # }
  
  ## Necessary bits ##
  integrand <- match.fun(integrand)
  dots <- list(...)
  XX <- cbind(X, A)
  pp <- ncol(X)
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
                                     X = xx[ , 1:pp]))
                 do.call(score_calc, args = args)
               })
  
  ## Reshape list into matrix ##
  out <- matrix(unlist(s.list, use.names = FALSE), 
                ncol = length(parameters), 
                byrow = TRUE,
                dimnames = list(gg, names(parameters)))
  
  out 
}

#-----------------------------------------------------------------------------#
# V Matrix
# 
# Computes the V matrix necessary for variance estimates. Used in 
# \code{\link{ipw_effect_calc}}. See web appendix of Perez et al. 2014 for more details.
# 
# @param scores the output of \code{\link{score_matrix}}
# @param point_estimates output of \code{\link{ipw_point_estimates}}
# @param allocation1 See details in \code{\link{ipw_effect_calc}}.
# @param trt.lvl1 See details in \code{\link{ipw_effect_calc}}.
# @param allocation2 See details in \code{\link{ipw_effect_calc}}.
# @param trt.lvl2 See details in \code{\link{ipw_effect_calc}}.
# @param effect_type See details in \code{\link{ipw_effect_calc}}.
# @param marginal See details in \code{\link{ipw_effect_calc}}.
# @return V matrix
# @export
# 
#-----------------------------------------------------------------------------#

V_matrix <- function(scores, 
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
      xx <- (hold_grp[ , ,a1,] - hold_oal[,,a1]) - (hold_grp[,, a2,] - hold_oal[,,a2])
    } else {
      #xx <- (hold_grp[ , a1, t1] - hold_oal[a1, t1]) - (hold_grp[, a2, t2] - hold_oal[a2, t2])
      xx <- (hold_grp[ , , a1, t1,] - hold_oal[, a1, t1,]) -
        (hold_grp[ , , a2, t2,] - hold_oal[, a2, t2,])
    }
  } 
  else if(effect_type == 'outcome'){
    if(marginal == TRUE){
      #xx <- hold_grp[ , a1] - hold_oal[a1]
      xx <- hold_grp[ , , a1,] - hold_oal[,a1,]
    } else {
      #xx <- hold_grp[ , a1, t1] - hold_oal[ a1, t1]
      xx <- hold_grp[ , , a1, t1,] - hold_oal[, a1, t1,]
    }
  }
  
  ee <- cbind(scores, xx)
  V <- crossprod(ee)/N
  V
}

#-----------------------------------------------------------------------------#
# V propensity_parameter
# 
# Computes the parameter necessary for propensity score estimates. Used in 
# \code{\link{ipw_interference}}. See web appendix of Perez et al. 2014 for more details.
# 

propensity_parameter <- function(formula,
                                 data,
                                 propensity_integrand = 'logit_integrand',
                                 loglihood_integrand = propensity_integrand,
                                 model_method = "glmer",
                                 model_options = list(family = stats::binomial(link = 'logit')),
                                 ...
)
{
  
  ## Necessary bits ##
  dots <- list(...)
  integrandFUN    <- match.fun(propensity_integrand)
  likelihoodFUN   <- match.fun(loglihood_integrand)
  oracle          <- model_method == 'oracle'
  cformula        <- Formula::Formula(formula)
  len_lhs         <- length(cformula)[1]
  len_rhs         <- length(cformula)[2]
  data            <- as.data.frame(data) # make sure data is data.frame not data.table, etc.
  
  ## For the sake of consistency downstream, reorder data frame by groups ##
  group_var <- attr(stats::terms(cformula, lhs = 0, rhs = len_rhs), 'term.labels')
  data <- data[order(data[ , group_var]), ]
  
  ## Parse out the formula into necessary pieces ##
  # Y <- Formula::model.part(cformula, data = data, lhs = 1, drop = TRUE)
  # A <- Formula::model.part(cformula, data = data, lhs = 2, drop = TRUE)
  # G <- Formula::model.part(cformula, data = data, rhs = len_rhs, drop = TRUE)
  
  # Used when there is 'participation' variable
  if(len_lhs > 2){
    B <- Formula::model.part(cformula, data = data, lhs = len_lhs, drop = TRUE)
  } else {
    B <- A
  }
  
  propensity_formula <- formula(stats::terms(cformula, lhs = len_lhs, rhs = -2))
  random.count <- length(lme4::findbars(propensity_formula))
  
  #trt_lvls     <- sort(unique(A))
  N            <- length(unique(G))
  #k            <- length(allocations)
  #l            <- length(trt_lvls)
  
  ## Warnings ##
  if(model_method == 'glm' & random.count > 0 ){
    stop('propensity_formula appears to include a random effect when "glm" was chosen \n 
         for parameter estimation. Set model_method to "glmer" to include a random effect')
  }
  
  if(propensity_integrand == "logit_integrand" & random.count > 1){
    stop('Logit integrand is designed to handle only 1 random effect.')
  }
  
  # if(min(allocations) < 0 | max(allocations) > 1){
  #   stop('Allocations must be between 0 and 1')
  # }
  # 
  # if(length(allocations) < 2){
  #   warning('At least 2 allocations must be specified in order to estimate indirect effects')
  # }
  
  if(N == 1){
    stop('The group variable must have at least 2 groups (more is better).')
  }
  
  # if(!(causal_estimation_options$variance_estimation %in% c('naive', 'robust'))){
  #   stop("The variance estimation method should be 'naive' or 'robust'")
  # }
  
  
  # Set outcome and covariates to numerical  
  data$H <- as.numeric(data$H)  
  X_names <- all.vars(propensity_formula[[3]])
  for (j in 1:(length(X_names)-1)){
    data[,colnames(data) == X_names[j]] <- 
      as.numeric(data[,colnames(data) == X_names[j]])
  }
  
  
  #### Compute Parameter Estimates ####
  
  estimation_args <- append(list(formula = propensity_formula, data = data), 
                            model_options)
  
  parameters <- list()
  
  if(model_method == "glmer"){
    propensity_model <- do.call(lme4::glmer, args = estimation_args)
    parameters$fixed_effects  <- lme4::getME(propensity_model, 'fixef')
    parameters$random_effects <- lme4::getME(propensity_model, 'theta')
    X <- lme4::getME(propensity_model, "X")
    
    if(sum(parameters$random_effects == 0) > 0){
      stop('At least one random effect was estimated as 0. This will lead to a
           non-invertible matrix if using robust variance estimation.')
    }
  } else if(model_method == "glm"){
    propensity_model <- do.call("glm", args = estimation_args)
    parameters$fixed_effects  <- stats::coef(propensity_model)
    parameters$random_effects <- NULL
    X <- stats::model.matrix(propensity_model)
  } else if(model_method == "oracle"){
    parameters$fixed_effects  <- model_options[[1]]
    parameters$random_effects <- model_options[[2]]
    X <- stats::model.matrix(propensity_formula, data)
    
    if(length(parameters$fixed_effects) != ncol(X)){
      stop('oracle fixed effects vector must have length of # of columns of X')
    }
  }
  
  results <- list(parameters, X)
  return(results)
}



