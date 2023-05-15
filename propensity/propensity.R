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

