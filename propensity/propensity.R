###################
# Wrap-up functions

ipw_propensity_variance <- function(parameters,
                                    allocations,
                                    causal_estimation_options = 
                                      list(variance_estimation = 'robust'),
                                    integrate_allocation = FALSE,
                                    H, X, A, G,
                                    effect_type = effect_type, #or contrast
                                    propensity_integrand = logit_integrand,
                                    ...){
  
  out <- list()
  dots <- list(...)
  ipw_args <- append(append(dots, causal_estimation_options),
                     list(propensity_integrand = logit_integrand, 
                          loglihood_integrand  = propensity_integrand,
                          allocations          = allocations,
                          parameters           = parameters,
                          runSilent            = TRUE, 
                          integrate_allocation = integrate_allocation,
                          H = H, X = X, A = A, B = A, G = G))
  ipw <- do.call(ipw_interference, args = ipw_args)
  out <- append(out, ipw)
  estimate_args <- list(obj = ipw,
                        variance_estimation = causal_estimation_options$variance_estimation,
                        causal_estimation_options$variance_estimation,
                        alpha1      = allocations[1],
                        trt.lvl1    = 1,
                        alpha2      = allocations[1],
                        trt.lvl2    = 0,
                        marginal    = FALSE,
                        effect_type = effect_type,
                        rescale.factor = 1,
                        conf.level = 0.95,
                        print = FALSE)
  
  est <- do.call(ipw_effect_calc, args = estimate_args)
  return(est)
}

#-----------------------------------------------------------------------------#
# Calculate IPW estimates
#
#  Computes either outcome or effect estimates from the object output by 
#  \code{\link{ipw_interference}}.  
#  
#  @details See \code{\link{direct_effect}}, \code{\link{indirect_effect}},
#  \code{\link{total_effect}}, and \code{\link{overall_effect}} for convenient
#  wrappers of \code{ipw_effect_calc} to compute common causal effects.
#  
#  This table summarizes the value that \code{ipw_effect_calc} returns.
#  \tabular{llc}{
#  Marginal \tab Effect_type    \tab Value returned \cr
#  FALSE    \tab 'outcome' \tab \eqn{\hat{Y}(trt.lvl1, alpha1)}{Yhat(trt.lvl1, alpha1)} \cr
#  TRUE     \tab 'outcome' \tab \eqn{\hat{Y}(alpha1)}{Yhat(alpha1)}  \cr
#  FALSE    \tab 'contrast' \tab 
#  \eqn{\hat{Y}(trt.lvl1, alpha1) - \hat{Y}(trt.lvl2, alpha2)}{Yhat(trt.lvl1, alpha1) - Yhat(trt.lvl2, alpha2)} \cr
#  TRUE    \tab 'contrast' \tab 
#  \eqn{\hat{Y}(alpha1) - \hat{Y}(alpha2)}{Yhat(alpha1) - Yhat(alpha2)} \cr
# }
#  
# @param obj the name of the object created by \code{\link{ipw_interference}}
# @param variance_estimation the variance estimation method.  See 
# \code{\link{interference}} for details
# @param alpha1 the allocation scheme for the outcome of interest or the first
# scheme in the constrast of interest. See details.
# @param trt.lvl1 the treatment level for the outcome of interest or the first
# treatment in the constrast of interest. If marginal = TRUE, this is ignored.
# @param alpha2 the second allocation scheme for the contrast of interest.
# Ignored if effect_type = 'outcome'.
# @param trt.lvl2  the second treatment in the constrast of interest. 
# If marginal = TRUE or effect_type = 'outcome', this is ignored.
# @param effect_type either 'contrast' or 'outcome'
# @param marginal TRUE or FALSE
# @param rescale.factor factor by which to rescale values. Defaults to 1.
# @param conf.level Confidence level for confidence intervals. Defaults to 0.95.
# @param print TRUE/FALSE. If TRUE, the point estimates and confidence interval
# are printed to the console. 
# @return A \code{data.frame} with 1 record and 4 variables: point (the point
#  estimate), variance (the variance estimate), ll (the lower bound of the 
#  confidence interval), and ul (the upper bound of the confidence interval).
# @export
#-----------------------------------------------------------------------------#


ipw_effect_calc <- function(obj, 
                            weights, #added parameter, w.matrix
                            variance_estimation,
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
      # pe          <- oal[a1, t1] - oal[a2, t2]
      # pe_grp_diff <- (grp[ , a1, t1] - oal[a1, t1]) - (grp[ , a2, t2] - oal[a2, t2])
      pe          <- oal[, a1, t1,] - oal[, a2, t2,]
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
    
    # V matrix
    V <- V_matrix(scores = obj$scores, 
                  point_estimates = obj$point_estimates, 
                  allocation1 = a1, allocation2 = a2, 
                  trt.lvl1 = t1, trt.lvl2 = t2, 
                  effect_type = effect_type, marginal = marginal)
    
    vdim <- dim(V)[1]
    
    V21 <- V[vdim, 1:(vdim - 1)] # Last row, up to last column
    V11 <- V[1:(vdim - 1), 1:(vdim - 1)] # up to last row, up to last column
    V22 <- V[vdim, vdim] # bottom right element
    
    ## Sandwich Variance Estimate ##
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
