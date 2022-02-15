# variance_estimation = "robust"
# effect_type = 'contrast'
# allocations = list(c(0.5,denominator_alphas),c(0.4,denominator_alphas))
# w.matrix <- wght_matrix(integrand, allocations, G, A, P)
# scores <- w.matrix
# point_estimates  <-  ipw_point_estimates(df$H, df$G, df$A, w.matrix)
# alphas   <- dimnames(scores)[[length(dim(scores))]]
# allocation1 <- alphas[1]
# allocation2 <- alphas[2]
# trt_lvls <- sort(unique(A))
# trt.lvl1 <- as.character(trt_lvls[1])
# trt.lvl2 <- as.character(trt_lvls[2])
# X <- NULL
# marginal <- TRUE

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
      xx <- (hold_grp[ , a1] - hold_oal[a1]) - (hold_grp[, a2] - hold_oal[a2])
    } else {
      xx <- (hold_grp[ , a1, t1] - hold_oal[a1, t1]) - (hold_grp[, a2, t2] - hold_oal[a2, t2])
    }
  } 
  else if(effect_type == 'outcome'){
    if(marginal == TRUE){
      xx <- hold_grp[ , a1] - hold_oal[a1]
    } else {
      xx <- hold_grp[  , a1, t1] - hold_oal[a1, t1]
    }
  }
  
  ee <- cbind(scores, xx)
  V <- crossprod(ee)/N
  V
}


ipw_effect_calc <- function(weights, 
                            point_estimates, 
                            effect_type, 
                            marginal,
                            allocation1, 
                            allocation2 = NA,
                            trt.lvl1 = 0, 
                            trt.lvl2 = 1, 
                            rescale.factor = 1,
                            conf.level = 0.95,
                            print = FALSE){  
  
  ## Necessary bits ##
  N  <- dim(weights)[1] 
  p  <- dim(weights)[2] 
  k  <- length(allocations)
  l  <- dim(point_estimates$outcomes$overall)[2]
  a1 <- as.character(allocation1)
  a2 <- as.character(allocation2)
  t1 <- as.character(trt.lvl1)
  t2 <- as.character(trt.lvl2)
  
  fff <- ifelse(marginal == TRUE, 'marginal_outcomes', 'outcomes')
  
  oal  <- point_estimates[[fff]]$overall
  grp  <- point_estimates[[fff]]$groups
  
  if(effect_type == 'contrast'){
    if(marginal == TRUE){
      pe          <- oal[a1] - oal[a2]
      pe_grp_diff <- (grp[ , a1] - oal[a1]) - (grp[, a2] - oal[a2])
      # if(variance_estimation == 'robust'){
      #   U_pe_grp    <- Ugrp[ , , a1] - Ugrp[ , , a2]
      # }
    } else {
      pe          <- oal[a1, t1] - oal[a2, t2]
      pe_grp_diff <- (grp[ , a1, t1] - oal[a1, t1]) - (grp[ , a2, t2] - oal[a2, t2])
      # if(variance_estimation == 'robust'){
      #   U_pe_grp    <- Ugrp[ , , a1, t1] - Ugrp[ , , a2, t2]
      # }
    }
  } else {
    if(marginal == TRUE){
      pe          <- oal[a1] 
      pe_grp_diff <- (grp[ , a1] - oal[a1])
      # if(variance_estimation == 'robust'){
      #   U_pe_grp    <- Ugrp[ , , a1]
      # }
    } else {
      pe          <- ifelse(k == 1, oal[t1], oal[a1, t1]) 
      pe_grp_diff <- ifelse(k == 1, (grp[ , t1] - oal[t1]),
                            (grp[ , a1, t1] - oal[a1, t1]))
      # if(variance_estimation == 'robust'){
      #   U_pe_grp    <- Ugrp[ , , a1, t1]
      # }
    }
  }
  
  ave <- (1/(N^2)) * (sum((pe_grp_diff)^2, na.rm = T)) * rescale.factor^2
  
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
                    conf.low = pe - me, conf.high = pe + me)
  return(out)
}

ipw_effect_calc <- Vectorize(ipw_effect_calc, 
                             vectorize.args = c("allocation1", 'allocation2', 'trt.lvl1', 'trt.lvl2', 
                                                'marginal', 'effect_type'),
                             SIMPLIFY = TRUE)
#ipw_effect_calc(weights, point_estimates, allocation1, trt.lvl1, allocation2, trt.lvl2, 'contrast', TRUE)
