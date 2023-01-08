source("ipw_point_estimate_tests.R")
source("regression_variance.R")

w.matrix = wght_matrix(plain_integrand, allocations, G, A, P)
fff = 'outcomes'
allocation1 = allocations[1]
allocation2 = NA
# t1 = 0 # default of outcome variance is control group
# t2 = 1

# TODO: check this
t1 = 0
t2 = 1
a1 <- as.character(allocation1)
a2 <- as.character(allocation2)
t1 <- as.character(t1)
t2 <- as.character(t2)

point_estimates_n <- ipw_point_estimates_mixed_test4(H_M, G, A, w.matrix, 
                                             neighinfo = neighinfo, x1 = x1, 
                                             X_type = X_type,  Con_type = "neigh")
point_estimates_g <- ipw_point_estimates_mixed_test4(H, G, A, w.matrix, 
                                              X = X, x0 = x0, 
                                              X_type = X_type, Con_type = "group")

#reg_coef <- point_estimates$outcomes$overall_coefH
reg_coef <- point_estimates_g$outcomes$overall_coefG
#H = H_M

x0 <- as.matrix(c(0.1, "M", 0.2))
x1 <- x0
ipw_regression_variance(H, w.matrix, point_estimates_g, effect_type ='contrast', 
                        marginal = FALSE, allocation1 = allocations[1], allocation2 = allocations[1], 
                        X = X, x0 = x0, X_type = X_type)

x1_num <- c(0.1)
ipw_regression_variance_neigh(H_M, w.matrix, point_estimates_n, effect_type ='contrast', 
                              marginal = FALSE, allocation1 = allocations[1], allocation2 = allocations[1], 
                              neighinfo = neighinfo, x1_num = x1_num)
point_estimates_n$outcomes$overall
