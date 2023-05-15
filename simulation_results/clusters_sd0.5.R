source("simulation_results/results_utils.R")

###### 1. influencer outcome model ######
#### 1.1 cat vars ####
# inf_model_cat_var
results = readRDS("simulation_results/cluster_results/inf_model_cat_var.RDS")
results_nocon = as.data.frame(matrix(unlist(lapply(results, function(l) l$nocon)), nrow = 500, byrow = TRUE))
colnames(results_nocon) <- colnames(results[[1]]$nocon)
result_stats(results_nocon, 1 + 1 * 0.5 + 2 * 0.5)
# 2.530 0.030 0.319 0.307 0.948

results_inf = as.data.frame(matrix(unlist(lapply(results, function(l) l$inf)), nrow = 500, byrow = TRUE))
colnames(results_inf) <- colnames(results[[1]]$inf)
result_stats(results_inf, 1 + 1 + 2 * 1)
# 4.048 0.048 0.513 0.589 0.944

results_sp = as.data.frame(matrix(unlist(lapply(results, function(l) l$sp)), nrow = 500, byrow = TRUE))
colnames(results_sp) <- colnames(results[[1]]$sp)
result_stats(results_sp,  1 + 1 * 0.5 + 2 * 0.5)
# 2.545 0.045 0.344 0.330 0.930

results_mixed = as.data.frame(matrix(unlist(lapply(results, function(l) l$mixed)), nrow = 500, byrow = TRUE))
colnames(results_mixed) <- colnames(results[[1]]$mixed)
result_stats(results_mixed, 1 + 1 + 2 * 1)
# 4.172 0.172 0.576 0.657 0.942

inf_cat = rbind(result_stats(results_nocon, 1 + 1 * 0.5 + 2 * 0.5),
                result_stats(results_inf, 1 + 1 + 2 * 1),
                result_stats(results_sp,  1 + 1 * 0.5 + 2 * 0.5),
                result_stats(results_mixed, 1 + 1 + 2 * 1))

#### 1.2 num vars ####
# inf_model_num_var
results = readRDS("simulation_results/cluster_results/sd = 0.5/inf_model_num_var(sd 0.5).RDS")
results_nocon = as.data.frame(matrix(unlist(lapply(results, function(l) l$nocon)), nrow = 500, byrow = TRUE))
colnames(results_nocon) <- colnames(results[[1]]$nocon)
result_stats(results_nocon, 1 + 1 * 0.5 + 2 * 0.5)
# 2.521 0.021 0.310 0.307 0.950

results_inf = as.data.frame(matrix(unlist(lapply(results, function(l) l$inf)), nrow = 500, byrow = TRUE))
colnames(results_inf) <- colnames(results[[1]]$inf)
result_stats(results_inf, 1 + 1 * 1 + 2 * 1)
# 4.002 0.002 0.293 0.292 0.944

results_sp = as.data.frame(matrix(unlist(lapply(results, function(l) l$sp)), nrow = 500, byrow = TRUE))
colnames(results_sp) <- colnames(results[[1]]$sp)
result_stats(results_sp,  1 + 1 * 0.5 + 2 * 0.5)
# 2.533 0.033 1.593 1.572 0.938

results_mixed = as.data.frame(matrix(unlist(lapply(results, function(l) l$mixed)), nrow = 500, byrow = TRUE))
colnames(results_mixed) <- colnames(results[[1]]$mixed)
result_stats(results_mixed,  1 + 1 * 1 + 2 * 1)
# 4.125 0.125 2.720 2.594 0.932

inf_num = rbind(result_stats(results_nocon, 1 + 1 * 0.5 + 2 * 0.5),
                result_stats(results_inf,  1 + 1 * 1 + 2 * 1),
                result_stats(results_sp,  1 + 1 * 0.5 + 2 * 0.5),
                result_stats(results_mixed, 1 + 1 * 1 + 2 * 1))

#### 1.3 both vars ####
results =  readRDS("simulation_results/cluster_results/sd = 0.5/inf_model_both_var(sd 0.5).RDS")
results_nocon = as.data.frame(matrix(unlist(lapply(results, function(l) l$nocon)), nrow = 500, byrow = TRUE))
colnames(results_nocon) <- colnames(results[[1]]$nocon)
result_stats(results_nocon, 1 + 1 * 0.5 + 2 * 0.5)
# SD = 0.5, 2.495 -0.005  0.313  0.303  0.934

results_inf = as.data.frame(matrix(unlist(lapply(results, function(l) l$inf)), nrow = 500, byrow = TRUE))
colnames(results_inf) <- colnames(results[[1]]$inf)
result_stats(results_inf, 1 + 1 + 2 * 1)
# SD = 2, 3.985 -0.015  0.538  0.275  0.680
# SD = 0.5,  4.006 0.006 0.349 0.308 0.904

results_sp = as.data.frame(matrix(unlist(lapply(results, function(l) l$sp)), nrow = 500, byrow = TRUE))
colnames(results_sp) <- colnames(results[[1]]$sp)
result_stats(results_sp,  1 + 1 * 0.5 + 2 * 0.5)
# SD = 0.5, 2.570 0.070 0.856 0.821 0.940

results_mixed = as.data.frame(matrix(unlist(lapply(results, function(l) l$mixed)), nrow = 500, byrow = TRUE))
colnames(results_mixed) <- colnames(results[[1]]$mixed)
result_stats(results_mixed, 1 + 1 + 2 * 1)
# SD = 2,  3.994 -0.006  0.866  0.647  0.844
# SD = 0.5, 4.123 0.123 1.685 1.531 0.926

inf_both = rbind(result_stats(results_nocon, 1 + 1 * 0.5 + 2 * 0.5),
                 result_stats(results_inf,  1 + 1 + 2 * 1),
                 result_stats(results_sp,  1 + 1 * 0.5 + 2 * 0.5),
                 result_stats(results_mixed, 1 + 1 + 2 * 1))

inf_table = as.data.frame(rbind(rep(c("p.est", "bias", "MC SE", "analytical SE", "CP"), 3),
                                cbind(inf_cat, inf_num, inf_both)))

rownames(inf_table) <- c("", "Avg", "Het Infl", "Het Sp", "Het Mixed")

colnames(inf_table) <- c("2 binary vars", rep(" ", 4),
                         "2 cont vars", rep(" ", 4),
                         "1 binary & 1 cont vars", rep(" ", 4))

##### 2. spillover outcome model ######
#### 2.1 cat vars ####
# sp_model_cat_var
results = readRDS("simulation_results/cluster_results/sp_model_cat_var.RDS")
results_nocon = as.data.frame(matrix(unlist(lapply(results, function(l) l$nocon)), nrow = 500, byrow = TRUE))
colnames(results_nocon) <- colnames(results[[1]]$nocon)
result_stats(results_nocon, 1 + 1 * 0.5 + 2 * 0.5)
# 2.554 0.054 0.337 0.349 0.956

results_inf = as.data.frame(matrix(unlist(lapply(results, function(l) l$inf)), nrow = 500, byrow = TRUE))
colnames(results_inf) <- colnames(results[[1]]$inf)
result_stats(results_inf, 1 + 1 * 0.5 + 2 * 0.5)
# 2.560 0.060 0.445 0.518 0.944

results_sp = as.data.frame(matrix(unlist(lapply(results, function(l) l$sp)), nrow = 500, byrow = TRUE))
colnames(results_sp) <- colnames(results[[1]]$sp)
result_stats(results_sp,  1 + 1 * 1 + 2 * 1)
# 4.096 0.096 0.541 0.527 0.950

results_mixed = as.data.frame(matrix(unlist(lapply(results, function(l) l$mixed)), nrow = 500, byrow = TRUE))
colnames(results_mixed) <- colnames(results[[1]]$mixed)
result_stats(results_mixed, 1 + 1 + 2 * 1)
# 4.242 0.242 0.647 0.735 0.956

sp_cat = rbind(result_stats(results_nocon, 1 + 1 * 0.5 + 2 * 0.5),
               result_stats(results_inf,  1 + 1 * 0.5 + 2 * 0.5),
               result_stats(results_sp,  1 + 1 * 1 + 2 * 1),
               result_stats(results_mixed, 1 + 1 * 1 + 2 * 1))

#### 2.2 num vars ####
# sp_model_num_var
results = readRDS("simulation_results/cluster_results/sd = 0.5/sp_model_num_var(sd 0.5).RDS")
results_nocon = as.data.frame(matrix(unlist(lapply(results, function(l) l$nocon)), nrow = 500, byrow = TRUE))
colnames(results_nocon) <- colnames(results[[1]]$nocon)
result_stats(results_nocon, 1 + 1 * 0.5 + 2 * 0.5)
# 2.529 0.029 0.357 0.343 0.932

results_inf = as.data.frame(matrix(unlist(lapply(results, function(l) l$inf)), nrow = 500, byrow = TRUE))
colnames(results_inf) <- colnames(results[[1]]$inf)
result_stats(results_inf, 1 + 1 * 0.5 + 2 * 0.5)
# 2.497 -0.003  0.294  0.285  0.948

results_sp = as.data.frame(matrix(unlist(lapply(results, function(l) l$sp)), nrow = 500, byrow = TRUE))
colnames(results_sp) <- colnames(results[[1]]$sp)
result_stats(results_sp,  1 + 1 * 1 + 2 * 1)
# 3.980 -0.020  1.336  1.447  0.960

results_mixed = as.data.frame(matrix(unlist(lapply(results, function(l) l$mixed)), nrow = 500, byrow = TRUE))
colnames(results_mixed) <- colnames(results[[1]]$mixed)
result_stats(results_mixed,  1 + 1 * 1 + 2 * 1)
# 4.013 0.013 1.998 2.123 0.960

sp_num = rbind(result_stats(results_nocon, 1 + 1 * 0.5 + 2 * 0.5),
               result_stats(results_inf,  1 + 1 * 0.5 + 2 * 0.5),
               result_stats(results_sp,  1 + 1 * 1 + 2 * 1),
               result_stats(results_mixed, 1 + 1 * 1 + 2 * 1))

#### 2.3 both vars ####
# sp_model_both_var
results = readRDS("simulation_results/cluster_results/sd = 2/sp_model_both_var(sd 2).RDS")
results_nocon = as.data.frame(matrix(unlist(lapply(results, function(l) l$nocon)), nrow = 500, byrow = TRUE))
colnames(results_nocon) <- colnames(results[[1]]$nocon)
result_stats(results_nocon, 1 + 1 * 0.5 + 2 * 0.5)
# SD = 1, 2.516 0.016 0.421 0.420 0.934
# SD = 0.5, 2.537 0.037 0.332 0.345 0.954

results_inf = as.data.frame(matrix(unlist(lapply(results, function(l) l$inf)), nrow = 500, byrow = TRUE))
colnames(results_inf) <- colnames(results[[1]]$inf)
result_stats(results_inf, 1 + 1 * 0.5 + 2 * 0.5)
# SD = 2, 2.539 0.039 0.556 0.458 0.896
# SD = 1, 2.487 -0.013  0.369  0.316  0.906
# SD = 0.5, 2.489 -0.011  0.304  0.299  0.938

results_sp = as.data.frame(matrix(unlist(lapply(results, function(l) l$sp)), nrow = 500, byrow = TRUE))
colnames(results_sp) <- colnames(results[[1]]$sp)
result_stats(results_sp,  1 + 1 * 1 + 2 * 1)
# SD = 1, 4.018 0.018 0.506 0.545 0.962
# SD = 0.5, 4.038 0.038 0.822 0.906 0.968

results_mixed = as.data.frame(matrix(unlist(lapply(results, function(l) l$mixed)), nrow = 500, byrow = TRUE))
colnames(results_mixed) <- colnames(results[[1]]$mixed)
result_stats(results_mixed, 1 + 1 + 2 * 1)
# SD = 2, 4.001 0.001 0.487 0.465 0.932
# SD = 1, 4.010 0.010 0.694 0.683 0.940
# SD = 0.5, 4.086 0.086 1.506 1.444 0.944

sp_both = rbind(result_stats(results_nocon, 1 + 1 * 0.5 + 2 * 0.5),
                result_stats(results_inf,  1 + 1 * 0.5 + 2 * 0.5),
                result_stats(results_sp,  1 + 1 * 1 + 2 * 1),
                result_stats(results_mixed, 1 + 1 + 2 * 1))

sp_table = as.data.frame(rbind(rep(c("p.est", "bias", "MC SE", "analytical SE", "CP"), 3),
                               cbind(sp_cat, sp_num, sp_both)))

rownames(sp_table) <- c("", "Avg", "Het Infl", "Het Sp", "Het Mixed")

colnames(sp_table) <- c("2 binary vars", rep(" ", 4),
                        "2 cont vars", rep(" ", 4),
                        "1 binary & 1 cont vars", rep(" ", 4))


##### 3. mixed outcome model ######
#### 3.1 cat vars ####
# mixed_model_cat_var
results = readRDS("simulation_results/cluster_results/mixed_model_cat_var.RDS")
results_nocon = as.data.frame(matrix(unlist(lapply(results, function(l) l$nocon)), nrow = 500, byrow = TRUE))
colnames(results_nocon) <- colnames(results[[1]]$nocon)
result_stats(results_nocon, 1 + 1 * 0.5 + 2 * 0.5 + 3 * 0.5 + 4 * 0.5)
# 6.067 0.067 0.857 0.810 0.936

results_inf = as.data.frame(matrix(unlist(lapply(results, function(l) l$inf)), nrow = 500, byrow = TRUE))
colnames(results_inf) <- colnames(results[[1]]$inf)
result_stats(results_inf, 1 + 1 * 1 + 2 * 1 + 3 * 0.5 + 4 * 0.5)
# 7.567 0.067 1.229 1.622 0.946

results_sp = as.data.frame(matrix(unlist(lapply(results, function(l) l$sp)), nrow = 500, byrow = TRUE))
colnames(results_sp) <- colnames(results[[1]]$sp)
result_stats(results_sp,  1 + 1 * 0.5 + 2 * 0.5 + 3 * 1 + 4 * 1)
# 9.616 0.116 1.306 1.235 0.934

results_mixed = as.data.frame(matrix(unlist(lapply(results, function(l) l$mixed)), nrow = 500, byrow = TRUE))
colnames(results_mixed) <- colnames(results[[1]]$mixed)
result_stats(results_mixed, 1 + 1 * 1 + 2 * 1 + 3 * 1 + 4 * 1)
# 11.451  0.451  1.736  2.360  0.944

mixed_cat = rbind( result_stats(results_nocon, 1 + 1 * 0.5 + 2 * 0.5 + 3 * 0.5 + 4 * 0.5),
                   result_stats(results_inf,  1 + 1 * 1 + 2 * 1 + 3 * 0.5 + 4 * 0.5),
                   result_stats(results_sp,  1 + 1 * 0.5 + 2 * 0.5 + 3 * 1 + 4 * 1),
                   result_stats(results_mixed, 1 + 1 * 1 + 2 * 1 + 3 * 1 + 4 * 1))

#### 3.2 num vars ####
# mixed_model_num_var
results = readRDS("simulation_results/cluster_results/sd = 0.5/mixed_model_num_var(sd 0.5).RDS")
results_nocon = as.data.frame(matrix(unlist(lapply(results, function(l) l$nocon)), nrow = 500, byrow = TRUE))
colnames(results_nocon) <- colnames(results[[1]]$nocon)
result_stats(results_nocon, 1 + 1 * 0.5 + 2 * 0.5 + 3 * 0.5 + 4 * 0.5)
# 5.992 -0.008  0.783  0.798  0.948

results_inf = as.data.frame(matrix(unlist(lapply(results, function(l) l$inf)), nrow = 500, byrow = TRUE))
colnames(results_inf) <- colnames(results[[1]]$inf)
result_stats(results_inf, 1 + 1 * 1 + 2 * 1 + 3 * 0.5 + 4 * 0.5)
# 7.517 0.017 0.701 0.727 0.964

results_sp = as.data.frame(matrix(unlist(lapply(results, function(l) l$sp)), nrow = 500, byrow = TRUE))
colnames(results_sp) <- colnames(results[[1]]$sp)
result_stats(results_sp,  1 + 1 * 0.5 + 2 * 0.5 + 3 * 1 + 4 * 1)
# 9.671 0.171 3.448 3.490 0.960

results_mixed = as.data.frame(matrix(unlist(lapply(results, function(l) l$mixed)), nrow = 500, byrow = TRUE))
colnames(results_mixed) <- colnames(results[[1]]$mixed)
result_stats(results_mixed,  1 + 1 * 1 + 2 * 1 + 3 * 1 + 4 * 1)
# 11.110  0.110  5.342  5.273  0.966

mixed_num = rbind( result_stats(results_nocon, 1 + 1 * 0.5 + 2 * 0.5 + 3 * 0.5 + 4 * 0.5),
                   result_stats(results_inf,  1 + 1 * 1 + 2 * 1 + 3 * 0.5 + 4 * 0.5),
                   result_stats(results_sp,  1 + 1 * 0.5 + 2 * 0.5 + 3 * 1 + 4 * 1),
                   result_stats(results_mixed, 1 + 1 * 1 + 2 * 1 + 3 * 1 + 4 * 1))

#### 3.3 both vars ####
# mixed_model_both_var
results = readRDS("simulation_results/cluster_results/sd = 0.5/mixed_model_both_var(sd 0.5).RDS")
results_nocon = as.data.frame(matrix(unlist(lapply(results, function(l) l$nocon)), nrow = 500, byrow = TRUE))
colnames(results_nocon) <- colnames(results[[1]]$nocon)
result_stats(results_nocon, 1 + 1 * 0.5 + 2 * 0.5 + 3 * 0.5 + 4 * 0.5)
# SD = 1, 6.062 0.062 1.013 0.941 0.924
# SD = 0.5, 6.151 0.151 0.897 0.821 0.924

results_inf = as.data.frame(matrix(unlist(lapply(results, function(l) l$inf)), nrow = 500, byrow = TRUE))
colnames(results_inf) <- colnames(results[[1]]$inf)
result_stats(results_inf,  1 + 1 + 2 * 1 + 3 * 0.5 + 4 * 0.5)
# SD = 2, 7.540 0.040 1.401 1.002 0.830
# SD = 1, 7.528 0.028 0.930 0.752 0.886
# SD = 0.5, 7.427 -0.073  0.827  0.771  0.936

results_sp = as.data.frame(matrix(unlist(lapply(results, function(l) l$sp)), nrow = 500, byrow = TRUE))
colnames(results_sp) <- colnames(results[[1]]$sp)
result_stats(results_sp,  1 + 1 * 0.5 + 2 * 0.5 + 3 * 1 + 4 * 1)
# SD = 1, 9.520 0.020 1.312 1.361 0.946
# SD = 0.5,  9.380 -0.120  2.175  2.273  0.948

results_mixed = as.data.frame(matrix(unlist(lapply(results, function(l) l$mixed)), nrow = 500, byrow = TRUE))
colnames(results_mixed) <- colnames(results[[1]]$mixed)
result_stats(results_mixed,  1 + 1 * 1 + 2 * 1 + 3 * 1 + 4 * 1)
# SD = 2, 11.046  0.046  1.461  1.268  0.908
# SD = 1, 10.914 -0.086  1.876  1.833  0.926
# SD = 0.5, 10.766 -0.234  3.900  3.764  0.934

mixed_both = rbind( result_stats(results_nocon, 1 + 1 * 0.5 + 2 * 0.5 + 3 * 0.5 + 4 * 0.5),
                    result_stats(results_inf,  1 + 1 * 1 + 2 * 1 + 3 * 0.5 + 4 * 0.5),
                    result_stats(results_sp,  1 + 1 * 0.5 + 2 * 0.5 + 3 * 1 + 4 * 1),
                    result_stats(results_mixed, 1 + 1 * 1 + 2 * 1 + 3 * 1 + 4 * 1))

#####

mixed_table = as.data.frame(rbind(rep(c("p.est", "bias", "MC SE", "analytical SE", "CP"), 3),
                                  cbind(mixed_cat, mixed_num, mixed_both)))

rownames(mixed_table) <- c("", "Avg", "Het Infl", "Het Sp", "Het Mixed")

colnames(mixed_table) <- c("2 binary vars", rep(" ", 4),
                           "2 cont vars", rep(" ", 4),
                           "1 binary & 1 cont vars", rep(" ", 4))

table_names <- c(rep(c("Influencer Outcome Model", rep("", 4)), 3),
                 rep(c("Spillover Outcome Model", rep("", 4)), 3),
                 rep(c("Mixed Outcome Model", rep("", 4)), 3))

table <- as.data.frame(rbind(table_names, 
                             cbind(inf_table, sp_table, mixed_table)))

library(gridExtra)
png("sd0.5.png", height = 50*nrow(table), width = 150*ncol(table))
grid.table(table)
dev.off()




