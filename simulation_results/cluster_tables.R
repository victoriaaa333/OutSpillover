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
results = rbind(readRDS("simulation_results/cluster_results/threshold_x_1/inf_model_num_var(x_1).RDS"),
                readRDS("simulation_results/cluster_results/threshold_x_1/inf_model_num_var(x_1)(1).RDS"))
results_nocon = as.data.frame(matrix(unlist(lapply(results, function(l) l$nocon)), nrow = 500, byrow = TRUE))
colnames(results_nocon) <- colnames(results[[1]]$nocon)
result_stats(results_nocon, 1 + 1 * 0.5 + 2 * 0.5)
# 2.543 0.043 0.318 0.314 0.940
# 2.539 0.039 0.328 0.313 0.934

results_inf = as.data.frame(matrix(unlist(lapply(results, function(l) l$inf)), nrow = 500, byrow = TRUE))
colnames(results_inf) <- colnames(results[[1]]$inf)
#result_stats(results_inf, 1 + 1 * 0.1 + 2 * 0.2)
# 1.504 0.004 0.221 0.299 0.984
result_stats(results_inf, 1 + 1 * 1 + 2 * 1)
# 3.998 -0.002  0.309  0.260  0.890

results_sp = as.data.frame(matrix(unlist(lapply(results, function(l) l$sp)), nrow = 500, byrow = TRUE))
colnames(results_sp) <- colnames(results[[1]]$sp)
result_stats(results_sp,  1 + 1 * 0.5 + 2 * 0.5)
# 2.516 0.016 0.822 0.849 0.952
# 2.544 0.044 1.025 1.000 0.940

results_mixed = as.data.frame(matrix(unlist(lapply(results, function(l) l$mixed)), nrow = 500, byrow = TRUE))
colnames(results_mixed) <- colnames(results[[1]]$mixed)
# result_stats(results_mixed,  1 + 1*0.1 + 2 * 0.2)
# 1.553 0.053 0.856 0.875 0.936
result_stats(results_mixed,  1 + 1 * 1 + 2 * 1)
# 4.002 0.002 1.374 1.263 0.926

inf_num = rbind(result_stats(results_nocon, 1 + 1 * 0.5 + 2 * 0.5),
            result_stats(results_inf,  1 + 1 * 1 + 2 * 1),
            result_stats(results_sp,  1 + 1 * 0.5 + 2 * 0.5),
            result_stats(results_mixed, 1 + 1 * 1 + 2 * 1))

#### 1.3 both vars ####
# inf_model_both_var
results = readRDS("simulation_results/cluster_results/threshold_x_1/inf_model_both_var(x_1).RDS")
results_nocon = as.data.frame(matrix(unlist(lapply(results, function(l) l$nocon)), nrow = 500, byrow = TRUE))
colnames(results_nocon) <- colnames(results[[1]]$nocon)
result_stats(results_nocon, 1 + 1 * 0.5 + 2 * 0.5)
#  2.525 0.025 0.314 0.306 0.932
#  2.532 0.032 0.300 0.306 0.956
#  SD  = 1, 2.516 0.016 0.301 0.310 0.952

results_inf = as.data.frame(matrix(unlist(lapply(results, function(l) l$inf)), nrow = 500, byrow = TRUE))
colnames(results_inf) <- colnames(results[[1]]$inf)
# result_stats(results_inf, 1 + 1 + 2 * 0.1)
# 2.241 0.041 0.695 0.677 0.962
result_stats(results_inf, 1 + 1 + 2 * 1)
# 3.991 -0.009  0.854  0.866  0.954
# SD = 1, 4.003 0.003 0.380 0.277 0.824

results_sp = as.data.frame(matrix(unlist(lapply(results, function(l) l$sp)), nrow = 500, byrow = TRUE))
colnames(results_sp) <- colnames(results[[1]]$sp)
result_stats(results_sp,  1 + 1 * 0.5 + 2 * 0.5)
# 2.581 0.081 2.864 3.049 0.956
# 2.480 -0.020  3.632  3.644  0.946
#  2.512 0.012 0.579 0.528 0.920

results_mixed = as.data.frame(matrix(unlist(lapply(results, function(l) l$mixed)), nrow = 500, byrow = TRUE))
colnames(results_mixed) <- colnames(results[[1]]$mixed)
# result_stats(results_mixed, 1 + 1 + 2 * 0.1)
# 2.216  0.016 14.338 13.660  0.912
result_stats(results_mixed, 1 + 1 + 2 * 1)
# 4.513  0.513 22.386 21.826  0.946
# 4.007 0.007 1.005 0.833 0.892

inf_both = rbind(result_stats(results_nocon, 1 + 1 * 0.5 + 2 * 0.5),
            result_stats(results_inf,  1 + 1 + 2 * 1),
            result_stats(results_sp,  1 + 1 * 0.5 + 2 * 0.5),
            result_stats(results_mixed, 1 + 1 + 2 * 1))

#####

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
results = readRDS("simulation_results/cluster_results/sp_model_num_var(threshold 1).RDS")
results_nocon = as.data.frame(matrix(unlist(lapply(results, function(l) l$nocon)), nrow = 500, byrow = TRUE))
colnames(results_nocon) <- colnames(results[[1]]$nocon)
result_stats(results_nocon, 1 + 1 * 0.5 + 2 * 0.5)
# 2.524 0.024 0.452 0.441 0.934
# 2.546 0.046 0.465 0.442 0.928

results_inf = as.data.frame(matrix(unlist(lapply(results, function(l) l$inf)), nrow = 500, byrow = TRUE))
colnames(results_inf) <- colnames(results[[1]]$inf)
result_stats(results_inf, 1 + 1 * 0.5 + 2 * 0.5)
# 2.512 0.012 0.272 0.265 0.950
# 2.524 0.024 0.290 0.287 0.942

results_sp = as.data.frame(matrix(unlist(lapply(results, function(l) l$sp)), nrow = 500, byrow = TRUE))
colnames(results_sp) <- colnames(results[[1]]$sp)
# result_stats(results_sp,  1 + 1 * 0.1 + 2 * 0.2)
# 1.497 -0.003  0.459  0.504  0.964
result_stats(results_sp,  1 + 1 * 1 + 2 * 1)
# 4.015 0.015 0.759 0.785 0.944

results_mixed = as.data.frame(matrix(unlist(lapply(results, function(l) l$mixed)), nrow = 500, byrow = TRUE))
colnames(results_mixed) <- colnames(results[[1]]$mixed)
# result_stats(results_mixed,  1 + 1 * 0.1 + 2 * 0.2)
# 1.506 0.006 0.504 0.540 0.970
result_stats(results_mixed,  1 + 1 * 1 + 2 * 1)
#  4.051 0.051 0.855 0.892 0.958

sp_num = rbind(result_stats(results_nocon, 1 + 1 * 0.5 + 2 * 0.5),
                result_stats(results_inf,  1 + 1 * 0.5 + 2 * 0.5),
                result_stats(results_sp,  1 + 1 * 1 + 2 * 1),
                result_stats(results_mixed, 1 + 1 * 1 + 2 * 1))

#### 2.3 both vars ####
# sp_model_both_var
results = readRDS("simulation_results/cluster_results/sp_model_both_var(threshold 1).RDS")
results_nocon = as.data.frame(matrix(unlist(lapply(results, function(l) l$nocon)), nrow = 500, byrow = TRUE))
colnames(results_nocon) <- colnames(results[[1]]$nocon)
result_stats(results_nocon, 1 + 1 * 0.5 + 2 * 0.5)
# 2.555 0.055 0.319 0.317 0.958
# 2.517 0.017 0.419 0.420 0.950

results_inf = as.data.frame(matrix(unlist(lapply(results, function(l) l$inf)), nrow = 500, byrow = TRUE))
colnames(results_inf) <- colnames(results[[1]]$inf)
result_stats(results_inf, 1 + 1 * 0.5 + 2 * 0.5)
# 2.508 0.008 0.720 0.695 0.938
# 2.497 -0.003  0.350  0.317  0.896

results_sp = as.data.frame(matrix(unlist(lapply(results, function(l) l$sp)), nrow = 500, byrow = TRUE))
colnames(results_sp) <- colnames(results[[1]]$sp)
# result_stats(results_sp,  1 + 1 * 1 + 2 * 0.1)
# 2.167 -0.033  3.180  3.366  0.964
result_stats(results_sp,  1 + 1 * 1 + 2 * 1)
# 4.018 0.018 0.492 0.546 0.962

results_mixed = as.data.frame(matrix(unlist(lapply(results, function(l) l$mixed)), nrow = 500, byrow = TRUE))
colnames(results_mixed) <- colnames(results[[1]]$mixed)
# result_stats(results_mixed, 1 + 1 * 1 + 2 * 0.1)
# 2.188 -0.012 15.784 15.238  0.946
result_stats(results_mixed, 1 + 1 * 1 + 2 * 1)
# 4.017 0.017 0.716 0.681 0.914

sp_both = rbind(result_stats(results_nocon, 1 + 1 * 0.5 + 2 * 0.5),
                 result_stats(results_inf,  1 + 1 * 0.5 + 2 * 0.5),
                 result_stats(results_sp,  1 + 1 * 1 + 2 * 1),
                 result_stats(results_mixed, 1 + 1 + 2 * 1))

#####

sp_table = as.data.frame(rbind(rep(c("p.est", "bias", "MC SE", "analytical SE", "CP"), 3),
                                cbind(sp_cat, sp_num, sp_both)))

rownames(sp_table) <- c("", "Avg", "Het Infl", "Het Sp", "Het Mixed")

colnames(sp_table) <- c("2 binary vars", rep(" ", 4),
                         "2 cont vars", rep(" ", 4),
                         "1 binary & 1 cont vars", rep(" ", 4))

###### 3. mixed outcome model ######
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
results = readRDS("simulation_results/cluster_results/mixed_model_num_var(threshold 1).RDS")
results_nocon = as.data.frame(matrix(unlist(lapply(results, function(l) l$nocon)), nrow = 500, byrow = TRUE))
colnames(results_nocon) <- colnames(results[[1]]$nocon)
result_stats(results_nocon, 1 + 1 * 0.5 + 2 * 0.5 + 3 * 0.5 + 4 * 0.5)
# 6.083 0.083 0.968 1.010 0.952
# 6.016 0.016 1.038 1.009 0.928

results_inf = as.data.frame(matrix(unlist(lapply(results, function(l) l$inf)), nrow = 500, byrow = TRUE))
colnames(results_inf) <- colnames(results[[1]]$inf)
result_stats(results_inf, 1 + 1 * 0.1 + 2 * 0.2 + 3 * 0.5 + 4 * 0.5)
# 5.012 0.012 0.687 0.631 0.942
result_stats(results_inf, 1 + 1 * 1 + 2 * 1 + 3 * 0.5 + 4 * 0.5)
# 7.519 0.019 0.771 0.725 0.944

results_sp = as.data.frame(matrix(unlist(lapply(results, function(l) l$sp)), nrow = 500, byrow = TRUE))
colnames(results_sp) <- colnames(results[[1]]$sp)
# result_stats(results_sp,  1 + 1 * 0.5 + 2 * 0.5 + 3 * 0.1 + 4 * 0.2)
# 3.685 0.085 1.322 1.296 0.940
result_stats(results_sp,  1 + 1 * 0.5 + 2 * 0.5 + 3 * 1 + 4 * 1)
# 9.673 0.173 1.877 2.019 0.958

results_mixed = as.data.frame(matrix(unlist(lapply(results, function(l) l$mixed)), nrow = 500, byrow = TRUE))
colnames(results_mixed) <- colnames(results[[1]]$mixed)
# result_stats(results_mixed,  1 + 1 * 0.1 + 2 * 0.2 + 3 * 0.1 + 4 * 0.2)
# 2.724 0.124 1.388 1.379 0.960
result_stats(results_mixed,  1 + 1 * 1 + 2 * 1 + 3 * 1 + 4 * 1)
# 11.131  0.131  2.318  2.403  0.948

mixed_num = rbind( result_stats(results_nocon, 1 + 1 * 0.5 + 2 * 0.5 + 3 * 0.5 + 4 * 0.5),
                   result_stats(results_inf,  1 + 1 * 1 + 2 * 1 + 3 * 0.5 + 4 * 0.5),
                   result_stats(results_sp,  1 + 1 * 0.5 + 2 * 0.5 + 3 * 1 + 4 * 1),
                   result_stats(results_mixed, 1 + 1 * 1 + 2 * 1 + 3 * 1 + 4 * 1))

#### 3.3 both vars ####
# mixed_model_both_var
# results = rbind( readRDS("simulation_results/cluster_results/mixed_model_both_var.RDS"),
#                  readRDS("simulation_results/cluster_results/mixed_model_both_var copy.RDS"))
results = readRDS("simulation_results/cluster_results/mixed_model_both_var(threshold 1).RDS")
results_nocon = as.data.frame(matrix(unlist(lapply(results, function(l) l$nocon)), nrow = 500, byrow = TRUE))
colnames(results_nocon) <- colnames(results[[1]]$nocon)
result_stats(results_nocon, 1 + 1 * 0.5 + 2 * 0.5 + 3 * 0.5 + 4 * 0.5)
# 6.043 0.043 0.766 0.761 0.941
# 6.088 0.088 0.994 0.949 0.922

results_inf = as.data.frame(matrix(unlist(lapply(results, function(l) l$inf)), nrow = 500, byrow = TRUE))
colnames(results_inf) <- colnames(results[[1]]$inf)
# result_stats(results_inf, 1 + 1 + 2 * 0.1 + 3 * 0.5 + 4 * 0.5)
# 5.698 -0.002  1.824  1.845  0.940
result_stats(results_inf, 1 + 1 + 2 * 1 + 3 * 0.5 + 4 * 0.5)
# 7.569 0.069 0.937 0.758 0.892

results_sp = as.data.frame(matrix(unlist(lapply(results, function(l) l$sp)), nrow = 500, byrow = TRUE))
colnames(results_sp) <- colnames(results[[1]]$sp)
# result_stats(results_sp,  1 + 1 * 0.5 + 2 * 0.5 + 3 * 1 + 4 * 0.1)
# 6.125 0.225 7.929 8.547 0.958
result_stats(results_sp,  1 + 1 * 0.5 + 2 * 0.5 + 3 * 1 + 4 * 1)
# 9.583 0.083 1.303 1.376 0.956

results_mixed = as.data.frame(matrix(unlist(lapply(results, function(l) l$mixed)), nrow = 500, byrow = TRUE))
colnames(results_mixed) <- colnames(results[[1]]$mixed)
# result_stats(results_mixed, 1 + 1 * 1 + 2 * 0.1 + 3 * 1 + 4 * 0.1)
# 6.606 1.006 40.481 38.786  0.936
result_stats(results_mixed, 1 + 1 * 1 + 2 * 1 + 3 * 1 + 4 * 1)
# 11.050  0.050  2.003  1.850  0.904

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
#colnames(table) <- table_names

library(gridExtra)
png("test(1).png", height = 50*nrow(table), width = 150*ncol(table))
grid.table(table)
dev.off()
