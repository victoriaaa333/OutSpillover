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
results =  readRDS("simulation_results/cluster_results/sd = 0.5/inf_model_both_var(sd 0.5).RDS")
results_nocon = as.data.frame(matrix(unlist(lapply(results, function(l) l$nocon)), nrow = 500, byrow = TRUE))
colnames(results_nocon) <- colnames(results[[1]]$nocon)
result_stats(results_nocon, 1 + 1 * 0.5 + 2 * 0.5)
# SD = 1, 2.531 0.031 0.317 0.312 0.950
# SD = 0.5, 2.495 -0.005  0.313  0.303  0.934

results_inf = as.data.frame(matrix(unlist(lapply(results, function(l) l$inf)), nrow = 500, byrow = TRUE))
colnames(results_inf) <- colnames(results[[1]]$inf)
result_stats(results_inf, 1 + 1 + 2 * 1)
# SD = 1, 3.980 -0.020  0.369  0.274  0.860
# SD = 0.5,  4.006 0.006 0.349 0.308 0.904

results_sp = as.data.frame(matrix(unlist(lapply(results, function(l) l$sp)), nrow = 500, byrow = TRUE))
colnames(results_sp) <- colnames(results[[1]]$sp)
result_stats(results_sp,  1 + 1 * 0.5 + 2 * 0.5)
# SD = 1, 2.517 0.017 0.569 0.525 0.922
# SD = 0.5, 2.570 0.070 0.856 0.821 0.940

results_mixed = as.data.frame(matrix(unlist(lapply(results, function(l) l$mixed)), nrow = 500, byrow = TRUE))
colnames(results_mixed) <- colnames(results[[1]]$mixed)
result_stats(results_mixed, 1 + 1 + 2 * 1)
# SD = 1, 3.977 -0.023  0.945  0.814  0.892
# SD = 0.5, 4.123 0.123 1.685 1.531 0.926

##### 2. spillover outcome model ######
# sp_model_both_var
results = readRDS("~/Downloads/sp_model_both_var(sd 0.5).RDS")
results_nocon = as.data.frame(matrix(unlist(lapply(results, function(l) l$nocon)), nrow = 500, byrow = TRUE))
colnames(results_nocon) <- colnames(results[[1]]$nocon)
result_stats(results_nocon, 1 + 1 * 0.5 + 2 * 0.5)
# SD = 1, 2.516 0.016 0.421 0.420 0.934
# SD = 0.5, 2.537 0.037 0.332 0.345 0.954

results_inf = as.data.frame(matrix(unlist(lapply(results, function(l) l$inf)), nrow = 500, byrow = TRUE))
colnames(results_inf) <- colnames(results[[1]]$inf)
result_stats(results_inf, 1 + 1 * 0.5 + 2 * 0.5)
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
# SD = 1, 4.010 0.010 0.694 0.683 0.940
# SD = 0.5, 4.086 0.086 1.506 1.444 0.944

##### 1. mixed outcome model ######
# mixed_model_both_var
results = readRDS("~/Downloads/mixed_model_both_var(sd 0.5).RDS")
results_nocon = as.data.frame(matrix(unlist(lapply(results, function(l) l$nocon)), nrow = 500, byrow = TRUE))
colnames(results_nocon) <- colnames(results[[1]]$nocon)
result_stats(results_nocon, 1 + 1 * 0.5 + 2 * 0.5 + 3 * 0.5 + 4 * 0.5)
# SD = 1, 6.062 0.062 1.013 0.941 0.924
# SD = 0.5, 6.151 0.151 0.897 0.821 0.924

results_inf = as.data.frame(matrix(unlist(lapply(results, function(l) l$inf)), nrow = 500, byrow = TRUE))
colnames(results_inf) <- colnames(results[[1]]$inf)
result_stats(results_inf,  1 + 1 + 2 * 1 + 3 * 0.5 + 4 * 0.5)
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
# SD = 1, 10.914 -0.086  1.876  1.833  0.926
# SD = 0.5, 10.766 -0.234  3.900  3.764  0.934
