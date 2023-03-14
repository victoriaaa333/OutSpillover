source("results/results_utils.R")

# Avg, Het infl, Het Sp and Het infl+Sp

#### Spillover outcome model ####

# 1. one binary covariate

# 2. one continuous covariate

# 3. one binary and one continuous

#### Influencer outcome model ####

# 1. one binary covariate
inf_binary_nocon =  rbind(readRDS("inf_results/small_coef/inf_binary_nocon.RDS"),
                        readRDS("inf_results/small_coef/inf_binary_nocon (1).RDS"))
result_stats(inf_binary_nocon, 3.5) # 3 + 1*0.5

inf_binary_inf =  rbind(readRDS("inf_results/small_coef/inf_binary_inf.RDS"),
                      readRDS("inf_results/small_coef/inf_binary_inf (1).RDS"))
result_stats(inf_binary_inf, 3 + 1*1)

inf_binary_sp =  rbind(readRDS("inf_results/small_coef/inf_binary_sp.RDS"),
                     readRDS("inf_results/small_coef/inf_binary_sp (1).RDS"))
result_stats(inf_binary_sp, 3.5)


# 2. one continuous covariate
inf_cont_nocon =  rbind(readRDS("inf_results/small_coef/inf_cont_nocon.RDS"),
                         readRDS("inf_results/small_coef/inf_cont_nocon (1).RDS"))
result_stats(inf_cont_nocon, 3.5) # 3 + 1*0.5

inf_cont_inf =  rbind(readRDS("inf_results/small_coef/inf_cont_inf.RDS"),
                       readRDS("inf_results/small_coef/inf_cont_inf (1).RDS"))
result_stats(inf_cont_inf, 3 + 1*0.1)

inf_cont_sp =  rbind(readRDS("inf_results/small_coef/inf_cont_sp.RDS"),
                      readRDS("inf_results/small_coef/inf_cont_sp (1).RDS"))
result_stats(inf_cont_sp, 3.5)

# 3. one binary and one continuous
# coef: 3, 1, 2
inf_mixed_nocon =  rbind(readRDS("inf_results/small_coef/inf_mixed_nocon.RDS"), 
                         readRDS("inf_results/small_coef/inf_mixed_nocon(1).RDS"), 
                         readRDS("inf_results/small_coef/inf_mixed_nocon(2).RDS"))
result_stats(inf_mixed_nocon, 4.5) # 3 + 1*0.5 + 2*0.5
# 4.458 -0.042  0.549  0.543  0.932

inf_mixed_inf =  rbind(readRDS("inf_results/small_coef/inf_mixed_inf.RDS"), 
                       readRDS("inf_results/small_coef/inf_mixed_inf(1).RDS"), 
                       readRDS("inf_results/small_coef/inf_mixed_inf(2).RDS"))
result_stats(inf_mixed_inf, 3 + 1*1 + 2*0.1)
# 4.254 0.054 0.456 0.399 0.921

inf_mixed_sp =  rbind(readRDS("inf_results/small_coef/inf_mixed_sp.RDS"), 
                       readRDS("inf_results/small_coef/inf_mixed_sp(1).RDS"), 
                       readRDS("inf_results/small_coef/inf_mixed_sp(2).RDS"))
result_stats(inf_mixed_sp, 3 + 1*0.5 + 2*0.5)
#  4.571 0.071 0.675 0.682 0.947
