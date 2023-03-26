source("results/results_utils.R")

#### Spillover outcome model ####
# nsim = 500
# 1. regular HT estimator
spillover_plain = readRDS("table_results/spillover_plain.RDS")
result_stats(spillover_plain, 13)

# 2. 2 binary variables
spillover_binary1 = readRDS("table_results/spillover_binary1.RDS")
result_stats(spillover_binary1, 21)

spillover_binary3 = readRDS("table_results/spillover_binary3.RDS")
result_stats(spillover_binary3, 13)

# 3. 2 continuous variables
spillover_cont1_2var = readRDS("table_results/spillover_cont1_2var.RDS")
result_stats(spillover_cont1_2var, 7.5)

spillover_cont2_2var = readRDS("table_results/spillover_cont2_2var.RDS")
result_stats(spillover_cont2_2var, 13)

# 4. 1 binary variable
spillover_binary1_1var = readRDS("table_results/spillover_binary1_1var.RDS")
result_stats(spillover_binary1_1var, 12)

spillover_binary3_1var = readRDS("table_results/spillover_binary3_1var.RDS")
result_stats(spillover_binary3_1var, 8.5)

# 5. 1 continuous variable
spillover_cont1_1var = readRDS("table_results/spillover_cont1_1var.RDS")
result_stats(spillover_cont1_1var, 5.7)

spillover_cont2_1var = readRDS("table_results/spillover_cont2_1var.RDS")
result_stats(spillover_cont2_1var, 8.5)

# 6. 1 categorical and 1 numeral variables
spillover_mixed1 =  readRDS("table_results/spillover_mixed1.RDS")
result_stats(spillover_mixed1, 14.7)

spillover_mixed2 = readRDS("table_results/spillover_mixed2.RDS")
result_stats(spillover_mixed2, 13)

#### Influencer outcome model ####
# 1. regular HT estimator
inf_plain = readRDS("table_results/inf_plain.RDS")
result_stats(inf_plain, 13)

# 2. 2 binary variables
inf_binary1 = readRDS("table_results/inf_binary1.RDS")
result_stats(inf_binary1, 13)

inf_binary3 = readRDS("table_results/inf_binary3.RDS")
result_stats(inf_binary3, 21)

# 3. 2 continuous variables
inf_cont1_2var = readRDS("table_results/inf_cont1.RDS")
result_stats(inf_cont1_2var, 13)

inf_cont2_2var = readRDS("table_results/inf_cont2.RDS")
result_stats(inf_cont2_2var, 7.5)

# 4. 1 binary variable
inf_binary1_1var = readRDS("table_results/inf_binary1_1var.RDS")
result_stats(inf_binary1_1var, 8.5)

inf_binary3_1var = readRDS("table_results/inf_binary3_1var.RDS")
result_stats(inf_binary3_1var, 12)


# 5. 1 continuous variable
inf_cont1_1var = readRDS("table_results/inf_cont1_1var.RDS")
result_stats(inf_cont1_1var, 8.5)

inf_cont2_1var = readRDS("table_results/inf_cont2_1var.RDS")
result_stats(inf_cont2_1var, 5.7)

# 6. 1 categorical and 1 numeral variables
inf_mixed1 = readRDS("table_results/inf_mixed1.RDS")
result_stats(inf_mixed1, 13)

inf_mixed2 = readRDS("table_results/inf_mixed2.RDS")
result_stats(inf_mixed2, 14.7)

### Mixed outcome model (point estimate only for now) ###
mixed_mixed = readRDS("table_results/mixed_mixed.RDS")
bias = -mean(mixed_mixed)-26.7
round(bias,3)

# first row is influencer effect, second is spillover
spillover_table = rbind(rep(c("p.est", "bias", "MC SE", "analytical SE", "CP"),3),
                        c(rep(NA, 5), rep(NA,5), result_stats(spillover_plain, 13)),
                        c(result_stats(spillover_binary3_1var, 8.5), 
                          result_stats(spillover_cont2_1var, 8.5),
                          result_stats(spillover_mixed2, 13)),
                        c(result_stats(spillover_binary1_1var, 12),
                          result_stats(spillover_cont1_1var, 5.7),
                          result_stats(spillover_mixed1, 14.7)))

inf_table = rbind(rep(c("p.est", "bias", "MC SE", "analytical SE", "CP"),3),
                  c(rep(NA, 5), rep(NA,5), result_stats(inf_plain, 13)),
                  c(result_stats(inf_binary3_1var, 12),
                    result_stats(inf_cont2_1var, 5.7),
                    result_stats(inf_mixed2, 14.7)),
                  c(result_stats(inf_binary1_1var, 8.5),
                    result_stats(inf_cont1_1var, 8.5),
                    result_stats(inf_mixed1, 13)))

#colnames(spillover_table) <- rep(c("p.est", "bias", "MC SE", "analytical SE", "CP"),3)
rownames(spillover_table) <- c(" ", "Avg", "Het Infl", "Het Sp")

#colnames(inf_table) <- rep(c("p.est", "bias", "MC SE", "analytical SE", "CP"),3)
rownames(inf_table) <- c(" ", "Avg", "Het Infl", "Het Sp")

colnames(spillover_table) <- c("binary", rep(" ", 4),
                               "continuous", rep(" ", 4),
                               "mixed", rep(" ", 4))
colnames(inf_table) <- c("binary", rep(" ", 4),
                               "continuous", rep(" ", 4),
                               "mixed", rep(" ", 4))

print(spillover_table, quote=FALSE)
print(inf_table, quote=FALSE)

