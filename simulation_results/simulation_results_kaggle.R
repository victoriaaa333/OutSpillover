# function for true cp calculation
coverage_prob <- function(dataset, true_effect){
  #CI_lower <- dataset$estimate - 1.96 * dataset$std.error
  #CI_upper <- dataset$estimate + 1.96 * dataset$std.error
  coverage <- intersect(which(dataset$conf.low <= true_effect), 
                        which(dataset$conf.high >= true_effect))
  return(length(coverage)/length(dataset$estimate))
}

coverage_prob0 <- function(dataset, true_effect, std){
  coverage <- intersect(which(dataset$estimate - 1.96 * std
                              <= true_effect), 
                        which(dataset$estimate + 1.96 * std
                              >= true_effect))
  return(length(coverage)/length(dataset$estimate))
}

#### Spillover Model Simulations ####
## No-condition ##
spillover_plain = readRDS("spillover_results/spillover_plain.RDS")
sd(spillover_plain$estimate) # 2.022351
mean(spillover_plain$std.error) # 1.972029
mean(spillover_plain$estimate) # -13.15786, bias = 0.15786
coverage_prob(spillover_plain, -13) #0.934
#coverage_prob0(spillover_plain, -13, sd(spillover_plain$estimate)) # 0.945

## Binary varibles ##
spillover_binary1 =  readRDS("spillover_results/spillover_binary1.RDS")
sd(spillover_binary1$estimate) # 5.43128
mean(spillover_binary1$std.error) # 5.196488
mean(spillover_binary1$estimate) # -21.03914, bias = -0.03914
coverage_prob(spillover_binary1, -21) #0.901
coverage_prob0(spillover_binary1, -21, sd(spillover_binary1$estimate)) # 0.956
coverage_prob0(spillover_binary1, -21, mean(spillover_binary1$std.error)) # 0.948

# Binary variable under group conditions, overall
spillover_binary2 =  readRDS("spillover_results/spillover_binary2.RDS")
sd(spillover_binary2$estimate) # 3.371333
mean(spillover_binary2$std.error) # 4.587584
mean(spillover_binary2$estimate) # -13.3068, bias = -0.3068
coverage_prob(spillover_binary2, -13) #0.951

# Binary variable under group conditions, w/i groups
spillover_binary3 =  readRDS("spillover_results/spillover_binary3.RDS")
sd(spillover_binary3$estimate) # 4.701511
mean(spillover_binary3$std.error) # 4.500119
mean(spillover_binary3$estimate) # -12.9894, bias = 0.0106
coverage_prob(spillover_binary3, -13) #0.878
# We should use this version when no regression involve?
coverage_prob0(spillover_binary3, -13, sd(spillover_binary3$estimate)) #0.962

## Continuous variables ##
spillover_cont1 = readRDS("spillover_results/spillover_cont1.RDS")
sd(spillover_cont1$estimate) # 5.070049
mean(spillover_cont1$std.error) # 4.899444
mean(spillover_cont1$estimate) # -7.613574, bias = -0.1357
coverage_prob(spillover_cont1, -7.5) #0.928

# Continuous variables under group conditions, overall
spillover_cont2 = readRDS("spillover_results/spillover_cont2.RDS")
sd(spillover_cont2$estimate) # 2.752394
mean(spillover_cont2$std.error) # 2.581217
mean(spillover_cont2$estimate) # -13.24717, bias = -0.24717
coverage_prob(spillover_cont2, -13) #0.922

# No group conditions because it's a regression estimator

## Mixed variables ##
spillover_mixed1 = readRDS("spillover_results/spillover_mixed1.RDS")
sd(spillover_mixed1$estimate) # 3.782947
mean(spillover_mixed1$std.error) # 3.557223
mean(spillover_mixed1$estimate) # -7.816454, bias = -0.316454
coverage_prob(spillover_mixed1, -7.5) #0.924

# Mixed variables under group conditions, overall
spillover_mixed2 = readRDS("spillover_results/spillover_mixed2.RDS")
sd(spillover_mixed2$estimate) # 3.852762
mean(spillover_mixed2$std.error) # 3.253912
mean(spillover_mixed2$estimate) # -13.22409, bias = -0.22409
coverage_prob(spillover_mixed2, -13) #0.883

#### Influencer Model Simulations ####
## No-condition ##
inf_plain = readRDS("inf_results/inf_plain.RDS")
sd(inf_plain$estimate) # 3.309021
mean(inf_plain$std.error) # 3.150794
mean(inf_plain$estimate) # -13.23253, bias = -0.23253
coverage_prob(inf_plain, -13) #0.907

inf_plain = readRDS("inf_results/inf_plain_cluster200.RDS")
sd(inf_plain$estimate) # 1.655949
mean(inf_plain$std.error) # 1.612916
mean(inf_plain$estimate) # -13.13284, bias = -0.23253
coverage_prob(inf_plain, -13) # 0.946

## Binary varibles ##
inf_binary1 =  readRDS("inf_results/inf_binary1.RDS")
sd(inf_binary1$estimate) # 3.424664
mean(inf_binary1$std.error) # 3.262785
mean(inf_binary1$estimate) # -12.99295, bias = 0.001
coverage_prob(inf_binary1, -13) #0.909

# Binary variable under group conditions, overall
inf_binary2 =  readRDS("inf_results/inf_binary2.RDS")
sd(inf_binary2$estimate) # 3.70819
mean(inf_binary2$std.error) # 5.136765
mean(inf_binary2$estimate) # -20.99706, bias = 0.003
coverage_prob(inf_binary2, -21) #0.974

# Binary variable under group conditions, w/i groups
inf_binary3 =  readRDS("inf_results/inf_binary3.RDS")
sd(inf_binary3$estimate) # 5.211216
mean(inf_binary3$std.error) # 5.021085
mean(inf_binary3$estimate) # -20.99063, bias = 0.01
coverage_prob(inf_binary3, -21) #0.917

## Continuous variables ##
inf_cont1 = readRDS("inf_results/inf_cont1.RDS")
sd(inf_cont1$estimate) # 7.82767
mean(inf_cont1$std.error) # 6.812539
mean(inf_cont1$estimate) # -13.368, bias = -0.368
coverage_prob(inf_cont1, -13) #0.885

# Continuous variables under group conditions, overall
inf_cont2 = readRDS("inf_results/inf_cont2.RDS")
sd(inf_cont2$estimate) # 2.391088
mean(inf_cont2$std.error) # 2.154077
mean(inf_cont2$estimate) # -7.903719, bias = -0.4037
coverage_prob(inf_cont2, -7.5) #0.935

# No group conditions because it's a regression estimator

## Mixed variables ##
inf_mixed1 = readRDS("inf_results/inf_mixed1(2).RDS")
sd(inf_mixed1$estimate) # 5.732741
mean(inf_mixed1$std.error) # 5.079389
mean(inf_mixed1$estimate) # -13.3825, bias = -0.3825
coverage_prob(inf_mixed1, -13) #0.904

# Mixed variables under group conditions, overall
inf_mixed2 = readRDS("inf_results/inf_mixed2(2).RDS")
sd(inf_mixed2$estimate) # 3.493663
mean(inf_mixed2$std.error) # 2.454987
mean(inf_mixed2$estimate) # -7.984672, bias = -0.4846
coverage_prob(inf_mixed2, -7.5) #0.844

## TODO:
# 1. bootstrapped variance for each 1000 simulation
# 2. increase the number of clusters to 200.
# 3. write down the table

#### Increase number of clusters to 200 ####
#### Spillover Model Simulations ####
## No-condition ##
# nsim = 600
spillover_plain = rbind(readRDS("spillover_results/num_of_clusters_200/spillover_plain.RDS"),
                        readRDS("spillover_results/num_of_clusters_200/spillover_plain(1).RDS"))

sd(spillover_plain$estimate) # 2.168015
mean(spillover_plain$std.error) # 2.01858
mean(spillover_plain$estimate) # -13.21375, bias = 0.21375
coverage_prob(spillover_plain, -13) #0.9266667

## Binary varibles ##
spillover_binary1 =  rbind(readRDS("spillover_results/num_of_clusters_200/spillover_binary1.RDS"),
                           readRDS("spillover_results/num_of_clusters_200/spillover_binary1(1).RDS"))
sd(spillover_binary1$estimate) # 2.778729
mean(spillover_binary1$std.error) # 2.702652
mean(spillover_binary1$estimate) # -21.16166, bias = 0.16166
coverage_prob(spillover_binary1, -21) #0.94

# Binary variable under group conditions, overall
# spillover_binary2 =  rbind(readRDS("spillover_results/num_of_clusters_200/spillover_binary2.RDS"),
#                            readRDS("spillover_results/num_of_clusters_200/spillover_binary2(1).RDS"))
# sd(spillover_binary2$estimate) # 3.371333
# mean(spillover_binary2$std.error) # 4.587584
# mean(spillover_binary2$estimate) # -13.3068, bias = -0.3068
# coverage_prob(spillover_binary2, -13) #0.951

# Binary variable under group conditions, w/i groups
spillover_binary3 =  rbind(readRDS("spillover_results/num_of_clusters_200/spillover_binary3.RDS"),
                           readRDS("spillover_results/num_of_clusters_200/spillover_binary3(1).RDS"))
sd(spillover_binary3$estimate) # 2.501396
mean(spillover_binary3$std.error) # 2.753614
mean(spillover_binary3$estimate) # -13.18557, bias = 0.0106
coverage_prob(spillover_binary3, -13) #0.924
# We should use this version when no regression involve?
coverage_prob0(spillover_binary3, -13, sd(spillover_binary3$estimate)) #0.95

## Continuous variables ##
spillover_cont1 = readRDS("spillover_results/num_of_clusters_200/spillover_cont1.RDS")
sd(spillover_cont1$estimate) # 2.677921
mean(spillover_cont1$std.error) # 2.628133
mean(spillover_cont1$estimate) # -7.720052, bias = -0.220052
coverage_prob(spillover_cont1, -7.5) #0.944

# Continuous variables under group conditions, overall
spillover_cont2 = readRDS("spillover_results/num_of_clusters_200/spillover_cont2.RDS")
sd(spillover_cont2$estimate) # 1.454335
mean(spillover_cont2$std.error) # 1.367367
mean(spillover_cont2$estimate) # -13.11965, bias = -0.211965
coverage_prob(spillover_cont2, -13) #0.938
# No group conditions because it's a regression estimator

## Mixed variables (2 num covariates in interaction) ##
spillover_mixed1 = readRDS("spillover_results/num_of_clusters_200/spillover_mixed1.RDS")
sd(spillover_mixed1$estimate) # 1.861457
mean(spillover_mixed1$std.error) # 1.893961
mean(spillover_mixed1$estimate) # -7.486774, bias = -0.014
coverage_prob(spillover_mixed1, -7.5) #0.95

# Mixed variables under group conditions, overall
spillover_mixed2 = readRDS("spillover_results/num_of_clusters_200/spillover_mixed2.RDS")
sd(spillover_mixed2$estimate) # 1.977831
mean(spillover_mixed2$std.error) # 1.708772
mean(spillover_mixed2$estimate) # -13.12521, bias = -0.12521
coverage_prob(spillover_mixed2, -13) #0.9

## Mixed variables (1 num and 1 cat covariate in interaction) ##
spillover_mixed1 = rbind(readRDS("spillover_results/1cat1num/spillover_mixed1.RDS"),
                         readRDS("spillover_results/1cat1num/spillover_mixed1(1).RDS"))
sd(spillover_mixed1$estimate) # 1.913753
mean(spillover_mixed1$std.error) # 2.037661
mean(spillover_mixed1$estimate) # -14.85719, bias = -
coverage_prob(spillover_mixed1, -14.7) #0.955

# Mixed variables under group conditions, overall
spillover_mixed2 = rbind(readRDS("spillover_results/1cat1num/spillover_mixed2.RDS"),
                         readRDS("spillover_results/1cat1num/spillover_mixed2(1).RDS"),
                         readRDS("spillover_results/1cat1num/spillover_mixed2(2).RDS"))
sd(spillover_mixed2$estimate) # 1.562126
mean(spillover_mixed2$std.error) # 1.359359
mean(spillover_mixed2$estimate) # -13.04766, bias = -
coverage_prob(spillover_mixed2, -13) #0.9188889

#### Influencer Model Simulations ####
# All the cases below are from cluster 200, sample size 100, nsim = 200
inf_table = c()
## No-condition ##
inf_plain = readRDS("inf_results/num_of_clusters_200/inf_plain.RDS")
sd(inf_plain$estimate) # 1.655949
mean(inf_plain$std.error) # 1.612916
mean(inf_plain$estimate) # -13.13284, bias = -0.13284
coverage_prob(inf_plain, -13) #0.946

## Binary varibles ##
inf_binary1 =  rbind(readRDS("inf_results/num_of_clusters_200/inf_binary1.RDS"),
                     readRDS("inf_results/num_of_clusters_200/inf_binary1(1).RDS"))
sd(inf_binary1$estimate) # 1.701592
mean(inf_binary1$std.error) # 1.712654
mean(inf_binary1$estimate) # -13.16731, bias = 0.16731
coverage_prob(inf_binary1, -13) #0.948

# Binary variable under group conditions, overall
# inf_binary2 =  rbind(readRDS("inf_results/num_of_clusters_200/inf_binary2.RDS"),
#                      readRDS("inf_results/num_of_clusters_200/inf_binary2(1).RDS"))
# sd(inf_binary2$estimate) # 1.768858
# mean(inf_binary2$std.error) # 3.439254
# mean(inf_binary2$estimate) # -20.9142, bias = -0.0948
# coverage_prob(inf_binary2, -21) #0.982

# Binary variable under group conditions, w/i groups
inf_binary3 =  rbind(readRDS("inf_results/num_of_clusters_200/inf_binary3.RDS"),
                     readRDS("inf_results/num_of_clusters_200/inf_binary3(1).RDS"))
sd(inf_binary3$estimate) # 2.54659
mean(inf_binary3$std.error) # 3.427047
mean(inf_binary3$estimate) # -21.19562, bias = -0.06475
coverage_prob(inf_binary3, -21) #0.946

## Continuous variables ##
inf_cont1 = rbind(readRDS("inf_results/num_of_clusters_200/inf_cont1.RDS"),
                  readRDS("inf_results/num_of_clusters_200/inf_cont1(1).RDS"))
sd(inf_cont1$estimate) # 4.061942
mean(inf_cont1$std.error) # 3.790154
mean(inf_cont1$estimate) # -13.27164, bias = 
coverage_prob(inf_cont1, -13) #0.922

# Continuous variables under group conditions, overall
inf_cont2 = rbind(readRDS("inf_results/num_of_clusters_200/inf_cont2.RDS"),
                  readRDS("inf_results/num_of_clusters_200/inf_cont2(1).RDS"))
sd(inf_cont2$estimate) # 1.119922
mean(inf_cont2$std.error) # 1.140617
mean(inf_cont2$estimate) # -7.59409, bias = -
coverage_prob(inf_cont2, -7.5) #0.954
# No group conditions because it's a regression estimator



#### Mixed variables (noc = 200, sample size = 200) ####
inf_mixed1 = rbind(readRDS("inf_results/num_of_clusters_200/inf_mixed1.RDS"),
                   readRDS("inf_results/num_of_clusters_200/inf_mixed1(2).RDS"))
sd(inf_mixed1$estimate) # 2.978993
mean(inf_mixed1$std.error) # 2.819687
mean(inf_mixed1$estimate) # -13.09401, bias = -0.09401
coverage_prob(inf_mixed1, -13) #0.93

# Mixed variables under group conditions, overall
inf_mixed2 = rbind(readRDS("inf_results/mixed/inf_mixed2(new var_outcome).RDS"),
                   readRDS("inf_results/mixed/inf_mixed2(new var_outcome)1.RDS"))
sd(inf_mixed2$estimate) # 1.792545
mean(inf_mixed2$std.error) # 1.311842
mean(inf_mixed2$estimate) # -7.449184, bias = 0.050816
coverage_prob(inf_mixed2, -7.5) #.8475



#### Mixed variables (1 cat 1 num) ####
inf_mixed1 = rbind(readRDS("inf_results/1cat1num/inf_mixed1.RDS"),
                   readRDS("inf_results/1cat1num/inf_mixed1(1).RDS"))
sd(inf_mixed1$estimate) # 2.220289
mean(inf_mixed1$std.error) # 2.107123
mean(inf_mixed1$estimate) # -13.30752, bias = 
coverage_prob(inf_mixed1, -13) #.95

inf_mixed2 = rbind(readRDS("inf_results/1cat1num/inf_mixed2.RDS"),
                   readRDS("inf_results/1cat1num/inf_mixed2(1).RDS"))
sd(inf_mixed2$estimate) # 1.517736
mean(inf_mixed2$std.error) # 1.368048
mean(inf_mixed2$estimate) # -14.8041, bias = 
# point estimate = 5 + 7*0.1 + 9 = 14.7
coverage_prob(inf_mixed2, -14.7) #.922

spillover_mixed1 = rbind(readRDS("spillover_results/1cat1num/spillover_mixed1.RDS"),
                   readRDS("spillover_results/1cat1num/spillover_mixed1(1).RDS"))
sd(spillover_mixed1$estimate) # 2.141295
mean(spillover_mixed1$std.error) # 2.209299
mean(spillover_mixed1$estimate) # -14.85719, bias = 
coverage_prob(spillover_mixed1, -14.7) #0.955

spillover_mixed2 = rbind(readRDS("spillover_results/1cat1num/spillover_mixed2.RDS"),
                         readRDS("spillover_results/1cat1num/spillover_mixed2(1).RDS"))
sd(spillover_mixed2$estimate) # 1.589753
mean(spillover_mixed2$std.error) # 1.365737
mean(spillover_mixed2$estimate) # -13.02027, bias = 
coverage_prob(spillover_mixed2, -13) #0.9116667

#### making tables ####
inf_table = c(-(mean(inf_plain$estimate)+13), sd(inf_plain$estimate), 
              mean(inf_plain$std.error), coverage_prob(inf_plain, -13))
inf_table = rbind(inf_table, c(-(mean(inf_binary1$estimate)+13), sd(inf_binary1$estimate), 
                               mean(inf_binary1$std.error), coverage_prob(inf_binary1, -13)))
inf_table = rbind(inf_table, c(-(mean(inf_binary3$estimate)+21), sd(inf_binary3$estimate), 
                               mean(inf_binary3$std.error), coverage_prob(inf_binary3, -21)))
inf_table = rbind(inf_table, c(-(mean(inf_cont1$estimate)+13), sd(inf_cont1$estimate), 
                               mean(inf_cont1$std.error), coverage_prob(inf_cont1, -13)))
inf_table = rbind(inf_table, c(-(mean(inf_cont2$estimate)+7.5), sd(inf_cont2$estimate), 
                               mean(inf_cont2$std.error), coverage_prob(inf_cont2, -7.5)))
inf_table = rbind(inf_table, c(-(mean(inf_mixed1$estimate)+13), sd(inf_mixed1$estimate), 
                               mean(inf_mixed1$std.error), coverage_prob(inf_mixed1, -13)))
inf_table = rbind(inf_table, c(-(mean(inf_mixed2$estimate)+7.5), sd(inf_mixed2$estimate), 
                               mean(inf_mixed2$std.error), coverage_prob(inf_mixed2, -7.5)))
colnames(inf_table) <- c("bias", "MC SE", "analytical SE", "CP")
rownames(inf_table) <- c("No-Con", "Binary, Het Sp", "Binary, Het Inf",
                         "Cont, Het Sp", "Cont, Het Inf",
                         "Mixed, Het Sp", "Mixed, Het Inf")
inf_table <- round(inf_table, 3)
#inf_mixed1 = rbind(readRDS("inf_results/mixed/inf_mixed1(new var_effect).RDS"))

spillover_table = c(-(mean(spillover_plain$estimate)+13), sd(spillover_plain$estimate), 
              mean(spillover_plain$std.error), coverage_prob(spillover_plain, -13))
spillover_table = rbind(spillover_table, c(-(mean(spillover_binary1$estimate)+21), sd(spillover_binary1$estimate), 
                               mean(spillover_binary1$std.error), coverage_prob(spillover_binary1, -21)))
spillover_table = rbind(spillover_table, c(-(mean(spillover_binary3$estimate)+13), sd(spillover_binary3$estimate), 
                               mean(spillover_binary3$std.error), coverage_prob(spillover_binary3, -13)))
spillover_table = rbind(spillover_table, c(-(mean(spillover_cont1$estimate)+7.5), sd(spillover_cont1$estimate), 
                               mean(spillover_cont1$std.error), coverage_prob(spillover_cont1, -7.5)))
spillover_table = rbind(spillover_table, c(-(mean(spillover_cont2$estimate)+13), sd(spillover_cont2$estimate), 
                               mean(spillover_cont2$std.error), coverage_prob(spillover_cont2, -13)))
spillover_table = rbind(spillover_table, c(-(mean(spillover_mixed1$estimate)+7.5), sd(spillover_mixed1$estimate), 
                               mean(spillover_mixed1$std.error), coverage_prob(spillover_mixed1, -7.5)))
spillover_table = rbind(spillover_table, c(-(mean(spillover_mixed2$estimate)+13), sd(spillover_mixed2$estimate), 
                               mean(spillover_mixed2$std.error), coverage_prob(spillover_mixed2, -13)))
colnames(spillover_table) <- c("bias", "MC SE", "analytical SE", "CP")
rownames(spillover_table) <- c("No-Con", "Binary, Het Sp", "Binary, Het Inf",
                         "Cont, Het Sp", "Cont, Het Inf",
                         "Mixed, Het Sp", "Mixed, Het Inf")

library(ftExtra)
spillover_table <- round(spillover_table, 3)
result_table <- as.data.frame(cbind(spillover_table, inf_table))
colnames(result_table)[1:4] <- paste("spillover", colnames(result_table)[1:4], sep = "_")
colnames(result_table)[5:8] <- paste("influencer", colnames(result_table)[5:8], sep = "_")
ft <- flextable(result_table)
ft %>% separate_header(sep = "_")
ft %>% span_header()

