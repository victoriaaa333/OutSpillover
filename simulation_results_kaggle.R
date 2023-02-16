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
sd(spillover_plain$estimate) # 4.459234
mean(spillover_plain$std.error) # 4.218611
mean(spillover_plain$estimate) # -12.70476, bias = 0.29524
coverage_prob(spillover_plain, -13) #0.867
coverage_prob0(spillover_plain, -13, sd(spillover_plain$estimate)) # 0.945

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

### Spillover Model Simulations ###
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

## Mixed variables ##
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

### Influencer Model Simulations ###
## No-condition ##
inf_plain = readRDS("inf_results/num_of_clusters_200/inf_plain.RDS")
sd(inf_plain$estimate) # 1.655949
mean(inf_plain$std.error) # 1.612916
mean(inf_plain$estimate) # -13.13284, bias = -0.13284
coverage_prob(inf_plain, -13) #0.946

## Binary varibles ##
inf_binary1 =  readRDS("inf_results/num_of_clusters_200/inf_binary1.RDS")
sd(inf_binary1$estimate) # 1.649183
mean(inf_binary1$std.error) # 1.701642
mean(inf_binary1$estimate) # -13.0517, bias = -0.0517
coverage_prob(inf_binary1, -13) #0.955

# Binary variable under group conditions, overall
inf_binary2 =  readRDS("inf_results/num_of_clusters_200/inf_binary2.RDS")
sd(inf_binary2$estimate) # 1.875211
mean(inf_binary2$std.error) # 3.209203
mean(inf_binary2$estimate) # -21.0948, bias = -0.0948
coverage_prob(inf_binary2, -21) #0.974

# Binary variable under group conditions, w/i groups
inf_binary3 =  readRDS("inf_results/num_of_clusters_200/inf_binary3.RDS")
sd(inf_binary3$estimate) # 2.434932
mean(inf_binary3$std.error) # 3.196304
mean(inf_binary3$estimate) # -21.06475, bias = -0.06475
coverage_prob(inf_binary3, -21) #0.96

## Continuous variables ##
inf_cont1 = readRDS("inf_results/num_of_clusters_200/inf_cont1.RDS")
sd(inf_cont1$estimate) # 4.062581
mean(inf_cont1$std.error) # 3.745038
mean(inf_cont1$estimate) # -12.76634, bias = 0.23366
coverage_prob(inf_cont1, -13) #0.91

# Continuous variables under group conditions, overall
inf_cont2 = readRDS("inf_results/num_of_clusters_200/inf_cont2.RDS")
sd(inf_cont2$estimate) # 1.08374
mean(inf_cont2$std.error) # 1.132037
mean(inf_cont2$estimate) # -7.55599, bias = -0.05599
coverage_prob(inf_cont2, -7.5) #0.96

# No group conditions because it's a regression estimator

#### Mixed variables (noc = 200, sample size = 200) ####
inf_mixed1 = rbind(readRDS("inf_results/num_of_clusters_200/inf_mixed1(sample size = 200).RDS"),
                   readRDS("inf_results/num_of_clusters_200/inf_mixed1(sample size = 200)2.RDS"))
sd(inf_mixed1$estimate) # 2.978993
mean(inf_mixed1$std.error) # 2.819687
mean(inf_mixed1$estimate) # -13.09401, bias = -0.09401
coverage_prob(inf_mixed1, -13) #0.93

# Mixed variables under group conditions, overall
inf_mixed2 = rbind(readRDS("inf_results/num_of_clusters_200/inf_mixed2(sample size = 200).RDS"),
                   readRDS("inf_results/num_of_clusters_200/inf_mixed2(sample size = 200)2.RDS"))
sd(inf_mixed2$estimate) # 1.792545
mean(inf_mixed2$std.error) # 1.311842
mean(inf_mixed2$estimate) # -7.449184, bias = 0.050816
coverage_prob(inf_mixed2, -7.5) #.8475


