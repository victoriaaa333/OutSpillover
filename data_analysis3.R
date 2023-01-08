# 1.05.23: 
# use test5 for point estimates (average across groups)

source("propensity.R")
source("integrand.R")
source("weight_matrix.R")
source("utils.R")
source("bootstrap_variance.R")
source("effects.R")
source("m_variance.R")
source("regression_variance.R")
source("regression_utils.R")
source("ipw_point_estimate_tests.R")
library(igraph)
library(lme4)
library(dplyr)

### 1. formatting dataset for input to our function
library(ggplot2)
theme_set(theme_bw())
options(repr.plot.width = 6)
options(repr.plot.height = 4)

library(icsw)
library(Matrix)

cai <- read.csv("cai_data/cai.all.csv")
head(cai)
colnames(cai)

trt <- ifelse(cai$delay == 0 & cai$intensive == 1, 1, 0)
cai$trt <- trt
prob_A = round(sum(cai$trt)/length(cai$trt), 3) # 0.2187091

# 1.1 necessary bits
numerator_alpha = 0.1 # hypothetical one
denominator_alpha = prob_A # assume only one-stage
P = 1

X<- cbind(cai$agpop, cai$age, cai$male, cai$educ,
          cai$literacy, cai$educ_good, cai$ricearea_2010,
          cai$disaster_loss, cai$disaster_yes, cai$insurance_repay,
          cai$disaster_prob, cai$understanding, cai$risk_averse, cai$n.peers)

X_type <- c("N", "N", "C", "C",
            "C", "C", "N",
            "N", "C", "C",
            "N", "N", "N", "N")
Y <- cai$takeup_survey
grp <- cai$village
cai <- cai %>% mutate(group_no = as.integer(as.factor(village)))
group_no <- cai$group_no

# 1.2 calculate H
load("cai_data/cai.adjacency.RData") #, header=F
unit_names <- cai$id
H <- rep(NULL, length(unit_names))
noinfo_neigh <- setdiff(colnames(A), unit_names) # units that are neighs but no information collected
filter_ind <- intersect(which(cai$delay == 1), which(cai$intensive == 0)) 
neigh_inds <- c()
filter_neigh_inds <- c()

for (i in 1:length(unit_names)) {
  neigh_names <- setdiff(as.vector(names(which(A[i,] == 1))),
                         noinfo_neigh)
  neigh_ind <- which(cai$id %in% neigh_names)
  neigh_inds <- c(neigh_inds, length(neigh_ind))
  filter_neigh_ind <- intersect(filter_ind, neigh_ind)
  filter_neigh_inds <- c(filter_neigh_inds, length(filter_neigh_ind))
  #H[i] <- ifelse(length(neigh_ind) > 0, mean(Y[neigh_ind]), 0)  
  H[i] <- ifelse(length(filter_neigh_ind) > 0, 
                 mean(Y[filter_neigh_ind]), NA)  
}

table(neigh_inds[trt == 1])
# 0   1   2   3   4   5 
# 53  54 108 222 340 226 

# mean(filter_neigh_inds[trt == 1]/neigh_inds[trt == 1], na.rm = TRUE)
# 0.2660526

table(filter_neigh_inds[trt == 1])
# 0   1   2   3   4   5 
# 403 352 188  55   3   2 
# 600/ 1003

table(neigh_inds)
# 0    1    2    3    4    5 
# 244  246  561 1007 1394 1134 

length(which(neigh_inds>0)) # 4342
length(intersect(which(trt == 1), which(neigh_inds>0))) # 950
# 0.2187932

# 1.3. derive the point estimates
allocations = list(c(numerator_alpha, denominator_alpha))
w.matrix = wght_matrix(plain_integrand, allocations, group_no, trt, P)
ps = ipw_point_estimates_mixed_test5(H, group_no, trt, w.matrix, 
                                     Con_type = "No-Con")
ipw_m_variance_groups(w.matrix, ps, effect_type ='contrast',
               marginal = FALSE, allocation1 = allocations[1], 
               allocation2 = allocations[1])

# numerator alpha = 0.2
# estimate std.error  conf.low conf.high
# -1.828399 0.1648542 -2.151507  -1.50529

# 0.1
# estimate  std.error   conf.low   conf.high
# -0.1920246 0.05535768 -0.3005237 -0.08352558

# 1. try to derive significant influencer effect
## 1. educ_good
X_con = cbind(cai$educ_good)
colnames(X_con) <- c("educ") 
x0_con = as.matrix(c(1))
X_type_con = as.matrix(c("C"))
point_estimates1 = ipw_point_estimates_mixed_test5(H, group_no, trt, w.matrix,
                                                   X = X_con, x0 = x0_con, 
                                                   X_type = X_type_con, Con_type = "group")
point_estimates1$outcomes$overall
ipw_m_variance_groups(w.matrix, point_estimates1, effect_type ='contrast',
               marginal = FALSE, allocation1 = allocations[1], 
               allocation2 = allocations[1])
# -2.066097 0.2668384 -2.589091 -1.543103
# -0.2016906 0.05197853 -0.3035667 -0.0998146

X_con = cbind(cai$educ_good)
colnames(X_con) <- c("educ") 
x0_con = as.matrix(c(0))
X_type_con = as.matrix(c("C"))
point_estimates1 = ipw_point_estimates_mixed_test5(H, group_no, trt, w.matrix,
                                                   X = X_con, x0 = x0_con, 
                                                   X_type = X_type_con, Con_type = "group")
point_estimates1$outcomes$overall
ipw_m_variance_groups(w.matrix, point_estimates1, effect_type ='contrast',
               marginal = FALSE, allocation1 = allocations[1], 
               allocation2 = allocations[1])
# -1.68993 0.1804809 -2.043666 -1.336194
# -0.1604209 0.06700564 -0.2917496 -0.02909228

## 2. gender
X_con = cbind(cai$male)
colnames(X_con) <- c("male") 
x0_con = as.matrix(c(1))
X_type_con = as.matrix(c("C"))
point_estimates1 = ipw_point_estimates_mixed_test5(H, group_no, trt, w.matrix,
                                                   X = X_con, x0 = x0_con, 
                                                   X_type = X_type_con, Con_type = "group")
point_estimates1$outcomes$overall
ipw_m_variance_groups(w.matrix, point_estimates1, effect_type ='contrast',
               marginal = FALSE, allocation1 = allocations[1], 
               allocation2 = allocations[1])
# -1.798278 0.1697622 -2.131005  -1.46555
# -0.1851188 0.05370878 -0.2903861 -0.07985151

X_con = cbind(cai$male)
colnames(X_con) <- c("male") 
x0_con = as.matrix(c(0))
X_type_con = as.matrix(c("C"))
point_estimates1 = ipw_point_estimates_mixed_test5(H, group_no, trt, w.matrix,
                                                   X = X_con, x0 = x0_con, 
                                                   X_type = X_type_con, Con_type = "group")
point_estimates1$outcomes$overall
ipw_m_variance_groups(w.matrix, point_estimates1, effect_type ='contrast',
               marginal = FALSE, allocation1 = allocations[1], 
               allocation2 = allocations[1])
# -2.101454 0.1480771 -2.39168 -1.811228
# -0.1720672 0.0463294 -0.2628712 -0.08126329

# 2. spillover effect for the two variables
## 1. educ_good
X_cat = as.matrix(cbind(cai$educ_good))
Cat_ind0 <- ifelse(X_cat == 0, 1, NA)
Cat_ind1 <- ifelse(X_cat == 1, 1, NA)
H_cat0_educ <- rep(NA, length(unit_names))
H_cat1_educ <- rep(NA, length(unit_names))

for (i in 1:length(unit_names)) {
  neigh_names <- setdiff(as.vector(names(which(A[i,] == 1))),
                         noinfo_neigh)
  neigh_ind <- which(cai$id %in% neigh_names)
  filter_neigh_ind <- intersect(filter_ind, neigh_ind)
  
  H_cat0_educ[i] <- ifelse(length(filter_neigh_ind) > 0,
                           mean(Y[filter_neigh_ind] * Cat_ind0[filter_neigh_ind], na.rm = TRUE), NA)
  H_cat1_educ[i] <- ifelse(length(filter_neigh_ind) > 0,
                           mean(Y[filter_neigh_ind] * Cat_ind1[filter_neigh_ind], na.rm = TRUE), NA)
  
  # H_cat0[i] <- ifelse(length(neigh_ind) > 0, 
  #                     mean(Y[neigh_ind] * Cat_ind0[neigh_ind], na.rm = TRUE), NA)  
  # H_cat1[i] <- ifelse(length(neigh_ind) > 0, 
  #                     mean(Y[neigh_ind] * Cat_ind1[neigh_ind], na.rm = TRUE), NA)  
}
point_estimates0 = ipw_point_estimates_mixed_test5(H_cat0_educ, group_no, trt, w.matrix)
point_estimates0$outcomes$overall
# 0.5134731 2.149417
# 0.01997627 0.1557994

ipw_m_variance_groups(w.matrix, point_estimates0, effect_type ='contrast',
               marginal = FALSE, allocation1 = allocations[1], 
               allocation2 = allocations[1])
# -1.635944 0.1807552 -1.990218  -1.28167
# -0.1358231 0.03907048 -0.2123998 -0.05924638

point_estimates1 = ipw_point_estimates_mixed_test5(H_cat1_educ, group_no, trt, w.matrix)
point_estimates1$outcomes$overall
# 0.5704437 2.622685
# 0.02204371 0.3038441
ipw_m_variance_groups(w.matrix, point_estimates1, effect_type ='contrast',
               marginal = FALSE, allocation1 = allocations[1], 
               allocation2 = allocations[1])
# -2.052241 0.2085514 -2.460995 -1.643488
# -0.2818004 0.1212181 -0.5193835 -0.0442173

## 2. gender
X_cat = as.matrix(cbind(cai$male))
Cat_ind0 <- ifelse(X_cat == 0, 1, NA)
Cat_ind1 <- ifelse(X_cat == 1, 1, NA)
H_cat0_gender <- rep(NA, length(unit_names))
H_cat1_gender <- rep(NA, length(unit_names))

for (i in 1:length(unit_names)) {
  neigh_names <- setdiff(as.vector(names(which(A[i,] == 1))),
                         noinfo_neigh)
  neigh_ind <- which(cai$id %in% neigh_names)
  filter_neigh_ind <- intersect(filter_ind, neigh_ind)
  
  H_cat0_gender[i] <- ifelse(length(filter_neigh_ind) > 0,
                             mean(Y[filter_neigh_ind] * Cat_ind0[filter_neigh_ind], na.rm = TRUE), NA)
  H_cat1_gender[i] <- ifelse(length(filter_neigh_ind) > 0,
                             mean(Y[filter_neigh_ind] * Cat_ind1[filter_neigh_ind], na.rm = TRUE), NA)
}

point_estimates0 = ipw_point_estimates_mixed_test5(H_cat0_gender, group_no, trt, w.matrix)
point_estimates0$outcomes$overall
# 0.443711 2.181702
# 0.006458385 0.09666967

ipw_m_variance_groups(w.matrix, point_estimates0, effect_type ='contrast',
               marginal = FALSE, allocation1 = allocations[1], 
               allocation2 = allocations[1])
# -1.737991 0.1545809 -2.040964 -1.435018
# -0.09021129 0.02146719 -0.1322862 -0.04813637

point_estimates1 = ipw_point_estimates_mixed_test5(H_cat1_gender, group_no, trt, w.matrix)
point_estimates1$outcomes$overall
# 0.529175 2.299631
# 0.02363512 0.2077956

ipw_m_variance_groups(w.matrix, point_estimates1, effect_type ='contrast',
               marginal = FALSE, allocation1 = allocations[1], 
               allocation2 = allocations[1])
# -1.770456 0.1746595 -2.112783  -1.42813
# -0.1841605 0.05402237 -0.2900424 -0.07827863

# 3. mixed influencer
## 1. educ_good
X_con = cbind(cai$educ_good)
colnames(X_con) <- c("educ") 
x0_con = as.matrix(c(1))
X_type_con = as.matrix(c("C"))
point_estimates1 = ipw_point_estimates_mixed_test5(H_cat0_gender, group_no, trt, w.matrix,
                                                   X = X_con, x0 = x0_con, 
                                                   X_type = X_type_con, Con_type = "group")
point_estimates1$outcomes$overall
ipw_m_variance_groups(w.matrix, point_estimates1, effect_type ='contrast',
               marginal = FALSE, allocation1 = allocations[1], 
               allocation2 = allocations[1])
# -1.763593 0.1218333 -2.002382 -1.524804
# -0.04977337 0.002989961 -0.05563359 -0.04391316

X_con = cbind(cai$educ_good)
colnames(X_con) <- c("educ") 
x0_con = as.matrix(c(0))
X_type_con = as.matrix(c("C"))
point_estimates1 = ipw_point_estimates_mixed_test5(H_cat1_gender, group_no, trt, w.matrix,
                                                   X = X_con, x0 = x0_con, 
                                                   X_type = X_type_con, Con_type = "group")
point_estimates1$outcomes$overall
ipw_m_variance_groups(w.matrix, point_estimates1, effect_type ='contrast',
               marginal = FALSE, allocation1 = allocations[1], 
               allocation2 = allocations[1])
# -1.626557 0.1871573 -1.993378 -1.259735
# -0.1517001 0.06716609 -0.2833432 -0.02005695

## 2. gender
X_con = cbind(cai$male)
colnames(X_con) <- c("male") 
x0_con = as.matrix(c(1))
X_type_con = as.matrix(c("C"))
point_estimates1 = ipw_point_estimates_mixed_test5(H_cat1_educ, group_no, trt, w.matrix,
                                                   X = X_con, x0 = x0_con, 
                                                   X_type = X_type_con, Con_type = "group")
point_estimates1$outcomes$overall
ipw_m_variance_groups(w.matrix, point_estimates1, effect_type ='contrast',
               marginal = FALSE, allocation1 = allocations[1], 
               allocation2 = allocations[1])
# -1.935911 0.2224576 -2.371919 -1.499902
# -0.271737 0.1198843 -0.5067059 -0.03676812

X_con = cbind(cai$male)
colnames(X_con) <- c("male") 
x0_con = as.matrix(c(1))
X_type_con = as.matrix(c("C"))
point_estimates1 = ipw_point_estimates_mixed_test5(H_cat0_educ, group_no, trt, w.matrix,
                                                   X = X_con, x0 = x0_con, 
                                                   X_type = X_type_con, Con_type = "group")
point_estimates1$outcomes$overall
ipw_m_variance_groups(w.matrix, point_estimates1, effect_type ='contrast',
               marginal = FALSE, allocation1 = allocations[1], 
               allocation2 = allocations[1])
# -1.637701 0.1826856 -1.995758 -1.279644
# -0.1358382 0.0392694 -0.2128049 -0.05887163

X_con = cbind(cai$male)
colnames(X_con) <- c("male") 
x0_con = as.matrix(c(0))
X_type_con = as.matrix(c("C"))
point_estimates1 = ipw_point_estimates_mixed_test5(H_cat1_educ, group_no, trt, w.matrix,
                                                   X = X_con, x0 = x0_con, 
                                                   X_type = X_type_con, Con_type = "group")
point_estimates1$outcomes$overall
ipw_m_variance_groups(w.matrix, point_estimates1, effect_type ='contrast',
               marginal = FALSE, allocation1 = allocations[1], 
               allocation2 = allocations[1])
#-2.115864 0.1246165 -2.360108  -1.87162
#-0.2029613 0.04058908 -0.2825145 -0.1234082

X_con = cbind(cai$male)
colnames(X_con) <- c("male") 
x0_con = as.matrix(c(0))
X_type_con = as.matrix(c("C"))
point_estimates1 = ipw_point_estimates_mixed_test5(H_cat0_educ, group_no, trt, w.matrix,
                                                   X = X_con, x0 = x0_con, 
                                                   X_type = X_type_con, Con_type = "group")
point_estimates1$outcomes$overall
ipw_m_variance_groups(w.matrix, point_estimates1, effect_type ='contrast',
               marginal = FALSE, allocation1 = allocations[1], 
               allocation2 = allocations[1])
# -1.364449  0.139413 -1.637693 -1.091204
# 4.391586e-05 0.002541246 -0.004936835 0.005024667

# Female has significant influence over educated groups compared to less-educated groups
#(CI: [-2.360108,  -1.87162], [-1.637693, -1.091204] not overlapping)
