# 11.29.22: data analysis with different alphas
# the first line of result is for numerator alpha = 0.2, the 2nd line is for 0.1

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
numerator_alpha = 0.2 # hypothetical one
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
ps = ipw_point_estimates_mixed_test4(H, group_no, trt, w.matrix)
ipw_m_variance(w.matrix, ps, effect_type ='contrast',
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
point_estimates1 = ipw_point_estimates_mixed_test4(H, group_no, trt, w.matrix,
                                                   X = X_con, x0 = x0_con, 
                                                   X_type = X_type_con, Con_type = "group")
point_estimates1$outcomes$overall
ipw_m_variance(w.matrix, point_estimates1, effect_type ='contrast',
               marginal = FALSE, allocation1 = allocations[1], 
               allocation2 = allocations[1])
# -0.06322612 0.3861603 -0.8200865 0.6936343
# -0.1051372 0.05255516 -0.2081434 -0.00213096

X_con = cbind(cai$educ_good)
colnames(X_con) <- c("educ") 
x0_con = as.matrix(c(0))
X_type_con = as.matrix(c("C"))
point_estimates1 = ipw_point_estimates_mixed_test4(H, group_no, trt, w.matrix,
                                                   X = X_con, x0 = x0_con, 
                                                   X_type = X_type_con, Con_type = "group")
point_estimates1$outcomes$overall
ipw_m_variance(w.matrix, point_estimates1, effect_type ='contrast',
               marginal = FALSE, allocation1 = allocations[1], 
               allocation2 = allocations[1])
# -0.01875443 0.3033071 -0.6132255 0.5757166
# -0.08439429 0.03397693 -0.1509878 -0.01780073

## 2. gender
X_con = cbind(cai$male)
colnames(X_con) <- c("male") 
x0_con = as.matrix(c(1))
X_type_con = as.matrix(c("C"))
point_estimates1 = ipw_point_estimates_mixed_test4(H, group_no, trt, w.matrix,
                                                   X = X_con, x0 = x0_con, 
                                                   X_type = X_type_con, Con_type = "group")
point_estimates1$outcomes$overall
ipw_m_variance(w.matrix, point_estimates1, effect_type ='contrast',
               marginal = FALSE, allocation1 = allocations[1], 
               allocation2 = allocations[1])
# -0.03280604 0.3084413 -0.6373398 0.5717277
# -0.09196189 0.02511788 -0.141192 -0.04273175

X_con = cbind(cai$male)
colnames(X_con) <- c("male") 
x0_con = as.matrix(c(0))
X_type_con = as.matrix(c("C"))
point_estimates1 = ipw_point_estimates_mixed_test4(H, group_no, trt, w.matrix,
                                                   X = X_con, x0 = x0_con, 
                                                   X_type = X_type_con, Con_type = "group")
point_estimates1$outcomes$overall
ipw_m_variance(w.matrix, point_estimates1, effect_type ='contrast',
               marginal = FALSE, allocation1 = allocations[1], 
               allocation2 = allocations[1])
# -0.0584406 0.2349916 -0.5190156 0.4021344
# -0.2117205 0.04259431 -0.2952038 -0.1282372

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
point_estimates0 = ipw_point_estimates_mixed_test4(H_cat0_educ, group_no, trt, w.matrix)
point_estimates0$outcomes$overall
# 0.5134731 2.149417
# 0.4089973 0.4381527

ipw_m_variance(w.matrix, point_estimates0, effect_type ='contrast',
               marginal = FALSE, allocation1 = allocations[1], 
               allocation2 = allocations[1])
# -1.635944 0.1807552 -1.990218  -1.28167
# -0.02915537 0.02703428 -0.08214159 0.02383085

point_estimates1 = ipw_point_estimates_mixed_test4(H_cat1_educ, group_no, trt, w.matrix)
point_estimates1$outcomes$overall
# 0.5704437 2.622685
# 0.4198014 0.5775893
ipw_m_variance(w.matrix, point_estimates1, effect_type ='contrast',
               marginal = FALSE, allocation1 = allocations[1], 
               allocation2 = allocations[1])
# -2.052241 0.2085514 -2.460995 -1.643488
# -0.1577879 0.03246266 -0.2214136 -0.09416229

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

point_estimates0 = ipw_point_estimates_mixed_test4(H_cat0_gender, group_no, trt, w.matrix)
point_estimates0$outcomes$overall
# 0.443711 2.181702
# 0.1030425 0.2993597
ipw_m_variance(w.matrix, point_estimates0, effect_type ='contrast',
               marginal = FALSE, allocation1 = allocations[1], 
               allocation2 = allocations[1])
# -1.737991 0.1545809 -2.040964 -1.435018
# -0.1963172 0.02944118 -0.2540209 -0.1386135

point_estimates1 = ipw_point_estimates_mixed_test4(H_cat1_gender, group_no, trt, w.matrix)
point_estimates1$outcomes$overall
# 0.529175 2.299631
# 0.4410313 0.4825769
ipw_m_variance(w.matrix, point_estimates1, effect_type ='contrast',
               marginal = FALSE, allocation1 = allocations[1], 
               allocation2 = allocations[1])
# -1.770456 0.1746595 -2.112783  -1.42813
# -0.04154555 0.02384584 -0.08828254 0.005191436

# 3. mixed influencer
## 1. educ_good
X_con = cbind(cai$educ_good)
colnames(X_con) <- c("educ") 
x0_con = as.matrix(c(1))
X_type_con = as.matrix(c("C"))
point_estimates1 = ipw_point_estimates_mixed_test4(H_cat0_gender, group_no, trt, w.matrix,
                                                   X = X_con, x0 = x0_con, 
                                                   X_type = X_type_con, Con_type = "group")
point_estimates1$outcomes$overall
ipw_m_variance(w.matrix, point_estimates1, effect_type ='contrast',
               marginal = FALSE, allocation1 = allocations[1], 
               allocation2 = allocations[1])
# -0.1207019 0.1175057 -0.3510089 0.1096051
# -0.1477642 0.04249263 -0.2310482 -0.06448017

X_con = cbind(cai$educ_good)
colnames(X_con) <- c("educ") 
x0_con = as.matrix(c(0))
X_type_con = as.matrix(c("C"))
point_estimates1 = ipw_point_estimates_mixed_test4(H_cat1_gender, group_no, trt, w.matrix,
                                                   X = X_con, x0 = x0_con, 
                                                   X_type = X_type_con, Con_type = "group")
point_estimates1$outcomes$overall
ipw_m_variance(w.matrix, point_estimates1, effect_type ='contrast',
               marginal = FALSE, allocation1 = allocations[1], 
               allocation2 = allocations[1])
# -0.003911905 0.3017426 -0.5953166 0.5874928
# -0.2345463 0.03482465 -0.3028014 -0.1662913

## 2. gender
X_con = cbind(cai$male)
colnames(X_con) <- c("male") 
x0_con = as.matrix(c(1))
X_type_con = as.matrix(c("C"))
point_estimates1 = ipw_point_estimates_mixed_test4(H_cat1_educ, group_no, trt, w.matrix,
                                                   X = X_con, x0 = x0_con, 
                                                   X_type = X_type_con, Con_type = "group")
point_estimates1$outcomes$overall
ipw_m_variance(w.matrix, point_estimates1, effect_type ='contrast',
               marginal = FALSE, allocation1 = allocations[1], 
               allocation2 = allocations[1])
# -0.04419716 0.3530058 -0.7360759 0.6476816
# -0.1660223 0.04008617 -0.2445898 -0.08745487

X_con = cbind(cai$male)
colnames(X_con) <- c("male") 
x0_con = as.matrix(c(1))
X_type_con = as.matrix(c("C"))
point_estimates1 = ipw_point_estimates_mixed_test4(H_cat0_educ, group_no, trt, w.matrix,
                                                   X = X_con, x0 = x0_con, 
                                                   X_type = X_type_con, Con_type = "group")
point_estimates1$outcomes$overall
ipw_m_variance(w.matrix, point_estimates1, effect_type ='contrast',
               marginal = FALSE, allocation1 = allocations[1], 
               allocation2 = allocations[1])
# -0.02845127 0.2974453 -0.6114333 0.5545308
# -0.02598186 0.02840162 -0.08164801 0.02968428
# male has significant influence over educated groups

X_con = cbind(cai$male)
colnames(X_con) <- c("male") 
x0_con = as.matrix(c(0))
X_type_con = as.matrix(c("C"))
point_estimates1 = ipw_point_estimates_mixed_test4(H_cat1_educ, group_no, trt, w.matrix,
                                                   X = X_con, x0 = x0_con, 
                                                   X_type = X_type_con, Con_type = "group")
point_estimates1$outcomes$overall
ipw_m_variance(w.matrix, point_estimates1, effect_type ='contrast',
               marginal = FALSE, allocation1 = allocations[1], 
               allocation2 = allocations[1])
#-0.09771414 0.2049657 -0.4994396 0.3040113
#-0.1034478 0.04500944 -0.1916647 -0.0152309

X_con = cbind(cai$male)
colnames(X_con) <- c("male") 
x0_con = as.matrix(c(0))
X_type_con = as.matrix(c("C"))
point_estimates1 = ipw_point_estimates_mixed_test4(H_cat0_educ, group_no, trt, w.matrix,
                                                   X = X_con, x0 = x0_con, 
                                                   X_type = X_type_con, Con_type = "group")
point_estimates1$outcomes$overall
ipw_m_variance(w.matrix, point_estimates1, effect_type ='contrast',
               marginal = FALSE, allocation1 = allocations[1], 
               allocation2 = allocations[1])
# 0.06687498 0.1824786 -0.2907765 0.4245265
# 0.07858636 0.04160822 -0.002964244  0.160137

# male influencer -0.09196189 0.02511788 -0.141192 -0.04273175
# female influencer -0.2117205 0.04259431 -0.2952038 -0.1282372
  

