source("propensity.R")
source("integrand.R")
source("weight_matrix.R")
source("utils.R")
source("bootstrap_variance.R")
source("effects.R")
source("m_variance.R")
source("regression_variance.R")
source("regression_utils.R")
#source("ipw_point_estimates_mixed.R")
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

for (i in 1:length(unit_names)) {
  neigh_names <- setdiff(as.vector(names(which(A[i,] == 1))),
                         noinfo_neigh)
  neigh_ind <- which(cai$id %in% neigh_names)
  H[i] <- ifelse(length(neigh_ind) > 0, mean(Y[neigh_ind]), 0)  
}

# 2. derive the point estimates
allocations = list(c(numerator_alpha, denominator_alpha))
w.matrix = wght_matrix(plain_integrand, allocations, group_no, trt, P)
#ipw_point_estimates_mixed2(H, group_no, trt, w.matrix)$outcomes$overall
ipw_point_estimates_mixed_test4(H, group_no, trt, w.matrix)$outcomes$overall

# 3. conditional on X
X_con = cbind(cai$agpop, cai$male, cai$educ)
colnames(X_con) <- c("agpop", "male", "educ")
x0_con = as.matrix(c(4, 1, 1))
X_type_con = as.matrix(c("N", "C", "C"))
ipw_point_estimates_mixed_test4(H, group_no, trt, w.matrix,
                                X = X_con, x0 = x0_con, 
                                X_type = X_type_con, Con_type = "group")$outcomes$overall

#######
# conditional within groups (influencer effect)
# case 1. conditional on 
# "education higher than primary school"
X_con = cbind(cai$educ_good)
colnames(X_con) <- c("educ") #, "educ"
x0_con = as.matrix(c(1))
X_type_con = as.matrix(c("C"))
point_estimates1 = ipw_point_estimates_mixed_test4(H, group_no, trt, w.matrix,
                                X = X_con, x0 = x0_con, 
                                X_type = X_type_con, Con_type = "group")
# 0.3 0.3525113 0.2330265
# 0.2 0.4656164 0.432609
ipw_m_variance(w.matrix, point_estimates1, effect_type ='contrast',
               marginal = FALSE, allocation1 = allocations[1], 
               allocation2 = allocations[1])
# 0.03300741 0.04025884 -0.04589847 0.1119133

x0_con = as.matrix(c(0))
X_type_con = as.matrix(c("C"))
point_estimates0 = ipw_point_estimates_mixed_test4(H, group_no, trt, w.matrix,
                                X = X_con, x0 = x0_con, 
                                X_type = X_type_con, Con_type = "group")
ipw_m_variance(w.matrix, point_estimates0, effect_type ='contrast',
               marginal = FALSE, allocation1 = allocations[1], 
               allocation2 = allocations[1])
#-0.003161112 0.01605324 -0.03462489 0.02830267

# 0.3 0.3088135 0.399043
# 0.2 0.4598913 0.4630524

# case 2. conditional on the gender of head
X_con = cbind(cai$male)
colnames(X_con) <- c("gender") #, "educ"
x0_con = as.matrix(c(1))
X_type_con = as.matrix(c("C"))
ipw_point_estimates_mixed_test4(H, group_no, trt, w.matrix,
                                X = X_con, x0 = x0_con, 
                                X_type = X_type_con, Con_type = "group")$outcomes$overall
# 0.2 0.4594254 0.4532883

x0_con = as.matrix(c(0))
X_type_con = as.matrix(c("C"))
ipw_point_estimates_mixed_test4(H, group_no, trt, w.matrix,
                                X = X_con, x0 = x0_con, 
                                X_type = X_type_con, Con_type = "group")$outcomes$overall
# 0.2 0.3993299 0.4179841


# case 3. conditional on disaster yes
X_con = cbind(cai$disaster_yes)
colnames(X_con) <- c("disaster") 
x0_con = as.matrix(c(1))
X_type_con = as.matrix(c("C"))
point_estimates1 = ipw_point_estimates_mixed_test4(H, group_no, trt, w.matrix,
                                X = X_con, x0 = x0_con, 
                                X_type = X_type_con, Con_type = "group")
# 0.2 0.4457274 0.4242831

ipw_m_variance(w.matrix, point_estimates1, effect_type ='contrast',
               marginal = FALSE, allocation1 = allocations[1], 
               allocation2 = allocations[1])
# 0.02144424 0.01662555 -0.01114125 0.05402972

x0_con = as.matrix(c(0))
X_type_con = as.matrix(c("C"))
point_estimates0 = ipw_point_estimates_mixed_test4(H, group_no, trt, w.matrix,
                                X = X_con, x0 = x0_con, 
                                X_type = X_type_con, Con_type = "group")
# 0.2 0.4851587 0.4939962
point_estimates0$outcomes$overall

ipw_m_variance(w.matrix, point_estimates0, effect_type ='contrast',
                   marginal = FALSE, allocation1 = allocations[1], 
                   allocation2 = allocations[1])
# -0.008837482 0.0286597 -0.06500947 0.0473345

# case 5. conditional on 
# continuous variables
X_con = cbind(cai$disaster_prob, cai$age)
colnames(X_con) <- c("disaster_prob" ,"age") #, "educ"
x0_con = as.matrix(c(50, 50))
X_type_con = as.matrix(c("N", "N"))
point_estimates1 = ipw_point_estimates_mixed_test4(H, group_no, trt, w.matrix,
                                                   X = X_con, x0 = x0_con, 
                                                   X_type = X_type_con, Con_type = "group")
point_estimates1$outcomes$overall
point_estimates1$outcomes$overall_coefG

# case 6. conditional on 
# "education higher than primary school"
X_con = cbind(cai$ricearea_2010)
colnames(X_con) <- c("ricearea_2010") #, "educ"
x0_con = as.matrix(c(15))
X_type_con = as.matrix(c("N"))
point_estimates1 = ipw_point_estimates_mixed_test4(H, group_no, trt, w.matrix,
                                                   X = X_con, x0 = x0_con, 
                                                   X_type = X_type_con, Con_type = "group")
point_estimates1$outcomes$overall
point_estimates1$outcomes$overall_coefG

ipw_regression_variance(H, w.matrix, point_estimates1, trt, effect_type ='contrast',
                        marginal = FALSE, allocation1 = allocations[1], allocation2 = allocations[1],
                        X = X_con, x0 = x0_con, X_type = X_type_con)
# 0.003270367 0.01277966 -0.02177731 0.02831804

# , , 0
# 
# c(0.2, 0.219)
# [1,]  0.4647449676
# [2,] -0.0001896032
# 
# , , 1
# 
# c(0.2, 0.219)
# [1,]  0.4653806662
# [2,] -0.0008484958



# higher rice production will influence ppl not to get insurance?



#######
# conditional within neighbors
X_cat = as.matrix(cbind(cai$educ_good))
Cat_ind <- ifelse(X_cat == 1, 1, NA)
H_cat <- rep(NA, length(unit_names))
for (i in 1:length(unit_names)) {
  neigh_names <- setdiff(as.vector(names(which(A[i,] == 1))),
                         noinfo_neigh)
  neigh_ind <- which(cai$id %in% neigh_names)
  
  H_cat[i] <- ifelse(length(neigh_ind) > 0, 
                 mean(Y[neigh_ind] * Cat_ind[neigh_ind], na.rm = TRUE), NA)  
}
ipw_point_estimates_mixed_test4(H_cat, group_no, trt, w.matrix)$outcomes$overall
# 0.6168022 2.570761

Cat_ind <- ifelse(X_cat == 0, 1, NA)
H_cat <- rep(NA, length(unit_names))
for (i in 1:length(unit_names)) {
  neigh_names <- setdiff(as.vector(names(which(A[i,] == 1))),
                         noinfo_neigh)
  neigh_ind <- which(cai$id %in% neigh_names)
  
  H_cat[i] <- ifelse(length(neigh_ind) > 0, 
                     mean(Y[neigh_ind] * Cat_ind[neigh_ind], na.rm = TRUE), NA)  
}
ipw_point_estimates_mixed_test(H_cat, group_no, trt, w.matrix)$outcomes$overall
# 0.5208146 2.077891




# ipw_point_estimates_mixed_test4(H_cat, group_no, trt, w.matrix,
#                                 X = X_con, x0 = x0_con, 
#                                 X_type = X_type_con, Con_type = "group")$outcomes$overall








# ipw_point_estimates_mixed2(H, group_no, trt, w.matrix, 
#           X = X_con, x0 = x0_con, X_type = X_type_con)$outcomes$overall

###

#? 1. how to deal with missing values



#######
# Per Laura's comment
# 1. "education higher than primary school"

# influencer effect
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
# 0.03300741 0.04025884 -0.04589847 0.1119133

x0_con = as.matrix(c(0))
X_type_con = as.matrix(c("C"))
point_estimates0 = ipw_point_estimates_mixed_test4(H, group_no, trt, w.matrix,
                                                   X = X_con, x0 = x0_con, 
                                                   X_type = X_type_con, Con_type = "group")
point_estimates0$outcomes$overall
ipw_m_variance(w.matrix, point_estimates0, effect_type ='contrast',
               marginal = FALSE, allocation1 = allocations[1], 
               allocation2 = allocations[1])
# -0.01128633 0.01721438 -0.0450259 0.02245324

# spillover effect
X_cat = as.matrix(cbind(cai$educ_good))
Cat_ind0 <- ifelse(X_cat == 0, 1, NA)
Cat_ind1 <- ifelse(X_cat == 1, 1, NA)
H_cat0 <- rep(NA, length(unit_names))
H_cat1 <- rep(NA, length(unit_names))

for (i in 1:length(unit_names)) {
  neigh_names <- setdiff(as.vector(names(which(A[i,] == 1))),
                         noinfo_neigh)
  neigh_ind <- which(cai$id %in% neigh_names)
  
  H_cat0[i] <- ifelse(length(neigh_ind) > 0, 
                     mean(Y[neigh_ind] * Cat_ind0[neigh_ind], na.rm = TRUE), NA)  
  H_cat1[i] <- ifelse(length(neigh_ind) > 0, 
                     mean(Y[neigh_ind] * Cat_ind1[neigh_ind], na.rm = TRUE), NA)  
  }
point_estimates0 = ipw_point_estimates_mixed_test4(H_cat0, group_no, trt, w.matrix)
point_estimates0$outcomes$overall
# 0.5208146 2.077891
ipw_m_variance(w.matrix, point_estimates0, effect_type ='contrast',
               marginal = FALSE, allocation1 = allocations[1], 
               allocation2 = allocations[1])
# -1.557076 0.1180226 -1.788396 -1.325756

point_estimates1 = ipw_point_estimates_mixed_test4(H_cat1, group_no, trt, w.matrix)
point_estimates1$outcomes$overall
# 0.6168022 2.570761
ipw_m_variance(w.matrix, point_estimates1, effect_type ='contrast',
               marginal = FALSE, allocation1 = allocations[1], 
               allocation2 = allocations[1])
# -1.953959 0.1494305 -2.246838 -1.661081


# 2. "education higher than primary school"

# spillover effect
X_cat = as.matrix(cbind(cai$male))
Cat_ind0 <- ifelse(X_cat == 0, 1, NA)
Cat_ind1 <- ifelse(X_cat == 1, 1, NA)
H_cat0 <- rep(NA, length(unit_names))
H_cat1 <- rep(NA, length(unit_names))

for (i in 1:length(unit_names)) {
  neigh_names <- setdiff(as.vector(names(which(A[i,] == 1))),
                         noinfo_neigh)
  neigh_ind <- which(cai$id %in% neigh_names)
  
  H_cat0[i] <- ifelse(length(neigh_ind) > 0, 
                      mean(Y[neigh_ind] * Cat_ind0[neigh_ind], na.rm = TRUE), NA)  
  H_cat1[i] <- ifelse(length(neigh_ind) > 0, 
                      mean(Y[neigh_ind] * Cat_ind1[neigh_ind], na.rm = TRUE), NA)  
}
point_estimates0 = ipw_point_estimates_mixed_test4(H_cat0, group_no, trt, w.matrix)
point_estimates0$outcomes$overall
# 0.4964561 2.305188
ipw_m_variance(w.matrix, point_estimates0, effect_type ='contrast',
               marginal = FALSE, allocation1 = allocations[1], 
               allocation2 = allocations[1])
# -1.808732 0.1577534 -2.117923 -1.499541

point_estimates1 = ipw_point_estimates_mixed_test4(H_cat1, group_no, trt, w.matrix)
point_estimates1$outcomes$overall
# 0.5579043 2.241763
ipw_m_variance(w.matrix, point_estimates1, effect_type ='contrast',
               marginal = FALSE, allocation1 = allocations[1], 
               allocation2 = allocations[1])
# -1.683858 0.1228697 -1.924678 -1.443038



# 3.  influence effect for educated people on female/male

# influencer effect of (high) education level on male
X_con = cbind(cai$educ_good)
colnames(X_con) <- c("educ") 
x0_con = as.matrix(c(1))
X_type_con = as.matrix(c("C"))
point_estimates1 = ipw_point_estimates_mixed_test4(H_cat1, group_no, trt, w.matrix,
                                                   X = X_con, x0 = x0_con, 
                                                   X_type = X_type_con, Con_type = "group")
point_estimates1$outcomes$overall
ipw_m_variance(w.matrix, point_estimates1, effect_type ='contrast',
               marginal = FALSE, allocation1 = allocations[1], 
               allocation2 = allocations[1])
# 0.03345483 0.0431369 -0.05109193 0.1180016

point_estimates0 = ipw_point_estimates_mixed_test4(H_cat0, group_no, trt, w.matrix,
                                                   X = X_con, x0 = x0_con, 
                                                   X_type = X_type_con, Con_type = "group")
point_estimates0$outcomes$overall
ipw_m_variance(w.matrix, point_estimates0, effect_type ='contrast',
               marginal = FALSE, allocation1 = allocations[1], 
               allocation2 = allocations[1])
# -0.01397612 0.04640187 -0.1049221 0.07696988

# influencer effect of (low) education level on male
x0_con = as.matrix(c(0))
X_type_con = as.matrix(c("C"))
point_estimates1 = ipw_point_estimates_mixed_test4(H_cat1, group_no, trt, w.matrix,
                                                   X = X_con, x0 = x0_con, 
                                                   X_type = X_type_con, Con_type = "group")
point_estimates1$outcomes$overall
ipw_m_variance(w.matrix, point_estimates0, effect_type ='contrast',
               marginal = FALSE, allocation1 = allocations[1], 
               allocation2 = allocations[1])
# -0.001604815 0.01788073 -0.03665039 0.03344076


# influencer effect of (low) education level on male
x0_con = as.matrix(c(0))
X_type_con = as.matrix(c("C"))
point_estimates0 = ipw_point_estimates_mixed_test4(H_cat0, group_no, trt, w.matrix,
                                                   X = X_con, x0 = x0_con, 
                                                   X_type = X_type_con, Con_type = "group")
point_estimates0$outcomes$overall
ipw_m_variance(w.matrix, point_estimates0, effect_type ='contrast',
               marginal = FALSE, allocation1 = allocations[1], 
               allocation2 = allocations[1])
# -0.06069786 0.03160707 -0.1226466 0.001250867

# older people rather than younger




