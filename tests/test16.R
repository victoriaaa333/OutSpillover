# Check unbiasedness for weighted regression for conditional H 
# (heterogenous spillover) without including extra terms, 
# and run simulations with larger sample to see if you omit extra term what happens. 


# Seems like the problem is because the number of 2nd-order
# neighbors is different for different clusters
# TODO: why this case?
hist(neigh2_treated[A == 1])
hist(neigh2_treated[A == 0])

# neighX is not correlated with treated 2nd-order neighbors, as expected
cor(neigh2_treated, neighX) #-0.01491738 -0.01212384
cor(neigh2_treated[1:100], neighX[1:100]) #first group: 0.03749419

# maybe should try regression within each group without additional terms?
