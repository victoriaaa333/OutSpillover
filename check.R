HAX <- cbind(ind_est, A, X)
influencer_cond = apply(as.matrix(HAX[,3:dim(HAX)[2]]), 1, 
                        function(x) ifelse(prod(x == "M"), 1, NA))
mean(as.numeric(ind_est[ 1:100]) * influencer_cond[1:100], na.rm = TRUE)
group_means(ind_est, A, G, X, x0 = "M", a = NA)[1:10]
