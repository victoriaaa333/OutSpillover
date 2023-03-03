coverage_prob <- function(dataset, true_effect){
  coverage <- intersect(which(dataset$conf.low <= true_effect), 
                        which(dataset$conf.high >= true_effect))
  return(length(coverage)/length(dataset$estimate))
}

result_stats <- function(result, true_effect){
  p.est = -mean(result$estimate)
  bias = p.est - true_effect
  mc.se = sd(result$estimate)
  ana.se = mean(result$std.error)
  cp = coverage_prob(result, -true_effect)
  return(round(c(p.est, bias, mc.se, ana.se, cp),3))
}