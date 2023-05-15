args = commandArgs(trailingOnly=TRUE)

print(args)
print(args[1])
print(typeof(args[1]))

filename = paste(
  "cluster_results/inf_model_both_var(sd=1, w/boot, run=", 
  args[1], 
  ").RDS", sep='')

print(filename)