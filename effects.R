group_direct_effect <- function(object)
{
  out <- object$outcomes$groups[,,2] - object$outcomes$groups[,,2]
  return(out)
}

population_direct_effect <- function(object)
{
  out <- object$outcomes$overall[,2] - object$outcomes$overall[,1]
  return(out)
}

