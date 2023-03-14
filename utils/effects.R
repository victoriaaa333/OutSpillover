group_direct_effect <- function(object)
{ if (!is.null(dim(object$outcomes$overall)[1])){
  out <- object$outcomes$groups[,,2] - object$outcomes$groups[,,1]
} else{
  out <- object$outcomes$groups[,2] - object$outcomes$groups[,1]
}
  return(out)
}

population_direct_effect <- function(object)
{ if (!is.null(dim(object$outcomes$overall)[1])){
  out <- object$outcomes$overall[,2] - object$outcomes$overall[,1]
}else{
  out <- as.vector(object$outcomes$overall[2] - object$outcomes$overall[1])
}
  return(out)
}

