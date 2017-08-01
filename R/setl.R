setl <- function(n_levels, family){
  lambda <- NULL
  if(family$family=="gaussian"){
    if(n_levels>=200) lambda <- 50
    if(n_levels<200 & n_levels>=100) lambda <- 20 
    if(n_levels<100 & n_levels>=25)  lambda <- 5
  }
  if(family$family=="binomial"){
    if(n_levels>=100) lambda <- 5
    if(n_levels<100 & n_levels>=25) lambda <- 1 
    if(n_levels<25  & n_levels>=15) lambda <- 0.01
  }
  return(lambda)
}