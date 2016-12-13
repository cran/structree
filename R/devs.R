devs <-
function(y,mu,family){
  
  family <- family$family
  
  if(family=="gaussian"){
    dev <- (y-mu)^2
  }
  if(family=="Gamma"){
    dev <- -2*(log(y/mu)+((y-mu)/mu))
  }
  if(family=="inverse.gaussian"){
    dev <- (y-mu)^2/(mu^2*y)
  }
  if(family=="binomial"){
    dev <- -2*log(1-abs(y-mu))
  }
  if(family=="poisson"){
    dev <-  2*(y*log(y/mu)-(y-mu))
  }
  return(dev)
}
