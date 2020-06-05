#' @export

structree.default <-
function(formula,
                               data,
                               family=gaussian,
                               stop_criterion=c("AIC","BIC","CV","pvalue"),
                               splits_max=NULL,
                               fold=5,
                               alpha=0.05,
                               grid_value=NULL,
                               min_border=NULL,
                               ridge=FALSE,
                               constant_covs=FALSE,
                               trace=TRUE,
                               plot=TRUE,
                               k=10,
                               weights=NULL,
                               offset=NULL,
                               ...){
  
  # check input 
  if(missing(formula)){
    stop("Argument formula is missing with no default.")
  }
  if(missing(data)){
    stop("Argument data is missing with no default.")
  }
  
  # check family 
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  if (family$family=="Gamma") family <- Gamma(link="log")
  
  stop_criterion <- match.arg(stop_criterion)
  
  # predefinition 
  comp <- specification(formula,data)  
  y            <- comp$y
  DM_kov       <- comp$DM_kov
  
  if(comp$fixef){
    secondlevel  <- comp$secondlevel
    slope        <- comp$slope
  } else{
    which_kat    <- comp$which_kat
    which_smooth <- comp$which_smooth
  }
  
  if(comp$fixef){
    if(stop_criterion == "CV"){
      stop("Cross-Validation not implemented for repeated measurements!")
    }
    if((family$family=="Gamma" & !is.null(ridge))){
      stop("Ridge estimation is only implemented for families: gaussian, binomial and poisson!")
    }
    lambda <- NULL 
    if(nrow(DM_kov)<1000){
      n_levels <- nlevels(DM_kov[,secondlevel])
      lambda   <- setl(n_levels,family)
    }
    output <- structree_fixef(y,DM_kov,secondlevel,slope,family,stop_criterion,splits_max,
                              alpha,grid_value,min_border,ridge,constant_covs,trace,weights,offset,lambda)
    
    attr(output,"type")        <- "fixedEff"
    attr(output,"secondlevel") <- secondlevel
    attr(output,"slope")       <- slope
  } else{
    output <- structree_cat(y,DM_kov,which_kat,which_smooth,family,stop_criterion,splits_max,fold,alpha,trace,plot,k,weights,offset)

    attr(output,"type")         <- "catPred"  
    attr(output,"which_kat")    <- which_kat
    attr(output,"which_smooth") <- which_smooth
    if(!is.null(which_smooth)){    
      attr(output,"k") <- k 
    }
  }
  
  output$call <- match.call() 
  
  # define class 
  class(output) <- "structree"
  
  return(output)
  
}
