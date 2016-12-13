specification <-
function(formula,data){
  
  formula <- paste(formula)
  name_y <- formula[2]
  
  ### y ###  
  if(!(name_y %in% names(data))){
    stop("Response y undefined", call.=FALSE)
  }
  y <- data[,name_y]
  
  ### DM_kov ### 
  x <- gsub(" ","",formula[3])
  x <- unlist(strsplit(x,"\\+"))

  # wrong function 
  fun_true <- c("tr","s")
  h0 <- "(.+)([(])(.+)"
  ind0 <- grep(h0,x)
  fun  <- unique(sub(h0,"\\1",x[ind0]))
  
  if(!all(fun %in% fun_true)){
    wrong <- paste(fun[which(!(fun%in%fun_true))],collapse="and")
    stop(paste("Function",wrong,"undefined!"), call.=FALSE)
  }
  
  fixef <- FALSE
  
  # tr()
  h1 <- "(tr[(])(.+)([)])"
  ind1 <- grepl(h1,x)
  if(!any(ind1)){
     stop("No variable to fuse!",call.=FALSE)
  } else {
     raneff <- grepl("^([^\\|]+?)\\|([^\\|]+)$",x[ind1])
     if(!all(raneff)){
       if(any(raneff)){
         stop("diftree can either handle categorical predictors or fixed effects, not both!",call.=FALSE)
       } else{
         which_kat <- sub(h1,"\\2",x[ind1])
         secondlevel <- NULL
         slope       <- NULL
       }
     } else{
       fixef <- TRUE
       h01 <- strsplit(sub(h1,"\\2",x[ind1]),"\\|")
       secondlevel <- unique(sapply(h01,function(x){x[[2]]}))
       if(length(secondlevel)>1){
         stop("Second-level unit not uniquely defined!",call.=FALSE)
       } else{
         slope <- sapply(h01,function(x){x[[1]]})
         if(!("1" %in% slope)){
           warning("Unit-specific intercept was not specified, but incorporated in the model!")
         }
         slope <- slope[slope!="1"]
         if(length(slope)==0){
           slope <- NULL
         }
       }
       which_kat <- NULL 
     }
  }
  
  # s() - which_smooth 
  h2 <- "(s[(])(.+)([)])"
  ind2 <- grepl(h2,x)
  if(!any(ind2)){
     which_smooth <- NULL
  } else {
     if(fixef){
      which_smooth <- NULL
      warning("Smooth estimates not implemented for repeated measurements!")
     } else{
      which_smooth <- sub(h2,"\\2",x[ind2])
     }
  }
  
  # others 
  which_linear <- x[!(ind1|ind2)]
  
  if(fixef){
    all_vars <- c(secondlevel,slope,which_smooth,which_linear)
  } else{
    all_vars <- c(which_kat,which_smooth,which_linear)
  }
  if(any(!(all_vars %in% names(data)))){
    wrong <- paste(all_vars[which(!(all_vars %in% names(data)))],collapse="and")
    stop(paste("Covariate",wrong,"undefined!"),call.=FALSE)
  }
  DM_kov <- data[,all_vars]
  
  return(list("y"=y,
              "DM_kov"=DM_kov,
              "fixef"=fixef,
              "which_kat"=which_kat,
              "secondlevel"=secondlevel,
              "slope"=slope, 
              "which_smooth"=which_smooth))
}
