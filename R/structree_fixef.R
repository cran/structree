structree_fixef <-
function(y,
                            DM_kov,
                            secondlevel,
                            slope,
                            family,
                            stop_criterion,
                            splits_max,
                            alpha,
                            grid_value,
                            min_border,
                            ridge,
                            constant_covs,
                            trace,
                            weights,
                            offset,
                            lambda){
  
  # global parameters 
  n            <- length(y)
  n_gamma      <- ncol(DM_kov)-length(slope)-1
  which_gamma  <- which(!(names(DM_kov) %in% c(secondlevel,slope)))
  which_slope  <- which(names(DM_kov) %in% slope)
  which_second <- which(names(DM_kov)==secondlevel)
  n_levels     <- nlevels(DM_kov[,secondlevel])
  if(is.null(weights)){ weights <- rep(1,n)}
  if(is.null(offset)) { offset  <- rep(0,n)}
  
  if(family$family=="gaussian"){model <- "linear"}
  if(family$family=="binomial"){model <- "logistic"}
  if(family$family=="poisson"){model <- "poisson"}
  
  # modify DM_kov 
  DM_kov_mod <- DM_kov 
  DM_kov_mod[,secondlevel] <- factor(DM_kov[,secondlevel],labels=1:n_levels)
  colnames(DM_kov_mod) <- paste("x",1:ncol(DM_kov),sep="")
  
  # big design  
  design_fm <- sapply(1:(n_levels-1), function(j) ifelse(DM_kov_mod[,which_second] ==j,1,0))
  if(!is.null(slope)){
    for(i in 1:length(slope)){
      design_fm <- cbind(design_fm,design_fm*DM_kov_mod[,which_slope[i]])
    }
  }
  n1 <- rep(paste0("p",c(0,which_slope)),each=(n_levels-1))
  n2 <- rep(1:(n_levels-1),length=ncol(design_fm))
  colnames(design_fm) <- paste0(n1,n2)
  
  # fixed-effects model
  if(!constant_covs){
    daten_ml <- as.data.frame(cbind(y,DM_kov[,c(which_slope,which_gamma),drop=FALSE],design_fm))  
    if(!ridge){
      est_ml <- glm(y~.,data=daten_ml,family=family,weights=weights,offset=offset)
    } else{
      h1 <- paste(names(daten_ml)[-1],collapse="+")
      h2 <- formula(paste("~",h1,"-1",sep=""))
      est_ml <- penalized(y,penalized=h2,lambda2=0.01,model=model,data=daten_ml,trace=FALSE)
    }
  } else{
    data_mixed <- cbind(y=y,DM_kov)
    form1 <- paste0(names(DM_kov)[which_gamma],collapse="+")
    if(is.null(slope)){
      form2 <- paste0("(1|",names(DM_kov)[which_second],")")
    } else{
      form2 <- NULL
      for(i in 1:length(slope)){
        temp  <- paste0("(",names(DM_kov[which_slope])[i],"|",names(DM_kov)[which_second],")")
        form2 <- paste0(c(form2,temp),collapse="+")
      }
    }
    form_mixed <- formula(paste0("y~",form1,"+",form2))
    est_mixed  <- lmer(form_mixed,data=data_mixed,weights=weights,offset=offset)
  }
  
  # get order 
  order <- lapply(0:length(slope),function(j){
    if(!constant_covs){
      return(order(c(tail(coefficients(est_ml),ncol(design_fm))[(1:(n_levels-1))+(n_levels-1)*j],0)))
    } else{
      return(order(ranef(est_mixed)[[secondlevel]][,j+1]))
    }
  })
  names(order)[1] <- secondlevel
  names(order)[-1] <- slope
  
  # restrict size of boarders 
  if(!is.null(min_border) & is.null(grid_value)){
    splits_grid <- lapply(0:length(slope), function(j){
      min_border:(n_levels-min_border)
    })
    names(splits_grid)[1] <- secondlevel
    names(splits_grid)[-1] <- slope
    n_splits_grid <- sapply(1:length(order),function(j) length(splits_grid[[j]]))
  }
  
  # compute only on a grid 
  if(!is.null(grid_value)){
    splits_grid <- lapply(0:length(slope),function(j){
      if(!constant_covs){
        vals     <- c(tail(coefficients(est_ml),ncol(design_fm))[(1:(n_levels-1))+(n_levels-1)*j],0)
      } else{
        vals     <- ranef(est_mixed)[[secondlevel]][,j]
      }
      vals_dif <- sort(vals)[-1]-sort(vals)[-n_levels]
      as.numeric(which(vals_dif>=grid_value[j+1]))
    })
    names(splits_grid)[1]  <- secondlevel
    names(splits_grid)[-1] <- slope
    n_splits_grid          <- sapply(1:length(order),function(j) length(splits_grid[[j]]))
    
    if(!is.null(min_border)){
      for(j in 1:length(order)){
        splits_grid[[i]] <- splits_grid[[j]][splits_grid[[j]]%in%(min_border:(n_levels-min_border))]
      }
      n_splits_grid <- sapply(1:length(order),function(j) length(splits_grid[[j]]))
    }
    
    if(is.null(splits_max)){
      splits_max    <- n_splits_grid
    } else{
      if(splits_max>n_splits_grid){
        splits_max  <- n_splits_grid
      }
    }
  }
  
  # build design 
  design_list <- lapply(1:length(order),function(i){                   
    sapply(2:n_levels,function(j) { ifelse(DM_kov_mod[,which_second] %in% order[[i]][j:n_levels],1,0)})
  })
  if((!is.null(min_border) | !is.null(grid_value))){
    for(i in 1:length(design_list)){
      design_list[[i]] <- design_list[[i]][,splits_grid[[i]]]
    }
  }
  if(!is.null(slope)){
    for(i in 1:length(slope)){
      design_list[[i+1]] <- design_list[[i+1]]*DM_kov_mod[,which_slope[i]]
    }
  }
  design_tree <- do.call(cbind,design_list)                
  n_s <- ncol(design_tree)
  if(is.null(splits_max)){
    splits_max <- n_s
  } else{
    if(splits_max>n_s){
      splits_max <- n_s
    }
  }
  
  if((is.null(min_border) & is.null(grid_value))){
    v <- unlist(lapply(1:length(order),function(j) order[[j]][-length(order[[j]])]))        
    w <- rep(paste("s",c(0,which_slope),sep=""),each=(n_levels-1))
  } else{
    v <- unlist(lapply(1:length(order),function(j) order[[j]][splits_grid[[j]]]))    
    w <- rep(paste("s",c(0,which_slope),sep=""),(n_splits_grid-1))
  }
  colnames(design_tree) <- paste(w,v,sep="")
  
  ####################################################################
  
  # function to estimate tree 
  tree <- function(dat,splits_max,we,off,silent=FALSE){
    
    # initalize
    ps <- c()
    splits_evtl <- 1:n_s    
    mod_potential <- mod_potential_up <- list()  
    splits <- c()            
    count <- 1             
    
    
    # null model 
    if(n_gamma>0){
      help0  <- paste("x",c(which_slope,which_gamma),sep="",collapse="+")          
      help1 <- formula(paste("y~",help0,sep=""))
      mod0 <-  mod_potential[[count]] <-  mod_potential_up[[count]] <- glm(help1,data=dat,family=family,weights=we,offset=off)       
    } else{ 
      help00 <- paste0("x",which_slope,collapse="+")
      help10 <- formula(paste("y~",help00,sep=""))
      mod0 <-  mod_potential[[count]] <-  mod_potential_up[[count]] <- glm(help10,data=dat,family=family,weights=we,offset=off)       
    }
    
    
    # successive estimation 
    while(count<=splits_max){                                
      
      pvalues <- rep(NA,n_s)                                        
      
      for(i in splits_evtl){                                       
        
        help2 <- paste(w,v,sep="")
        help3 <- formula(paste("~.+",help2[i],sep=""))              
        
        mod <- update(mod0,help3)                                   
        
        dev <-  deviance(mod0)-deviance(mod)
        pvalues[i] <- pchisq(dev,1,lower.tail=FALSE,log.p=TRUE)    
        
      }
      ps[count] <- min(pvalues,na.rm=TRUE)
      splits[count] <- which.min(pvalues)                           
      
      help4 <- formula(paste("~.+",help2[splits[count]],sep=""))     
      mod_potential_up[[count+1]] <- mod0 <- update(mod0,help4) 
      if(!is.null(lambda)){
        h3 <- formula(paste0("~",help0))
        h4 <- paste(names(coefficients(mod0))[-(1:(n_gamma+1))],collapse="+")
        h5 <- formula(paste0("~",h4))
        mod_potential[[count+1]] <- penalized(y,penalized=h5,unpenalized=h3,
                                              lambda2=lambda,model=model,
                                              data=dat,trace=FALSE)
      } else{
        mod_potential[[count+1]] <- mod_potential_up[[count+1]]
      }
      splits_evtl <- splits_evtl[-which(splits_evtl==splits[count])]
      count <- count + 1   
      if(trace & !silent){
        cat(paste("Split",count-1,"\n",sep=" "))
      }
    }
    return(list(splits,mod_potential_up,mod_potential,ps))
  }
  
  # estimation 
  dat <- as.data.frame(cbind(y,design_tree,DM_kov_mod[,which_slope,drop=FALSE],DM_kov_mod[,which_gamma,drop=FALSE])) 
  
  estimation    <- tree(dat,splits_max,weights,offset)
  mod_up        <- estimation[[2]]
  mod_potential <- estimation[[3]]
  splits        <- estimation[[1]] 
  
  if(!constant_covs){
    if(stop_criterion=="AIC"){
      AIC <- sapply(1:length(mod_up),function(j) AIC(mod_up[[j]])) 
      which_min <- which.min(AIC) 
      tune_values <- AIC
    }
    if(stop_criterion=="BIC"){
      logLiks <- sapply(1:length(mod_up), function(j) logLik(mod_up[[j]]))
      dfs <- sapply(1:length(mod_up), function(j) length(coefficients(mod_up[[j]])))
      BIC <- -2*logLiks+log(n)*dfs
      which_min <- which.min(BIC) 
      tune_values <- BIC
    } 
    if(stop_criterion=="pvalue"){
      dev_ml <- 0 
      for(i in 1:n){
        dev_ml <- dev_ml + devs(y[i],fitted(est_ml)[i],family)
      }
      ps  <- c() 
      sig <- TRUE
      sig_count <- 1 
      while(sig){
        lrstat <- deviance(mod_up[[sig_count]])-dev_ml
        dfs    <- n_levels*(1+length(slope))-sig_count
        p      <- ps[sig] <- pchisq(lrstat,dfs,lower.tail=FALSE)
        proof  <- (p<alpha)
        if(!proof | splits_max==(sig_count+1)){
          sig <- FALSE
          which_min <- sig_count
        } else{ 
          sig_count <- sig_count+1 
        }
      }
      tune_values <- ps 
    }
  } else{
    X <- as.matrix(dat[,-1])
    ps  <- c() 
    sig <- TRUE
    sig_count <- 1 
    while(sig){
      score       <- t(X)%*%(y-fitted(mod_up[[sig_count]]))
      Fhochminus1 <- solve(t(X)%*%X)
      scorestat   <- t(score)%*%Fhochminus1%*%score
      dfs         <- n_levels-sig_count
      p           <- ps[sig] <- pchisq(scorestat,dfs,lower.tail=FALSE)
      proof       <- (p<alpha)
      if(!proof){
        sig <- FALSE
        which_min <- sig_count
      } else{ 
        sig_count <- sig_count+1 
      }
    }
    tune_values <- ps
  }
  mod_opt <- mod_potential[[which_min]]
  
  #############################################################
  
  # save all splits 
  if( (is.null(min_border) & is.null(grid_value)) ){
    help5 <- c(0,(n_levels-1)+n_levels*(0:length(slope)))           
  } else{
    help5 <- c(0,cumsum(n_splits_grid-1)) 
  }
  
  splits_done <- lapply(1:length(order), function(j){
    help6 <- splits[splits %in% (help5[j]+1):help5[j+1]] 
    return(v[help6]) 
  })   
  names(splits_done)[1]  <- paste0("int|",secondlevel)
  names(splits_done)[-1] <- paste0(slope,"|",secondlevel)
  
  
  # build partitions 
  partitions <- lapply(1:length(order), function(j) {   
    
    splits_aktuell <- splits_done[[j]]           
    
    if(length(splits_aktuell)==0){                          
      fusion <- matrix(1,nrow=n_levels,ncol=1)  
      rownames(fusion) <- levels(DM_kov[,secondlevel])
      colnames(fusion) <- "Partition 0"
      return(fusion)
    } 
    fusion <- matrix(1,nrow=n_levels,ncol=length(splits_aktuell)+1) 
    rownames(fusion) <- levels(DM_kov[,secondlevel])
    colnames(fusion) <- paste("Partition",0:length(splits_aktuell))  
    
    for(i in 1:length(splits_aktuell)){                             
      help7 <- which(order[[j]]==splits_aktuell[i])+1              
      help8 <- order[[j]][help7:n_levels]                   
      fusion[help8,((i+1):ncol(fusion))] <- fusion[help8,((i+1):ncol(fusion))]+1  
    }
    return(fusion)
  })
  names(partitions) <- names(splits_done)

  
  # save coefficients 
  beta_hat <- lapply(1:length(order),function(j){
    
    which_splits <- which(splits %in% (help5[j]+1):help5[j+1])                    
    
    beta <- matrix(0,nrow=n_levels,ncol=length(mod_potential))
    
    if(length(which_splits)==0){ return(beta) } 
    
    
    
    for(i in 1:(length(mod_potential)-1)){
      
      splits_aktuell <- splits_done[[j]][which(which_splits <=i)]                 
      
      help9 <- paste0(rep(paste0("s",c(0,which_slope)[j]),length(splits_aktuell)),splits_aktuell) 
      alpha <- coefficients(mod_potential[[i+1]])[help9]                                    
      alpha_sort <- alpha[colnames(dat)]
      alpha_sort <- alpha_sort[!is.na(alpha_sort)]
      koeffizienten <- c(0,cumsum(alpha_sort))
      
      partition <- partitions[[j]][,length(splits_aktuell)+1]
      beta[,i+1] <- koeffizienten[partition]
      rownames(beta) <- levels(DM_kov[,secondlevel])
    }
    return(beta)
  })
  names(beta_hat) <- names(splits_done)
  
  # correct betas 
  cor_factors <- lapply(1:length(order), function(j){
    sapply(1:length(mod_potential), function(k) {coefficients(mod_potential[[k]])[j]})
  })
  for(j in 1:length(order)){
    for(i in 1:ncol(beta_hat[[j]])){
      beta_hat[[j]][,i] <- beta_hat[[j]][,i]+cor_factors[[j]][i]
    }
  }
  
   
  if(which_min>1){
    
    # save splits until stopp
    splits <- splits[1:(which_min-1)]                
  
    splits_done <- lapply(1:length(order), function(j) {
      which_splits <- which(splits %in% (help5[j]+1):help5[j+1]) 
      return(v[splits][which_splits])
    })
    names(splits_done) <- names(beta_hat)
  

    # save optimal partition 
    partitions_opt <- sapply(1:length(order), function(j) {
      splits_aktuell  <- splits_done[[j]]
      return(length(splits_aktuell))
    })
    
  } else{
    partitions_opt <- rep(0,length(order))
  }  
  
  
  # save coefficients of groups 
  coefs_group <- lapply(1:length(order),function(j) { 
    koeffizienten <- unique(beta_hat[[j]][,which_min][order[[j]]])
    names(koeffizienten) <- paste0("beta",c(0,which_slope)[j],1:length(koeffizienten))
    return(koeffizienten)
  })
  names(coefs_group) <- names(beta_hat)
  
  names(order) <- names(beta_hat)
  
  # save group ID 
  group_ID <- sapply(1:length(order),function(i){       
    gruppen_vektor <- rep(0,n)  
    for(j in 1:n_levels){    
      gruppen_vektor[DM_kov_mod[,which_second]==j] <- partitions[[i]][j,(partitions_opt[i]+1)] 
    }
    return(gruppen_vektor)
  })
  names(group_ID) <- names(beta_hat)


  #save coefficients of optimal model 
  coefs_all <- coefs_group
    
  if(n_gamma>0){
    for(i in 1:n_gamma){ 
      which_coefs <- grep(paste0("x",which_gamma[i]),names(coefficients(mod_opt)))
      coefs_all[[colnames(DM_kov)[which_gamma][i]]] <- coefficients(mod_opt)[which_coefs]
    }
  }
    
  coefs_end <- unlist(lapply(1:length(coefs_all), function(j) coefs_all[[j]])) 
  
  for(i in 1:n_gamma){
    po <- grep(paste0("x",which_gamma[i]),names(coefs_end))
    if(is.factor(DM_kov[,which_gamma[i]])){
      wh <- paste0(colnames(DM_kov)[which_gamma[i]],levels(DM_kov[,which_gamma[i]])[-1])
    } else{
      wh <- colnames(DM_kov[which_gamma[i]])
    }
    names(coefs_end)[po] <- wh 
  }

  n_groups <- sapply(1:length(order),function(j) length(coefs_group[[j]])) 
  help11 <- unlist(sapply(1:length(order), function(j) seq(1,n_groups[j])))
  help12 <- paste0(rep(names(beta_hat),n_groups),help11)
  help13 <- unlist(sapply(1:length(order),function(j) names(coefs_group[[j]])))
  help14 <- which(coefs_end %in% coefs_end[help13])
  names(coefs_end)[help14] <- help12
  
  
  # list to return 
  to_return <- list("coefs_end"=coefs_end,
                    "partitions"=partitions,
                    "beta_hat"=beta_hat,
                    "which_opt"=which_min,
                    "opts"=partitions_opt,
                    "order"=order,
                    "tune_values"=tune_values,
                    "group_ID"=group_ID,
                    "coefs_group"=coefs_group,
                    "y"=y,
                    "DM_kov"=DM_kov,
                    "model"=mod_opt)

  # return
  return(to_return)
  
}
