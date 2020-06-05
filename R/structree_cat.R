structree_cat <-
function(y,                    
                           DM_kov,              
                           which_kat,
                           which_smooth,
                           family,
                           stop_criterion,
                           splits_max,
                           fold,
                           alpha,
                           trace,
                           plot,
                           k,
                           weights,
                           offset){
  
  
  # global parameters
  n            <- length(y)                            
  n_fak        <- length(which_kat)                    
  which_kat    <- which(names(DM_kov) %in% which_kat)  
  if(!is.null(which_smooth)){
    which_smooth <- which(names(DM_kov) %in% which_smooth)
  }
  n_gamma      <- ncol(DM_kov)-n_fak                  
  n_levels     <- sapply(which_kat,function(j) length(levels(DM_kov[,j]))-1)
  if(is.null(weights)){ weights <- rep(1,n)}
  if(is.null(offset)) { offset  <- rep(0,n)}
  
  # check which_kat 
  is_factor <- which(sapply(which_kat,function(j) is.factor(DM_kov[,j]))!=rep(TRUE,n_fak))
  if(length(is_factor)==1){
    stop("One variable to fuse is not a factor")
  }
  if(length(is_factor)>1){
    stop("Variables to fuse are no factors")
  }
  which_levels   <- lapply(which_kat,function(j) levels(DM_kov[,j]))
  
  # modify DM_kov 
  DM_kov_mod <- DM_kov 
  
  for(i in which_kat){
    if(is.ordered(DM_kov[,i])){
      DM_kov_mod[,i] <- ordered(DM_kov[,i],labels=paste(0:n_levels[which(which_kat==i)]))       # use labels 0,...,k 
    } else{ DM_kov_mod[,i] <- factor(DM_kov[,i],labels=paste(0:n_levels[which(which_kat==i)]))  
    }
  }
  colnames(DM_kov_mod) <- paste("x",1:ncol(DM_kov),sep="") 
  
  # ML estimation 
  all_ordered <- sapply(1:n_fak, function(j) is.ordered(DM_kov_mod[,which_kat[j]]))
  if(sum(all_ordered)<n_fak){
    daten_ml <- data.frame(y=y,DM_kov_mod)   
    if(is.null(which_smooth)){
      est_ml   <- glm(y~.,data=daten_ml,family=family,weights=weights,offset=offset)
    } else{
      lp_smooth <- paste0("s(x",which_smooth,",bs='cr',k=",k,")",collapse="+")
      lp_linear <- paste0("x",(1:ncol(DM_kov_mod))[-c(which_smooth)],collapse="+")
      lp        <- formula(paste0("y~",lp_linear,"+",lp_smooth))
      est_ml    <- gam(lp,data=daten_ml,family=family,weights=weights,offset=offset)
    }
  }
  
  # get order 
  order <- lapply(1:n_fak,function(i){  
    
    if(is.ordered(DM_kov_mod[,which_kat[i]])){
      as.numeric(levels(DM_kov_mod[,which_kat[i]]))        
    } else { which_coefs <- coef(est_ml)[paste("x",which_kat[i],as.numeric(levels(DM_kov_mod[,which_kat[i]])[-1]),sep="")] 
    as.numeric(order(c(0,which_coefs)))-1}
    
  })
  names(order) <- paste("x",which_kat,sep="")
  
  # build design 
  n_s <- sum(n_levels)    
  if(is.null(splits_max)){
    splits_max <- n_s
  } else{
    if(splits_max>n_s){
      splits_max <- n_s
    }
  }

  design_list <- lapply(1:n_fak,function(i){               
    sapply(2:(n_levels[i]+1),function(k) {
      sapply(1:n,function(j) {
        ifelse(DM_kov_mod[j,which_kat[i]] %in% order[[i]][k:(n_levels[i]+1)],1,0)
      })
    })
  })
  design_tree <- do.call(cbind,design_list)                
  
  v <- unlist(lapply(1:n_fak,function(j) order[[j]][-length(order[[j]])]))        
  w <- rep(paste("s",which_kat,sep=""),(n_levels))
  colnames(design_tree) <- paste(w,v,sep="")
  
  #########################################################################################

  # function to estimate tree 
  tree <- function(dat,splits_max,we,off,silent=FALSE){
    
    # initals 
    ps <- c()
    splits_evtl   <- 1:n_s     
    mod_potential <- list()  
    splits <- c()           
    count <- 1         
    dat$we  <- we
    dat$off <- off
    
    # null model 
    if(n_gamma>0){
      if(is.null(which_smooth)){
        help0  <- paste0("x",(1:ncol(DM_kov_mod))[-which_kat],collapse="+")          
        help1  <- formula(paste0("y~",help0))
        mod0   <- mod_potential[[count]] <- glm(help1,data=dat,family=family,weights=we,offset=off)   
      } else{
        if(n_gamma>length(which_smooth)){
          help0  <- paste0("x",(1:ncol(DM_kov_mod))[-c(which_kat,which_smooth)],"+",collapse="+") 
        } else { 
          help0 <- NULL 
        }
        help1 <- formula(paste0("y~",help0,lp_smooth))
        mod0  <- mod_potential[[count]] <- gam(help1,data=dat,family=family,weights=we,offset=off)
      }
    } else{ 
      mod0 <- mod_potential[[count]] <- glm(y~1,data=dat,family=family,weights=we,offset=off)                                       
    }
    
    # sucessive estimtation 
    while(count<=splits_max){                               
      
      pvalues <- rep(NA,n_s)                                        
      
      for(i in splits_evtl){                                       
        
        help2 <- paste(w,v,sep="")
        help3 <- formula(paste("~.+",help2[i],sep=""))               
        # hinzu 
        
        mod <- update(mod0,help3)                                 
        
        dev <-  deviance(mod0)-deviance(mod)
        pvalues[i] <- pchisq(dev,1,lower.tail=FALSE,log.p=TRUE)                
        
      }
      ps[count] <- min(pvalues,na.rm=TRUE)
      splits[count] <- which.min(pvalues)                          
      
      help4 <- formula(paste("~.+",help2[splits[count]],sep=""))     
      mod_potential[[count+1]] <- mod0 <- update(mod0,help4)          
      splits_evtl <- splits_evtl[-which(splits_evtl==splits[count])] 
      count <- count + 1   
      if(trace & !silent){
        cat(paste("Split",count-1,"\n",sep=" "))
      }
      
    }
    return(list(splits,mod_potential,ps))
  }
  
  # estimation 
  dat           <- as.data.frame(cbind(y,design_tree,DM_kov_mod[,-which_kat,drop=FALSE])) 
  schaetzung    <- tree(dat,splits_max,weights,offset)
  mod_potential <- schaetzung[[2]]
  splits        <- schaetzung[[1]] 
  ps            <- schaetzung[[3]]
  ps            <- exp(ps)
  
  # determine optimal model 
  if(stop_criterion=="pvalue"){
    sig <- TRUE
    sig_count <- 1 
    while(sig){
      p <- ps[sig_count+1]
      alpha_adj <- alpha/(n_s-sig_count)  # adjust p-values
      proof <- (p < alpha_adj)
      if(!proof | splits_max==(sig_count+1)){
        sig <- FALSE
        which_min <- sig_count
      } else{ 
        sig_count <- sig_count+1 
      }
    }
    tune_values <- ps[1:which_min]
  }
  if(stop_criterion=="AIC"){
    AIC <- sapply(1:length(mod_potential),function(j) AIC(mod_potential[[j]])) 
    which_min <- which.min(AIC) 
    tune_values <- AIC
  }
  if(stop_criterion=="BIC"){
    logLiks <- sapply(1:length(mod_potential), function(j) logLik(mod_potential[[j]]))
    dfs <- sapply(1:length(mod_potential), function(j) length(coef(mod_potential[[j]])))
    BIC <- -2*logLiks+log(n)*dfs
    which_min <- which.min(BIC) 
    tune_values <- BIC
  } 
  if(stop_criterion=="CV"){
    
    t_index <- split(sample(1:n), rep(1:fold, length = n)) 
    l_index <- lapply(t_index, function(i) setdiff(1:n, i))
    
    cv <- matrix(NA,nrow=fold,ncol=length(mod_potential)) 
    
    for(k in 1:fold){
      dat_temp <- as.data.frame(cbind(y=y[l_index[[k]]],design_tree[l_index[[k]],],DM_kov_mod[l_index[[k]],-which_kat, drop=FALSE])) 
      weights_temp <- weights[l_index[[k]]]
      offset_temp  <- offset[l_index[[k]]]
      
      schaetzung_temp    <- tree(dat_temp,splits_max,weights_temp,offset_temp,silent=TRUE)
      mod_potential_temp <- schaetzung_temp[[2]]
      new_data           <- as.data.frame(cbind(y=y[t_index[[k]]],design_tree[t_index[[k]],],DM_kov_mod[t_index[[k]],-which_kat, drop=FALSE]))
      new_data$we        <- weights[t_index[[k]]]
      new_data$off       <- offset[t_index[[k]]]
      mu_hat             <- sapply(1:length(mod_potential_temp), function(j) predict(mod_potential_temp[[j]],newdata=new_data,type="response"))
      dif <- matrix(NA,nrow=length(t_index[[k]]),ncol=length(mod_potential_temp))
      for(i in 1:length(t_index[[k]])){
        for(m in 1:length(mod_potential)){
          dif[i,m] <- devs(y[t_index[[k]]][i],mu_hat[i,m],family=family)
        }
      }
      cv[k,1:length(mod_potential_temp)] <- apply(dif,2,sum) 
      if(trace){
        cat(paste("CV: Fold",k,"\n",sep=" "))
      }
    }
    CV <- apply(cv,2,sum)
    which_min <- which.min(CV)
    tune_values <- CV
  }
  mod_opt <- mod_potential[[which_min]]
  
  #########################################################################################
  
  # save all splits
  splits_done <- lapply(1:n_fak, function(j) c())    
  names(splits_done) <- colnames(DM_kov)[which_kat]  
  
  help5 <- c(0,cumsum(n_levels))                      
  
  for(j in 1:n_fak){
    help6 <- splits[splits %in% (help5[j]+1):help5[j+1]] 
    splits_done[[j]] <- v[help6] 
  } 
  
  # build partitions   
  partitions <- lapply(1:n_fak, function(j) {   
    
    splits_aktuell <- splits_done[[j]]           
    
    if(length(splits_aktuell)==0){                         
      fusion <- matrix(1,nrow=n_levels[j]+1,ncol=1) 
      rownames(fusion) <- which_levels[[j]]  
      colnames(fusion) <- "Partition 0"
      return(fusion)
    } 
    fusion <- matrix(1,nrow=n_levels[j]+1,ncol=length(splits_aktuell)+1) 
    rownames(fusion) <- which_levels[[j]]                                 
    colnames(fusion) <- paste("Partition",0:length(splits_aktuell))  
    
    for(i in 1:length(splits_aktuell)){                               
      help7 <- which(order[[j]]==splits_aktuell[i])+1               
      help8 <- order[[j]][help7:(n_levels[j]+1)]                    
      fusion[help8+1,((i+1):ncol(fusion))] <- fusion[help8+1,((i+1):ncol(fusion))]+1 
    }
    return(fusion)
  })
  names(partitions) <- colnames(DM_kov)[which_kat]
  
  # save coefficients 
  beta_hat <- lapply(1:n_fak,function(j){
    
    which_splits <- which(splits %in% (help5[j]+1):help5[j+1])                   
    
    beta <- matrix(0,nrow=n_levels[j]+1,ncol=length(mod_potential))
    
    if(length(which_splits)==0){ return(beta) } 
    
    for(i in 1:length(mod_potential)-1){
      
      splits_aktuell <- splits_done[[j]][which(which_splits <=i)]                 
      
      help9 <- paste0(rep(paste0("s",which_kat[j]),length(splits_aktuell)),splits_aktuell) 
      alpha <- coef(mod_potential[[i+1]])[help9]                                      
      alpha_sort <- alpha[colnames(dat)]
      alpha_sort <- alpha_sort[!is.na(alpha_sort)]
      koeffizienten <- c(0,cumsum(alpha_sort))
      
      partition <- partitions[[j]][,length(splits_aktuell)+1]
      beta[,i+1] <- koeffizienten[partition]
      rownames(beta) <- which_levels[[j]]
    }
    return(beta)
  })
  names(beta_hat) <- colnames(DM_kov)[which_kat] 
  
  # correct intercepts 
  coefs_0 <- sapply(1:n_fak, function(j) beta_hat[[j]][1,])
  if(ncol(coefs_0)==1){
    cor_factor <- as.vector(coefs_0)[-1]
  } else{
    cor_factor <- apply(coefs_0[-1,],1,sum)
  }
  intercepts <- sapply(1:length(mod_potential), function(j) (coef(mod_potential[[j]])[1]+cor_factor[j]))
  
  # correct betas 
  for(j in 1:n_fak){
    for(i in 1:ncol(beta_hat[[j]])){
      beta_hat[[j]][,i] <- beta_hat[[j]][,i]-beta_hat[[j]][1,i] 
    }
  }
  
  if(which_min>1){
  
    # save splits until stopp 
    splits <- splits[1:(which_min-1)]                 
  
    splits_done <- lapply(1:n_fak, function(j) {
      which_splits <- which(splits %in% (help5[j]+1):help5[j+1]) 
      return(v[splits][which_splits])
    })
    names(splits_done) <- colnames(DM_kov)[which_kat]   
  
  
    # save optimal partition 
    partitions_opt <- sapply(1:n_fak, function(j) {
      splits_aktuell  <- splits_done[[j]]
      return(length(splits_aktuell))
    })
    
  } else{
    partitions_opt <- rep(0,n_fak)
  }
  
  
  # save coefficients of groups 
  coefs_group <- lapply(1:n_fak,function(j) { 
    koeffizienten <- unique(beta_hat[[j]][-1,which_min][order[[j]]])
    names(koeffizienten) <- paste0("beta",which_kat[j],1:length(koeffizienten))
    return(koeffizienten)
  })
  names(coefs_group) <- colnames(DM_kov)[which_kat]  
  
  
  # save group ID 
  group_ID <- sapply(1:n_fak,function(i){       
    gruppen_vektor <- rep(0,n)  
    for(j in 0:n_levels[i]){    
      gruppen_vektor[DM_kov_mod[,which_kat[i]]==j] <- partitions[[i]][(j+1),(partitions_opt[i]+1)] 
    }
    return(gruppen_vektor)
  })
  colnames(group_ID) <- colnames(DM_kov)[which_kat] 
  
  names(order) <- colnames(DM_kov)[which_kat] 
  for(i in 1:length(order)){
    order[[i]] <- order[[i]]+1
  }
  
  # save coefficients of optimal model 
  coefs_all <- coefs_group
  
  if(n_gamma>0){
    if(is.null(which_smooth)){
      for(i in 1:n_gamma){ 
        coefs_all[[colnames(DM_kov)[-which_kat][i]]] <- coef(mod_opt)[1+i]
      }
    } else{
      if(n_gamma>length(which_smooth)){
        for(i in 1:(n_gamma-length(which_smooth))){ 
          coefs_all[[colnames(DM_kov)[-c(which_kat,which_smooth)][i]]] <- coef(mod_opt)[1+i]
        }
      }
    }
  }

  # sort
  if(is.null(which_smooth)){
    coefs_end <- unlist(lapply(colnames(DM_kov), function(j) coefs_all[[j]]))
  } else{
    coefs_end <- unlist(lapply(colnames(DM_kov)[-which_smooth], function(j) coefs_all[[j]])) 
  }
  coefs_end <- c(intercepts[which_min],coefs_end) 
  
  # name
  if(n_gamma>0){
    if(is.null(which_smooth)){
      names(coefs_end)[(length(coefs_end)-n_gamma+1):length(coefs_end)]<- colnames(DM_kov)[-which_kat]
    } else{
      if(n_gamma>length(which_smooth)){
        names(coefs_end)[(length(coefs_end)-n_gamma+length(which_smooth)+1):length(coefs_end)] <- colnames(DM_kov)[-c(which_kat,which_smooth)] 
      }
    }
  }
    
  n_groups <- sapply(1:n_fak,function(j) length(coefs_group[[j]]))  
  help11 <- unlist(sapply(1:n_fak, function(j) seq(1,n_groups[j])))
  help12 <- paste0(rep(colnames(DM_kov)[which_kat],n_groups),help11)
  help13 <- unlist(sapply(1:n_fak,function(j) names(coefs_group[[j]])))
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

  # plot smooth estimate
  if(!is.null(which_smooth) & plot==TRUE){
    plot(mod_opt)
  }
  
  # return 
  return(to_return)
  
}
