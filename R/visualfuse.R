visualfuse <-
function(x,
                       var,
                       main,
                       cex.txt,
                       cex.main){

  partition_org <- x$partitions[[var]]
                            
    
  partition <- partition_org[,x$opts[var]+1]
  
  n_cluster <- length(unique(partition))
        
  which_levels <- names(partition)
  
  cluster <- c()
  
  for(j in 1:n_cluster){
    cluster[j] <- paste(which_levels[which(partition==j)],collapse=", ")
  }
  coefs_org <- format(round(x$coefs_group[[var]],3),nsmall=3)
  coefs <- rep(NA,n_cluster)
  if(length(which(partition==partition[1]))==1){
    coefs[partition[1]] <- 0 
    coefs[is.na(coefs)] <- coefs_org 
  } else{
    coefs <- coefs_org 
  }
  
  plot.new()
  plot.window(ylim=c(1,n_cluster+3),xlim=c(0,5))
  box()
  y_values <- seq(n_cluster+1,1)
  for(k in 1:n_cluster){
    text(0.1,y_values[k],k,cex=cex.txt)
    text(2,y_values[k],cluster[k],cex=cex.txt)
    text(4.5,y_values[k],coefs[k],cex=cex.txt)
    text(2,n_cluster+2.5,"partition",cex=1.2*cex.txt)
    text(4.5,n_cluster+2.5,"coefficients",cex=1.2*cex.txt)
  }
  if(!is.null(main)){
    title(main,cex.main=cex.main)
  }
  invisible(x)
}
