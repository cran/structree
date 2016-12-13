ptree <-
function(x,
                  var,
                  xlab,
                  main,
                  lwd,
                  cex.txt,
                  cex.axis,
                  cex.lab,
                  cex.main){
  
 partitions <- x$partitions
 partition_org <- partitions[[var]]
 
 n_splits <- ncol(partition_org)-1
 n_cluster <- length(unique(partition_org[,ncol(partition_org)]))
 n_levels <- nrow(partition_org)
 
 x_grid <-  seq(1,(n_splits+5))
 y_grid <-  seq(1,(n_cluster+1))
 
 
 
 plot.new()
 plot.window(ylim=c(y_grid[1],y_grid[length(y_grid)]),xlim=c(x_grid[1],x_grid[length(x_grid)]))
 box()
 
 partition <- as.matrix(partition_org[x$order[[var]],])
 
   hilfspunkte <- list() 
   split_points <- rep(0,n_splits)
   
   hilfspunkte[[n_splits+2]] <- matrix(0,nrow=n_cluster,ncol=2)
   hilfspunkte[[n_splits+2]][,1] <- rep(n_splits+2,n_cluster)
   hilfspunkte[[n_splits+2]][,2] <- n_cluster:1
   
 
   if(n_splits>1){ 
     seq <- n_splits:2
     
     for(i in seq){
       where_split <- partition[min(which(partition[,i]!=partition[,i+1])),i]
       split_kats  <- c(where_split,where_split+1)
       
       which_fused <- unique(partition[which(partition[,i+1] %in% split_kats),(n_splits+1)])
       split_points[i] <- rownames(partition)[max(which(partition[,i+1]==where_split))]
       
       
       hilfspunkte[[i+1]] <- matrix(0,nrow=n_cluster,ncol=2)
       hilfspunkte[[i+1]][,1] <- rep(n_splits-which(seq==i)+2,n_cluster)
       hilfspunkte[[i+1]][-which_fused,2] <- hilfspunkte[[i+2]][-which_fused,2]
       hilfspunkte[[i+1]][which_fused,2]  <- sum(hilfspunkte[[i+2]][which_fused,2])/length(which_fused)
     }
   }
   
   hilfspunkte[[2]] <- matrix(2,nrow=n_cluster,ncol=2)
   hilfspunkte[[2]][,2] <- sum(hilfspunkte[[3]][,2])/n_cluster
   split_points[1] <- rownames(partition)[max(which(partition[,2]==1))]
   
   hilfspunkte[[1]] <- matrix(1,nrow=n_cluster,ncol=2)
   hilfspunkte[[1]][,2] <- hilfspunkte[[2]][,2]
   
   for(j in length(hilfspunkte):2){
    for(i in 1:n_cluster){
      lines(c(hilfspunkte[[j-1]][i,1],hilfspunkte[[j]][i,1]),c(hilfspunkte[[j-1]][i,2],hilfspunkte[[j]][i,2]),
            lwd=lwd)
    }
   }
   
   if(!is.null(main)){
    title(main,cex.main=cex.main)
   }
   
   text((n_splits+3)/2,y_grid[length(y_grid)],"partition",cex=cex.txt)
   if(attr(x,"type")=="catPred"){
     text(n_splits+4,y_grid[length(y_grid)],"category",cex=cex.txt)
   } else{
     text(n_splits+4,y_grid[length(y_grid)],"unit",cex=cex.txt)
   }
   text(rep(n_splits+4,n_cluster),y_grid[n_cluster:1],rownames(partition),cex.txt)
   
   partition_opt <- x$opts[var]
   lines(c(partition_opt+2,partition_opt+2),c(0.5,n_cluster+0.5),lty="dashed",lwd=lwd)
   
   
   punkte_b <- list()
   for(i in (n_splits+2):3){
     punkte_b[[i-2]] <- hilfspunkte[[i]][which(hilfspunkte[[i-1]][,2]!=hilfspunkte[[i]][,2]),]
   }
   help <- matrix(NA,nrow=n_splits,ncol=2)
   help[,1] <- 3:(n_splits+2)
   help[,2] <- sapply(1:n_splits,function(j) sum(punkte_b[[j]][,2])/nrow(punkte_b[[j]]))
   
   text_function <- function(i,cex.txt){
     text(help[i,1],help[i,2],paste("<=",split_points[i]),cex=cex.txt)
   }
   do.call(text_function,list(c(1:n_splits),cex.txt=cex.txt))
   axis(1,1:n_splits,at=3:(n_splits+2),cex.axis=cex.axis)
   if(!is.null(xlab)){
    mtext(xlab,line=3,side=1,at=(n_splits+5)/2,cex=cex.lab)
   }
}
