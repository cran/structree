ppaths <-
function(x,
                   var,
                   xlab,
                   ylab,
                   main,
                   lwd,
                   cex.axis,
                   cex.txt,
                   cex.lab,
                   cex.main){
  
  betas <- x$beta_hat

    beta <- betas[[var]] 
    n_splits <- ncol(beta)-1
    plot(beta[1,],type="l",ylim=c(min(beta),max(beta)),ylab=ylab,xlab=xlab,main=main,xlim=c(0,1.4*n_splits),axes=FALSE,lwd=lwd,cex.lab=cex.lab,
         cex.main=cex.main)
    box()
    for(i in 2:nrow(beta)){
      lines(beta[i,],lwd=lwd)
    }
    text(rep(1.2*n_splits,nrow(beta)),beta[,(n_splits+1)],rownames(beta),cex=cex.txt)
    abline(v=x$which_opt,lty="dashed",lwd=lwd)
    axis(1,round(seq(0,n_splits,by=n_splits/4),0),at=round(seq(1,n_splits+1,by=n_splits/4),0),cex.axis=cex.axis)
    axis(2,cex.axis=cex.axis)
  
  invisible(x)  
}
