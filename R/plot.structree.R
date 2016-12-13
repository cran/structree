#' Plotting Results of Tree-Structured Clustering 
#' 
#' @description
#' Takes a fitted \code{structree} object and plots the results of the tree component of the model.  
#' 
#' @param x Object of class \code{\link[structree]{structree}}.
#' @param select Elements of the tree component that are plotted; 
#' if \code{select} is not specified, by default all components are pictured in one plot. 
#' @param paths If true, the coefficient paths are plotted. 
#' @param result If true, the resulting partition is displayed.
#' @param ask If true, each element chosen by \code{select} is plotted seperately. 
#' @param xlab Label of x-axis.
#' @param ylab Label of y-axis.
#' @param main Title of the plot.
#' @param lwd Linewidth.
#' @param cex.txt Size of the text. 
#' @param cex.axis Size of the axis.
#' @param cex.lab Size of the labels. 
#' @param cex.main Size of title.
#' @param ... Further arguments passed to or from other methods.
#'
#' @details 
#' By default the function pictures the estimated trees against all splits. 
#' If \code{select=NULL} the trees for all the predictors will be plotted. 
#' 
#' @author Moritz Berger <moritz.berger@stat.uni-muenchen.de> \cr \url{http://www.statistik.lmu.de/~mberger/}
#' 
#' @seealso \code{\link[structree]{structree}}
#'
#' @references 
#' Tutz, Gerhard and Berger, Moritz (2015): Tree-Structured Modelling of Categorical Predictors in Regression, 
#' Cornell University Library, arXiv: 1504:04700.
#' 
#' Berger, Moritz and Tutz, Gerhard (2015): Tree-Structured Clustering in Fixed Effects Models, 
#' Cornell University Library, arXiv: 1512.05169. 
#' 
#' @examples 
#' data(rent)
#' 
#' \dontrun{
#' mod <- structree(nmqm~tr(bez)+tr(bj)+tr(rooms)+badkach,data=rent,
#'                  family=gaussian,stop_criterion="CV")
#' 
#' plot(mod, paths=TRUE)
#' }
#' 
#' 
#' @method plot structree
#' @export
#' @importFrom graphics abline axis box layout lines mtext par plot plot.new plot.window text title



plot.structree <-
function(x,
                           select=NULL,
                           paths=FALSE,
                           result=FALSE,
                           ask=FALSE,
                           xlab=NULL,
                           ylab=NULL,
                           main=NULL,
                           lwd=1,
                           cex.txt=1,
                           cex.axis=1,
                           cex.lab=1,
                           cex.main=1,
                           ...){
  
  if(missing(x)){
    stop("Argument x is missing with no default.")
  }
  
  if(x$which_opt==1){
    stop("No split performed, no plot available!")
  }
  
  n_fak <- length(x$partitions)
  
  if(is.null(select)){
    select <- 1:n_fak 
  }
  if(!all(select %in% 1:n_fak)){
    stop("Selection not valid!")
  } 

  if(ask){
    par(ask=TRUE)
  } else{
    rows <- floor(sqrt(length(select)))
    cols <- ceiling(length(select)/rows)
    layout(matrix(1:(rows*cols),nrow=rows,byrow=TRUE))
  }
  
  if(is.null(main)){
    main=names(x$beta_hat)
  }
  if(is.null(xlab)){
    xlab="number of splits"
  }
  if(is.null(ylab)){
    ylab="coefficients"
  }
  
  for(i in select){
    if(paths==TRUE){
      ppaths(x,i,xlab,ylab,main[i],lwd,cex.axis,cex.txt,cex.lab,cex.main)
    } else {
      if(result==TRUE){
      visualfuse(x,i,main[i],cex.txt,cex.main)
      } else{
        ptree(x,i,xlab,main[i],lwd,cex.txt,cex.axis,cex.lab,cex.main)
      }
    }
  }
  
  if(ask){
    par(ask=FALSE)
  } else{
    layout(matrix(1))
  }
  invisible(x)
}
