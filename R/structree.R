#' Tree-Structured Clustering 
#' 
#' @description
#' Fusion of categories of ordinal or nominal predictors or fusion of measurement units by tree-structured clustering. 
#' 
#' @param formula Object of class \code{\link{formula}}: a symbolic description of the model to be fitted. See detail. 
#' @param data Data.frame of class \code{\link{data.frame}} containing the variables of the model. 
#' @param family a description of the error distribution and link function to be used in the model. 
#' This can be a character string naming a family function, a family function or the result of a call to a family function. 
#' See \code{\link{family}} for details of family functions. 
#' @param stop_criterion Criterion to determine the optimal number of splits in the tree component of the model; 
#' one out of \code{"AIC"}, \code{"BIC"}, \code{"CV"} and \code{"pvalue"}.
#' @param splits_max Maximal number of splits in the tree component. 
#' @param fold Number of folds; only for stop criterion \code{"CV"}.
#' @param alpha Significance level; only for stop criterion \code{"pvalue"}.
#' @param grid_value An optional parameter; \code{grid_value} is a scalar giving the minimal distance between 
#' two adjacent observation units that are used as candidates for splitting; only for repeated measurements.
#' @param min_border An optional parameter; \code{min_border} is a integer giving the minimal size of the outer 
#' nodes of the tree; only for repeated measurements.
#' @param ridge If true, a small ridge penalty is added to obtain the order of measurement units; only for repeated measurements. 
#' @param constant_covs Must be set to true, if constant covariates are available; only for repeated measurments 
#' (currently only available for Gaussian response).  
#' @param trace If true, information about the estimation progress is printed.
#' @param plot If true, the smooth components of the model are plottet; only for categorical predictors.  
#' @param k Dimension of the B-spline basis that is used to fit smooth components. For details see \code{\link[mgcv]{s}};
#' only for categorical predictors.  
#' @param weights An optional vector of prior weights to be used in the fitting process; see also \code{\link{glm}}.
#' @param offset An a priori known component to be included in the linear predictor during fitting; see also \code{\link{glm}}.
#' @param x,object Object of class \code{"structree"}.
#' @param ... Further arguments passed to or from other methods.
#' 
#' @details 
#' A typical \link{formula} has the form \code{response ~ predictors}, where \code{response} is the name of the response variable 
#' and \code{predictors} is a series of terms that specify the predictor of the model. 
#' 
#' For an ordinal or nominal predictors z one has to enter \code{tr(x)} into the formula. 
#' 
#' For smooth components x one has to enter \code{s(x)} into the formula; currently not implemented for repeated measurements. 
#' 
#' For fixed effects z of observation units u one has to enter \code{tr(z|u)} into the formula. 
#' An unit-specific intercept is specified by \code{tr(1|u)}.
#' 
#' The framework only allows for categorical predictors or observations units in the tree component, but not both. 
#' All other predictors with a linear term are entered as usual by \code{x1+...+xp}. 
#' 
#' 
#' @return 
#' Object of class \code{"structree"}. 
#' An object of class \code{"structree"} is a list containing the following components:
#' 
#' \item{coefs_end}{all coefficients of the estimated model}
#' \item{partitions}{list of matrices containing the partitions of the predictors in the tree component including all iterations}
#' \item{beta_hat}{list of matrices with the fitted coefficients in the tree component including all iterations}
#' \item{which_opt}{number of the optimal model (total number of splits-1)} 
#' \item{opts}{number of splits per predictor in the tree component}
#' \item{order}{list of ordered split-points of the predictors in the tree component}
#' \item{tune_values}{value of the stopping criterion that determine the optimal model}
#' \item{group_ID}{list of the group IDs for each observations}
#' \item{coefs_group}{list of coefficients of the estimated model}
#' \item{y}{Response vector}
#' \item{DM_kov}{Design matrix}
#' 
#' @author Moritz Berger <Moritz.Berger@imbie.uni-bonn.de> \cr \url{http://www.imbie.uni-bonn.de/personen/dr-moritz-berger/}
#' 
#' @seealso \code{\link[structree]{plot.structree}}
#'
#' @references 
#' Tutz, Gerhard and Berger, Moritz (2018): Tree-structured modelling of categorical predictors in regression, 
#' Advances in Data Analysis and Classification 12(3), 737-758.
#' 
#' Berger, Moritz and Tutz, Gerhard (2018): Tree-structured clustering in fixed effects models, 
#' Journal of Computational and Graphical Statistics 27(2), 380-392.
#' 
#' @examples 
#' data(rent)
#' 
#' \dontrun{
#' mod <- structree(nmqm~tr(bez)+tr(bj)+tr(rooms)+badkach0,data=rent,
#'                  family=gaussian,stop_criterion="CV")
#' 
#' print(mod)
#' coef(mod)
#' }
#' 
#' @export
#' @import mgcv
#' @import lme4
#' @importFrom penalized penalized coefficients 
#' @importFrom stats Gamma coef deviance fitted formula gaussian glm logLik pchisq predict update
#' @importFrom utils tail


structree <-
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
  UseMethod("structree")
}

#' @rdname structree 
#' @method print structree
#' @export

print.structree <-
  function(x,
           ...){
    
    type <- attr(x,"type")
    
    if(type=="catPred"){
      zv <- attr(x,"which_kat")
      xv <- attr(x,"which_smooth")
      xv <- c(names(x$DM_kov)[!(names(x$DM_kov) %in% c(xv,zv))],xv)
    } else{
      slu <- attr(x,"secondlevel")
      zv   <- c("Intercept",attr(x,"slope"))
      xv   <- names(x$DM_kov)[!(names(x$DM_kov) %in% c(slu,zv))]
    }
    x$opts <- x$opts+1
    opts        <- x$opts
    names(opts) <- zv
    
    if(type=="fixedEff"){
      cat("Tree Structured Clustering of observation units:\n")
    } else{
      cat("Tree Structured Clustering of categorical predictors:\n")
    }  
    cat("\n")
    cat("Call:",paste(strwrap(paste(deparse(x$call)),width=100),collapse="\n"),"\n")
    cat("\n")
    if(type=="catPred"){
      cat("Variables in tree component:",paste0(zv,collapse=", "),"\n")
      cat("Variables in parametric component:",paste0(xv,collapse=", "),"\n")
    } else{
      cat("Second-level unit:",slu,"\n")
      cat("Unit specific effects for:",paste0(zv,collapse=", "),"\n")
      cat("Fixed effects for:",paste0(xv,collapse=", "),"\n")
    }
    cat("Number of Splits:",x$which_opt-1,"\n")
    cat("\n")
    cat("Number of Groups:\n")
    print(opts)
    
    invisible(x)
  }

#' @rdname structree
#' @method coef structree
#' @export

coef.structree <-
  function(object,
           ...){
    
    if(is.null(attr(object,"which_smooth"))){
      coefficients <- object$coefs_end
    } else{
      coefficients <- object$coefs_end
      k <- attr(object,"k")
      m <- summary(object$model)$m
      if(length(k)==1){
        k <- rep(k,m)
      }
      coefs_s <- tail(object$model$coefficients,sum(k-1))
      which_smooth <- attr(object,"which_smooth")
      v <- unlist(lapply(1:m,function(j) seq(1,k[j]-1)))
      names(coefs_s) <- paste0(rep(paste0("s(",which_smooth,")."),k-1),v)
      coefficients <- c(coefficients,coefs_s)
    }
    return(coefficients)
}