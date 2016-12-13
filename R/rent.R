#' @name rent
#' @title Munich Rent Data  
#' @description The data set is part of the Munich rent index in 2003. It is available from the data archive 
#' of the Department of Statistics at the University of Munich (http://www.statistik.lmu.de/service/datenarchiv). 
#' @docType data
#' @usage data(rent)
#' @format  A data frame containing 2053 observations on 11 variables:
#' 
#'  \bold{nmqm}   net rent per square meter (metric)
#'  
#'  \bold{wfl}  floor space (metric)
#'
#'  \bold{rooms}  number of rooms (ordinal)
#' 
#'  \bold{bj}  year of construction (ordinal)
#' 
#'  \bold{bez}  residential area (norminal)
#'  
#'  \bold{ww0} hot water supply (1: no, 0: yes)
#'  
#'  \bold{zh0} central heating (1: no, 0: yes)
#'  
#'  \bold{badkach0} tiled bathroom (1: no, 0: yes)
#'  
#'  \bold{badextra} supplementary equipment in bathroom (1: yes, 0: no)
#'  
#'  \bold{kueche} well equipped kitchen (1: yes, 0: no)
#'  
#'  \bold{quality} quality of residential area (ordinal)
#' @references 
#' Fahrmeir, L. and Kuenstler, R. and Pigeot, I. and Tutz, G. (2004): Statistik: der Weg zur Datenanalyse. 5. Auflage, Springer, Berlin.
#' @examples 
#' data(rent)
#' 
#' y <- rent$nmqm
#' X <- rent[,-1]
#' 
#' boxplot(y)
#' summary(X)
#'  
NULL
