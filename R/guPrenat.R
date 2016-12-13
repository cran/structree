#' @name guPrenat
#' @title Prenatal Care in Guatemala 
#' @description A data set derived from the National Survey of Maternal and Child Health in Guatemala in 1987. 
#' The data contains observations of children that were born in the 5-year period before the survey. 
#' @docType data
#' @usage data(guPrenat)
#' @format A data frame containing 1211 observations on 9 variables:
#' 
#' \bold{cluster}  community (nominal)
#'  
#' \bold{prenat}  prenatal care (0: traditional, 1: modern)
#'
#' \bold{motherAge}  mother 25 years or older (0: no, 1: yes)
#' 
#' \bold{indig}  mother's ethnicity (nominal)
#' 
#' \bold{momEd}  mother's level of education (nominal)
#'  
#' \bold{husEd} husband's level of education (nominal)
#'  
#' \bold{husEmpl} husband's employment status (nominal)
#'  
#' \bold{toilet} modern toilet in house (0: no, 1: yes)
#'  
#' \bold{TV} frequency of TV usage (nominal)
#'  
#' @references 
#' Rodriguez, Germa'n and Goldman, Noreen (1995), "Improved estimation procedures for multilevel
#' models with binary response: a case-study", Journal of the Royal Statistical Society, Series A, 164,
#' 339-355.
#' 
#' Douglas Bates and Martin Maechler and Ben Bolker (2014). mlmRev: Examples from Multilevel Modelling Software Review. R
#' package version 1.0-6. https://CRAN.R-project.org/package=mlmRev
#' @examples 
#' data(guPrenat)
#' 
#' y <- guPrenat$prenat
#' community <- guPrenat$cluster
#' 
#' table(y)
#' hist(table(community))
#'  
NULL
