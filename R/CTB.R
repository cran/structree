#' @name CTB 
#' @title Achievement Test from CTB/McGraw-Hill 
#' @description The data set contains results of an achievement test that measures different objectives and subskills of subjects in mathematics and science. 
#' Inter alia, the students had to respond to 56 multiple-choice items (31 mathematics, 25 science). 
#' For the original description, see Section 5.6 of Chapter 5 in De Boeck and Wilson (2004).   
#' @docType data
#' @usage data(CTB)
#' @format A data frame containing 1211 observations on 9 variables:
#' 
#' \bold{score}  number of correctly solved items (metric)
#'  
#' \bold{school}  school ID (nominal)
#'
#' \bold{size}  number of students in the school, in hundreds (metric)
#' 
#' \bold{bachelor}  transformed and standardized percentage of adults with BA degree or higher in area with school zip code (metric)
#' 
#' \bold{born}  transformed and standardized percentage of adults in the school area who were born in the state where they now reside (metric)
#'  
#' \bold{mortgage} transformed and standardized median of the monthly mortgage in the school area (metric)
#'  
#' \bold{language} transformed and standardized percentage of foreign language households in the school area (metric)
#'  
#' \bold{type} type of school (1: catholic, 2: private, 3: public)
#'  
#' \bold{gender}  gender (0: male, 1: female)
#' @references
#' De Boeck, P. and M. Wilson (2004). Explanatory item response models: A generalized linear and nonlinear approach. Springer Verlag. 
#' @examples 
#' data(CTB)
#' 
#' y <- CTB$score
#' x <- CTB$gender
#' 
#' hist(y)
#' table(x)
#'  
NULL
