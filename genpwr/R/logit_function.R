#' Logit Function
#'
#' Calculates the logit of a specified value. 
#'
#' @param x a number between 0 and 1.  
#' 
#' @return The logit of x. 
#'
#' @export
#'
logit<-
  function (x, min = 0, max = 1) 
{
  p <- (x - min)/(max - min)
  log(p/(1 - p))
}