#' Logit Function
#'
#' Calculates the logit of a specified value. 
#'
#' @param x a number between 0 and 1.  
#' @param min minimum
#' @param max maximum
#' 
#' @return The logit of x. 
#'
#' @examples
#' logit(0.5)
#'
#' @export
#'
logit<-
  function (x, min = 0, max = 1) 
{
  p <- (x - min)/(max - min)
  if(p/(1-p) < 0){
  	return(NA)
  }else{return(log(p/(1 - p)))}
  
}
