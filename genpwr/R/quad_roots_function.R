#' Function to Solve Quadratic Equations
#'
#' Finds the positive root of a quadratic equation ax^2 + bx + c . 
#'
#' @param a the coefficient for x^2
#' @param b the coefficient for x
#' @param c the constant
#' 
#' @return The positive root of the quadratic equation ax^2 + bx + c 
#'
#' @examples 
#' pw<-power.calc(N=c(1000,2000), Case.Rate=c(0.5), k=NULL, MAF=seq(0.05, 0.1, 0.01), OR=c(3,4),Alpha=c(0.05), True.Model='All', Test.Model='All')
#'
#' @export
#'
quad_roots<-function(a,b,c){
  c(((-b-sqrt(b^2-4*a*c))/(2*a)),((-b+sqrt(b^2-4*a*c))/(2*a)))
}
