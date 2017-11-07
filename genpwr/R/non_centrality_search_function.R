####################################################################
#Search for Non-Centrality Parameter for SS calcualtion for 2Df Test
####################################################################


#' Function to Determine Non-Centrality Parameter of the Chi-squared distribution
#'
#' This function is set to 0 and solved for x, the non-centrality parameter to determine the sample size in \code{\link{ss.calc}}
#'
#' @param x the non-centrality parameter 
#' @param power the desired power
#' @param stat the detectable likelihood ratio test statistic
#' @param Alpha the desired type 1 error rate 
#' @param df the degrees of freedom for the likelihood ratio test
#' 
#' @return numeric value of the function
#'
#' @export
#'
ncp.search<-
  function(x, power, stat, Alpha, df){
  (1-power)-pchisq(qchisq(1-Alpha, df, ncp=0), df=1, ncp = x) }
