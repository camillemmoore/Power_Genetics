

####################################################################################
# odds ratios for different models
####################################################################################

#' Odds ratio calculation
#'
#' Calculates odds ratio for given parameters. Used by the function odds_ratio_function.
#'
#' @param a The probability of a case given homogeneity for the major allele
#' @param b The probability of a case given heterozygosity 
#' @param c The probability of a case given homogeneity for the minor allele
#' @param d The probability of a control given homogeneity for the major allele
#' @param e The probability of a control given heterozygosity 
#' @param f The probability of a control given homogeneity for the minor allele#' 
#' @param mod The model to be used (eg "Dominant", "Recessive", etc)
#' @param risk_allele Is this allele a risk allele? use T or F
#' @return Odds ratio
#'
#' @examples
#' or_calc(a = 0.3649185, b = 0.12797197, c = 0.007109554, d= 0.4450815, e= 0.05202803, f = 0.002890446, mod = "Dominant", risk_allele = T)
#'
#' @export
#'
or_calc <- function(a,b,c,d,e,f, mod, risk_allele)
{
  if(mod == "Recessive"){
    if(risk_allele) res <- (c(c*e/(b*f))[c(c*e/(b*f)) >= 1])
    if(!risk_allele) res <- (c(c*e/(b*f))[c(c*e/(b*f)) <= 1])
  }
  if(mod == "Dominant"){
    if(risk_allele) res <- (c(b*d/(a*e))[c(b*d/(a*e)) >= 1])
    if(!risk_allele) res <- (c(b*d/(a*e))[c(b*d/(a*e)) <= 1])
  }
  if(mod == "Additive"){
    if(risk_allele) res <- (c(c(b*d/(a*e))[b*d/(a*e) >= 1]))
    if(!risk_allele) res <- (c(c(b*d/(a*e))[b*d/(a*e) <= 1]))
  }
  if(length(res) == 0) res <- NA
  return(res)
}