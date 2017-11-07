###########################################################################
#Functions to calculate the likelihood for each testing model
#Input MLE's and a 2x3 table
###########################################################################


#' Function to Calculate Additive Log Likelihood
#'
#' Calculates the log likelihood for a given set of coefficients under an additive genetic model. 
#'
#' @param beta Vector of logistic regression coefficients.
#' @param t A 2x3 table of joint probabilities of disease and genotype.  Rows = case vs. control, columns=genotype. 
#' 
#' @return The log likelihood. 
#'
#' @export
#'
additive.ll<-function(beta, t){
  beta0 <- beta[1]
  beta1 <-beta[2]
  
  N_AA_case <- t[1,1]
  N_AB_case <- t[1,2]
  N_BB_case <- t[1,3]
  N_AA_control <- t[2,1]
  N_AB_control <- t[2,2]
  N_BB_control <- t[2,3]
  
  ll<- N_AA_case*log(exp(beta0)/(1+exp(beta0))) + 
    N_AB_case*log(exp(beta0+beta1)/(1+exp(beta0+beta1)))+
    N_BB_case*log(exp(beta0+2*beta1)/(1+exp(beta0+2*beta1)))+
    N_AA_control*log(1/(1+exp(beta0))) + 
    N_AB_control*log(1/(1+exp(beta0+beta1)))+
    N_BB_control*log(1/(1+exp(beta0+2*beta1)))
  
  return(ll)
}

#' Function to Calculate Dominant Log Likelihood
#'
#' Calculates the log likelihood for a given set of coefficients under a dominant genetic model. 
#'
#' @param beta Vector of logistic regression coefficients.
#' @param t A 2x3 table of joint probabilities of disease and genotype.  Rows = case vs. control, columns=genotype. 
#' 
#' @return The log likelihood. 
#'
#' @export
#'
dominant.ll<-function(beta, t){
  beta0 <- beta[1]
  beta1 <- beta[2]
  
  N_AA_case <- t[1,1]
  N_AB_case <- t[1,2]
  N_BB_case <- t[1,3]
  N_AA_control <- t[2,1]
  N_AB_control <- t[2,2]
  N_BB_control <- t[2,3]
  
  ll <- N_AA_case*log(exp(beta0)/(1+exp(beta0))) + 
    N_AB_case*log(exp(beta0+beta1)/(1+exp(beta0+beta1)))+
    N_BB_case*log(exp(beta0+beta1)/(1+exp(beta0+beta1)))+
    N_AA_control*log(1/(1+exp(beta0))) + 
    N_AB_control*log(1/(1+exp(beta0+beta1)))+
    N_BB_control*log(1/(1+exp(beta0+beta1)))
  
  return(ll)
}

#' Function to Calculate Recessive Log Likelihood
#'
#' Calculates the log likelihood for a given set of coefficients under a recessive genetic model. 
#'
#' @param beta Vector of logistic regression coefficients.
#' @param t A 2x3 table of joint probabilities of disease and genotype.  Rows = case vs. control, columns=genotype. 
#' 
#' @return The log likelihood. 
#'
#' @export
#'
recessive.ll<-function(beta, t){
  beta0 <- beta[1]
  beta1 <- beta[2]
  
  N_AA_case <- t[1,1]
  N_AB_case <- t[1,2]
  N_BB_case <- t[1,3]
  N_AA_control <- t[2,1]
  N_AB_control <- t[2,2]
  N_BB_control <- t[2,3]
  
  ll <- N_AA_case*log(exp(beta0)/(1+exp(beta0))) + 
    N_AB_case*log(exp(beta0)/(1+exp(beta0)))+
    N_BB_case*log(exp(beta0+beta1)/(1+exp(beta0+beta1)))+
    N_AA_control*log(1/(1+exp(beta0))) + 
    N_AB_control*log(1/(1+exp(beta0)))+
    N_BB_control*log(1/(1+exp(beta0+beta1)))
  
  return(ll)
}

#' Function to Calculate 2df Log Likelihood
#'
#' Calculates the log likelihood for a given set of coefficients under an unspecificed/2df genetic model. 
#'
#' @param beta Vector of logistic regression coefficients.
#' @param t A 2x3 table of joint probabilities of disease and genotype.  Rows = case vs. control, columns=genotype. 
#' 
#' @return The log likelihood. 
#'
#' @export
#'
df2.ll<-function(beta, t){
  beta0 <- beta[1]
  beta1 <- beta[2]
  beta2 <- beta[3]
  
  N_AA_case <- t[1,1]
  N_AB_case <- t[1,2]
  N_BB_case <- t[1,3]
  N_AA_control <- t[2,1]
  N_AB_control <- t[2,2]
  N_BB_control <- t[2,3]
  
  ll <- N_AA_case*log(exp(beta0)/(1+exp(beta0))) + 
    N_AB_case*log(exp(beta0+beta1)/(1+exp(beta0+beta1)))+
    N_BB_case*log(exp(beta0+beta2)/(1+exp(beta0+beta2)))+
    N_AA_control*log(1/(1+exp(beta0))) + 
    N_AB_control*log(1/(1+exp(beta0+beta1)))+
    N_BB_control*log(1/(1+exp(beta0+beta2)))
  
  return(ll)
}

#' Function to Calculate Null Log Likelihood
#'
#' Calculates the log likelihood for a given set of coefficients under the null. 
#'
#' @param t A 2x3 table of joint probabilities of disease and genotype.  Rows = case vs. control, columns=genotype. 
#' 
#' @return The log likelihood. 
#'
#' @export
#'
null.ll<-function(t){
  Case.Rate<-sum(t[1,])/sum(t)
  beta0 = logit(Case.Rate)
  ll<- Case.Rate*log(exp(beta0)/(1+exp(beta0))) + (1-Case.Rate)*log(1/(1+exp(beta0))) 
  return(ll)
}

#' Function to Calculate Log Likelihood
#'
#' Convenience function to calculate the log likeilhood of a specified model.
#'
#' @param beta Vector of logistic regression coefficients.
#' @param t A 2x3 table of joint probabilities of disease and genotype.  Rows = case vs. control, columns=genotype. 
#' @param model The genetic model in the logisitic regression: "Dominant", "Additive", "Recessive", "2df" or "null"
#' 
#' @return The log likelihood. 
#'
#' @export
#'
calc.like<-function(beta, t, model){
  if(model=='Dominant'){ll <- dominant.ll(beta,t)}
  if(model=='Additive'){ll <- additive.ll(beta,t)}
  if(model=='Recessive'){ll <- recessive.ll(beta,t)}
  if(model=='2df'){ll <- df2.ll(beta, t)}
  if(model=='null'){ll <- null.ll(t)}
  return(ll)
}