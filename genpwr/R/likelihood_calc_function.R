###########################################################################
#Functions to calculate the likelihood for each testing model
#Input MLE's and a 2x3 table
###########################################################################

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

null.ll<-function(t){
  Case.Rate<-sum(t[1,])/sum(t)
  beta0 = logit(Case.Rate)
  ll<- Case.Rate*log(exp(beta0)/(1+exp(beta0))) + (1-Case.Rate)*log(1/(1+exp(beta0))) 
  return(ll)
}

#Convenience function to calculate all of above likelihoods
calc.like<-function(beta, t, model){
  if(model=='Dominant'){ll <- dominant.ll(beta,t)}
  if(model=='Additive'){ll <- additive.ll(beta,t)}
  if(model=='Recessive'){ll <- recessive.ll(beta,t)}
  if(model=='2df'){ll <- df2.ll(beta, t)}
  if(model=='null'){ll <- null.ll(t)}
  return(ll)
}