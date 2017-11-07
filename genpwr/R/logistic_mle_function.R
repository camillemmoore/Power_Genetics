####################################################################################
#Function to calculate MLE's for logistic models
#Specify a table and the testing model (Dominant, Recessive, Additive, 2df)
####################################################################################

logistic.mles<-function(t,model){
  
  Case.Rate <- sum(t[1,])/sum(t)
  N_AA_case <- t[1,1]
  N_AB_case <- t[1,2]
  N_BB_case <- t[1,3]
  
  N_AA <- sum(t[,1])
  N_AB <- sum(t[,2])
  N_BB <- sum(t[,3])
  
  if (model=='null'){#Null MOdel
    beta = logit(Case.Rate)
  }
  
  #Dominant
  if (model=='Dominant'){
    beta0 = logit(N_AA_case/N_AA)
    beta1 = logit((N_AB_case + N_BB_case)/(N_AB + N_BB)) - beta0
    beta = c(beta0, beta1)
  }
  
  #Recessive
  if (model=='Recessive'){
    beta0 = logit((N_AA_case+N_AB_case)/(N_AA+N_AB))
    beta1 = logit((N_BB_case)/(N_BB)) - beta0
    beta = c(beta0, beta1)
  }
  
  #2DF
  if (model=='2df'){
    beta0 = logit(N_AA_case/N_AA)
    beta1 = logit((N_AB_case)/(N_AB)) - beta0
    beta2 = logit((N_BB_case)/(N_BB)) - beta0
    beta = c(beta0, beta1, beta2)
  }
  
  #Additive Model - need to use an optimizer to solve
  if (model=='Additive'){
    beta<-optim(c(0,0),function(x) -additive.ll(x,t), control=c(abstol = 0.00001))$par
  }
  return(beta)
}
