###########################################################################
#Functions to calculate the likelihood for each testing model
#For Logisitic Models: Input MLE's and a 2x3 table
#For Linear Models: Input MLE's, vector of effect sizes, minor allele frequency, and SD of Y given X
###########################################################################


#' Function to Calculate Additive Log Likelihood for a Logistic Regression Model
#'
#' Calculates the log likelihood for a given set of logistic regression coefficients under an additive genetic model.
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

#' Function to Calculate Dominant Log Likelihood for a Logistic Regression Model
#'
#' Calculates the log likelihood for a given set of logistic regression coefficients under a dominant genetic model.
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

#' Function to Calculate Recessive Log Likelihood for a Logistic Regression Model
#'
#' Calculates the log likelihood for a given set of logistic regression coefficients under a recessive genetic model.
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

#' Function to Calculate 2df Log Likelihood for a Logistic Regression Model
#'
#' Calculates the log likelihood for a given set of logistic regression coefficients under an unspecificed/2df genetic model.
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

#' Function to Calculate Null Log Likelihood for a Logistic Regression Model
#'
#' Calculates the log likelihood for a given set of logistic regression coefficients under the null.
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

#' Function to Calculate Log Likelihood for a Logistic Regression Model
#'
#' Convenience function to calculate the log likelihood of a specified model.
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

##################################################################
#Functions to calculate expected likelihoods for linear models
##################################################################

#' Function to Calculate Expected Log Likelihood for a Single Genotype
#'
#' Calculates the expected log likelihood for a single genotype given the true and estimated mean and standard deviation for the outcome.
#'
#' @param mean_truth Mean of the outcome given X(predictors/genotype) under the true model.
#' @param mean_model Mean of the outcome given X(predictors/genotype) under the test model.
#' @param sd_y_x_model The standard deviation of Y (the outcome) given X (predictors/genotype) under the test model.
#' @param sd_y_x_truth The standard deviation of Y given X (predictors/genotype) given genotype under the true model.
#'
#' @return The log likelihood.
#'
#' @export
#'
expected.linear.ll <- function(mean_truth, mean_model, sd_y_x_truth, sd_y_x_model){
  -0.5*log(2*pi*sd_y_x_model*sd_y_x_model)-
  (sd_y_x_truth*sd_y_x_truth + mean_truth^2)/(2*sd_y_x_model*sd_y_x_model)+
  mean_model*mean_truth/(sd_y_x_model*sd_y_x_model)-
  mean_model^2/(2*sd_y_x_model*sd_y_x_model)
  }

#' Function to Calculate Expected Null Log Likelihood for a Linear Regression Model
#'
#' Calculates the expected log likelihood for a given set of linear regression coefficients under the null.
#'
#' @param beta Vector of linear regression coefficients.
#' @param m Minor allele frequency.
#' @param es Vector of effect sizes with two elements, (mean AB - mean AA) and (mean BB - mean AA).
#' @param sd_y_x_model The standard deviation of Y (the outcome) given X (predictors/genotype) under the test model.
#' @param sd_y_x_truth The standard deviation of Y given X (predictors/genotype) given genotype under the true model.
#'
#' @return The log likelihood.
#'
#' @export
#'
null.ll.linear<-function(beta, m, es, sd_y_x_model, sd_y_x_truth){
  ll<- ((1-m)^2)*expected.linear.ll(mean_truth=0, mean_model=beta, sd_y_x_truth, sd_y_x_model)+
      2*m*(1-m)*expected.linear.ll(mean_truth=es[1], mean_model=beta, sd_y_x_truth, sd_y_x_model)+
    (m^2)*expected.linear.ll(mean_truth=es[2], mean_model=beta, sd_y_x_truth, sd_y_x_model)

      return(ll)
}


#' Function to Calculate Additive Log Likelihood for a Linear Regression Model
#'
#' Calculates the log likelihood for a given set of linear regression coefficients under an additive genetic model.
#'
#' @param beta Vector of linear regression coefficients.
#' @param m Minor allele frequency.
#' @param es Vector of effect sizes with two elements, (mean AB - mean AA) and (mean BB - mean AA).
#' @param sd_y_x_model The standard deviation of Y (the outcome) given X (predictors/genotype) under the test model.
#' @param sd_y_x_truth The standard deviation of Y given X (predictors/genotype) given genotype under the true model.
#'
#' @return The log likelihood.
#'
#' @export
#'
additive.ll.linear<-function(beta, m, es, sd_y_x_model, sd_y_x_truth){
  beta0 = beta[1]
  beta1 = beta[2]

  ll<- ((1-m)^2)*expected.linear.ll(mean_truth=0, mean_model=beta0, sd_y_x_truth, sd_y_x_model)+
  2*m*(1-m)*expected.linear.ll(mean_truth=es[1], mean_model=beta0+beta1, sd_y_x_truth, sd_y_x_model)+
    (m^2)*expected.linear.ll(mean_truth=es[2], mean_model=beta0+beta1*2, sd_y_x_truth, sd_y_x_model)

  return(ll)
}

#' Function to Calculate Dominant Log Likelihood for a Linear Regression Model
#'
#' Calculates the log likelihood for a given set of linear regression coefficients under a dominant genetic model.
#'
#' @param beta Vector of linear regression coefficients.
#' @param m Minor allele frequency.
#' @param es Vector of effect sizes with two elements, (mean AB - mean AA) and (mean BB - mean AA).
#' @param sd_y_x_model The standard deviation of Y (the outcome) given X (predictors/genotype) under the test model.
#' @param sd_y_x_truth The standard deviation of Y given X (predictors/genotype) given genotype under the true model.
#'
#' @return The log likelihood.
#'
#' @export
#'
dominant.ll.linear<-function(beta, m, es, sd_y_x_model, sd_y_x_truth){
  beta0 = beta[1]
  beta1 = beta[2]

  ll<- ((1-m)^2)*expected.linear.ll(mean_truth=0, mean_model=beta0, sd_y_x_truth, sd_y_x_model)+
  2*m*(1-m)*expected.linear.ll(mean_truth=es[1], mean_model=beta0+beta1, sd_y_x_truth, sd_y_x_model)+
    (m^2)*expected.linear.ll(mean_truth=es[2], mean_model=beta0+beta1, sd_y_x_truth, sd_y_x_model)

  return(ll)
}

#' Function to Calculate Recessive Log Likelihood for a Linear Regression Model
#'
#' Calculates the log likelihood for a given set of linear regression coefficients under a recessive genetic model.
#'
#' @param beta Vector of linear regression coefficients.
#' @param m Minor allele frequency.
#' @param es Vector of effect sizes with two elements, (mean AB - mean AA) and (mean BB - mean AA).
#' @param sd_y_x_model The standard deviation of Y (the outcome) given X (predictors/genotype) under the test model.
#' @param sd_y_x_truth The standard deviation of Y given X (predictors/genotype) given genotype under the true model.
#'
#' @return The log likelihood.
#'
#' @export
#'
recessive.ll.linear<-function(beta, m, es, sd_y_x_model, sd_y_x_truth){
  beta0 = beta[1]
  beta1 = beta[2]
  ll<- ((1-m)^2)*expected.linear.ll(mean_truth=0, mean_model=beta0, sd_y_x_truth, sd_y_x_model)+
  2*m*(1-m)*expected.linear.ll(mean_truth=es[1], mean_model=beta0, sd_y_x_truth, sd_y_x_model)+
    (m^2)*expected.linear.ll(mean_truth=es[2], mean_model=beta0+beta1, sd_y_x_truth, sd_y_x_model)

  return(ll)
}

#' Function to Calculate 2 Degree of Freedom Log Likelihood for a Linear Regression Model
#'
#' Calculates the log likelihood for a given set of linear regression coefficients under a the 2df model.
#'
#' @param beta Vector of linear regression coefficients.
#' @param m Minor allele frequency.
#' @param es Vector of effect sizes with two elements, (mean AB - mean AA) and (mean BB - mean AA).
#' @param sd_y_x_model The standard deviation of Y (the outcome) given X (predictors/genotype) under the test model.
#' @param sd_y_x_truth The standard deviation of Y given X (predictors/genotype) given genotype under the true model.
#'
#' @return The log likelihood.
#'
#' @export
#'
df2.ll.linear<-function(beta, m, es, sd_y_x_model, sd_y_x_truth){
  beta0 = beta[1]
  beta1 = beta[2]
  beta2 = beta[3]
  ll<- ((1-m)^2)*expected.linear.ll(mean_truth=0, mean_model=beta0, sd_y_x_truth, sd_y_x_model)+
  2*m*(1-m)*expected.linear.ll(mean_truth=es[1], mean_model=beta0+beta1, sd_y_x_truth, sd_y_x_model)+
    (m^2)*expected.linear.ll(mean_truth=es[2], mean_model=beta0+beta2, sd_y_x_truth, sd_y_x_model)
  return(ll)
}

#' Function to Calculate Log Likelihood for a Linear Regression Model
#'
#' Convenience function to calculate the log likelihood of a specified model.
#'
#' @param beta Vector of linear regression coefficients.
#' @param m Minor allele frequency.
#' @param es_ab effect size for mean AB - mean AA 
#' @param es_bb effect size for mean BB - mean AA 
#' @param sd_y_x_model The standard deviation of Y (the outcome) given X (predictors/genotype) under the test model.
#' @param sd_y_x_truth The standard deviation of Y given X (predictors/genotype) given genotype under the true model.
#' @param model The genetic model in the linear regression: "Dominant", "Additive", "Recessive", "2df" or "null"
#'
#' @return The log likelihood.
#'
#' @export
#'
calc.like.linear<-function(beta, m, es_ab, es_bb, sd_y_x_model, sd_y_x_truth, model){
  es <- c(es_ab, es_bb)
  if(model=='Dominant'){ll <- dominant.ll.linear(beta, m, es, sd_y_x_model, sd_y_x_truth)}
  if(model=='Additive'){ll <- additive.ll.linear(beta, m, es, sd_y_x_model, sd_y_x_truth)}
  if(model=='Recessive'){ll <- recessive.ll.linear(beta, m, es, sd_y_x_model, sd_y_x_truth)}
  if(model=='2df'){ll <- df2.ll.linear(beta, m, es, sd_y_x_model, sd_y_x_truth)}
  if(model=='null'){ll <- null.ll.linear(beta, m, es, sd_y_x_model, sd_y_x_truth)}
  return(ll)
}

#' Function to return log likelihood function for specified model type
#'
#' Convenience function to return log likelihood function for specified model type
#'
#' @param model The genetic model in the linear regression: "Dominant", "Additive", "Recessive", "2df" or "null"
#'
#' @return Log likelihood function for specified model type
#'
#' @export
#'
ll.linear.selector <- function(model){
  if(model=='Dominant'){res <- dominant.ll.linear}
  if(model=='Additive'){res <- additive.ll.linear}
  if(model=='Recessive'){res <- recessive.ll.linear}
  if(model=='2df'){res <- df2.ll.linear}
  if(model=='null'){res <- null.ll.linear}
  return(res)
}
