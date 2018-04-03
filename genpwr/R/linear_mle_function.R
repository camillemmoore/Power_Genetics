# MLE Functions for Normal/Continuous Outcomes
#' Function to calculate MLE's for linear models
#'
#' Finds the maximum likelihood estiamtes for a given MAF under the specified genetic model and effect size.
#'
#' @param m minor allele frequency
#' @param es vector of effect sizes with two elements, (mean AB - mean AA) and (mean BB - mean AA)
#' @param model The assumed genetic model(s) used in testing: 'Dominant', 'Additive', 'Recessive', '2df'
#'
#' @return A vector of linear regression model coefficients.
#'
#' @export
#'
linear.mles<-
  function(m,es,model){

    if (model=='null'){#Null Model
      beta0 = 2*m*(1-m)*es[1]+m*m*es[2]
      beta = beta0
    }

    #Dominant
    if (model=='Dominant'){
      beta0 = 0
      beta1 = (2*m*(1-m)*es[1]+m*m*es[2])/(1-(1-m)^2)
      beta = c(beta0, beta1)
    }

    #Recessive
    if (model=='Recessive'){
      beta0 = 2*m*(1-m)*es[1]/(1-m^2)
      beta1 = es[2]
      beta = c(beta0, beta1)
    }

    #2DF
    if (model=='2df'){
      beta0 = 0
      beta1 = es[1]
      beta2 = es[2]
      beta = c(beta0, beta1, beta2)
    }

    #Additive Model
    if (model=='Additive'){
      expected_x = 2*m*(1-m)*1+ m*m*2
      expected_x2 = 2*m*(1-m)*1^2+ m*m*2^2
      expected_y = 2*m*(1-m)*es[1]+ m*m*es[2]
      expected_xy = 2*m*(1-m)*1*es[1]+ m*m*2*es[2]
      cov_xy = expected_xy-expected_x*expected_y
      var_x = expected_x2 - expected_x^2
      beta1 = cov_xy/var_x
      beta0 = expected_y - beta1*expected_x
      beta = c(beta0, beta1)
    }
    return(beta)
  }
