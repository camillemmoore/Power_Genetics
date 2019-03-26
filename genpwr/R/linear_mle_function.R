# MLE Functions for Normal/Continuous Outcomes
#' Function to calculate MLE's for linear models
#'
#' Finds the maximum likelihood estimates for a given MAF under the specified genetic model and effect size.
#'
#' @param m minor allele frequency
#' @param es_ab effect size for mean AB - mean AA 
#' @param es_bb effect size for mean BB - mean AA 
#' @param model The assumed genetic model(s) used in testing: 'Dominant', 'Additive', 'Recessive', '2df'
#'
#' @return A vector of linear regression model coefficients.
#'
#' @examples
#' linear.mles(m = 0.1, es_ab = 0, es_bb = 3, model = "Dominant")
#'
#' @export
#'
linear.mles<-
	function(m,es_ab, es_bb,model)
{

	es = c(es_ab, es_bb)

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
		beta1 = es[2]-beta0
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

#' Functions to Calculate Residual SD for Normal/Continuous Outcomes
#' Function to calculate the standard deviation of y given x for linear models
#'
#' @param m minor allele frequency
#' @param es_ab effect size for mean AB - mean AA 
#' @param es_bb effect size for mean BB - mean AA
#' @param sd_y the standard deviation of y in the overall population.
#' @param model The assumed genetic model(s) used in testing: 'Dominant', 'Additive', 'Recessive', '2df'
#'
#' @return A vector of linear regression model coefficients.
#'
#' @examples 
#' linear.sds(m = 0.1, es_ab = 0, es_bb = 3, sd_y = 1, model = "Dominant")
#'
#' @export
#'
linear.sds<-function(m, es_ab, es_bb, sd_y, model)
{

	es = c(es_ab, es_bb)

		if (model=='null'){#Null Model
			sd_y_x = sd_y
		}

		#Dominant
		if (model=='Dominant'){
			beta1 = (2*m*(1-m)*es[1]+m*m*es[2])/(1-(1-m)^2)
			expected_x = 2*m*(1-m)*1+ m*m*1
			expected_x2 = 2*m*(1-m)*1^2+ m*m*1^2
			var_x = expected_x2 - expected_x^2
			sd_y_x = sqrt(sd_y^2 - (beta1^2)*var_x)
		}

		#Recessive
		if (model=='Recessive'){
			beta0 = 2*m*(1-m)*es[1]/(1-m^2)
			beta1 = es[2]-beta0
			expected_x = m*m*1
			expected_x2 = m*m*1^2
			var_x = expected_x2 - expected_x^2
			sd_y_x = sqrt(sd_y^2 - (beta1^2)*var_x)
		}

		#2DF
		if (model=='2df'){
			expected_y = 2*m*(1-m)*es[1]+ m*m*es[2]
			temp = (expected_y^2)*((1-m)^2) +
						 ((es[1]-expected_y)^2)*2*(1-m)*m+
						 ((es[2]-expected_y)^2)*m*m
			sd_y_x = sqrt(sd_y^2 - temp)
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
			sd_y_x = sqrt(sd_y^2 - (beta1^2)*var_x)
		}
		return(sd_y_x)
}
