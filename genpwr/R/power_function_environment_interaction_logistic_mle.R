#' Function to calculate MLE's for logistic models with logistic environment interaction
#'
#' Finds the maximum likelihood estimates for a given 2x3 table under the specified genetic model. 
#'
#' @param t A 2x6 table of the joint probabilities of disease, genotype, and environment. Rows are case vs. control and columns are genotypes. 
#' @param N Sample size
#' @param power Power
#' @param Alpha Alpha
#' @param mod Test model
#' 
#' @return A vector of logistic regression model coefficients.  
#'
#' @export
#'
ll.ge.logistic <- function(t, N = NULL, power = NULL, Alpha, mod){
	if(all(c(is.null(N), is.null(power)))) stop("must specify either N or power")
	if(!any(c(is.null(N), is.null(power)))) stop("must specify either N or power, not both")

	prob_AA_control_e0 <- t[2, 1]
	prob_AB_control_e0 <- t[2, 2]
	prob_BB_control_e0 <- t[2, 3]
	prob_AA_case_e0 <- t[1, 1]
	prob_AB_case_e0 <- t[1, 2]
	prob_BB_case_e0 <- t[1, 3]
	prob_AA_control_e <- t[2, 4]
	prob_AB_control_e <- t[2, 5]
	prob_BB_control_e <- t[2, 6]
	prob_AA_case_e <- t[1, 4]
	prob_AB_case_e <- t[1, 5]
	prob_BB_case_e <- t[1, 6]

	if(mod == "Dominant"){
		ll_fun <- function(beta){
			prob_AA_case_e0*log(exp(beta[1])/(1 + exp(beta[1])))+
			prob_AB_case_e0*log(exp(beta[1] + beta[2])/(1 + exp(beta[1] + beta[2])))+
			prob_BB_case_e0*log(exp(beta[1] + beta[2])/(1 + exp(beta[1] + beta[2])))+
			prob_AA_case_e*log(exp(beta[1] + beta[3])/(1 + exp(beta[1] + beta[3])))+
			prob_AB_case_e*log(exp(beta[1] + beta[2] + beta[3] + beta[4])/(1 + exp(beta[1] + beta[2] + beta[3] + beta[4])))+
			prob_BB_case_e*log(exp(beta[1] + beta[2] + beta[3] + beta[4])/(1 + exp(beta[1] + beta[2] + beta[3] + beta[4])))+
			prob_AA_control_e0*log(1/(1 + exp(beta[1])))+
			prob_AB_control_e0*log(1/(1 + exp(beta[1] + beta[2])))+
			prob_BB_control_e0*log(1/(1 + exp(beta[1] + beta[2])))+
			prob_AA_control_e*log(1/(1 + exp(beta[1] + beta[3])))+
			prob_AB_control_e*log(1/(1 + exp(beta[1] + beta[2] + beta[3] + beta[4])))+
			prob_BB_control_e*log(1/(1 + exp(beta[1] + beta[2] + beta[3] + beta[4])))
		}

		ll_g_e_fun <- function(beta){
			prob_AA_case_e0*log(exp(beta[1])/(1 + exp(beta[1])))+
			prob_AB_case_e0*log(exp(beta[1] + beta[2])/(1 + exp(beta[1] + beta[2])))+
			prob_BB_case_e0*log(exp(beta[1] + beta[2])/(1 + exp(beta[1] + beta[2])))+
			prob_AA_case_e*log(exp(beta[1] + beta[3])/(1 + exp(beta[1] + beta[3])))+
			prob_AB_case_e*log(exp(beta[1] + beta[2] + beta[3])/(1 + exp(beta[1] + beta[2] + beta[3])))+
			prob_BB_case_e*log(exp(beta[1] + beta[2] + beta[3])/(1 + exp(beta[1] + beta[2] + beta[3])))+
			prob_AA_control_e0*log(1/(1 + exp(beta[1])))+
			prob_AB_control_e0*log(1/(1 + exp(beta[1] + beta[2])))+
			prob_BB_control_e0*log(1/(1 + exp(beta[1] + beta[2])))+
			prob_AA_control_e*log(1/(1 + exp(beta[1] + beta[3])))+
			prob_AB_control_e*log(1/(1 + exp(beta[1] + beta[2] + beta[3])))+
			prob_BB_control_e*log(1/(1 + exp(beta[1] + beta[2] + beta[3])))
		}

		ll_e_fun <- function(beta){
			prob_AA_case_e0*log(exp(beta[1])/(1 + exp(beta[1])))+
			prob_AB_case_e0*log(exp(beta[1])/(1 + exp(beta[1])))+
			prob_BB_case_e0*log(exp(beta[1])/(1 + exp(beta[1])))+
			prob_AA_case_e*log(exp(beta[1] + beta[2])/(1 + exp(beta[1] + beta[2])))+
			prob_AB_case_e*log(exp(beta[1] + beta[2])/(1 + exp(beta[1] + beta[2])))+
			prob_BB_case_e*log(exp(beta[1] + beta[2])/(1 + exp(beta[1] + beta[2])))+
			prob_AA_control_e0*log(1/(1 + exp(beta[1])))+
			prob_AB_control_e0*log(1/(1 + exp(beta[1])))+
			prob_BB_control_e0*log(1/(1 + exp(beta[1])))+
			prob_AA_control_e*log(1/(1 + exp(beta[1] + beta[2])))+
			prob_AB_control_e*log(1/(1 + exp(beta[1] + beta[2])))+
			prob_BB_control_e*log(1/(1 + exp(beta[1] + beta[2])))
		}

		ll_g_fun <- function(beta){
			prob_AA_case_e0*log(exp(beta[1])/(1 + exp(beta[1])))+
			prob_AB_case_e0*log(exp(beta[1] + beta[2])/(1 + exp(beta[1] + beta[2])))+
			prob_BB_case_e0*log(exp(beta[1] + beta[2])/(1 + exp(beta[1] + beta[2])))+
			prob_AA_case_e*log(exp(beta[1])/(1 + exp(beta[1])))+
			prob_AB_case_e*log(exp(beta[1] + beta[2])/(1 + exp(beta[1] + beta[2])))+
			prob_BB_case_e*log(exp(beta[1] + beta[2])/(1 + exp(beta[1] + beta[2])))+
			prob_AA_control_e0*log(1/(1 + exp(beta[1])))+
			prob_AB_control_e0*log(1/(1 + exp(beta[1] + beta[2])))+
			prob_BB_control_e0*log(1/(1 + exp(beta[1] + beta[2])))+
			prob_AA_control_e*log(1/(1 + exp(beta[1])))+
			prob_AB_control_e*log(1/(1 + exp(beta[1] + beta[2])))+
			prob_BB_control_e*log(1/(1 + exp(beta[1] + beta[2])))
		}
	}

	if(mod == "Recessive"){
		ll_fun <- function(beta){
			prob_AA_case_e0*log(exp(beta[1])/(1 + exp(beta[1])))+
			prob_AB_case_e0*log(exp(beta[1])/(1 + exp(beta[1])))+
			prob_BB_case_e0*log(exp(beta[1] + beta[2])/(1 + exp(beta[1] + beta[2])))+
			prob_AA_case_e*log(exp(beta[1] + beta[3])/(1 + exp(beta[1] + beta[3])))+
			prob_AB_case_e*log(exp(beta[1] + beta[3])/(1 + exp(beta[1] + beta[3])))+
			prob_BB_case_e*log(exp(beta[1] + beta[2] + beta[3] + beta[4])/(1 + exp(beta[1] + beta[2] + beta[3] + beta[4])))+
			prob_AA_control_e0*log(1/(1 + exp(beta[1])))+
			prob_AB_control_e0*log(1/(1 + exp(beta[1])))+
			prob_BB_control_e0*log(1/(1 + exp(beta[1] + beta[2])))+
			prob_AA_control_e*log(1/(1 + exp(beta[1] + beta[3])))+
			prob_AB_control_e*log(1/(1 + exp(beta[1] + beta[3])))+
			prob_BB_control_e*log(1/(1 + exp(beta[1] + beta[2] + beta[3] + beta[4])))
		}

		ll_g_e_fun <- function(beta){
			prob_AA_case_e0*log(exp(beta[1])/(1 + exp(beta[1])))+
			prob_AB_case_e0*log(exp(beta[1])/(1 + exp(beta[1])))+
			prob_BB_case_e0*log(exp(beta[1] + beta[2])/(1 + exp(beta[1] + beta[2])))+
			prob_AA_case_e*log(exp(beta[1] + beta[3])/(1 + exp(beta[1] + beta[3])))+
			prob_AB_case_e*log(exp(beta[1] + beta[3])/(1 + exp(beta[1] + beta[3])))+
			prob_BB_case_e*log(exp(beta[1] + beta[2] + beta[3])/(1 + exp(beta[1] + beta[2] + beta[3])))+
			prob_AA_control_e0*log(1/(1 + exp(beta[1])))+
			prob_AB_control_e0*log(1/(1 + exp(beta[1])))+
			prob_BB_control_e0*log(1/(1 + exp(beta[1] + beta[2])))+
			prob_AA_control_e*log(1/(1 + exp(beta[1] + beta[3])))+
			prob_AB_control_e*log(1/(1 + exp(beta[1] + beta[3])))+
			prob_BB_control_e*log(1/(1 + exp(beta[1] + beta[2] + beta[3])))
		}

		ll_e_fun <- function(beta){
			prob_AA_case_e0*log(exp(beta[1])/(1 + exp(beta[1])))+
			prob_AB_case_e0*log(exp(beta[1])/(1 + exp(beta[1])))+
			prob_BB_case_e0*log(exp(beta[1])/(1 + exp(beta[1])))+
			prob_AA_case_e*log(exp(beta[1] + beta[2])/(1 + exp(beta[1] + beta[2])))+
			prob_AB_case_e*log(exp(beta[1] + beta[2])/(1 + exp(beta[1] + beta[2])))+
			prob_BB_case_e*log(exp(beta[1] + beta[2])/(1 + exp(beta[1] + beta[2])))+
			prob_AA_control_e0*log(1/(1 + exp(beta[1])))+
			prob_AB_control_e0*log(1/(1 + exp(beta[1])))+
			prob_BB_control_e0*log(1/(1 + exp(beta[1])))+
			prob_AA_control_e*log(1/(1 + exp(beta[1] + beta[2])))+
			prob_AB_control_e*log(1/(1 + exp(beta[1] + beta[2])))+
			prob_BB_control_e*log(1/(1 + exp(beta[1] + beta[2])))
		}

		ll_g_fun <- function(beta){
			prob_AA_case_e0*log(exp(beta[1])/(1 + exp(beta[1])))+
			prob_AB_case_e0*log(exp(beta[1])/(1 + exp(beta[1])))+
			prob_BB_case_e0*log(exp(beta[1] + beta[2])/(1 + exp(beta[1] + beta[2])))+
			prob_AA_case_e*log(exp(beta[1])/(1 + exp(beta[1])))+
			prob_AB_case_e*log(exp(beta[1])/(1 + exp(beta[1])))+
			prob_BB_case_e*log(exp(beta[1] + beta[2])/(1 + exp(beta[1] + beta[2])))+
			prob_AA_control_e0*log(1/(1 + exp(beta[1])))+
			prob_AB_control_e0*log(1/(1 + exp(beta[1])))+
			prob_BB_control_e0*log(1/(1 + exp(beta[1] + beta[2])))+
			prob_AA_control_e*log(1/(1 + exp(beta[1])))+
			prob_AB_control_e*log(1/(1 + exp(beta[1])))+
			prob_BB_control_e*log(1/(1 + exp(beta[1] + beta[2])))
		}
	}

	if(mod == "Additive"){
		ll_fun <- function(beta){
			prob_AA_case_e0*log(exp(beta[1])/(1 + exp(beta[1])))+
			prob_AB_case_e0*log(exp(beta[1] + beta[2])/(1 + exp(beta[1] + beta[2])))+
			prob_BB_case_e0*log(exp(beta[1] + 2*beta[2])/(1 + exp(beta[1] + 2*beta[2])))+
			prob_AA_case_e*log(exp(beta[1] + beta[3])/(1 + exp(beta[1] + beta[3])))+
			prob_AB_case_e*log(exp(beta[1] + beta[2] + beta[3] + beta[4])/(1 + exp(beta[1] + beta[2] + beta[3] + beta[4])))+
			prob_BB_case_e*log(exp(beta[1] + 2*beta[2] + beta[3] + 2*beta[4])/(1 + exp(beta[1] + 2*beta[2] + beta[3] + 2*beta[4])))+
			prob_AA_control_e0*log(1/(1 + exp(beta[1])))+
			prob_AB_control_e0*log(1/(1 + exp(beta[1] + beta[2])))+
			prob_BB_control_e0*log(1/(1 + exp(beta[1] + 2*beta[2])))+
			prob_AA_control_e*log(1/(1 + exp(beta[1] + beta[3])))+
			prob_AB_control_e*log(1/(1 + exp(beta[1] + beta[2] + beta[3] + beta[4])))+
			prob_BB_control_e*log(1/(1 + exp(beta[1] + 2*beta[2] + beta[3] + 2*beta[4])))
		}

		ll_g_e_fun <- function(beta){
			prob_AA_case_e0*log(exp(beta[1])/(1 + exp(beta[1])))+
			prob_AB_case_e0*log(exp(beta[1] + beta[2])/(1 + exp(beta[1] + beta[2])))+
			prob_BB_case_e0*log(exp(beta[1] + 2*beta[2])/(1 + exp(beta[1] + 2*beta[2])))+
			prob_AA_case_e*log(exp(beta[1] + beta[3])/(1 + exp(beta[1] + beta[3])))+
			prob_AB_case_e*log(exp(beta[1] + beta[2] + beta[3])/(1 + exp(beta[1] + beta[2] + beta[3])))+
			prob_BB_case_e*log(exp(beta[1] + 2*beta[2] + beta[3])/(1 + exp(beta[1] + 2*beta[2] + beta[3])))+
			prob_AA_control_e0*log(1/(1 + exp(beta[1])))+
			prob_AB_control_e0*log(1/(1 + exp(beta[1] + beta[2])))+
			prob_BB_control_e0*log(1/(1 + exp(beta[1] + 2*beta[2])))+
			prob_AA_control_e*log(1/(1 + exp(beta[1] + beta[3])))+
			prob_AB_control_e*log(1/(1 + exp(beta[1] + beta[2] + beta[3])))+
			prob_BB_control_e*log(1/(1 + exp(beta[1] + 2*beta[2] + beta[3])))
		}

		ll_e_fun <- function(beta){
			prob_AA_case_e0*log(exp(beta[1])/(1 + exp(beta[1])))+
			prob_AB_case_e0*log(exp(beta[1])/(1 + exp(beta[1])))+
			prob_BB_case_e0*log(exp(beta[1])/(1 + exp(beta[1])))+
			prob_AA_case_e*log(exp(beta[1] + beta[2])/(1 + exp(beta[1] + beta[2])))+
			prob_AB_case_e*log(exp(beta[1] + beta[2])/(1 + exp(beta[1] + beta[2])))+
			prob_BB_case_e*log(exp(beta[1] + beta[2])/(1 + exp(beta[1] + beta[2])))+
			prob_AA_control_e0*log(1/(1 + exp(beta[1])))+
			prob_AB_control_e0*log(1/(1 + exp(beta[1])))+
			prob_BB_control_e0*log(1/(1 + exp(beta[1])))+
			prob_AA_control_e*log(1/(1 + exp(beta[1] + beta[2])))+
			prob_AB_control_e*log(1/(1 + exp(beta[1] + beta[2])))+
			prob_BB_control_e*log(1/(1 + exp(beta[1] + beta[2])))
		}

		ll_g_fun <- function(beta){
			prob_AA_case_e0*log(exp(beta[1])/(1 + exp(beta[1])))+
			prob_AB_case_e0*log(exp(beta[1] + beta[2])/(1 + exp(beta[1] + beta[2])))+
			prob_BB_case_e0*log(exp(beta[1] + 2*beta[2])/(1 + exp(beta[1] + 2*beta[2])))+
			prob_AA_case_e*log(exp(beta[1])/(1 + exp(beta[1])))+
			prob_AB_case_e*log(exp(beta[1] + beta[2])/(1 + exp(beta[1] + beta[2])))+
			prob_BB_case_e*log(exp(beta[1] + 2*beta[2])/(1 + exp(beta[1] + 2*beta[2])))+
			prob_AA_control_e0*log(1/(1 + exp(beta[1])))+
			prob_AB_control_e0*log(1/(1 + exp(beta[1] + beta[2])))+
			prob_BB_control_e0*log(1/(1 + exp(beta[1] + 2*beta[2])))+
			prob_AA_control_e*log(1/(1 + exp(beta[1])))+
			prob_AB_control_e*log(1/(1 + exp(beta[1] + beta[2])))+
			prob_BB_control_e*log(1/(1 + exp(beta[1] + 2*beta[2])))
		}
	}

	if(mod == "2df"){
		# 1 = beta0
		# 2 = AB
		# 3 = BB
		# 4 = E
		# 5 = AB/E
		# 6 = BB/E
		ll_fun <- function(beta){
			prob_AA_case_e0*log(exp(beta[1])/(1 + exp(beta[1])))+
			prob_AB_case_e0*log(exp(beta[1] + beta[2])/(1 + exp(beta[1] + beta[2])))+
			prob_BB_case_e0*log(exp(beta[1] + beta[3])/(1 + exp(beta[1] + beta[3])))+
			prob_AA_case_e*log(exp(beta[1] + beta[4])/(1 + exp(beta[1] + beta[4])))+
			prob_AB_case_e*log(exp(beta[1] + beta[2] + beta[4] + beta[5])/(1 + exp(beta[1] + beta[2] + beta[4] + beta[5])))+
			prob_BB_case_e*log(exp(beta[1] + beta[3] + beta[4] + beta[6])/(1 + exp(beta[1] + beta[3] + beta[4] + beta[6])))+
			prob_AA_control_e0*log(1/(1 + exp(beta[1])))+
			prob_AB_control_e0*log(1/(1 + exp(beta[1] + beta[2])))+
			prob_BB_control_e0*log(1/(1 + exp(beta[1] + beta[3])))+
			prob_AA_control_e*log(1/(1 + exp(beta[1] + beta[4])))+
			prob_AB_control_e*log(1/(1 + exp(beta[1] + beta[2] + beta[4] + beta[5])))+
			prob_BB_control_e*log(1/(1 + exp(beta[1] + beta[3] + beta[4] + beta[6])))
		}

		ll_g_e_fun <- function(beta){
			prob_AA_case_e0*log(exp(beta[1])/(1 + exp(beta[1])))+
			prob_AB_case_e0*log(exp(beta[1] + beta[2])/(1 + exp(beta[1] + beta[2])))+
			prob_BB_case_e0*log(exp(beta[1] + beta[3])/(1 + exp(beta[1] + beta[3])))+
			prob_AA_case_e*log(exp(beta[1] + beta[4])/(1 + exp(beta[1] + beta[4])))+
			prob_AB_case_e*log(exp(beta[1] + beta[2] + beta[4])/(1 + exp(beta[1] + beta[2] + beta[4])))+
			prob_BB_case_e*log(exp(beta[1] + beta[3] + beta[4])/(1 + exp(beta[1] + beta[3] + beta[4])))+
			prob_AA_control_e0*log(1/(1 + exp(beta[1])))+
			prob_AB_control_e0*log(1/(1 + exp(beta[1] + beta[2])))+
			prob_BB_control_e0*log(1/(1 + exp(beta[1] + beta[3])))+
			prob_AA_control_e*log(1/(1 + exp(beta[1] + beta[4])))+
			prob_AB_control_e*log(1/(1 + exp(beta[1] + beta[2] + beta[4])))+
			prob_BB_control_e*log(1/(1 + exp(beta[1] + beta[3] + beta[4])))
		}

		ll_e_fun <- function(beta){
			prob_AA_case_e0*log(exp(beta[1])/(1 + exp(beta[1])))+
			prob_AB_case_e0*log(exp(beta[1])/(1 + exp(beta[1])))+
			prob_BB_case_e0*log(exp(beta[1])/(1 + exp(beta[1])))+
			prob_AA_case_e*log(exp(beta[1] + beta[2])/(1 + exp(beta[1] + beta[2])))+
			prob_AB_case_e*log(exp(beta[1] + beta[2])/(1 + exp(beta[1] + beta[2])))+
			prob_BB_case_e*log(exp(beta[1] + beta[2])/(1 + exp(beta[1] + beta[2])))+
			prob_AA_control_e0*log(1/(1 + exp(beta[1])))+
			prob_AB_control_e0*log(1/(1 + exp(beta[1])))+
			prob_BB_control_e0*log(1/(1 + exp(beta[1])))+
			prob_AA_control_e*log(1/(1 + exp(beta[1] + beta[2])))+
			prob_AB_control_e*log(1/(1 + exp(beta[1] + beta[2])))+
			prob_BB_control_e*log(1/(1 + exp(beta[1] + beta[2])))
		}

		ll_g_fun <- function(beta){
			prob_AA_case_e0*log(exp(beta[1])/(1 + exp(beta[1])))+
			prob_AB_case_e0*log(exp(beta[1] + beta[2])/(1 + exp(beta[1] + beta[2])))+
			prob_BB_case_e0*log(exp(beta[1] + beta[3])/(1 + exp(beta[1] + beta[3])))+
			prob_AA_case_e*log(exp(beta[1])/(1 + exp(beta[1])))+
			prob_AB_case_e*log(exp(beta[1] + beta[2])/(1 + exp(beta[1] + beta[2])))+
			prob_BB_case_e*log(exp(beta[1] + beta[3])/(1 + exp(beta[1] + beta[3])))+
			prob_AA_control_e0*log(1/(1 + exp(beta[1])))+
			prob_AB_control_e0*log(1/(1 + exp(beta[1] + beta[2])))+
			prob_BB_control_e0*log(1/(1 + exp(beta[1] + beta[3])))+
			prob_AA_control_e*log(1/(1 + exp(beta[1])))+
			prob_AB_control_e*log(1/(1 + exp(beta[1] + beta[2])))+
			prob_BB_control_e*log(1/(1 + exp(beta[1] + beta[3])))
		}
	}

	if(mod == "2df"){
		beta <- optim(c(logit(prob_AA_case_e0/(prob_AA_control_e0 + prob_AA_case_e0)),0,0,
			logit(prob_AA_case_e/((prob_AA_control_e)+(prob_AA_case_e))),0,0), function(x) -ll_fun(x), control=c(abstol = 1e-5))$par
		# ll0 <- ll_fun(beta)
		# diff <- 1
		# while(diff > 1e-6){
		# 	beta <- optim(beta, function(x) -ll_fun(x), control=c(abstol = 1e-5))$par
		# 	ll0 <- c(ll0, ll_fun(beta))
		# 	diff <- ll0[length(ll0)] - ll0[length(ll0) - 1]
		# }
		# # beta <- optim(beta, function(x) -ll_fun(x), control = c(abstol = 1e-5))$par
		# ll <- -optim(c(0,0,0,0,0,0), function(x) -ll_fun(x), control=c(abstol = 1e-5))$value
		# beta_g_e <- optim(c(0,0,0,0), function(x) -ll_g_e_fun(x), control=c(abstol = 1e-5))$par
		ll_g_e <- -optim(c(0,0,0,0), function(x) -ll_g_e_fun(x), control=c(abstol = 1e-5))$value
		# beta_e <- optim(c(0,0), function(x) -ll_e_fun(x), control=c(abstol = 1e-5))$par
		ll_e <- -optim(c(0,0), function(x) -ll_e_fun(x), control=c(abstol = 1e-5))$value
		# beta_g <- optim(c(0,0,0), function(x) -ll_g_fun(x), control=c(abstol = 1e-5))$par
		ll_g <- -optim(c(0,0,0), function(x) -ll_g_fun(x), control=c(abstol = 1e-5))$value
	}else{
		if(mod == "Dominant"){
			beta0 <- logit(prob_AA_case_e0/(prob_AA_control_e0 + prob_AA_case_e0))
			betaG <- logit((prob_AB_case_e0+prob_BB_case_e0)/(prob_AB_control_e0+prob_BB_control_e0+prob_AB_case_e0+prob_BB_case_e0))
			betaE <- logit(prob_AA_case_e/((prob_AA_control_e)+(prob_AA_case_e)))
			betaGE <- logit((prob_AB_case_e+prob_BB_case_e)/((prob_AB_control_e+prob_BB_control_e)+(prob_AB_case_e+prob_BB_case_e)))
			beta <- c(beta0, betaG - beta0, betaE - beta0, betaGE + beta0 - betaG - betaE)
		}else if(mod == "Recessive"){
			beta0 <- logit((prob_AA_case_e0+prob_AB_case_e0)/(prob_AA_control_e0+prob_AB_control_e0+prob_AA_case_e0+prob_AB_case_e0))
			betaG <- logit(prob_BB_case_e0/(prob_BB_control_e0+prob_BB_case_e0))
			betaE <- logit((prob_AA_case_e+prob_AB_case_e)/(prob_AA_control_e+prob_AB_control_e+prob_AA_case_e+prob_AB_case_e))
			betaGE <- logit((prob_BB_case_e)/(prob_BB_control_e+prob_BB_case_e))
			beta <- c(beta0, betaG - beta0, betaE - beta0, betaGE + beta0 - betaG - betaE)
		}else if(mod == "Additive"){
			beta <- optim(c(logit(prob_AA_case_e0/(prob_AA_control_e0 + prob_AA_case_e0)),0,logit(prob_AA_case_e/((prob_AA_control_e)+(prob_AA_case_e))),0), 
				function(x) -ll_fun(x), control=c(abstol = 1e-5))$par
		}
		# betacalc <- beta
		# beta <- optim(beta, function(x) -ll_fun(x), control=c(abstol = 1e-5))$par
		# if(!all(betacalc == beta)) stop("beta different from betacalc")
		# ll <- -optim(betacalc, function(x) -ll_fun(x), control=c(abstol = 1e-5))$value
		# beta_g_e <- optim(c(0,0,0), function(x) -ll_g_e_fun(x), control=c(abstol = 1e-5))$par
		ll_g_e <- -optim(c(0,0,0), function(x) -ll_g_e_fun(x), control=c(abstol = 1e-5))$value
		# beta_e <- optim(c(0,0), function(x) -ll_e_fun(x), control=c(abstol = 1e-5))$par
		ll_e <- -optim(c(0,0), function(x) -ll_e_fun(x), control=c(abstol = 1e-5))$value
		# beta_g <- optim(c(0,0), function(x) -ll_g_fun(x), control=c(abstol = 1e-5))$par
		ll_g <- -optim(c(0,0), function(x) -ll_g_fun(x), control=c(abstol = 1e-5))$value
	}

	# optimization loop to make sure beta is optimized correctly
	diff <- 1
	ll0 <- ll_fun(beta)
	while(diff > 1e-6){
		beta <- optim(beta, function(x) -ll_fun(x), control=c(abstol = 1e-5))$par
		ll0 <- c(ll0, ll_fun(beta))
		diff <- ll0[length(ll0)] - ll0[length(ll0) - 1]
	}
	ll <- ll_fun(beta)
	# ll_g_e <- ll_g_e_fun(beta_g_e)
	# ll_e <- ll_e_fun(beta_e)
	# ll_g <- ll_g_fun(beta_g)

	Case.Rate <- sum(t[1,])	

	ll.null <- Case.Rate*log(exp(logit(sum(t[1,])/sum(t)))/(1+exp(logit(sum(t[1,])/sum(t))))) + (1-Case.Rate)*log(1/(1+exp(logit(sum(t[1,])/sum(t)))))
	# logit(sum(t[1,])/sum(t))
	if(is.null(power)){
		if(mod == "2df"){
			power_res <- 1-pchisq(qchisq(1-Alpha, df=2, ncp=0), df=2, ncp = N*stat) 
		}else{
			power_res <- pnorm(sqrt(N*2*(ll-ll_g_e)) - qnorm(1-Alpha/2)) + pnorm(-sqrt(N*2*(ll-ll_g_e)) - qnorm(1-Alpha/2))*1
		}
		return(power_res)
	}
	if(is.null(N)){
		stat <- 2*(as.numeric(ll-ll_g_e))
		if(mod=='2df'){
			ss <- uniroot(function(x) ncp.search(x, power, stat, Alpha[q], df=2),
							lower=0, upper=1000, extendInt = 'upX', tol=0.00001)$root/stat
		}else{
			ss <- (qnorm(1-Alpha/2)+qnorm(power))^2/(2*(ll-ll_g_e))
		}
		return(ss)
	}
}
