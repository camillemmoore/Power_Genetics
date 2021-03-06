#' Function to generate integrand for mle for cases
#'
#' Returns the standard deviation of y given x for linear models with linear environment interaction
#'
#' @param x1 "true" part of model
#' @param x2 "test" part of model
#'
#' @return a function to be used as the integrand for the mle
#'
#' @examples
#' integrand_funct_case(-1.462531 + 1*0.1823216,
#' 	-1.462531 + 1*0.1823216) 
#'
#' @export
#'
integrand_funct_case <- function(x1, x2){
	return(exp(x1) / (1 + exp(x1)) * log(exp(x2) / (1 + exp(x2))))
}

#' Function to generate integrand for mle for controls
#'
#' Returns the standard deviation of y given x for linear models with linear environment interaction
#'
#' @param x1 "true" part of model
#' @param x2 "test" part of model
#'
#' @return a function to be used as the integrand for the mle
#'
#' @examples
#' integrand_funct_control(-1.462531 + 1*0.1823216,
#' 	-1.462531 + 1*0.1823216) 
#'
#' @export
#'
integrand_funct_control <- function(x1, x2){
	# xposneg <- function(anx) if(anx >= 1) return(exp(-anx)/(1 + exp(-anx))) else return(1/(1 + exp(anx)))
	return(1 / (1 + exp(x1)) * log(1 / (1 + exp(x2))))
}

#' Function to output log likelihood for logistic outcome with linear environment variables
#'
#' Returns the standard deviation of y given x for linear models with linear environment interaction
#'
#' @param sd_e Standard deviation of the environmental variable
#' @param N desired sample size
#' @param MAF Vector of minor allele frequencies
#' @param power desired power
#' @param beta0 the beta0 coefficient in the logistic model
#' @param OR_G Vector of genetic odds ratios to detect
#' @param OR_E Vector of environmental odds ratios to detect
#' @param OR_GE Vector of genetic/environmental interaction odds ratios to detect
#' @param Alpha the desired type 1 error rate(s)
#' @param True.Model A vector specifying the true underlying genetic model(s): 'Dominant', 'Additive', 'Recessive' or 'All'
#' @param Test.Model A vector specifying the assumed genetic model(s) used in testing: 'Dominant', 'Additive', 'Recessive' or 'All'
#'
#' @return a function to be used as the integrand for the mle
#'
#' @examples
#' ll.ge.logistic.lin.envir(sd_e = 1, MAF = 0.2, N = 30, beta0 = -1.462531, OR_G = 1.1, 
#' 	OR_E = 1.2, OR_GE = 1.5, Alpha = 0.05, True.Model = "Dominant", Test.Model = "Dominant")
#' @export
#'
ll.ge.logistic.lin.envir <- function(sd_e, N = NULL, MAF, power = NULL, beta0, OR_G, OR_E, OR_GE, Alpha, True.Model, Test.Model)
{


	if(all(c(is.null(N), is.null(power)))) stop("must specify either N or power")
	if(!any(c(is.null(N), is.null(power)))) stop("must specify either N or power, not both")

	P_AA <- (1 - MAF)^2
	P_AB <- 2*(1-MAF)*MAF
	P_BB <- MAF^2

	beta_true <- c(beta0, log(OR_G), log(OR_E), log(OR_GE))

	mod_num <- function(amod){
		if(amod == "Dominant") mod_num <- c(1,1)
		if(amod == "Recessive") mod_num <- c(0,1)
		if(amod == "Additive") mod_num <- c(1,2)
		return(mod_num)
	}

	if(Test.Model != "2df"){
		mntr <- mod_num(True.Model)
		mnte <- mod_num(Test.Model)
		ll_fun <- function(beta){
			ll_fun00 <- function(xe){
				1 / sqrt(2 * pi * sd_e^2) * exp(-xe^2/(2 * sd_e^2)) * (
					P_AA * integrand_funct_case(beta_true[1] + xe*beta_true[3],
						beta[1] + xe*beta[3])+
					P_AB * integrand_funct_case(beta_true[1] + mntr[1]*beta_true[2] + xe*(beta_true[3]+mntr[1]*beta_true[4]),
						beta[1] + mnte[1]*beta[2] + xe*(beta[3]+mnte[1]*beta[4]))+
					P_BB * integrand_funct_case(beta_true[1] + mntr[2]*beta_true[2] + xe*(beta_true[3]+mntr[2]*beta_true[4]),
						beta[1] + mnte[2]*beta[2] + xe*(beta[3]+mnte[2]*beta[4]))+
					P_AA * integrand_funct_control(beta_true[1] + xe*beta_true[3],
						beta[1] + xe*beta[3])+
					P_AB * integrand_funct_control(beta_true[1] + mntr[1]*beta_true[2] + xe*(beta_true[3]+mntr[1]*beta_true[4]),
						beta[1] + mnte[1]*beta[2] + xe*(beta[3]+mnte[1]*beta[4]))+
					P_BB * integrand_funct_control(beta_true[1] + mntr[2]*beta_true[2] + xe*(beta_true[3]+mntr[2]*beta_true[4]),
						beta[1] + mnte[2]*beta[2] + xe*(beta[3]+mnte[2]*beta[4]))
					)
			}
			res <- tryCatch(integrate(ll_fun00, -30*sd_e, 30*sd_e)$value, error = function(e) "error", warning = function(w) "warning")
			zzz <- 1
			while(res == "error"){
				zzz <- 0.5 + zzz
				res <- tryCatch(integrate(ll_fun00, -30/zzz*sd_e, 30/zzz*sd_e)$value, error = function(e) "error", warning = function(w) "warning")
			}
			return(res)
		}
		ll_g_e_fun <- function(beta){
			ll_fun00 <- function(xe){
				1 / sqrt(2 * pi * sd_e^2) * exp(-xe^2/(2 * sd_e^2)) * (
					P_AA * integrand_funct_case(beta_true[1] + xe*beta_true[3],
						beta[1] + xe*beta[3])+
					P_AB * integrand_funct_case(beta_true[1] + mntr[1]*beta_true[2] + xe*(beta_true[3]+mntr[1]*beta_true[4]),
						beta[1] + mnte[1]*beta[2] + xe*(beta[3]))+
					P_BB * integrand_funct_case(beta_true[1] + mntr[2]*beta_true[2] + xe*(beta_true[3]+mntr[2]*beta_true[4]),
						beta[1] + mnte[2]*beta[2] + xe*beta[3])+ 
					P_AA * integrand_funct_control(beta_true[1] + xe*beta_true[3],
						beta[1] + xe*beta[3])+
					P_AB * integrand_funct_control(beta_true[1] + mntr[1]*beta_true[2] + xe*(beta_true[3]+mntr[1]*beta_true[4]),
						beta[1] + mnte[1]*beta[2] + xe*(beta[3]))+
					P_BB * integrand_funct_control(beta_true[1] + mntr[2]*beta_true[2] + xe*(beta_true[3]+mntr[2]*beta_true[4]),
						beta[1] + mnte[2]*beta[2] + xe*beta[3]) 
					)
			}
			res <- tryCatch(integrate(ll_fun00, -30*sd_e, 30*sd_e)$value, error = function(e) "error", warning = function(w) "warning")
			zzz <- 1
			while(res == "error"){
				zzz <- 0.5 + zzz
				res <- tryCatch(integrate(ll_fun00, -30/zzz*sd_e, 30/zzz*sd_e)$value, error = function(e) "error", warning = function(w) "warning")
			}
			return(res)
		}
	}


	if(Test.Model == "2df"){
		# 1 = beta0
		# 2 = AB
		# 3 = BB
		# 4 = E
		# 5 = AB/E
		# 6 = BB/E
		mntr <- mod_num(True.Model)
		ll_fun <- function(beta){
			ll_fun00 <- function(xe){
				1 / sqrt(2 * pi * sd_e^2) * exp(-xe^2/(2 * sd_e^2)) * (
					P_AA * integrand_funct_case(beta_true[1] + xe*beta_true[3],
						beta[1] + xe*beta[4])+
					P_AB * integrand_funct_case(beta_true[1] + mntr[1]*beta_true[2] + xe*(beta_true[3]+mntr[1]*beta_true[4]),
						beta[1] + beta[2] + xe*(beta[4]+beta[5]))+
					P_BB * integrand_funct_case(beta_true[1] + mntr[2]*beta_true[2] + xe*(beta_true[3]+mntr[2]*beta_true[4]),
						beta[1] + beta[3] + xe*(beta[4]+beta[6]))+
					P_AA * integrand_funct_control(beta_true[1] + xe*beta_true[3],
						beta[1] + xe*beta[4])+
					P_AB * integrand_funct_control(beta_true[1] + mntr[1]*beta_true[2] + xe*(beta_true[3]+mntr[1]*beta_true[4]),
						beta[1] + beta[2] + xe*(beta[4]+beta[5]))+
					P_BB * integrand_funct_control(beta_true[1] + mntr[2]*beta_true[2] + xe*(beta_true[3]+mntr[2]*beta_true[4]),
						beta[1] + beta[3] + xe*(beta[4]+beta[6]))
					)
			}
			res <- tryCatch(integrate(ll_fun00, -30*sd_e, 30*sd_e)$value, error = function(e) "error", warning = function(w) "warning")
			zzz <- 1
			while(res == "error"){
				zzz <- 0.5 + zzz
				res <- tryCatch(integrate(ll_fun00, -30/zzz*sd_e, 30/zzz*sd_e)$value, error = function(e) "error", warning = function(w) "warning")
			}
			return(res)
		}
		ll_g_e_fun <- function(beta){
			ll_fun00 <- function(xe){
				1 / sqrt(2 * pi * sd_e^2) * exp(-xe^2/(2 * sd_e^2)) * (
					P_AA * integrand_funct_case(beta_true[1] + xe*beta_true[3],
						beta[1] + xe*beta[4])+
					P_AB * integrand_funct_case(beta_true[1] + mntr[1]*beta_true[2] + xe*(beta_true[3]+mntr[1]*beta_true[4]),
						beta[1] + beta[2] + xe*(beta[4]))+
					P_BB * integrand_funct_case(beta_true[1] + mntr[2]*beta_true[2] + xe*(beta_true[3]+mntr[2]*beta_true[4]),
						beta[1] + beta[3] + xe*(beta[4]))+ 
					P_AA * integrand_funct_control(beta_true[1] + xe*beta_true[3],
						beta[1] + xe*beta[4])+
					P_AB * integrand_funct_control(beta_true[1] + mntr[1]*beta_true[2] + xe*(beta_true[3]+mntr[1]*beta_true[4]),
						beta[1] + beta[2] + xe*(beta[4]))+
					P_BB * integrand_funct_control(beta_true[1] + mntr[2]*beta_true[2] + xe*(beta_true[3]+mntr[2]*beta_true[4]),
						beta[1] + beta[3] + xe*(beta[4])) 
					)
			}
			res <- tryCatch(integrate(ll_fun00, -30*sd_e, 30*sd_e)$value, error = function(e) "error", warning = function(w) "warning")
			zzz <- 1
			while(res == "error"){
				zzz <- 0.5 + zzz
				res <- tryCatch(integrate(ll_fun00, -30/zzz*sd_e, 30/zzz*sd_e)$value, error = function(e) "error", warning = function(w) "warning")
			}
			return(res)
		}

		beta <- optim(c(0,0,0,0,0,0), function(x) -ll_fun(x), control=c(abstol = 1e-5))$par
		ll_g_e <- -optim(c(0,0,0,0), function(x) -ll_g_e_fun(x), control=c(abstol = 1e-5))$value
	}else{
		
		beta <- optim(beta_true, function(x) -ll_fun(x), control=c(abstol = 1e-5))$par #"L-BFGS-B"	
		ll_g_e <- -optim(c(0,0,0), function(x) -ll_g_e_fun(x), control=c(abstol = 1e-5))$value

	}

	# optimization loop to make sure beta is optimized correctly
	diff <- 1
	ll0 <- ll_fun(beta) #, ll_fun(beta_true))
	while(diff > 1e-6){
		beta <- optim(beta, function(x) -ll_fun(x), control=c(abstol = 1e-5))$par #"L-BFGS-B"
		ll0 <- c(ll0, ll_fun(beta))
		diff <- ll0[length(ll0)] - ll0[length(ll0) - 1]
	}
	ll <- ll_fun(beta)
	# ll_g_e <- ll_g_e_fun(beta_g_e)
	# ll_e <- ll_e_fun(beta_e)
	# ll_g <- ll_g_fun(beta_g)

	# Case.Rate <- sum(t[1,])	
	if(is.null(power)){
		stat <- 2*(as.numeric(ll-ll_g_e))
		if(Test.Model=='2df'){
			power_res <- 1-pchisq(qchisq(1-Alpha, df=2, ncp=0), df=2, ncp = N*stat)
		}else{
			power_res <- pnorm(sqrt(N*2*(ll-ll_g_e)) - qnorm(1-Alpha/2)) + pnorm(-sqrt(N*2*(ll-ll_g_e)) - qnorm(1-Alpha/2))*1
		}
		return(power_res)
	}
	if(is.null(N)){
		stat <- 2*(as.numeric(ll-ll_g_e))
		if(Test.Model=='2df'){
			ss <- uniroot(function(x) ncp.search(x, power, Alpha, df=2),
							lower=0, upper=1000, extendInt = 'upX', tol=0.00001)$root/stat
		}else{
		  ss <- uniroot(function(x) ncp.search(x, power, Alpha, df=1),
		                      lower=0, upper=1000, extendInt = 'upX', tol=0.00001)$root/stat
		}
		return(ss)
	}
}
