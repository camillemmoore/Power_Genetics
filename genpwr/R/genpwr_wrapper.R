#' Function to Calculate Power for Linear Models with logistic environment interaction
#'
#' Calculates the power to detect an difference in means/effect size/regression coefficient, at a given sample size, N, with type 1 error rate, Alpha
#'
#' @param calc What kind of calculation to perform? sample size ("ss"), power ("power"), or effect size ("es")
#' @param model Distribution of the outcome variable? ("logistic" or "linear")
#' @param ge.interaction If no environment interaction, should be NULL, otherwise should be "logistic" or "linear"
#' @param N Vector of the desired sample size(s)
#' @param Power Vector of the desired power(s)
#' @param Alpha the desired type 1 error rate(s)
#' @param MAF Vector of minor allele frequencies
#' @param sd_e Standard deviation of the environmental variable
#' @param sd_y Standard deviation of the outcome in the population (ignoring genotype). Either sd_y_x or sd_y must be specified.
#' @param OR Vector of genetic odds ratios to detect in absence of environmental odds ratios
#' @param OR_G Vector of genetic odds ratios to detect
#' @param OR_E Vector of environmental odds ratios to detect
#' @param OR_GE Vector of genetic/environmental interaction odds ratios to detect
#' @param ES Vector of effect sizes (difference in means) to detect. Either ES or R2 must be specified.
#' @param ES_G Vector of genetic effect sizes (difference in means) to detect. Either ES_G, ES_E, and ES_EG or R2_G, R2_E, and R2_EG must be specified.
#' @param ES_E Vector of environmental effect sizes (difference in means) to detect. Either ES_G, ES_E, and ES_EG or R2_G, R2_E, and R2_EG must be specified.
#' @param ES_GE Vector of genetic/environment interaction effect sizes (difference in means) to detect. Either ES_G, ES_E, and ES_EG or R2_G, R2_E, and R2_EG must be specified.
#' @param R2 Vector of R-squared values to detect. Either ES or R2 must be specified.
#' @param R2_G Vector of genetic R-squared values to detect. Either ES_G, ES_E, and ES_EG or R2_G, R2_E, and R2_EG must be specified.
#' @param R2_E Vector of environmental R-squared values to detect. Either ES_G, ES_E, and ES_EG or R2_G, R2_E, and R2_EG must be specified.
#' @param R2_GE Vector of genetic/environment interaction R-squared values Either ES_G, ES_E, and ES_EG or R2_G, R2_E, and R2_EG must be specified.
#' @param Case.Rate Standard deviation of the outcome in the population (ignoring genotype). Either Case.Rate_x or Case.Rate must be specified.
#' @param k Vector of the number of controls per case. Either k or Case.Rate must be specified.
#' @param P_e Vector of proportions of the population with exposure to the environmental effect
#' @param risk_allele Logical: If OR > 1, the allele is classified as a "risk allele"
#' @param True.Model A vector specifying the true underlying genetic model(s): 'Dominant', 'Additive', 'Recessive' or 'All'
#' @param Test.Model A vector specifying the assumed genetic model(s) used in testing: 'Dominant', 'Additive', 'Recessive' or 'All'
#'
#' @return A data frame including the power for all combinations of the specified parameters (Case.Rate, ES, Power, etc)
#'
#' @examples
#' pw <- genpwr.calc(calc = "power", model = "logistic", ge.interaction = "continuous",
#' 	N=c(30, 100), OR_G=c(1.1,2), OR_E=c(1.2, 1.4), OR_GE=c(1.5, 2), 
#' 	sd_e = c(1, 1.1), MAF=c(0.1, 0.2), Case.Rate = c(0.2, 0.3),Alpha=c(0.02, 0.05),
#' 	True.Model="All", Test.Model="All")
#' 
#'
#' @export
#'
genpwr.calc <- function(calc, model, ge.interaction = NULL,
	N = NULL, Power = NULL,
	MAF = NULL, Alpha = 0.05,
	P_e = NULL, sd_e = NULL,
	sd_y = NULL, Case.Rate = NULL, k = NULL,
	OR = NULL, OR_G = NULL, OR_E = NULL, OR_GE = NULL, risk_allele = TRUE,
	ES = NULL, ES_G = NULL, ES_E = NULL, ES_GE = NULL,
	R2 = NULL, R2_G = NULL, R2_E = NULL, R2_GE = NULL,
	True.Model = "All", Test.Model = "All")
{
	"%ni%" <- Negate("%in%")
	calc <- tolower(calc)
	if(calc %ni% c("ss", "n", "es", "pow", "power")) stop("'calc' must be one of: 'ss' (sample size), 'es' (effect size), 'power'")
	if(calc %in% c("n", "power")){
		dd <- c("ss", "pow")
		names(dd) <- c("n", "power")
		calc <- dd[calc]
	}
	wnulls <- function(x) print(sprintf("The following parameters were given a value but not used: %s", paste0(x, collapse = ", ")))
	if(model %ni% c("logistic", "linear")) stop("'model' must be either 'logistic' or 'linear'")
	if(!is.null(ge.interaction)){
		if(ge.interaction %ni% c("binary", "continuous"))
			stop("'ge.interaction' must be either 'binary' or 'continuous' if not NULL")
	}
	# if(!is.logical(environment.logistic)) stop("'environment.logistic' must be TRUE or FALSE")
	#foo <-  c("N", "Power", "MAF", "P_e", "sd_e", "sd_y", "Case.Rate", "k", "OR", "OR_G", "OR_E", "OR_GE", "ES", "ES_G", "ES_E", "ES_GE", "R2", "R2_G", "R2_E", "R2_GE")
	if(calc == "pow"){
		if(model == "logistic"){
			if(!is.null(ge.interaction)){
				if(ge.interaction == "binary"){
					tnulls <- c("Power", "sd_e", "sd_y", "OR", "ES", "ES_G", "ES_E", "ES_GE", "R2", "R2_G", "R2_E", "R2_GE")
					if(any(sapply(tnulls, function(tn) !is.null(get(tn))))) wnulls(tnulls[sapply(tnulls, function(tn) !is.null(get(tn)))])
					power_envir.calc(N=N, Case.Rate=Case.Rate, k=k, MAF=MAF, OR_G=OR_G, OR_E=OR_E, OR_GE=OR_GE, P_e=P_e,
						Alpha=Alpha, True.Model=True.Model, Test.Model=Test.Model)
				}else{
					tnulls <- c("Power", "P_e", "sd_y", "OR", "ES", "ES_G", "ES_E", "ES_GE", "R2", "R2_G", "R2_E", "R2_GE")
					if(any(sapply(tnulls, function(tn) !is.null(get(tn))))) wnulls(tnulls[sapply(tnulls, function(tn) !is.null(get(tn)))])
					power_linear_envir.calc.logistic_outcome(N=N, MAF=MAF, OR_G=OR_G, OR_E=OR_E, OR_GE=OR_GE, sd_e=sd_e, 
						Case.Rate=Case.Rate, k=k, Alpha=Alpha, True.Model=True.Model, Test.Model=Test.Model)
				}
			}else{
				tnulls <- c("Power", "P_e", "sd_e", "sd_y", "OR_G", "OR_E", "OR_GE", "ES", "ES_G", "ES_E", "ES_GE", "R2", "R2_G", "R2_E", "R2_GE")
				if(any(sapply(tnulls, function(tn) !is.null(get(tn))))) wnulls(tnulls[sapply(tnulls, function(tn) !is.null(get(tn)))])
				power.calc(N=N, Case.Rate=Case.Rate, k=k, MAF=MAF, OR=OR,
					Alpha=Alpha, True.Model=True.Model, Test.Model=Test.Model)
			}
		}else if(model == "linear"){
			if(!is.null(ge.interaction)){
				if(ge.interaction == "binary"){
					tnulls <- c("Power", "sd_e", "Case.Rate", "k", "OR", "OR_G", "OR_E", "OR_GE", "ES", "R2")
					if(any(sapply(tnulls, function(tn) !is.null(get(tn))))) wnulls(tnulls[sapply(tnulls, function(tn) !is.null(get(tn)))])
					power_envir.calc.linear_outcome(N=N, MAF=MAF, ES_G=ES_G, ES_E=ES_E, ES_GE=ES_GE, P_e=P_e, 
						R2_G=R2_G, R2_E=R2_E, R2_GE=R2_GE, sd_y=sd_y, Alpha=Alpha, True.Model=True.Model, Test.Model=Test.Model)
				}else{
					tnulls <- c("Power", "P_e", "Case.Rate", "k", "OR", "OR_G", "OR_E", "OR_GE", "ES", "R2")
					if(any(sapply(tnulls, function(tn) !is.null(get(tn))))) wnulls(tnulls[sapply(tnulls, function(tn) !is.null(get(tn)))])
					power_linear_envir.calc.linear_outcome(N=N, MAF=MAF, ES_G=ES_G, ES_E=ES_E, ES_GE=ES_GE, sd_e=sd_e, 
						R2_G=R2_G, R2_E=R2_E, R2_GE=R2_GE, sd_y=sd_y, Alpha=Alpha, True.Model=True.Model, Test.Model=Test.Model)
				}
			}else{
				tnulls <- c("Power", "P_e", "sd_e", "Case.Rate", "k", "OR", "OR_G", "OR_E", "OR_GE", "ES_G", "ES_E", "ES_GE", "R2_G", "R2_E", "R2_GE")
				if(any(sapply(tnulls, function(tn) !is.null(get(tn))))) wnulls(tnulls[sapply(tnulls, function(tn) !is.null(get(tn)))])
				power.calc.linear(N=N, MAF=MAF, ES=ES,R2=R2, sd_y=sd_y,
					Alpha=Alpha, True.Model=True.Model, Test.Model=Test.Model)
			}
		}
	}else if(calc == "ss"){
		if(model == "logistic"){
			if(!is.null(ge.interaction)){
				if(ge.interaction == "binary"){
					tnulls <- c("N", "sd_e", "sd_y", "OR", "ES", "ES_G", "ES_E", "ES_GE", "R2", "R2_G", "R2_E", "R2_GE")
					if(any(sapply(tnulls, function(tn) !is.null(get(tn))))) wnulls(tnulls[sapply(tnulls, function(tn) !is.null(get(tn)))])
					ss_envir.calc(power=Power, Case.Rate=Case.Rate, k=k, MAF=MAF, OR_G=OR_G, OR_E=OR_E, OR_GE=OR_GE, P_e=P_e,
						Alpha=Alpha, True.Model=True.Model, Test.Model=Test.Model)
				}else{
					tnulls <- c("N","P_e", "sd_y", "OR", "ES", "ES_G", "ES_E", "ES_GE", "R2", "R2_G", "R2_E", "R2_GE")
					if(any(sapply(tnulls, function(tn) !is.null(get(tn))))) wnulls(tnulls[sapply(tnulls, function(tn) !is.null(get(tn)))])
					ss_linear_envir.calc.logistic_outcome(power=Power, MAF=MAF, OR_G=OR_G, OR_E=OR_E, OR_GE=OR_GE, sd_e=sd_e, 
						Case.Rate=Case.Rate, k=k, Alpha=Alpha, True.Model=True.Model, Test.Model=Test.Model)
				}
			}else{
				tnulls <- c("N", "P_e", "sd_e", "sd_y", "OR_G", "OR_E", "OR_GE", "ES", "ES_G", "ES_E", "ES_GE", "R2", "R2_G", "R2_E", "R2_GE")
				if(any(sapply(tnulls, function(tn) !is.null(get(tn))))) wnulls(tnulls[sapply(tnulls, function(tn) !is.null(get(tn)))])
				ss.calc(power=Power, Case.Rate=Case.Rate, k=k, MAF=MAF, OR=OR,
					Alpha=Alpha, True.Model=True.Model, Test.Model=Test.Model)
			}
		}else if(model == "linear"){
			if(!is.null(ge.interaction)){
				if(ge.interaction == "binary"){
					tnulls <- c("N", "Power", "sd_e", "Case.Rate", "k", "OR", "OR_G", "OR_E", "OR_GE", "ES", "R2")
					if(any(sapply(tnulls, function(tn) !is.null(get(tn))))) wnulls(tnulls[sapply(tnulls, function(tn) !is.null(get(tn)))])
					ss_envir.calc.linear_outcome(pow=Power, MAF=MAF, ES_G=ES_G, ES_E=ES_E, ES_GE=ES_GE, P_e=P_e, 
						R2_G=R2_G, R2_E=R2_E, R2_GE=R2_GE, sd_y=sd_y, Alpha=Alpha, True.Model=True.Model, Test.Model=Test.Model)
				}else{
					tnulls <- c("N", "Power", "P_e", "Case.Rate", "k", "OR", "OR_G", "OR_E", "OR_GE", "ES", "R2")
					if(any(sapply(tnulls, function(tn) !is.null(get(tn))))) wnulls(tnulls[sapply(tnulls, function(tn) !is.null(get(tn)))])
					ss_linear_envir.calc.linear_outcome(pow=Power, MAF=MAF, ES_G=ES_G, ES_E=ES_E, ES_GE=ES_GE, sd_e=sd_e, 
						R2_G=R2_G, R2_E=R2_E, R2_GE=R2_GE, sd_y=sd_y, Alpha=Alpha, True.Model=True.Model, Test.Model=Test.Model)
				}
			}else{
				tnulls <- c("N", "P_e", "sd_e", "Case.Rate", "k", "OR", "OR_G", "OR_E", "OR_GE", "ES_G", "ES_E", "ES_GE", "R2_G", "R2_E", "R2_GE")
				if(any(sapply(tnulls, function(tn) !is.null(get(tn))))) wnulls(tnulls[sapply(tnulls, function(tn) !is.null(get(tn)))])
				ss.calc.linear(power=Power, MAF=MAF, ES=ES, R2=OR, sd_y=sd_y,
					Alpha=Alpha, True.Model=True.Model, Test.Model=Test.Model)
			}
		}
	}else if(calc == "es"){
		if(model == "logistic"){
			if(!is.null(ge.interaction)){
				# if(environment.logistic){
					stop("functionality for effect size calculation with environment interaction does not yet exist")
				# }else{
				# }
			}else{
				tnulls <- c("P_e", "sd_e", "sd_y", "OR", "OR_G", "OR_E", "OR_GE", "ES", "ES_G", "ES_E", "ES_GE", "R2", "R2_G", "R2_E", "R2_GE")
				if(any(sapply(tnulls, function(tn) !is.null(get(tn))))) wnulls(tnulls[sapply(tnulls, function(tn) !is.null(get(tn)))])
				odds_ratio_function(N=N, Case.Rate=Case.Rate, k=k, MAF=MAF, power=Power, risk_allele = risk_allele,
					Alpha=Alpha, True.Model=True.Model, Test.Model=Test.Model)
			}
		}else if(model == "linear"){
			if(!is.null(ge.interaction)){
				# if(environment.logistic){
					stop("functionality for effect size calculation with environment interaction does not yet exist")
				# }else{
				# }
			}else{
				tnulls <- c("P_e", "sd_e", "Case.Rate", "OR", "OR_G", "OR_E", "OR_GE", "ES", "ES_G", "ES_E", "ES_GE", "R2", "R2_G", "R2_E", "R2_GE")
				if(any(sapply(tnulls, function(tn) !is.null(get(tn))))) wnulls(tnulls[sapply(tnulls, function(tn) !is.null(get(tn)))])
				es.calc.linear(power=Power, N=N, MAF=MAF, sd_y=sd_y,
					Alpha=Alpha, True.Model=True.Model, Test.Model=Test.Model)
			}
		}
	}
}
