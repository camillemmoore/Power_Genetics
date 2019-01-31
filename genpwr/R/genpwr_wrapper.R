# Wrapper function

genpwr.calc <- function(calc, model, ge.interaction = F, environment.logistic = F,
	N = NULL, Power = NULL,
	MAF = NULL, Alpha = 0.05,
	P_e = NULL, sd_e = NULL,
	sd_y = NULL, Case.Rate = NULL, k = NULL,
	OR = NULL, OR_G = NULL, OR_E = NULL, OR_GE = NULL, risk_allele = T,
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
	if(model %ni% c("logistic", "linear")) stop("'model' must be either logistic or linear")
	if(!is.logical(ge.interaction)) stop("'ge.interaction' must be TRUE or FALSE")
	if(!is.logical(environment.logistic)) stop("'environment.logistic' must be TRUE or FALSE")
	#foo <-  c("N", "Power", "MAF", "P_e", "sd_e", "sd_y", "Case.Rate", "k", "OR", "OR_G", "OR_E", "OR_GE", "ES", "ES_G", "ES_E", "ES_GE", "R2", "R2_G", "R2_E", "R2_GE")
	if(calc == "pow"){
		if(model == "logistic"){
			if(ge.interaction){
				if(environment.logistic){
					tnulls <- c("Power", "sd_e", "sd_y", "OR", "ES", "ES_G", "ES_E", "ES_GE", "R2", "R2_G", "R2_E", "R2_GE")
					if(any(!is.null(tnulls))) wnulls(tnulls[!is.null(tnulls)])
					power_envir.calc(N=N, Case.Rate=Case.Rate, k=k, MAF=MAF, OR_G=OR_G, OR_E=OR_E, OR_GE=OR_GE, P_e=P_e,
						Alpha=Alpha, True.Model=True.Model, Test.Model=Test.Model)
				}else{
					tnulls <- c("Power", "P_e", "sd_y", "OR", "ES", "ES_G", "ES_E", "ES_GE", "R2", "R2_G", "R2_E", "R2_GE")
					if(any(!is.null(tnulls))) wnulls(tnulls[!is.null(tnulls)])
					power_linear_envir.calc.logistic_outcome(N=N, MAF=MAF, OR_G=OR_G, OR_E=OR_E, OR_GE=OR_GE, sd_e=sd_e, 
						Case.Rate=Case.Rate, k=k, Alpha=Alpha, True.Model=True.Model, Test.Model=Test.Model)
				}
			}else{
				tnulls <- c("Power", "P_e", "sd_e", "sd_y", "OR_G", "OR_E", "OR_GE", "ES", "ES_G", "ES_E", "ES_GE", "R2", "R2_G", "R2_E", "R2_GE")
				power.calc(N=N, Case.Rate=Case.Rate, k=k, MAF=MAF, OR=OR,
					Alpha=Alpha, True.Model=True.Model, Test.Model=Test.Model)
			}
		}else if(model == "linear"){
			if(ge.interaction){
				if(environment.logistic){
					tnulls <- c("Power", "sd_e", "Case.Rate", "k", "OR", "OR_G", "OR_E", "OR_GE", "ES", "R2")
					if(any(!is.null(tnulls))) wnulls(tnulls[!is.null(tnulls)])
					power_envir.calc.linear_outcome(N=N, MAF=MAF, ES_G=ES_G, ES_E=ES_E, ES_GE=ES_GE, P_e=P_e, 
						R2_G=R2_G, R2_E=R2_E, R2_GE=R2_GE, sd_y=sd_y, Alpha=Alpha, True.Model=True.Model, Test.Model=Test.Model)
				}else{
					tnulls <- c("Power", "P_e", "Case.Rate", "k", "OR", "OR_G", "OR_E", "OR_GE", "ES", "R2")
					if(any(!is.null(tnulls))) wnulls(tnulls[!is.null(tnulls)])
					power_linear_envir.calc.linear_outcome(N=N, MAF=MAF, ES_G=ES_G, ES_E=ES_E, ES_GE=ES_GE, sd_e=sd_e, 
						R2_G=R2_G, R2_E=R2_E, R2_GE=R2_GE, sd_y=sd_y, Alpha=Alpha, True.Model=True.Model, Test.Model=Test.Model)
				}
			}else{
				tnulls <- c("Power", "P_e", "sd_e", "Case.Rate", "k", "OR", "OR_G", "OR_E", "OR_GE", "ES_G", "ES_E", "ES_GE", "R2_G", "R2_E", "R2_GE")
				if(any(!is.null(tnulls))) wnulls(tnulls[!is.null(tnulls)])
				power.calc.linear(N=N, MAF=MAF, ES=ES,R2=R2, sd_y=sd_y,
					Alpha=Alpha, True.Model=True.Model, Test.Model=Test.Model)
			}
		}
	}else if(calc == "ss"){
		if(model == "logistic"){
			if(ge.interaction){
				if(environment.logistic){
					tnulls <- c("N", "Power", "sd_e", "sd_y", "OR", "ES", "ES_G", "ES_E", "ES_GE", "R2", "R2_G", "R2_E", "R2_GE")
					if(any(!is.null(tnulls))) wnulls(tnulls[!is.null(tnulls)])
					ss_envir.calc(power=Power, Case.Rate=Case.Rate, k=k, MAF=MAF, OR_G=OR_G, OR_E=OR_E, OR_GE=OR_GE, P_e=P_e,
						Alpha=Alpha, True.Model=True.Model, Test.Model=Test.Model)
				}else{
					tnulls <- c("N", "Power", "P_e", "sd_y", "OR", "ES", "ES_G", "ES_E", "ES_GE", "R2", "R2_G", "R2_E", "R2_GE")
					if(any(!is.null(tnulls))) wnulls(tnulls[!is.null(tnulls)])
					ss_linear_envir.calc.logistic_outcome(power=Power, MAF=MAF, OR_G=OR_G, OR_E=OR_E, OR_GE=OR_GE, sd_e=sd_e, 
						Case.Rate=Case.Rate, k=k, Alpha=Alpha, True.Model=True.Model, Test.Model=Test.Model)
				}
			}else{
				tnulls <- c("N", "Power", "P_e", "sd_e", "sd_y", "OR_G", "OR_E", "OR_GE", "ES", "ES_G", "ES_E", "ES_GE", "R2", "R2_G", "R2_E", "R2_GE")
				if(any(!is.null(tnulls))) wnulls(tnulls[!is.null(tnulls)])
				ss.calc(power=Power, Case.Rate=Case.Rate, k=k, MAF=MAF, OR=OR,
					Alpha=Alpha, True.Model=True.Model, Test.Model=Test.Model)
			}
		}else if(model == "linear"){
			if(ge.interaction){
				if(environment.logistic){
					tnulls <- c("N", "Power", "sd_e", "Case.Rate", "k", "OR", "OR_G", "OR_E", "OR_GE", "ES", "R2")
					if(any(!is.null(tnulls))) wnulls(tnulls[!is.null(tnulls)])
					ss_envir.calc.linear_outcome(pow=Power, MAF=MAF, ES_G=ES_G, ES_E=ES_E, ES_GE=ES_GE, P_e=P_e, 
						R2_G=R2_G, R2_E=R2_E, R2_GE=R2_GE, sd_y=sd_y, Alpha=Alpha, True.Model=True.Model, Test.Model=Test.Model)
				}else{
					tnulls <- c("N", "Power", "P_e", "Case.Rate", "k", "OR", "OR_G", "OR_E", "OR_GE", "ES", "R2")
					if(any(!is.null(tnulls))) wnulls(tnulls[!is.null(tnulls)])
					ss_linear_envir.calc.linear_outcome(pow=Power, MAF=MAF, ES_G=ES_G, ES_E=ES_E, ES_GE=ES_GE, sd_e=sd_e, 
						R2_G=R2_G, R2_E=R2_E, R2_GE=R2_GE, sd_y=sd_y, Alpha=Alpha, True.Model=True.Model, Test.Model=Test.Model)
				}
			}else{
				tnulls <- c("N", "P_e", "sd_e", "Case.Rate", "k", "OR", "OR_G", "OR_E", "OR_GE", "ES_G", "ES_E", "ES_GE", "R2_G", "R2_E", "R2_GE")
				if(any(!is.null(tnulls))) wnulls(tnulls[!is.null(tnulls)])
				ss.calc.linear(power=Power, MAF=MAF, ES=ES, R2=OR, sd_y=sd_y,
					Alpha=Alpha, True.Model=True.Model, Test.Model=Test.Model)
			}
		}
	}else if(calc == "es"){
		if(model == "logistic"){
			if(ge.interaction){
				# if(environment.logistic){
					stop("functionality for effect size calculation with environment interaction does not yet exist")
				# }else{
				# }
			}else{
				tnulls <- c("P_e", "sd_e", "sd_y", "OR", "OR_G", "OR_E", "OR_GE", "ES", "ES_G", "ES_E", "ES_GE", "R2", "R2_G", "R2_E", "R2_GE")
				if(any(!is.null(tnulls))) wnulls(tnulls[!is.null(tnulls)])
				odds_ratio_function(N=N, Case.Rate=Case.Rate, k=k, MAF=MAF, power=Power, risk_allele = risk_allele,
					Alpha=Alpha, True.Model=True.Model, Test.Model=Test.Model)
			}
		}else if(model == "linear"){
			if(ge.interaction){
				# if(environment.logistic){
					stop("functionality for effect size calculation with environment interaction does not yet exist")
				# }else{
				# }
			}else{
				tnulls <- c("P_e", "sd_e", "Case.Rate", "OR", "OR_G", "OR_E", "OR_GE", "ES", "ES_G", "ES_E", "ES_GE", "R2", "R2_G", "R2_E", "R2_GE")
				if(any(!is.null(tnulls))) wnulls(tnulls[!is.null(tnulls)])
				es.calc.linear(power=Power, N=N, MAF=MAF, sd_y=sd_y,
					Alpha=Alpha, True.Model=True.Model, Test.Model=Test.Model)
			}
		}
	}
}
