#' Function to Calculate Power for Linear Models with logistic environment interaction
#'
#' Calculates the power to detect an difference in means/effect size/regression coefficient, at a given sample size, N, with type 1 error rate, Alpha
#'
#' @param N Vector of the desired sample size(s)
#' @param Alpha the desired type 1 error rate(s)
#' @param MAF Vector of minor allele frequencies
#' @param sd_e Standard deviation of the environmental variable
#' @param OR_G Vector of genetic odds ratios to detect
#' @param OR_E Vector of environmental odds ratios to detect
#' @param OR_GE Vector of genetic/environmental interaction odds ratios to detect
#' @param Case.Rate Standard deviation of the outcome in the population (ignoring genotype). Either Case.Rate_x or Case.Rate must be specified.
#' @param k Vector of the number of controls per case. Either k or Case.Rate must be specified.
#' @param True.Model A vector specifying the true underlying genetic model(s): 'Dominant', 'Additive', 'Recessive' or 'All'
#' @param Test.Model A vector specifying the assumed genetic model(s) used in testing: 'Dominant', 'Additive', 'Recessive' or 'All'
#' @param compareQuanto For comparison with Quanto results - uses Quanto's formula to calculate results
#'
#' @return A data frame including the power for all combinations of the specified parameters (Case.Rate, ES, Power, etc)
#'
#' @examples
#' pw <- power_linear_envir.calc.logistic_outcome(N=c(30, 100), OR_G=c(1.1,2), OR_E=c(1.2, 1.4), OR_GE=c(1.5, 2), 
#' 	sd_e = c(1, 1.1), MAF=c(0.1, 0.2), Case.Rate = c(0.2, 0.3),Alpha=c(0.02, 0.05),
#' 	True.Model="All", Test.Model="All")
#' 
#'
#' @export
#'
power_linear_envir.calc.logistic_outcome <- function(N=NULL, MAF=NULL, OR_G=NULL, OR_E=NULL, OR_GE=NULL, sd_e=NULL, 
		Case.Rate=NULL, k=NULL, Alpha=0.05, True.Model='All', Test.Model='All', compareQuanto = 0)
{

	compareQuanto = 1 * compareQuanto
	library(MASS)

	############################################################################################################
	#Error Messages for insufficient sample size information, MAF, and case vs. control ratio
	############################################################################################################
	if(is.null(N)){
		stop("N, the total sample size, must be specified.")
	}

	if(is.null(MAF)){
		stop("MAF (minor allele frequency) must be specified.")
	}

	if(all(is.null(c(OR_G, OR_E, OR_GE)))) {
		stop(paste0("OR_G, OR_E, and OR_GE (odds ratios for gene, environment,", 
			"and gene-environment interaction) must be specified."))
	}

	if(is.null(k)==F & is.null(Case.Rate)==F){
		stop("Specify one of k, the number of controls per case, or Case.Rate, the proportion of cases in the study sample, not both.")
	}


	############################################################################################################
	#Error Messages for out of range values
	############################################################################################################
	if(sum(Case.Rate>=1)>0 | sum(Case.Rate<=0)>0){
		stop("Case.Rate must be greater than 0 and less than 1.")
	}

	if(sum(MAF>=1)>0 | sum(MAF<=0)>0){
		stop("MAF must be greater than 0 and less than 1.")
	}

	if(sum(N<=0)>0){
		stop("N must be greater than 0.")
	}

	if(sum(Alpha>=1)>0 | sum(Alpha<=0)>0){
		stop("Alpha must be greater than 0 and less than 1.")
	}

	if(any(Test.Model %in% c("Additive1", "Additive2")) | any(True.Model %in% c("Additive1", "Additive2"))) {
		print(paste0("For additive models, this function treats effect size as the difference between the homozygous for major",
					"allele and heterozygous, and 'Additive1' and 'Additive2' will be converted to 'Additive'"))
	}

	if(sum(!(Test.Model %in% c("Dominant", "Recessive", "Additive", "2df", "All")))>0){
		stop(paste("Invalid Test.Model:",
			   paste(Test.Model[!(Test.Model %in% c("Dominant", "Recessive", "Additive", "2df", "All"))], collapse=', ')))
	}

	if(length(setdiff(True.Model, c("Dominant", "Recessive", "Additive", "All")))>0){
		stop(paste("Invalid True.Model:",
			paste(setdiff(True.Model, c("Dominant", "Recessive", "Additive", "All")), collapse=', ')))
	}
	############################################################################################################
	#Create model vectors if model = 'All'
	############################################################################################################
	#Test model vector
	if('All' %in% Test.Model){Test.Model<-c("Dominant", "Recessive", "Additive", "2df")}

	#True model vector
	if('All' %in% True.Model){True.Model<-c("Dominant", "Recessive", "Additive")}

	############################################################################################################
	# calculate beta0 for different models
	############################################################################################################
	
	beta_fun_dist <- function(xe, beta0, mod, sd_e, P_AA, P_AB, P_BB, OR_G, OR_E, OR_GE){
		beta_g <- log(OR_G)
		beta_e <- log(OR_E)
		beta_ge <- log(OR_GE)
		1 / sqrt(2 * pi * sd_e^2) * exp(-xe^2/(2 * sd_e^2)) * (
			P_AA*(exp(beta0 + xe*beta_e)/(1 + exp(beta0 + xe*beta_e)))+
			P_AB*(exp(beta0 + beta_g + xe*(beta_e + beta_ge))/(1 + exp(beta0 + beta_g + xe*(beta_e + beta_ge))))+
			P_BB*(exp(beta0 + beta_g + xe*(beta_e + beta_ge))/(1 + exp(beta0 + beta_g + xe*(beta_e + beta_ge))))
			) 
	}

	# P_AA <- (1-MAF)^2
	# P_AB <- 2*MAF*(1-MAF)
	# P_BB <- MAF^2
	# beta0 <- c(-OR_G*(P_AB + P_BB), -OR_G*P_BB, -OR_G*(P_AB + 2*P_BB))
	# names(beta0) <- c("Dominant", "Recessive", "Additive")
	e.save.tab <- expand.grid(OR_G, OR_E, OR_GE, MAF, sd_e, Case.Rate, True.Model)
	names(e.save.tab) <- c("OR_G", "OR_E", "OR_GE", "MAF", "sd_e", "Case.Rate", "True.Model") #Alpha
	e.save.tab$P_AA <- (1-e.save.tab$MAF)^2
	e.save.tab$P_AB <- 2*e.save.tab$MAF*(1-e.save.tab$MAF)
	e.save.tab$P_BB <- e.save.tab$MAF^2
	e.save.tab$beta0 <- NA

	#calculating our beta0's
	for(iij in 1:nrow(e.save.tab)){
		beta_fun_dist2 <- function(xe, beta0) beta_fun_dist(xe = xe, beta0 = beta0, sd_e = e.save.tab[iij, "sd_e"],
			OR_G = e.save.tab[iij, "OR_G"], OR_E = e.save.tab[iij, "OR_E"], OR_GE = e.save.tab[iij, "OR_GE"], 
			mod = e.save.tab[iij, "mod_iij"], P_AA = e.save.tab[iij, "P_AA"], P_AB = e.save.tab[iij, "P_AB"], 
			P_BB = e.save.tab[iij, "P_BB"])
		beta_fun_dist_int <- function(beta0){
			# beta_fun_dist0 <- function(xe) beta_fun_dist(xe, beta0)
			integrate(function(xe) beta_fun_dist2(xe, beta0), 
				-30*e.save.tab[iij,"sd_e"], 30*e.save.tab[iij,"sd_e"])$value - e.save.tab[iij, "Case.Rate"]
		}
		e.save.tab$beta0[iij] <- nleqslv(x = 0, fn = beta_fun_dist_int)$x
	}

	# e.save.tab <- rbindAlpha(
	# 	cbind(e.save.tab, True.Model = "Dominant", 
	# 		beta0 = nleqslv(x = 0, fn = function(beta0){
	# 				integrate(function(xe) beta_fun_dist(xe, beta0), -Inf, Inf)$value - Case.Rate
	# 			})$x),
	# 	cbind(e.save.tab, True.Model = "Additive", beta0 = -e.save.tab$OR_G*(e.save.tab$P_AB + 2*e.save.tab$P_BB)),
	# 	cbind(e.save.tab, True.Model = "Recessive", beta0 = -e.save.tab$OR_G*e.save.tab$P_BB))

	
	############################################################################################################
	# Calculate Power for each scenario in e.save.tab under the specified testing model
	############################################################################################################
	power.tab <- NULL

	############################################################################################################
	#Loop over sample size
	############################################################################################################
	for (ss in N){

		################################################################################################
		#Loop over all of the testing models and calculate power for each ES, SD, and MAF scenario
		################################################################################################
		for (mod in Test.Model){
			temp <- NULL
			for(alpha0 in Alpha){
				temp.0 <- NULL

				pow <- mapply(function(x){ll.ge.logistic.lin.envir(
					sd_e = e.save.tab[x,"sd_e"],
					MAF = e.save.tab[x,"MAF"],
					N =ss,
					beta0 = e.save.tab[x,"beta0"],
					OR_G = e.save.tab[x,"OR_G"],
					OR_E = e.save.tab[x,"OR_E"],
					OR_GE = e.save.tab[x,"OR_GE"],
					Alpha = alpha0,
					True.Model = e.save.tab[x,"True.Model"],
					Test.Model = mod,
					compareQuanto = compareQuanto)}, seq(1:nrow(e.save.tab)))

				temp.0 <- cbind(temp.0, pow)
				colnames(temp.0) <- paste0("Power_at_Alpha_", alpha0)
				temp <- cbind(temp.0, temp)

				# ll.stat = 2*(ll.alt-ll.reduced)
				# if(length(Alpha)>1){pow <- t(pow)
				# rownames(pow) <- seq(1:nrow(pow))}

			}
			#Save the power calculations for each testing model in a final table for the sample size and case rate
			power.tab<-rbind(power.tab,cbind(data.frame(Test.Model=mod, True.Model = as.character(e.save.tab[, "True.Model"]),
							MAF = e.save.tab[, "MAF"], N = ss, sd_e = e.save.tab[, "sd_e"], OR_G = e.save.tab[, "OR_G"], 
							OR_E = e.save.tab[, "OR_E"], OR_GE = e.save.tab[, "OR_GE"], 
							Case.Rate = e.save.tab[, "Case.Rate"]), temp),row.names = NULL)

		}
	}
	colnames(power.tab)<-c("Test.Model", "True.Model", "MAF", "N_total", "SD_E", "OR_G", "OR_E", "OR_GE", 
				"Case.Rate", paste0("Power_at_Alpha_", Alpha))

	return(power.tab)
}
