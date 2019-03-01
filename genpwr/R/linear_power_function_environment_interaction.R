#' Function to Calculate Power for Linear Models with logistic environment interaction
#'
#' Calculates the power to detect an difference in means/effect size/regression coefficient, at a given sample size, N, with type 1 error rate, Alpha
#'
#' @param N Vector of the desired sample size(s)
#' @param Alpha the desired type 1 error rate(s)
#' @param MAF Vector of minor allele frequencies
#' @param sd_y Standard deviation of the outcome in the population (ignoring genotype). Either sd_y_x or sd_y must be specified.
#' @param ES_G Vector of genetic effect sizes (difference in means) to detect. Either ES_G, ES_E, and ES_EG or R2_G, R2_E, and R2_EG must be specified.
#' @param ES_E Vector of environmental effect sizes (difference in means) to detect. Either ES_G, ES_E, and ES_EG or R2_G, R2_E, and R2_EG must be specified.
#' @param ES_GE Vector of genetic/environment interaction effect sizes (difference in means) to detect. Either ES_G, ES_E, and ES_EG or R2_G, R2_E, and R2_EG must be specified.
#' @param R2_G Vector of genetic R-squared values to detect. Either ES_G, ES_E, and ES_EG or R2_G, R2_E, and R2_EG must be specified.
#' @param R2_E Vector of environmental R-squared values to detect. Either ES_G, ES_E, and ES_EG or R2_G, R2_E, and R2_EG must be specified.
#' @param R2_GE Vector of genetic/environment interaction R-squared valuesEither ES_G, ES_E, and ES_EG or R2_G, R2_E, and R2_EG must be specified.
#' @param P_e Vector of proportions of the population with exposure to the environmental effect
#' @param True.Model A vector specifying the true underlying genetic model(s): 'Dominant', 'Additive', 'Recessive' or 'All'
#' @param Test.Model A vector specifying the assumed genetic model(s) used in testing: 'Dominant', 'Additive', 'Recessive' or 'All'
#'
#' @return A data frame including the power for all combinations of the specified parameters (Case.Rate, ES, Power, etc)
#'
#' @examples
#' power_envir.calc.linear_outcome(N=c(100,200), ES_G = c(1.2, 1.5), ES_E = c(1.3, 2.1), 
#' 	ES_GE = c(0.5, 2), Alpha = 0.05, MAF = c(0.1, 0.2), P_e = c(0.2, 0.3), 
#' 	sd_y = 10, True.Model = "All", Test.Model = "All")
#'
#' @export
#'
power_envir.calc.linear_outcome <- function(N=NULL, MAF=NULL, ES_G=NULL, ES_E=NULL, ES_GE=NULL, P_e=NULL, 
		R2_G=NULL, R2_E=NULL, R2_GE=NULL, sd_y=NULL,Alpha=0.05, True.Model='All', Test.Model='All')
{
	#library(MASS)

	############################################################################################################
	#Error Messages for insufficient sample size information, MAF, and case vs. control ratio
	############################################################################################################
	if(is.null(N)){
		stop("N, the total sample size, must be specified.")
	}

	if(is.null(MAF)){
		stop("MAF (minor allele frequency) must be specified.")
	}

	if(all(is.null(c(ES_G, ES_E, ES_GE))) & all(is.null(c(R2_G, R2_E, R2_GE)))){
		stop(paste0("Either ES_G, ES_E, and ES_GE (detectable effect sizes for gene, environment,", 
			"and gene-environment interaction) or R2_G, R2_E, and R2_GE (detectable R-squared for ", 
			"gene, environment, and gene-environment interaction) must be specified."))
	}

	# print(c(ES_G, ES_E, ES_GE, R2_G, R2_E, R2_GE))

	if(any(!is.null(c(ES_G, ES_E, ES_GE))) & any(!is.null(c(R2_G, R2_E, R2_GE)))){
		stop(paste0("Specify either all of ES_G, ES_E, and ES_GE (detectable effect sizes for gene, environment,", 
			"and gene-environment interaction) or all of R2_G, R2_E, and R2_GE (detectable R-squared for ", 
			"gene, environment, and gene-environment interaction), not any combinations of both"))
	}

	if(any(!is.null(c(ES_G, ES_E, ES_GE))) & !all(!is.null(c(ES_G, ES_E, ES_GE)))){
		stop(paste0("Specify all of ES_G, ES_E, and ES_GE (detectable effect sizes for gene, environment,", 
			"and gene-environment interaction)"))
	}

	if(any(!is.null(c(R2_G, R2_E, R2_GE))) & !all(!is.null(c(R2_G, R2_E, R2_GE)))){
		stop(paste0("Specify all of R2_G, R2_E, and R2_GE (detectable R-squared for gene, environment,", 
			"and gene-environment interaction)"))
	}	

	if(is.null(sd_y)){
		stop("sd_y, the standard deviation of the outcome in the overall population, must be specified.")
	}
	
	if(is.null(c(P_e))){
		stop("P_e (Environmental Factor Prevalence) must be specified.")
	}
	############################################################################################################
	#Error Messages for out of range values
	############################################################################################################
	if(sum(sd_y<=0)>0){
		stop("sd_y must be greater than 0.")
	}

	if(sum(P_e>=1)>0 | sum(P_e<=0)>0){
		stop("P_e must be greater than 0 and less than 1.")
	}



	if(sum(R2_G>=1)>0 | sum(R2_E<=0)>0 | sum(R2_GE>=1)>0){
		stop("R2_G, R2_E, and R2_GE must be greater than 0 and less than 1.")
	}

	if(sum(sd_y<=0)>0){
		stop("sd_y must be greater than 0.")
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
	# Calculate variances to go between R2 and ES (genotype) to do effect size calculations
	############################################################################################################
	#Getting R squared  and "bar" values

	var_x_dom = (1^2)*(1-(1-MAF)^2)-(1*(1-(1-MAF)^2))^2
	var_x_add = (1^2)*(2*MAF*(1-MAF))+(2^2)*(MAF^2)-(1*(2*MAF*(1-MAF))+2*(MAF^2))^2
	var_x_rec = (1^2)*(MAF^2)-(1*(MAF^2))^2

	var_x <- data.frame(True.Model=c(rep('Dominant', length(var_x_dom)),
									rep('Additive', length(var_x_add)),
									rep('Recessive', length(var_x_rec))),
						MAF = rep(MAF, 3),
						var_x = c(var_x_dom, var_x_add, var_x_rec)
						)
	mu_g_dom = 1 - (1 - MAF)^2
	mu_g_add = 2*(1 - MAF) * MAF + 2 * MAF^2
	mu_g_rec = MAF^2
	mu_g <- data.frame(True.Model=c(rep('Dominant', length(MAF)),
									rep('Additive', length(MAF)),
									rep('Recessive', length(MAF))),
						MAF = rep(MAF, 3),
						mu_g = c(mu_g_dom, mu_g_add, mu_g_rec)
						)
	# mu_g <- list(Recessive = MAF^2, Dominant = 1 - (1 - MAF)^2, Additive = 2*(1 - MAF) * MAF + 2 * MAF^2)
	# var_g <- list(Recessive = MAF^2 * (2* MAF * (1 - MAF) + (1 - MAF)^2), 
	# 			Dominant = (1 - MAF)^2 * (2* MAF * (1 - MAF) + MAF^2), 
	# 			Additive = 2*MAF - 2 * MAF^2)

	############################################################################################################
	# Create a data.frame with all possible combinations of MAF, SD and effect size
	# Will calculate power for each of these scenarios
	############################################################################################################
	# Depending on if ES or R2 is calculated, calcuate the other effect size measures
	
	# betaG_bar <- ES_G + ES_GE * (P_e)
	# betaE_bar <- ES_E + ES_GE * (mu_g)

	# R2_G <- betaG_bar^2 * var_g / (sd_y)^2
	# R2_E <- betaE_bar^2 * P_e * (1 - P_e) / (sd_y)^2
	# R2_GE <- ES_GE^2 * (var_g * P_e * (1 - P_e)) / sd_y^2
	



	if(all(is.null(c(ES_G, ES_E, ES_GE)))){
		e.save.tab = expand.grid(True.Model, MAF, P_e, sd_y, R2_G, R2_E, R2_GE, stringsAsFactors = F)
		colnames(e.save.tab) <- c("True.Model", "MAF", "P_e", "sd_y", "R2_G", "R2_E", "R2_GE")
		e.save.tab <- merge(e.save.tab, var_x)
		e.save.tab <- merge(e.save.tab, mu_g)

		e.save.tab$ES_G_bar <- sqrt(e.save.tab$R2_G * e.save.tab$sd_y^2 / (e.save.tab$var_x))
		e.save.tab$ES_E_bar <- sqrt(e.save.tab$R2_E * e.save.tab$sd_y^2 / (e.save.tab$P_e * (1 - e.save.tab$P_e)))
		e.save.tab$ES_GE <- sqrt(e.save.tab$R2_GE * e.save.tab$sd_y^2 / (e.save.tab$var_x * P_e * (1 - P_e)))

		e.save.tab$ES_G <- e.save.tab$ES_G_bar - e.save.tab$ES_GE * e.save.tab$P_e
		e.save.tab$ES_E <- e.save.tab$ES_E_bar - e.save.tab$ES_GE * mu_g[mu_g$True.Model == e.save.tab$True.Model, "mu_g"]

		e.save.tab <- e.save.tab[,c("True.Model", "MAF", "P_e", "sd_y", "var_x", "ES_G", "ES_E", "ES_GE",
			"ES_G_bar", "ES_E_bar", "R2_G", "R2_E", "R2_GE")]
	}

	if(all(is.null(c(R2_G, R2_E, R2_GE)))){
		e.save.tab = expand.grid(True.Model, MAF, P_e, sd_y, ES_G, ES_E, ES_GE, stringsAsFactors = F)
		colnames(e.save.tab) <- c("True.Model", "MAF", "P_e", "sd_y", "ES_G", "ES_E", "ES_GE")
		e.save.tab <- merge(e.save.tab, var_x)
		e.save.tab <- merge(e.save.tab, mu_g)

		e.save.tab$ES_G_bar <- e.save.tab$ES_G + e.save.tab$ES_GE * e.save.tab$P_e
		e.save.tab$ES_E_bar <- e.save.tab$ES_E + e.save.tab$ES_GE * e.save.tab$mu_g

		e.save.tab$R2_G <- e.save.tab$ES_G_bar^2 * e.save.tab$var_x / e.save.tab$sd_y^2
		e.save.tab$R2_E <- e.save.tab$ES_E_bar^2 * e.save.tab$P_e * (1 - e.save.tab$P_e) / e.save.tab$sd_y^2
		e.save.tab$R2_GE <- e.save.tab$ES_GE^2 * e.save.tab$var_x * P_e * (1 - P_e) / e.save.tab$sd_y^2

		e.save.tab <- e.save.tab[,c("True.Model", "MAF", "P_e", "sd_y", "var_x", "ES_G", "ES_E", "ES_GE",
			"ES_G_bar", "ES_E_bar", "R2_G", "R2_E", "R2_GE")]


		if(any(apply(e.save.tab[, c("R2_G", "R2_E", "R2_GE")], 1, function(x) sum(x) >= 1))){
			excluded <- e.save.tab[apply(e.save.tab[, c("R2_G", "R2_E", "R2_GE")], 1, function(x) sum(x) >= 1), ]
			e.save.tab <- e.save.tab[!apply(e.save.tab[, c("R2_G", "R2_E", "R2_GE")], 1, function(x) sum(x) >= 1),]

			message("Some combinations of the specified ES and sd_y imply R2>1 and 0 or negative variance of y within a genotype.\n
				Power was not calculated for the following combinations: \n", paste(capture.output(print(excluded)), collapse = "\n"))
		}
		if(nrow(e.save.tab)==0){
			stop("All combinations of the specified ES and sd_y imply R2>1 and 0 or negative variance of y within a genotype.\n
				Power could not be calculated. Try using smaller ES and/or larger sd_y.")
		}
	}

	# For each scenario calculate the SD of Y give X for the true model
	e.save.tab$sd_y_x_true = mapply(function(x){ # MAF=    P_e=    ES_G=   ES_E=   ES_GE=
		linear.outcome.log.envir.interaction.sds(MAF = e.save.tab[x,'MAF'], P_e = e.save.tab[x,'P_e'], 
			sd_y = e.save.tab[x,'sd_y'], ES_G = e.save.tab[x,'ES_G'], ES_E = e.save.tab[x,'ES_E'], ES_GE = e.save.tab[x,'ES_GE'], 
			True.Model = e.save.tab[x, "True.Model"], mod = e.save.tab[x,"True.Model"])}, seq(1:nrow(e.save.tab)))
	# sd for no interaction for reduced model
	e.save.tab$sd_y_x_true_0int = mapply(function(x){ # MAF=    P_e=    ES_G=   ES_E=   ES_GE=
		linear.outcome.log.envir.interaction.sds(MAF = e.save.tab[x,'MAF'], P_e = e.save.tab[x,'P_e'], 
			sd_y = e.save.tab[x,'sd_y'], ES_G = e.save.tab[x,'ES_G_bar'], ES_E = e.save.tab[x,'ES_E_bar'], ES_GE = 0, 
			True.Model = e.save.tab[x, "True.Model"], mod = e.save.tab[x,"True.Model"])}, seq(1:nrow(e.save.tab)))


	############################################################################################################
	# Calculate Power for each scenario in e.save.tab under the specified testing model
	############################################################################################################
	power.tab <- NULL

	############################################################################################################
	#Loop over sample size
	############################################################################################################
	for (n in N){
		################################################################################################
		#Loop over all of the testing models and calculate power for each ES, SD, and MAF scenario
		################################################################################################
		for (mod in Test.Model){
			# Calculate SD of Y given X for each scenario, given the test model
			sd_y_x <- mapply(function(x){
				linear.outcome.log.envir.interaction.sds(MAF = e.save.tab[x,"MAF"], 
					sd_y = e.save.tab[x, "sd_y"], P_e = e.save.tab[x,"P_e"], ES_G = e.save.tab[x,"ES_G"], 
					ES_E = e.save.tab[x,"ES_E"], ES_GE = e.save.tab[x,"ES_GE"], mod = mod, True.Model = e.save.tab[x, "True.Model"])
				}, seq(1:nrow(e.save.tab)))
			# Calculate SD of Y given X for each scenario, given the test model with no effect size for the reduced model
			sd_y_x_0int <- mapply(function(x){linear.outcome.log.envir.interaction.sds(MAF = e.save.tab[x,'MAF'], 
				sd_y = e.save.tab[x, "sd_y"], P_e = e.save.tab[x,'P_e'], ES_G = e.save.tab[x,'ES_G_bar'], 
				ES_E = e.save.tab[x,'ES_E_bar'], ES_GE = 0, mod = mod, True.Model = e.save.tab[x, "True.Model"])}, seq(1:nrow(e.save.tab)))
			sd_y_x_0int_new <- mapply(function(x){linear.outcome.log.envir.interaction.sds(MAF = e.save.tab[x,'MAF'], 
				sd_y = e.save.tab[x, "sd_y"], P_e = e.save.tab[x,'P_e'], ES_G = e.save.tab[x,'ES_G'], reduced = T,
				ES_E = e.save.tab[x,'ES_E'], ES_GE = e.save.tab[x, "ES_GE"], mod = mod, True.Model = e.save.tab[x, "True.Model"])}, seq(1:nrow(e.save.tab)))
			if(!all(sapply(sd_y_x_0int-sd_y_x_0int_new, all.equal, 0))) stop("sd's don't match")

			ll.alt <- mapply(function(x){calc.like.linear.log.envir.interaction(
				linear.mles.log.envir.interaction(MAF = e.save.tab[x,"MAF"], P_e = e.save.tab[x, "P_e"], ES_G = e.save.tab[x,"ES_G"], 
					ES_E = e.save.tab[x,"ES_E"], ES_GE = e.save.tab[x,"ES_GE"], Test.Model = mod, True.Model = e.save.tab[x, "True.Model"]),
									MAF = e.save.tab[x,"MAF"],
									P_e = e.save.tab[x, "P_e"],
									ES_G = e.save.tab[x,"ES_G"],
									ES_E = e.save.tab[x,"ES_E"],
									ES_GE = e.save.tab[x,"ES_GE"],
									sd_y_x_truth = e.save.tab[x, "sd_y_x_true"],
									sd_y_x_model = sd_y_x[x],
									True.Model = e.save.tab[x, "True.Model"],
									Test.Model=mod)}, seq(1:nrow(e.save.tab)))
			
			# reduced is the same calculation, except with ES_GE equal to 0
			ll.reduced <- mapply(function(x){calc.like.linear.log.envir.interaction(
				linear.mles.log.envir.interaction(MAF = e.save.tab[x,"MAF"], P_e = e.save.tab[x, "P_e"], ES_G = e.save.tab[x,"ES_G_bar"], 
					ES_E = e.save.tab[x,"ES_E_bar"], ES_GE = 0, Test.Model = mod, True.Model = e.save.tab[x, "True.Model"]),
									MAF = e.save.tab[x,"MAF"],
									P_e = e.save.tab[x, "P_e"],
									ES_G = e.save.tab[x,"ES_G_bar"],
									ES_E = e.save.tab[x,"ES_E_bar"],
									ES_GE = 0,
									sd_y_x_truth = e.save.tab[x, "sd_y_x_true_0int"],
									sd_y_x_model = sd_y_x_0int[x],
									True.Model = e.save.tab[x, "True.Model"],
									Test.Model=mod)}, seq(1:nrow(e.save.tab)))
			ll.reduced_new <- mapply(function(x){calc.like.linear.log.envir.interaction(
				linear.mles.log.envir.interaction(MAF = e.save.tab[x,"MAF"], P_e = e.save.tab[x, "P_e"], ES_G = e.save.tab[x,"ES_G"], 
					ES_E = e.save.tab[x,"ES_E"], ES_GE = e.save.tab[x, "ES_GE"], Test.Model = mod, True.Model = e.save.tab[x, "True.Model"], reduced = T),
									reduced = T,
									MAF = e.save.tab[x,"MAF"],
									P_e = e.save.tab[x, "P_e"],
									ES_G = e.save.tab[x,"ES_G"],
									ES_E = e.save.tab[x,"ES_E"],
									ES_GE = e.save.tab[x, "ES_GE"],
									sd_y_x_truth = e.save.tab[x, "sd_y_x_true"],
									sd_y_x_model = sd_y_x_0int[x],
									True.Model = e.save.tab[x, "True.Model"],
									Test.Model=mod)}, seq(1:nrow(e.save.tab)))
			if(!all(sapply(ll.reduced_new-ll.reduced, all.equal, 0))) stop("ll's don't match")

			ll.stat = 2*(ll.alt-ll.reduced)

			#Calculate the power for the given sample size for a range of Alpha levels
			if(mod=='2df'){
				pow = mapply(function(stat) 1-pchisq(qchisq(1-Alpha, df=2, ncp=0), df=2, ncp = n*stat), ll.stat)
			}else{
				pow = mapply(function(stat) pnorm(sqrt(n*stat) - qnorm(1-Alpha/2))+pnorm(-sqrt(n*stat) - qnorm(1-Alpha/2))*1, ll.stat)
			}
			if(length(Alpha)>1){
				pow <- t(pow)
				rownames(pow) <- seq(1:nrow(pow))
			}

			#Save the power calculations for each testing model in a final table for the sample size and case rate
			power.tab<-rbind(power.tab,data.frame(Test.Model=mod, True.Model = as.character(e.save.tab[, "True.Model"]),
					MAF = e.save.tab[, "MAF"], N = n, P_e = e.save.tab[, "P_e"], ES_G = e.save.tab[, "ES_G"], 
					ES_E = e.save.tab[, "ES_E"], ES_GE = e.save.tab[, "ES_GE"], R2_G=e.save.tab[, "R2_G"], 
					R2_E=e.save.tab[, "R2_E"], R2_GE=e.save.tab[, "R2_GE"], SD=e.save.tab[, "sd_y"], pow),row.names = NULL)

		}
	}
	colnames(power.tab)<-c("Test.Model", "True.Model", "MAF", "N_total", "P_e", "ES_G", "ES_E", "ES_GE", 
				"R2_G", "R2_E", "R2_GE", "SD_Y", paste("Power_at_Alpha_", Alpha, sep=''))

	return(power.tab)
}
