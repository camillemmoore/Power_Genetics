#' Function to Calculate Power for Logistic Models with Environment Interaction
#'
#' Calculates the power to detect an difference in means/effect size/regression coefficient, at a given sample size, N, with type 1 error rate, Alpha
#'
#' @param power Vector of the desired power(s)
#' @param Case.Rate proportion of cases in the sample (cases/(cases + controls)). 
#' @param k Vector of the number of controls per case. Either k or Case.Rate must be specified.
#' @param MAF Vector of minor allele frequencies
#' @param OR_G Vector of genetic odds ratios to detect
#' @param OR_E Vector of environmental odds ratios to detect
#' @param OR_GE Vector of genetic/environmental interaction odds ratios to detect
#' @param P_e Vector of proportions of the population with exposure to the environmental effect
#' @param Alpha the desired type 1 error rate(s)
#' @param True.Model A vector specifying the true underlying genetic model(s): 'Dominant', 'Additive1', 'Additive2', 'Recessive' or 'All'
#' @param Test.Model A vector specifying the assumed genetic model(s) used in testing: 'Dominant', 'Additive', 'Recessive' or 'All'
#'
#' @return A data frame including the power for all combinations of the specified parameters (Case.Rate, ES, Power, etc)
#'
#' @examples
#' ssc <- ss_envir.calc(P_e = 0.2, MAF = 0.1, power = 0.6, Case.Rate = 0.5, Alpha = 0.05, 
#'		OR_G = 1.5, OR_E = 2, OR_GE = 1.8, Test.Model = "All", True.Model = "All")
#'
#'
#' @export
#'
ss_envir.calc <- 
	function(power=0.8, Case.Rate=NULL, k=NULL, MAF=NULL, OR_G=NULL, OR_E=NULL, OR_GE=NULL, P_e = NULL,
					Alpha=0.05, True.Model='All', Test.Model='All')
{

	############################################################################################################
	#Error Messages for insufficient sample size information, MAF, odds ratios, and case vs. control ratio
	############################################################################################################

	if(is.null(OR_E)){
		stop("OR_E, the odds ratio for the environmental exposure, must be specified.")
	}

	if(is.null(OR_G)){
		stop("OR_G, the odds ratio for the genetic mutation, must be specified.")
	}

	if(is.null(OR_GE)){
		stop("OR_GE, the odds ratio for the interaction between the genetic mutation and environmental exposure, must be specified.")
	}

	if(is.null(k) & is.null(Case.Rate)){
		stop("k, the number of controls per case, or Case.Rate, the proportion of cases in the study sample, must be specified.")
	}

	if(is.null(k) & is.null(Case.Rate)){
		stop("Specify one of k, the number of controls per case, or Case.Rate, the proportion of cases in the study sample, not both.")
	}

	if(is.null(MAF)){
		stop("MAF (minor allele frequency) must be specified.")
	}

	if(any(is.null(c(OR_G, OR_E, OR_GE)))){
		stop("OR (detectable odds ratio) must be specified.")
	}

	if(is.null(c(P_e))){
		stop("P_e (Environmental Factor Prevalence) must be specified.")
	}

	############################################################################################################
	#Error Messages for out of range values
	############################################################################################################

	if(sum(Case.Rate>=1)>0 | sum(Case.Rate<=0)>0){
		stop("R2 must be greater than 0 and less than 1.")
	}

	if(sum(MAF>=1)>0 | sum(MAF<=0)>0){
		stop("MAF must be greater than 0 and less than 1.")
	}

	if(sum(power>=1)>0 | sum(power<=0)>0){
    	stop("Power must be greater than 0 and less than 1.")
	}

	if(sum(k<=0)>0){
		stop("k must be greater than 0.")
	}

	if(sum(c(OR_E, OR_G, OR_GE)<=0)>0){
		stop("OR must be greater than 0.")
	}

	if(sum(Alpha>=1)>0 | sum(Alpha<=0)>0){
		stop("Alpha must be greater than 0 and less than 1.")
	}

	if(sum(!(Test.Model %in% c("Dominant", "Recessive", "Additive", "2df", "All")))>0){
		stop(paste("Invalid Test.Model:",
			   paste(Test.Model[!(Test.Model %in% c("Dominant", "Recessive", "Additive", "2df", "All"))], collapse=', ')))
	}

	if(sum(!(True.Model %in% c("Dominant", "Recessive", "Additive1", "Additive2", "Additive", "All")))>0){
		stop(paste("Invalid True.Model:",
			   paste(True.Model[!(True.Model %in% c("Dominant", "Recessive", "Additive1", "Additive2", "All"))], collapse=', ')))
	}

	############################################################################################################
	#Calculate needed sample size information from provided inputs
	############################################################################################################
	#Test model vector
	# if('All' %in% Test.Model){Test.Model<-c("Dominant", "Recessive", "Additive")}
	if('All' %in% Test.Model){Test.Model<-c("Dominant", "Recessive", "Additive", "2df")}

	#True model vector
	if('All' %in% True.Model){True.Model<-c("Dominant", "Recessive", "Additive")}

	#If k is provided calculate the Case.Rate
	if(is.null(Case.Rate)==T){Case.Rate = 1/(1+k)}

	#Find all possible combinations of N and case_rate
	power.tab <- expand.grid(power, Case.Rate)
	colnames(power.tab) <- c('power', 'Case.Rate')
	iter <- nrow(power.tab)
	
	final.ss.tab <- NULL

	for (zz in 1:iter){
		power <- power.tab[zz,"power"]
		Case.Rate <- power.tab[zz,"Case.Rate"]

		############################################################################################################
		#Use OR's and MAF's to calculate true distriubiton of genotypes and disease
		############################################################################################################

		##############################################################################
		#For each Environment Percentage
		##############################################################################

		pe.save.tab <- NULL
			
		for(P_e_ii in P_e){

			##############################################################################
			#For each OR
			##############################################################################
			o.save.tab <-NULL
			o.grid <- expand.grid(OR_G, OR_E, OR_GE)
			names(o.grid) <- c("OR_G", "OR_E", "OR_GE")
			for (o_ii in 1:nrow(o.grid)){

				##########################################################################
				#For each MAF
				##########################################################################

				m.save.tab<-NULL
				
				for (m in MAF){
					#Temporary place to save 2x3 tables for each true model with MAF=m and OR=o
					save.tab <- NULL

					#Proportion with each genotype. This is the same for all true.models.
					P_AA <- (1-m)^2
					P_AB <- 2*m*(1-m)
					P_BB <- m^2

					############################################################################
					#Create 2x6 tables of joint probabilities for each true model of interest
					############################################################################

					if("Dominant" %in% True.Model){
						#Dominant Model
						# m <<- m; P_e_ii<<-P_e_ii;Case.Rate<<-Case.Rate;o.grid<<-o.grid;o_ii<<-o_ii
						t <- dom.fun.t(MAF = m, P_e = P_e_ii, Case.Rate = Case.Rate,
							OR_G = o.grid$OR_G[o_ii], OR_E = o.grid$OR_E[o_ii], OR_GE = o.grid$OR_GE[o_ii])
						dom.tab <- data.frame(model = rep("Dominant", 2), table = t)
						save.tab<-rbind(save.tab, dom.tab)
					}
					
					if("Recessive" %in% True.Model){
						#Recessive Model
						# m <<- m; P_e_ii<<-P_e_ii;Case.Rate<<-Case.Rate;o.grid<<-o.grid;o_ii<<-o_ii
						t <- rec.fun.t(MAF = m, P_e = P_e_ii, Case.Rate = Case.Rate,
							OR_G = o.grid$OR_G[o_ii], OR_E = o.grid$OR_E[o_ii], OR_GE = o.grid$OR_GE[o_ii])
						rec.tab <- data.frame(model = rep("Recessive", 2), table = t)
						save.tab<-rbind(save.tab, rec.tab)
					}

					if("Additive" %in% True.Model){
						#Additive Model
						# m <<- m; P_e_ii<<-P_e_ii;Case.Rate<<-Case.Rate;o.grid<<-o.grid;o_ii<<-o_ii
						t <- add.fun.t(MAF = m, P_e = P_e_ii, Case.Rate = Case.Rate,
							OR_G = o.grid$OR_G[o_ii], OR_E = o.grid$OR_E[o_ii], OR_GE = o.grid$OR_GE[o_ii])
						add.tab <- data.frame(model = rep("Additive", 2), table = t)
						save.tab<-rbind(save.tab, add.tab)
					}

					m.save.tab<-rbind(m.save.tab,
						data.frame(True.Model = save.tab$model, MAF=m, P_e = P_e_ii, matrix(rep(as.numeric(o.grid[o_ii,]),nrow(save.tab)), ncol = 3, byrow = T),
							Disease.Status = rep(c("case", "control"),nrow(save.tab)/2),
							Geno.AA.e0 = save.tab[,2],Geno.AB.e0 = save.tab[,3], Geno.BB.e0 = save.tab[,4],
							Geno.AA.e = save.tab[,5],Geno.AB.e = save.tab[,6], Geno.BB.e = save.tab[,7]))

				}
				o.save.tab<-rbind(o.save.tab, m.save.tab)

			}
			pe.save.tab <- rbind(pe.save.tab, o.save.tab)
		}

		############################################################################################################
		#Calculate Power for each scenario under the specified testing model
		############################################################################################################


		################################################################################################
		#Loop over all of the testing models and calculate power for each OR and MAF scenario
		################################################################################################
		
		for (mod in Test.Model){
			temp<-NULL
			for(alpha0 in Alpha){
				temp.0<-NULL

				#Repeat calculation for each OR/MAF combination
				for (j in seq(1, nrow(pe.save.tab),2)){
					t<-pe.save.tab[j:(j+1),c("Geno.AA.e0", "Geno.AB.e0", "Geno.BB.e0", "Geno.AA.e", "Geno.AB.e", "Geno.BB.e")]

					mres <- ll.ge.logistic(t, power = power, Alpha = alpha0, mod = mod)
					# mres <<- mres
					# #Calculate the mreser for the given sample size for a range of Alpha levels
					# if(mod=='2df'){mres = 1-pchisq(qchisq(1-Alpha, df=2, ncp=0), df=2, ncp = N*stat)
					# }else{mres = pnorm(sqrt(N*stat) - qnorm(1-Alpha/2))+pnorm(-sqrt(N*stat) - qnorm(1-Alpha/2))*0}#

					temp.0<-rbind(temp.0, mres)
				}
				colnames(temp.0) <- paste0("SampleSize", gsub("pow", "", c("pow_GE")), "_at_Alpha_", alpha0)
				temp <- cbind(temp, temp.0)
				#Save the power calculations for each testing model in a final table for the sample size and case rate
			}
			ss.tab <- data.frame(Test.Model=mod, pe.save.tab[seq(1, nrow(pe.save.tab),2),1:6],
													Power = power, Case.Rate,temp)
			names(ss.tab)[5:7] <- c("OR_G", "OR_E", "OR_GE")
			# ss.tab<-rbind(ss.tab, ss.tab0)
			# colnames(ss.tab)[1:11]<-c("Test.Model", "True.Model", "MAF", "P_e", "OR_G", "OR_E", "OR_GE", "N", "N_cases", "N_controls", "Case.Rate")#, 
										# paste0("Power", gsub("pow", "", c("pow_GE", "pow_G", "pow_E")), "_at_Alpha_", Alpha))
			final.ss.tab<-rbind(final.ss.tab, ss.tab)
		}
	}
	final.ss.tab <- final.ss.tab[,c(1:4, 7, 5:6, 8:ncol(final.ss.tab))]
	return(final.ss.tab)

}
