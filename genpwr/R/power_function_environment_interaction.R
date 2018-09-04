
library(nleqslv)
source("~/drive/tf/camille/genpwr/Power_Genetics/genpwr/R/logit_function.R")
source("~/drive/tf/camille/genpwr/Power_Genetics/genpwr/R/zero_finder_nleqslv.R")
source("~/drive/tf/camille/genpwr/Power_Genetics/genpwr/R/power_function_environment_interaction_t_matrix_eqs.R")
source("~/drive/tf/camille/genpwr/Power_Genetics/genpwr/R/power_function_environment_interaction_logistic_mle.R")

# In Quanto, with 
# mode: dominant
# Environmental Prevalence:
#> P_e = 0.2
#> MAF = 0.1
#> P_AA <- P_AA <- (1 - MAF)^2; P_AB <- 2*(1-MAF)*MAF;P_BB <- MAF^2 # (P_AA = 0.81)
#> N = 200 #(100 each)
# disease percentage: 
#> Case.Rate = 0.5
#> Alpha <- alpha <- 0.05
#> OR_G = 1.5
#> OR_E = 2
#> OR_GE = 1.8
#> Test.Model = "Dominant";True.Model="Dominant"
#> k = NULL
#> compareQuanto = F
# results:
# interaction power: 0.0872
# gene power: 0.2673
# environment power: 0.5774

power_envir.calc <- 
	function(N=NULL, Case.Rate=NULL, k=NULL, MAF=NULL, OR_G=NULL, OR_E=NULL, OR_GE=NULL, P_e = NULL,
					Alpha=0.05, True.Model='All', Test.Model='All', compareQuanto = 0)
{
	if(is.logical(compareQuanto)) compareQuanto = 1 * compareQuanto

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

	if(is.null(N)){
		stop("N, the total sample size, must be specified.")}

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

	if(sum(N<=0)>0){
		stop("N must be greater than 0.")
	}

	if(sum(k<=0)>0){
		stop("k must be greater than 0.")
	}

	if(sum(c(OR_E, OR_G, OR_GE)<=0)>0){
		stop("OR_G, OR_E, and OR_GE must all be greater than 0.")
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
	sample.size.tab <- expand.grid(N, Case.Rate)
	colnames(sample.size.tab) <- c('N', 'Case.Rate')
	sample.size.tab$N_cases <- sample.size.tab$N*sample.size.tab$Case.Rate
	sample.size.tab$N_controls <- sample.size.tab$N - sample.size.tab$N_cases

	iter <- nrow(sample.size.tab)

	final.pow.tab <- NULL

	for (zz in 1:iter){
		N <- sample.size.tab[zz,'N']
		Case.Rate <- sample.size.tab[zz,'Case.Rate']
		N_cases <- sample.size.tab[zz,'N_cases']
		N_controls <- sample.size.tab[zz,'N_controls']

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
						data.frame(True.Model = save.tab$model, MAF=m, P_e = P_e_ii, matrix(rep(data.frame(o.grid[o_ii,]),nrow(save.tab)), ncol = 3, byrow = T),
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

					mres <- ll.ge.logistic(t, N = N, Alpha = alpha0, mod = mod, compareQuanto = compareQuanto)
					# mres <<- mres
					# #Calculate the mreser for the given sample size for a range of Alpha levels
					# if(mod=='2df'){mres = 1-pchisq(qchisq(1-Alpha, df=2, ncp=0), df=2, ncp = N*stat)
					# }else{mres = pnorm(sqrt(N*stat) - qnorm(1-Alpha/2))+pnorm(-sqrt(N*stat) - qnorm(1-Alpha/2))*0}#

					temp.0<-rbind(temp.0, mres)
				}
				colnames(temp.0) <- paste0("Power", gsub("pow", "", c("pow_GE", "pow_G", "pow_E")), "_at_Alpha_", alpha0)
				temp <- cbind(temp, temp.0)
				#Save the power calculations for each testing model in a final table for the sample size and case rate
			}
			power.tab <- data.frame(Test.Model=mod, pe.save.tab[seq(1, nrow(pe.save.tab),2),1:6],
													N, N_cases, N_controls, Case.Rate,temp)
			names(power.tab)[5:7] <- c("OR_G", "OR_E", "OR_GE")
			# power.tab<-rbind(power.tab, power.tab0)
			# colnames(power.tab)[1:11]<-c("Test.Model", "True.Model", "MAF", "P_e", "OR_G", "OR_E", "OR_GE", "N", "N_cases", "N_controls", "Case.Rate")#, 
										# paste0("Power", gsub("pow", "", c("pow_GE", "pow_G", "pow_E")), "_at_Alpha_", Alpha))
			final.pow.tab<-rbind(final.pow.tab, power.tab)
		}
	}
	final.pow.tab <- final.pow.tab[,c(1:4, 7, 5:6, 8:ncol(final.pow.tab))]
	return(final.pow.tab)

}



if(F){ #vestigial stuff to be deleted

	# q + a + b + c
	# r + w + x + y
	# all.equal(q + r, P_AA * (1 - P_e))
	# all.equal(a + w, P_AA * P_e)
	# all.equal(b + x, (1 - P_AA) * (1 - P_e))
	# all.equal(c + y, (1 - P_AA) * P_e)


	dominant.ll<-function(t)
	{

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

		beta <- optim(c(0,0,0,0), function(x) -ll_fun(x), control=c(abstol = 0.00001))$par
		beta_g_e <- optim(c(0,0,0), function(x) -ll_g_e_fun(x), control=c(abstol = 0.00001))$par
		beta_e <- optim(c(0,0), function(x) -ll_e_fun(x), control=c(abstol = 0.00001))$par
		beta_g <- optim(c(0,0), function(x) -ll_g_fun(x), control=c(abstol = 0.00001))$par

		ll <- ll_fun(beta)
		ll_g_e <- ll_g_e_fun(beta_g_e)
		ll_e <- ll_e_fun(beta_e)
		ll_g <- ll_g_fun(beta_g)

		Case.Rate <- sum(t[1,])	

		ll.null <- Case.Rate*log(exp(logit(sum(t[1,])/sum(t)))/(1+exp(logit(sum(t[1,])/sum(t))))) + (1-Case.Rate)*log(1/(1+exp(logit(sum(t[1,])/sum(t)))))
		# logit(sum(t[1,])/sum(t))
		pow_res <- numeric(3)
		pow_res[1] <- pnorm(sqrt(N*2*(ll-ll_g_e)) - qnorm(1-Alpha/2)) + pnorm(-sqrt(N*2*(ll-ll_g_e)) - qnorm(1-Alpha/2))*0
		pow_res[2] <- pnorm(sqrt(N*2*(ll_g-ll.null)) - qnorm(1-Alpha/2)) + pnorm(-sqrt(N*2*(ll_g-ll.null)) - qnorm(1-Alpha/2))*0
		# pnorm(sqrt(N*2*(ll_g_e-ll_e)) - qnorm(1-Alpha/2)) + pnorm(-sqrt(N*2*(ll_g_e-ll_e)) - qnorm(1-Alpha/2))*0
		pow_res[3] <- pnorm(sqrt(N*2*(ll_e-ll.null)) - qnorm(1-Alpha/2)) + pnorm(-sqrt(N*2*(ll_e-ll.null)) - qnorm(1-Alpha/2))*0
		return(pow_res)
	}

	pow <- dominant.ll(t)




	fun_ex <- function(x) {
		f <- numeric(3)
		f[1] <- sum(x) - 2 * x[1] + 1
		f[2] <- sum(x) + x[2] + 4
		f[3] <- sum(x) + 3 * x[3] - 2 
		f
	}

	x.start <- c(1,1,1)
	nleqslv(x.start, fun_ex)




	## logistic mles
	# NOTE: NOT SURE IF I NEED THIS
	logistic.mles<-function(t,model)
	{

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
		prob_BB_case_e <- t[1, 6] #generate dataset based on these probabilities and see if the betas match up
		Case.Rate <- sum(t[1,])

		if (model=='null'){#Null MOdel
			beta0 = logit(Case.Rate)
			ll<- Case.Rate*log(exp(beta0)/(1+exp(beta0))) + (1-Case.Rate)*log(1/(1+exp(beta0)))
		}

		#Dominant
		if (model=='Dominant'){
			beta0 <- logit(prob_AA_case_e0/(prob_AA_control_e0 + prob_AA_case_e0))
			betaG <- logit((prob_AB_case_e0+prob_BB_case_e0)/(prob_AB_control_e0+prob_BB_control_e0+prob_AB_case_e0+prob_BB_case_e0)) - beta0
			betaE <- logit(prob_AA_case_e/(prob_AA_control_e+prob_AA_case_e)) - beta0
			betaGE <- logit((prob_AB_case_e+prob_BB_case_e)/(prob_AB_control_e+prob_BB_control_e+prob_AB_case_e+prob_BB_case_e)) - beta0 - betaG - betaE
			beta <- c(beta0, betaG, betaE, betaGE)
		}

		#Recessive
		if (model=='Recessive'){
		}

		#2DF
		if (model=='2df'){
		}

		#Additive Model - need to use an optimizer to solve
		if (model=='Additive'){
		}
		return(beta)
	}
	# NOTE: NOT SURE IF I NEED THIS
	null.ll<-function(t){
		Case.Rate<-sum(t[1,])/sum(t)
		beta0 = logit(Case.Rate)
		ll<- Case.Rate*log(exp(beta0)/(1+exp(beta0))) + (1-Case.Rate)*log(1/(1+exp(beta0)))
		return(ll)
	}

	# a = (alpha omega epsilon - alpha epsilon + alpha - sqrt((-alpha omega epsilon + alpha epsilon - alpha + b omega - b + c omega - c + delta omega - delta - omega + 1)^2 - 4 (omega - 1) (-alpha delta omega epsilon + alpha omega epsilon + alpha (-b) omega epsilon - alpha c omega epsilon)) + b (-omega) + b - c omega + c - delta omega + delta + omega - 1)/(2 (omega - 1)) 

	# step 2:

	# (b ( -1 -alpha epsilon + alpha + delta - delta gamma + gamma) + (alpha - 1) (delta - 1) (epsilon - 1) gamma)/(b (-1 + gamma) - (-alpha epsilon + alpha + epsilon + 1) gamma) = (a ( -1 -alpha epsilon + alpha + delta - delta gamma + gamma) + alpha (delta - 1) epsilon gamma)/(a (-1 + gamma) - alpha epsilon gamma)

	if(F){
		# none of this is working
		# looking below, I think that works better
		# want to constrain x values to be greater than 0

		# constr_0_1 <- function(z) (atan(z) + pi/2) / pi
		constr_0_1 <- function(z) 1/(1 + exp(-z))
		# constr_0_1_inv <- function(y) -log(1/y - 1)
		fun2 <- function(z){
			x <- constr_0_1(z) # constrained to be greater than 1
			# x = constr_0_1(z) #function constrained to be between 0 and 1
			f <- fun(x)
			# f <- numeric(3)

			# f[1] <- (P_AA * (1 - P_e) - 1 + Case.Rate + sum(x)) * x[1]/((1 - Case.Rate - sum(x)) * (P_AA * P_e - x[1])) - OR_E
			# f[2] <- (P_AA * (1 - P_e) - 1 + Case.Rate + sum(x)) * x[2]/((1 - Case.Rate - sum(x)) * ((1 - P_AA) * (1 - P_e) - x[2])) - OR_G
			# f[3] <- (P_AA * (1 - P_e) - 1 + Case.Rate + sum(x)) * x[3]/((1 - Case.Rate - sum(x)) * ((1 - P_AA) * P_e - x[3])) - OR_GE#OR_E*OR_G*OR_GE
			# f
		}
		z.start <- c(-2.2,-2.2,-2.2)
		# solving gives you your z's:
		res_z <- nleqslv(z.start, fun2)#, control=list(trace=0,btol=.01,delta="cauchy"))
		res <- constr_0_1(res_z$x)
	}

}
