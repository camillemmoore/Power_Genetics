#' Function to output probability vector used in calculation of MLE's for linear outcome with linear environment interaction
#'
#' Returns probability vector used in calculation of MLE's for linear outcome with linear environment interaction
#'
#' @param mod type of model
#'
#' @return A probability vector to be used in MLE calculation for linear outcome with linear environment interaction
#'
#' @export
#'
p_vec_returner_lin_env <- function(MAF)
{
	p_AA0 <- (1 - MAF)^2
	p_AB0 <- 2*(1 - MAF)*MAF
	p_BB0 <- MAF^2
	p_AA1 <- (1 - MAF)^2
	p_AB1 <- 2*(1 - MAF)*MAF
	p_BB1 <- MAF^2
	p_vec <- c(p_AA0, p_AB0, p_BB0, p_AA1, p_AB1, p_BB1)
	return(p_vec)
}

X_mat_returner_lle <- function(mod)
{

	X_mat_dom <- rbind(
		c(1, 0, -1, 0),
		c(1, 1, -1, -1),
		c(1, 1, -1, -1),
		c(1, 0, 1, 0),
		c(1, 1, 1, 1),
		c(1, 1, 1, 1))

	X_mat_rec <- rbind(
		c(1, 0, -1, 0),
		c(1, 0, -1, 0),
		c(1, 1, -1, -1),
		c(1, 0, 1, 0),
		c(1, 0, 1, 0),
		c(1, 1, 1, 1))

	X_mat_add <- rbind(
		c(1, 0, -1, 0),
		c(1, 1, -1, -1),
		c(1, 2, -1, -2),
		c(1, 0, 1, 0),
		c(1, 1, 1, 1),
		c(1, 2, 1, 2))

	X_mat_2df <- rbind(
		c(1, 0, 0, -1, 0, 0),
		c(1, 1, 0, -1, -1, 0),
		c(1, 0, 1, -1, 0, -1),
		c(1, 0, 0, 1, 0, 0),
		c(1, 1, 0, 1, 1, 0),
		c(1, 0, 1, 1, 0, 1))

	X_mat_null <- rbind(
		c(1, 0, 0, 0),
		c(1, 0, 0, 0),
		c(1, 0, 0, 0),
		c(1, 0, 0, 0),
		c(1, 0, 0, 0),
		c(1, 0, 0, 0))

	X_mat_list <- list("Dominant" = X_mat_dom, "Recessive" = X_mat_rec, "Additive" = X_mat_add, "2df" = X_mat_2df, "null" = X_mat_null)

	return(X_mat_list[[mod]])
}


#' Function to calculate the standard deviation of y given x for linear models with linear environment interaction
#'
#' Returns the standard deviation of y given x for linear models with linear environment interaction
#'
#' @param MAF Minor allele Frequency
#' @param beta0 baseline value for the outcome
#' @param ES_G Genetic Effect size
#' @param ES_E Environment Effect size
#' @param ES_GE Environment x Genetic interaction Effect size
#' @param True.Model True Model
#' @param Test.Model Test Model
#'
#' @return The standard deviation of y given x for linear models with linear environment interaction
#'
#' @export
#'
linear.mles.lin.envir.interaction <- function(MAF, beta0, ES_G, ES_E, ES_GE, Test.Model, True.Model)
{

	# create elements we will use to calculate MLE
	p_vec <- p_vec_returner_lin_env(MAF)
	W_mat <- diag(p_vec)
	X_matF <- X_mat_returner_lle(Test.Model) # the MLE is calculated using the test model

	# the effect sizes are calculated from the true model
	ES_vec <- (X_mat_returner_lle(True.Model) %*% c(beta0, ES_G, ES_E, ES_GE))[,1]

	beta_hat <- as.vector(ginv(crossprod(X_matF, W_mat) %*% X_matF) %*% crossprod(X_matF, W_mat) %*% ES_vec)
	return(beta_hat)
}

#' Function to calculate the standard deviation of y given x for linear models with linear environment interaction 
#' for the reduced model without GxE interaction
#'
#' Returns the standard deviation of y given x for linear models with linear environment interaction
#'
#' @param MAF Minor allele Frequency
#' @param beta0 baseline value for the outcome
#' @param ES_G Genetic Effect size
#' @param ES_E Environment Effect size
#' @param ES_GE Environment x Genetic interaction Effect size
#' @param True.Model True Model
#' @param Test.Model Test Model
#'
#' @return The standard deviation of y given x for linear models with linear environment interaction
#'
#' @export
#'
linear.mles.lin.envir.interaction_reduced <- function(MAF, beta0, ES_G, ES_E, ES_GE, Test.Model, True.Model)
{
  
  # create elements we will use to calculate MLE
  p_vec <- p_vec_returner_lin_env(MAF)
  W_mat <- diag(p_vec)
  X_matF <- X_mat_returner_lle(Test.Model) # the MLE is calculated using the test model
  if(Test.Model=='2df'){X_matF<-X_matF[,1:4]}else{X_matF<-X_matF[,1:3]}
  
  # the effect sizes are calculated from the true model
  ES_vec <- (X_mat_returner_lle(True.Model) %*% c(beta0, ES_G, ES_E, ES_GE))[,1]
  
  beta_hat <- as.vector(ginv(crossprod(X_matF, W_mat) %*% X_matF) %*% crossprod(X_matF, W_mat) %*% ES_vec)
  return(beta_hat)
}



#' Function to calculate the standard deviation of y given x for linear models with linear environment interaction
#'
#' Returns the standard deviation of y given x for linear models with linear environment interaction
#'
#' @param MAF Minor allele Frequency
#' @param sd_e Standard deviation of linear environmental factor
#' @param beta0 baseline value for the outcome
#' @param ES_G Genetic Effect size
#' @param ES_E Environment Effect size
#' @param ES_GE Environment x Genetic interaction Effect size
#' @param mod Test model
#' @param True.Model True model
#' @param sd_y Standard deviation of y
#'
#' @return The standard deviation of y given x for linear models with linear environment interaction
#'
#' @export
#'
linear.outcome.lin.envir.interaction.sds <- function(MAF, sd_e, beta0, ES_G, ES_E, ES_GE, mod, True.Model, sd_y)
{

	p_vec <- p_vec_returner_lin_env(MAF)

	ES_test <- linear.mles.lin.envir.interaction(MAF = MAF, ES_G = ES_G, beta0 = beta0,
					ES_E = ES_E, ES_GE = ES_GE, Test.Model = mod, True.Model = True.Model)


	ES_vec <- (X_mat_returner_lle(mod) %*% ES_test)[,1]
	ES_vec_truemodel  <- (X_mat_returner_lle(True.Model) %*% c(beta0, ES_G, ES_E, ES_GE))[,1]

	gvec_true <- X_mat_returner_lle(True.Model)[1:3,2]
	gvec_test <- X_mat_returner_lle(mod)[1:3,2]
	if(mod == "2df") gvec_test <- X_mat_returner_lle("Dominant")[1:3,2]
	prob_vec <- c((1-MAF)^2, 2*(1-MAF)*MAF, MAF^2)

	# mod indices change if mod is 2df to select the proper genotype scenario since the coefficients for AB are independent of those for BB
	if(mod == "2df") mod_indices <- list(c(1,2,4,5), c(1,2,4,5), c(1,3,4,6)) else mod_indices <- list(c(1,2,3,4), c(1,2,3,4), c(1,2,3,4))

	sd_y_x <- sqrt(sd_y^2 + sum(sapply(1:3, function(x){
				-2 * prob_vec[x] * ((beta0 + ES_G*gvec_true[x])*(ES_test[mod_indices[[x]][1]] + ES_test[mod_indices[[x]][2]]*gvec_test[x]) + 
					(ES_E + ES_GE*gvec_true[x])*(ES_test[mod_indices[[x]][3]] + ES_test[mod_indices[[x]][4]]*gvec_test[x])*sd_e^2) + 
				prob_vec[x] * ((ES_test[mod_indices[[x]][1]] + ES_test[mod_indices[[x]][2]]*gvec_test[x])^2 + (ES_test[mod_indices[[x]][3]] + 
					ES_test[mod_indices[[x]][4]]*gvec_test[x])^2*sd_e^2)
				})))
	# this old way works if you aren't using the "mod_indices", which are necessary for 2df
	# sd_y_x <- sqrt(sd_y^2 + sum(sapply(1:3, function(x){
	# 			-2 * prob_vec[x] * ((beta0 + ES_G*gvec_true[x])*(ES_test[1] + ES_test[2]*gvec_test[x]) + 
	# 				(ES_E + ES_GE*gvec_true[x])*(ES_test[3] + ES_test[4]*gvec_test[x])*sd_e^2) + 
	# 			prob_vec[x] * ((ES_test[1] + ES_test[2]*gvec_test[x])^2 + (ES_test[3] + ES_test[4]*gvec_test[x])^2*sd_e^2)
	# 			})))

	return(sd_y_x)
}

#' Function to calculate the standard deviation of y given x for linear models with linear environment interaction
#'
#' Returns the standard deviation of y given x for linear models with linear environment interaction
#'
#' @param MAF Minor allele Frequency
#' @param sd_e Standard deviation of linear environmental factor
#' @param beta0 baseline value for the outcome
#' @param ES_G Genetic Effect size
#' @param ES_E Environment Effect size
#' @param ES_GE Environment x Genetic interaction Effect size
#' @param mod Test model
#' @param True.Model True model
#' @param sd_y Standard deviation of y
#'
#' @return The standard deviation of y given x for linear models with linear environment interaction
#'
#' @export
#'
linear.outcome.lin.envir.interaction.sds_reduced <- function(MAF, sd_e, beta0, ES_G, ES_E, ES_GE, mod, True.Model, sd_y)
{
  
  p_vec <- p_vec_returner_lin_env(MAF)
  
  ES_test <- linear.mles.lin.envir.interaction_reduced(MAF = MAF, ES_G = ES_G, beta0 = beta0,
                                                       ES_E = ES_E, ES_GE = ES_GE, Test.Model = mod, True.Model = True.Model)
  
  if(mod=='2df'){X_mat <- X_mat_returner_lle(mod)[,1:4]}else{X_mat <- X_mat_returner_lle(mod)[,1:3]}
  
  ES_vec <- (X_mat %*% ES_test)[,1]
  ES_vec_truemodel  <- (X_mat_returner_lle(True.Model) %*% c(beta0, ES_G, ES_E, ES_GE))[,1]
  
  gvec_true <- X_mat_returner_lle(True.Model)[1:3,2]
  gvec_test <- X_mat_returner_lle(mod)[1:3,2]
  if(mod == "2df") gvec_test <- X_mat_returner_lle("Dominant")[1:3,2] 
  prob_vec <- c((1-MAF)^2, 2*(1-MAF)*MAF, MAF^2)
  
  # mod indices change if mod is 2df to select the proper genotype scenario since the coefficients for AB are independent of those for BB
  if(mod == "2df") mod_indices <- list(c(1,2,4,5), c(1,2,4,5), c(1,3,4,6)) else mod_indices <- list(c(1,2,3,4), c(1,2,3,4), c(1,2,3,4))
  
  sd_y_x <- sqrt(sd_y^2 + sum(sapply(1:3, function(x){
    -2 * prob_vec[x] * ((beta0 + ES_G*gvec_true[x])*(ES_test[mod_indices[[x]][1]] + ES_test[mod_indices[[x]][2]]*gvec_test[x]) + 
                          (ES_E + ES_GE*gvec_true[x])*(ES_test[mod_indices[[x]][3]] + 0*gvec_test[x])*sd_e^2) + 
      prob_vec[x] * ((ES_test[mod_indices[[x]][1]] + ES_test[mod_indices[[x]][2]]*gvec_test[x])^2 + (ES_test[mod_indices[[x]][3]] + 
                                                                                                       0*gvec_test[x])^2*sd_e^2)
  })))
  
  return(sd_y_x)
}


#' Function to calculate the standard deviation of y given x for linear models with linear environment interaction
#'
#' Returns the standard deviation of y given x for linear models with linear environment interaction
#'
#' @param beta_hat Effect sizes from MLE
#' @param MAF Minor allele Frequency
#' @param sd_e Standard deviation of linear environmental factor
#' @param ES_G Genetic Effect size
#' @param ES_E Environment Effect size
#' @param ES_GE Environment x Genetic interaction Effect size
#' @param sd_y_x_truth Standard deviation of y for the true model
#' @param sd_y_x_test Standard deviation of y for the test model
#' @param sd_y Standard deviation of y
#' @param Test.model Test model
#' @param True.Model True model
#'
#' @return The standard deviation of y given x for linear models with linear environment interaction
#'
#' @export
#'
calc.like.linear.lin.envir.interaction <- function(beta_hat, MAF, sd_e, ES_G, ES_E, ES_GE, beta0,
	sd_y_x_truth, sd_y_x_model, Test.Model, True.Model)
{
	p_vec <- p_vec_returner_lin_env(MAF)

	betavec <- (X_mat_returner(Test.Model) %*% beta_hat)[,1]
	ES_vec <- (X_mat_returner(True.Model) %*% c(beta0, ES_G, ES_E, ES_GE))[,1]

	# want to have info on beta0 + betaG * G as well as betaE + betaGE * G, where G denotes the genotype coding based on the model type

	betamat <- cbind(betavec[1:3], betavec[4:6] - betavec[1:3])
	ES_mat <- cbind(ES_vec[1:3], ES_vec[4:6] - ES_vec[1:3])

	ll <- sum(p_vec[1:3] * apply(cbind(ES_mat, betamat), 1, function(x) 
		expected.linear.ll.lin.env_old(mean_truth = x[1:2], mean_model = x[3:4], sd_y_x_truth, sd_y_x_model, sd_e = sd_e)))
		# expected.linear.ll(mean_truth = x[1:2], mean_model = x[3:4], sd_y_x_truth, sd_y_x_model)))

	return(ll)
}

#' Function to Calculate Expected Log Likelihood for a Single Genotype with linear environment interaction
#'
#' Calculates the expected log likelihood for a single genotype with linear environment interaction 
#'  given the true and estimated mean and standard deviation for the outcome.
#'
#' @param sd_y_x_model The standard deviation of Y (the outcome) given X (predictors/genotype) under the test model.
#'
#' @return The log likelihood.
#'
#' @export
#'
expected.linear.ll.lin.env <- function(sd_y_x_model)
{
  -0.5*log(2*pi*sd_y_x_model^2) - 0.5
}
expected.linear.ll.lin.env_old <- function(mean_truth, mean_model, sd_y_x_truth, sd_y_x_model, sd_e)
{
	-0.5*log(2*pi*sd_y_x_model^2) -
		1/(2*sd_y_x_model^2) * (sd_y_x_truth^2 + (mean_truth[1] - mean_model[1])^2+
		sd_e^2 * (mean_truth[2] - mean_model[2])^2)
}
#### >> @param mean_truth Mean of the outcome given X(predictors/genotype) under the true model.
#### >> @param mean_model Mean of the outcome given X(predictors/genotype) under the test model.
#### >> @param sd_y_x_truth The standard deviation of Y given X (predictors/genotype) given genotype under the true model.
#### >> @param sd_e Standard deviation of linear environmental factor



# linear.outcome.lin.envir.interaction.sds.2 <- function(MAF, sd_e, ES_G, ES_E, ES_GE, sd_y, beta0, True.Model){
# 	P_AA <- (1-MAF)^2
# 	P_AB <- 2*MAF*(1-MAF)
# 	P_BB <- MAF^2
	
# 	sd_x_y_dom = sqrt(sd_y^2 - (P_AA*(beta0^2 + ES_E^2*sd_e^2) + (P_AB + P_BB)*((beta0 + ES_G)^2 + (ES_E + ES_GE)^2*sd_e^2)))
# 	sd_x_y_rec = sqrt(sd_y^2 - ((P_AA + P_AB)*(beta0^2 + ES_E^2*sd_e^2) + P_BB*((beta0 + ES_G)^2 + (ES_E + ES_GE)^2*sd_e^2)))
# 	sd_x_y_add = sqrt(sd_y^2 - (P_AA*(beta0^2 + ES_E^2*sd_e^2) + P_AB*((beta0 + ES_G)^2 + (ES_E + ES_GE)^2*sd_e^2) + P_BB*((beta0 + 2*ES_G)^2 + (ES_E + 2*ES_GE)^2*sd_e^2)))

# 	sd_x_y <- data.frame(True.Model=c(rep('Dominant', length(sd_x_y_dom)),
# 									rep('Additive', length(sd_x_y_add)),
# 									rep('Recessive', length(sd_x_y_rec))),
# 						MAF = rep(MAF, 3),
# 						sd_x_y = c(sd_x_y_dom, sd_x_y_add, sd_x_y_rec)
# 						)
# 	return(sd_x_y[sd_x_y$True.Model == True.Model, "sd_x_y"])
# }
#

# linear.outcome.lin.envir.interaction.sds <- function(MAF, sd_e, ES_G, ES_E, ES_GE, sd_y, beta0, True.Model){
# 	P_AA <- (1-MAF)^2
# 	P_AB <- 2*MAF*(1-MAF)
# 	P_BB <- MAF^2
	
# 	sd_x_y_dom = sqrt(sd_y^2 - (P_AA*(beta0^2 + ES_E^2*sd_e^2) + (P_AB + P_BB)*((beta0 + ES_G)^2 + (ES_E + ES_GE)^2*sd_e^2)))
# 	sd_x_y_rec = sqrt(sd_y^2 - ((P_AA + P_AB)*(beta0^2 + ES_E^2*sd_e^2) + P_BB*((beta0 + ES_G)^2 + (ES_E + ES_GE)^2*sd_e^2)))
# 	sd_x_y_add = sqrt(sd_y^2 - (P_AA*(beta0^2 + ES_E^2*sd_e^2) + P_AB*((beta0 + ES_G)^2 + (ES_E + ES_GE)^2*sd_e^2) + P_BB*((beta0 + 2*ES_G)^2 + (ES_E + 2*ES_GE)^2*sd_e^2)))

# 	sd_x_y <- data.frame(True.Model=c(rep('Dominant', length(sd_x_y_dom)),
# 									rep('Additive', length(sd_x_y_add)),
# 									rep('Recessive', length(sd_x_y_rec))),
# 						MAF = rep(MAF, 3),
# 						sd_x_y = c(sd_x_y_dom, sd_x_y_add, sd_x_y_rec)
# 						)
# 	return(sd_x_y[sd_x_y$True.Model == True.Model, "sd_x_y"])
# }