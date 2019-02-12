# MLE Functions for Normal/Continuous Outcomes with logistic environment interaction
#' Function to output X matrices used in calculation of MLE's for linear outcome with logistic environment interaction
#'
#' Returns X matrices used in calculation of MLE's for linear outcome with logistic environment interaction
#'
#' @param mod type of model
#'
#' @return A matrix to be used in MLE calculation for linear outcome with logistic environment interaction
#'
#' @export
#'
X_mat_returner <- function(mod, reduced = F)
{

	X_mat_dom <- rbind(
		c(1, 0, 0, 0),
		c(1, 1, 0, 0),
		c(1, 1, 0, 0),
		c(1, 0, 1, 0),
		c(1, 1, 1, 1),
		c(1, 1, 1, 1))

	X_mat_rec <- rbind(
		c(1, 0, 0, 0),
		c(1, 0, 0, 0),
		c(1, 1, 0, 0),
		c(1, 0, 1, 0),
		c(1, 0, 1, 0),
		c(1, 1, 1, 1))

	X_mat_add <- rbind(
		c(1, 0, 0, 0),
		c(1, 1, 0, 0),
		c(1, 2, 0, 0),
		c(1, 0, 1, 0),
		c(1, 1, 1, 1),
		c(1, 2, 1, 2))

	X_mat_2df <- rbind(
		c(1, 0, 0, 0, 0, 0),
		c(1, 1, 0, 0, 0, 0),
		c(1, 0, 1, 0, 0, 0),
		c(1, 0, 0, 1, 0, 0),
		c(1, 1, 0, 1, 1, 0),
		c(1, 0, 1, 1, 0, 1))
	if(reduced) X_mat_2df <- X_mat_2df[,1:(ncol(X_mat_2df) - 1)]

	X_mat_null <- rbind(
		c(1, 0, 0, 0),
		c(1, 0, 0, 0),
		c(1, 0, 0, 0),
		c(1, 0, 0, 0),
		c(1, 0, 0, 0),
		c(1, 0, 0, 0))

	X_mat_list <- list("Dominant" = X_mat_dom, "Recessive" = X_mat_rec, "Additive" = X_mat_add, "2df" = X_mat_2df, "null" = X_mat_null)

	if(reduced) X_mat_list <- lapply(X_mat_list, function(amat) amat <- amat[,1:(ncol(amat) - 1)])

	return(X_mat_list[[mod]])
}

#' Function to output probability vector used in calculation of MLE's for linear outcome with logistic environment interaction
#'
#' Returns probability vector used in calculation of MLE's for linear outcome with logistic environment interaction
#'
#' @param mod type of model
#' @param P_e Population prevalence of logistic environmental factor
#'
#' @return A probability vector to be used in MLE calculation for linear outcome with logistic environment interaction
#'
#' @export
#'
p_vec_returner <- function(MAF, P_e)
{
	p_AA0 <- (1 - MAF)^2 * (1 - P_e)
	p_AB0 <- 2*(1 - MAF)*MAF * (1 - P_e)
	p_BB0 <- MAF^2 *(1 - P_e)
	p_AA1 <- (1 - MAF)^2 * P_e
	p_AB1 <- 2*(1 - MAF)*MAF * P_e
	p_BB1 <- MAF^2 * P_e
	p_vec <- c(p_AA0, p_AB0, p_BB0, p_AA1, p_AB1, p_BB1)
	return(p_vec)
}

#' Function to calculate the standard deviation of y given x for linear models with logistic environment interaction
#'
#' Returns the standard deviation of y given x for linear models with logistic environment interaction
#'
#' @param MAF Minor allele Frequency
#' @param P_e Population prevalence of logistic environmental factor
#' @param ES_G Genetic Effect size
#' @param ES_E Environment Effect size
#' @param ES_GE Environment x Genetic interaction Effect size
#' @param mod model
#'
#' @return The standard deviation of y given x for linear models with logistic environment interaction
#'
#' @export
#'
linear.mles.log.envir.interaction <- function(MAF, P_e, ES_G, ES_E, ES_GE, Test.Model, True.Model, reduced = F)
{

	# create elements we will use to calculate MLE
	p_vec <- p_vec_returner(MAF, P_e)
	W_mat <- diag(p_vec)
	X_matF <- X_mat_returner(Test.Model, reduced = reduced) # the MLE is calculated using the test model

	# the effect sizes are calculated from the true model
	ES_vec <- (X_mat_returner(True.Model) %*% c(0, ES_G, ES_E, ES_GE))[,1]

	beta_hat <- as.vector(ginv(crossprod(X_matF, W_mat) %*% X_matF) %*% crossprod(X_matF, W_mat) %*% ES_vec)
	return(beta_hat)
}

#' Function to calculate the standard deviation of y given x for linear models with logistic environment interaction
#'
#' Returns the standard deviation of y given x for linear models with logistic environment interaction
#'
#' @param MAF Minor allele Frequency
#' @param P_e Population prevalence of logistic environmental factor
#' @param ES_G Genetic Effect size
#' @param ES_E Environment Effect size
#' @param ES_GE Environment x Genetic interaction Effect size
#' @param mod Test model
#' @param True.Model True model
#' @param sd_y Standard deviation of y
#'
#' @return The standard deviation of y given x for linear models with logistic environment interaction
#'
#' @export
#'
linear.outcome.log.envir.interaction.sds <- function(MAF, P_e, ES_G, ES_E, ES_GE, mod, True.Model, sd_y, reduced = F)
{

	p_vec <- p_vec_returner(MAF, P_e)

	ES_test <- linear.mles.log.envir.interaction(MAF = MAF, P_e = P_e, ES_G = ES_G, 
					ES_E = ES_E, ES_GE = ES_GE, Test.Model = mod, True.Model = True.Model, reduced = reduced)


	ES_vec <- (X_mat_returner(mod, reduced = reduced) %*% ES_test)[,1]
	ES_vec_truemodel  <- (X_mat_returner(True.Model) %*% c(0, ES_G, ES_E, ES_GE))[,1]

	y_bar <- p_vec_returner(MAF, P_e) * ES_vec_truemodel
	y_bar_sum <- sum(y_bar)

	sd_y_x <- sqrt(sd_y^2  - sum(p_vec * (ES_vec - y_bar_sum)^2))

	return(sd_y_x)
}

#' Function to calculate the standard deviation of y given x for linear models with logistic environment interaction
#'
#' Returns the standard deviation of y given x for linear models with logistic environment interaction
#'
#' @param beta_hat Effect sizes from MLE
#' @param MAF Minor allele Frequency
#' @param P_e Population prevalence of logistic environmental factor
#' @param ES_G Genetic Effect size
#' @param ES_E Environment Effect size
#' @param ES_GE Environment x Genetic interaction Effect size
#' @param sd_y_x_truth Standard deviation of y for the true model
#' @param sd_y_x_test Standard deviation of y for the test model
#' @param sd_y Standard deviation of y
#' @param Test.model Test model
#' @param True.Model True model
#'
#' @return The standard deviation of y given x for linear models with logistic environment interaction
#'
#' @export
#'
calc.like.linear.log.envir.interaction <- function(beta_hat, MAF, P_e, ES_G, ES_E, ES_GE, sd_y_x_truth, sd_y_x_model, Test.Model, True.Model, reduced  = F)
{
	p_vec <- p_vec_returner(MAF, P_e)

	betavec <- (X_mat_returner(Test.Model, reduced = reduced) %*% beta_hat)[,1]
	ES_vec <- (X_mat_returner(True.Model) %*% c(0, ES_G, ES_E, ES_GE))[,1]


	ll <- sum(p_vec * apply(cbind(ES_vec, betavec), 1, function(x) 
		expected.linear.ll(mean_truth = x[1], mean_model = x[2], sd_y_x_truth, sd_y_x_model)))

	return(ll)
}
