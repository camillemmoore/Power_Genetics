#' Function to Calculate t matrix for logistic outcome with binary environment interaction in dominant model
#'
#' Calculates the t matrix for logistic outcome with binary environment interaction in dominant model
#'
#' @param MAF Vector of minor allele frequencies
#' @param P_e Vector of proportions of the population with exposure to the environmental effect
#' @param Case.Rate proportion of cases in the sample (cases/(cases + controls)). 
#' @param OR_G Vector of genetic odds ratios to detect
#' @param OR_E Vector of environmental odds ratios to detect
#' @param OR_GE Vector of genetic/environmental interaction odds ratios to detect
#'
#' @return t matrix for all combinations of environment/outcome
#'
#' @examples
#' dom.fun.t(MAF = 0.1, P_e = 0.2, Case.Rate = 0.5,
#' 	OR_G = 1.5, OR_E = 2, OR_GE = 1.8)
#'
#' @export
#'
dom.fun.t <- function(MAF, P_e, OR_E, OR_G, OR_GE, Case.Rate)
{
	P_AA <- (1 - MAF)^2
	P_AB <- 2*MAF*(1 - MAF)
	P_BB <- MAF^2

	dom.fun <- function(x){
		f <- numeric(3)

		# f[1] <- (P_AA * (1 - P_e) - 1 + Case.Rate + sum(x)) * x[1]/((1 - Case.Rate - sum(x)) * (P_AA * P_e - x[1])) - OR_E
		# f[2] <- (P_AA * (1 - P_e) - 1 + Case.Rate + sum(x)) * x[2]/((1 - Case.Rate - sum(x)) * ((1 - P_AA) * (1 - P_e) - x[2])) - OR_G
		# f[3] <- (P_AA * (1 - P_e) - 1 + Case.Rate + sum(x)) * x[3]/((1 - Case.Rate - sum(x)) * ((1 - P_AA) * P_e - x[3])) - OR_E*OR_G*OR_GE
		f[1] <- ((1 - Case.Rate - sum(x)) * (P_AA * P_e - x[1]))/((P_AA * (1 - P_e) - 1 + Case.Rate + sum(x)) * x[1]) - OR_E
		f[2] <- ((1 - Case.Rate - sum(x)) * ((1 - P_AA) * (1 - P_e) - x[2]))/((P_AA * (1 - P_e) - 1 + Case.Rate + sum(x)) * x[2]) - OR_G
		f[3] <- ((1 - Case.Rate - sum(x)) * ((1 - P_AA) * P_e - x[3]))/((P_AA * (1 - P_e) - 1 + Case.Rate + sum(x)) * x[3]) - OR_E*OR_G*OR_GE
		f
	}
	tol = 0.6
	a<-b<-c<-x<-q<-r<-w<-y<- -1
	
	x.start.vals <- c(
		(1 - Case.Rate) * P_AA * (1 - P_e),
		(1 - Case.Rate) * (1 - P_AA) * (1 - P_e),
		(1 - Case.Rate) * (1 - P_AA) * P_e)

	while(any(c(a,b,c,x,q,r,w,y) < 0 | any(c(a,b,c,x,q,r,w,y) > 1))){
		tol <- max(c(tol*2/3, 0.2))
		res <- zero_finder_nleqslv(dom.fun, veclength = 3, x.start.vals = x.start.vals)

		a <- res[1]
		b <- res[2]
		c <- res[3]
		x <- (1 - P_AA) * (1 - P_e) - b
		q <- 1 - Case.Rate - a - b - c
		r <- P_AA * (1 - P_e) - 1 + Case.Rate + a + b + c
		w <- P_AA * P_e - a
		y <- (1 - P_AA)*P_e - c
	}

	# dominant
	prob_AA_control_e0 <- q
	prob_AB_control_e0 <- b * P_AB / (P_AB + P_BB)
	prob_BB_control_e0 <- b * P_BB / (P_AB + P_BB)
	prob_AA_case_e0 <- r
	prob_AB_case_e0 <- x * P_AB / (P_AB + P_BB)
	prob_BB_case_e0 <- x * P_BB / (P_AB + P_BB)
	prob_AA_control_e <- a
	prob_AB_control_e <- c * P_AB / (P_AB + P_BB)
	prob_BB_control_e <- c * P_BB / (P_AB + P_BB)
	prob_AA_case_e <- w
	prob_AB_case_e <- y * P_AB / (P_AB + P_BB)
	prob_BB_case_e <- y * P_BB / (P_AB + P_BB)

	t <- rbind(
			c(prob_AA_case_e0, prob_AB_case_e0, prob_BB_case_e0, prob_AA_case_e, prob_AB_case_e, prob_BB_case_e), 
			c(prob_AA_control_e0, prob_AB_control_e0, prob_BB_control_e0, prob_AA_control_e, prob_AB_control_e, prob_BB_control_e)
	)
	return(t)
}

#' Function to Calculate t matrix for logistic outcome with binary environment interaction in recessive model
#'
#' Calculates the t matrix for logistic outcome with binary environment interaction in recessive model
#'
#' @param MAF Vector of minor allele frequencies
#' @param P_e Vector of proportions of the population with exposure to the environmental effect
#' @param Case.Rate proportion of cases in the sample (cases/(cases + controls)). 
#' @param OR_G Vector of genetic odds ratios to detect
#' @param OR_E Vector of environmental odds ratios to detect
#' @param OR_GE Vector of genetic/environmental interaction odds ratios to detect
#'
#' @return t matrix for all combinations of environment/outcome
#'
#' @examples
#' rec.fun.t(MAF = 0.1, P_e = 0.2, Case.Rate = 0.5,
#' 	OR_G = 1.5, OR_E = 2, OR_GE = 1.8)
#'
#' @export
#'
rec.fun.t <- function(MAF, P_e, OR_E, OR_G, OR_GE, Case.Rate)
{
	P_AA <- (1 - MAF)^2
	P_AB <- 2*MAF*(1 - MAF)
	P_BB <- MAF^2

	rec.fun <- function(x){
		f <- numeric(3)

		f[1] <- ((1 - Case.Rate - sum(x)) * ((1 - P_BB) * P_e - x[1]))/(((1 - P_BB) * (1 - P_e) - 1 + Case.Rate + sum(x)) * x[1]) - OR_E
		f[2] <- ((1 - Case.Rate - sum(x)) * (P_BB * (1 - P_e) - x[2]))/(((1 - P_BB) * (1 - P_e) - 1 + Case.Rate + sum(x)) * x[2]) - OR_G
		f[3] <- ((1 - Case.Rate - sum(x)) * (P_BB * P_e - x[3]))/(((1 - P_BB) * (1 - P_e) - 1 + Case.Rate + sum(x)) * x[3]) - OR_E*OR_G*OR_GE
		f
	}
	tol = 0.6
	a<-b<-c<-x<-q<-r<-w<-y<- -1

	x.start.vals <- c(
		(1 - Case.Rate) * (1 - P_BB) * (1 - P_e),
		(1 - Case.Rate) * P_BB * (1 - P_e),
		(1 - Case.Rate) * P_BB * P_e)

	while(any(c(a,b,c,x,q,r,w,y) < 0 | any(c(a,b,c,x,q,r,w,y) > 1))){
		tol <- max(c(tol*2/3, 0.2))
		res <- zero_finder_nleqslv(rec.fun, veclength = 3, x.start.vals = x.start.vals)

		a <- res[1]
		b <- res[2]
		c <- res[3]
		x <- P_BB * (1 - P_e) - b
		q <- 1 - Case.Rate - a - b - c
		r <- (1 - P_BB) * (1 - P_e) - 1 + Case.Rate + a + b + c
		w <- (1 - P_BB) * P_e - a
		y <- P_BB*P_e - c
	}

	# recessive
	prob_AA_control_e0 <- q * P_AA / (P_AA + P_AB)
	prob_AB_control_e0 <- q * P_AB / (P_AA + P_AB)
	prob_BB_control_e0 <- b
	prob_AA_case_e0 <- r * P_AA / (P_AA + P_AB) 
	prob_AB_case_e0 <- r * P_AB / (P_AA + P_AB)
	prob_BB_case_e0 <- x
	prob_AA_control_e <- a * P_AA / (P_AA + P_AB)
	prob_AB_control_e <- a * P_AB / (P_AA + P_AB)
	prob_BB_control_e <- c
	prob_AA_case_e <- w * P_AA / (P_AA + P_AB)
	prob_AB_case_e <- w * P_AB / (P_AA + P_AB)
	prob_BB_case_e <- y

	t <- rbind(
			c(prob_AA_case_e0, prob_AB_case_e0, prob_BB_case_e0, prob_AA_case_e, prob_AB_case_e, prob_BB_case_e), 
			c(prob_AA_control_e0, prob_AB_control_e0, prob_BB_control_e0, prob_AA_control_e, prob_AB_control_e, prob_BB_control_e)
	)
	return(t)
}

#' Function to Calculate t matrix for logistic outcome with binary environment interaction in additive model
#'
#' Calculates the t matrix for logistic outcome with binary environment interaction in additive model
#'
#' @param MAF Vector of minor allele frequencies
#' @param P_e Vector of proportions of the population with exposure to the environmental effect
#' @param Case.Rate proportion of cases in the sample (cases/(cases + controls)). 
#' @param OR_G Vector of genetic odds ratios to detect
#' @param OR_E Vector of environmental odds ratios to detect
#' @param OR_GE Vector of genetic/environmental interaction odds ratios to detect
#'
#' @return t matrix for all combinations of environment/outcome
#'
#' @examples
#' add.fun.t(MAF = 0.1, P_e = 0.2, Case.Rate = 0.5,
#' 	OR_G = 1.5, OR_E = 2, OR_GE = 1.8)
#'
#' @export
#'
add.fun.t <- function(MAF, P_e, OR_E, OR_G, OR_GE, Case.Rate)
{
	P_AA <- (1 - MAF)^2
	P_AB <- 2*MAF*(1 - MAF)
	P_BB <- MAF^2

	add.fun <- function(x){
		f <- numeric(5)

		# f[2] <- (P_AA * P_e - x[1]) * x[2] / (x[1] * (P_AB * (1 - P_e) - x[2])) - OR_G / OR_E
		# f[3] <- ((P_AB * (1 - P_e) - x[2])) * x[3] / (x[2] * (P_AB * P_e - x[3])) - OR_E * OR_GE
		# f[4] <- (P_AB * P_e - x[3]) * x[4] / (x[3] * (P_BB * (1 - P_e) - x[4])) - OR_G / (OR_E * OR_GE)
		# f[5] <- (P_AB * P_e - x[3]) * x[5]/(x[3] * (P_BB * P_e - x[5])) - OR_G*OR_GE
		# f[5] <- (P_BB * (1 - P_e) - x[4]) * x[5]/(x[4] * (P_BB * P_e - x[5])) - OR_E * OR_GE^2
		
		# f[1] <- (P_AA * (1 - P_e) - 1 + Case.Rate + sum(x)) * x[1]/((1 - Case.Rate - sum(x)) * (P_AA * P_e - x[1])) - OR_E
		# f[2] <- (P_AA * (1 - P_e) - 1 + Case.Rate + sum(x)) * x[2]/((1 - Case.Rate - sum(x)) * (P_AB * (1 - P_e) - x[2])) - OR_G
		# f[3] <- (P_AA * (1 - P_e) - 1 + Case.Rate + sum(x)) * x[3]/((1 - Case.Rate - sum(x)) * (P_AB * P_e - x[3])) - OR_G*OR_E*OR_GE
		# f[4] <- (P_AA * (1 - P_e) - 1 + Case.Rate + sum(x)) * x[4]/((1 - Case.Rate - sum(x)) * (P_BB * (1 - P_e) - x[4])) - OR_G^2		
		# f[5] <- (P_AA * (1 - P_e) - 1 + Case.Rate + sum(x)) * x[5]/((1 - Case.Rate - sum(x)) * (P_BB * P_e - x[5])) - OR_G^2*OR_E*OR_GE^2
		f[1] <- ((1 - Case.Rate - sum(x)) * (P_AA * P_e - x[1])) / ((P_AA * (1 - P_e) - 1 + Case.Rate + sum(x)) * x[1]) - OR_E
		f[2] <- ((1 - Case.Rate - sum(x)) * (P_AB * (1 - P_e) - x[2])) / ((P_AA * (1 - P_e) - 1 + Case.Rate + sum(x)) * x[2]) - OR_G
		f[3] <- ((1 - Case.Rate - sum(x)) * (P_AB * P_e - x[3])) / ((P_AA * (1 - P_e) - 1 + Case.Rate + sum(x)) * x[3]) - OR_G*OR_E*OR_GE
		f[4] <- ((1 - Case.Rate - sum(x)) * (P_BB * (1 - P_e) - x[4])) / ((P_AA * (1 - P_e) - 1 + Case.Rate + sum(x)) * x[4]) - OR_G^2		
		f[5] <- ((1 - Case.Rate - sum(x)) * (P_BB * P_e - x[5])) / ((P_AA * (1 - P_e) - 1 + Case.Rate + sum(x)) * x[5]) - OR_G^2*OR_E*OR_GE^2

		
		
		f
	}
	tol = 0.6
	a<-b<-c<-d<-e<-q<-r<-t<-w<-x<-y<-z<- -1

	x.start.vals <- c(
		(1 - Case.Rate) * P_AA * P_e,
		(1 - Case.Rate) * P_AB * (1 - P_e),
		(1 - Case.Rate) * P_AB * P_e,
		(1 - Case.Rate) * P_BB * (1 - P_e),
		(1 - Case.Rate) * P_BB * P_e)

	while(any(c(a,b,c,d,e,q,r,t,w,x,y,z) < 0 | any(c(a,b,c,d,e,q,r,t,w,x,y,z) > 1))){
		tol <- max(c(tol*2/3, 0.2))
		res <- zero_finder_nleqslv(add.fun, veclength = 5, x.start.vals = x.start.vals)

		a <- res[1]
		b <- res[2]
		c <- res[3]
		d <- res[4]
		e <- res[5]
		z <- P_BB * P_e - e
		y <- P_BB * (1 - P_e) - d
		x <- P_AB * P_e - c
		w <- P_AB * (1 - P_e) - b
		t <- P_AA * P_e - a
		r <- P_AA * (1 - P_e) - 1 + Case.Rate + a + b + c + d + e
		q <- 1 - Case.Rate - a - b - c - d - e
	}

	# additive
	prob_AA_control_e0 <- q
	prob_AB_control_e0 <- b
	prob_BB_control_e0 <- d
	prob_AA_case_e0 <- r
	prob_AB_case_e0 <- w
	prob_BB_case_e0 <- y
	prob_AA_control_e <- a
	prob_AB_control_e <- c
	prob_BB_control_e <- e
	prob_AA_case_e <- t
	prob_AB_case_e <- x
	prob_BB_case_e <- z

	t <- rbind(
			c(prob_AA_case_e0, prob_AB_case_e0, prob_BB_case_e0, prob_AA_case_e, prob_AB_case_e, prob_BB_case_e), 
			c(prob_AA_control_e0, prob_AB_control_e0, prob_BB_control_e0, prob_AA_control_e, prob_AB_control_e, prob_BB_control_e)
	)
	return(t)
}
