#' Zero finder
#'
#' Finds zeros of multinomial functions using the nleqslv package
#'
#' @param afun The function to find zeros
#' @param veclength The dimension of the system of equations
#' @param tol The range within which to set start values for the function to use to find zeros
#' @param x.start.vals Optional user defined start values
#' @return Predicted zeros of the given equation
#'
#' @examples
#' afun <- function(x) {
#' 	y <- numeric(2)
#' 	y[1] <- x[1]^2 + x[2]^2 - 1
#' 	y[2] <- exp(x[1]-1) + x[2]^3 - 1
#' 	y
#' }
#'
#' zero_finder_nleqslv(afun)
#'
#' @export
#'
zero_finder_nleqslv <- function(afun, veclength, tol = 0.4, x.start.vals = NULL){
	library(nleqslv)
	# x.start gives example values to base x.start values on. If null, x.start values are randomly selected from runif
	conv <- reps <- 0
	while(conv == 0){
		reps <- reps + 1
		if(is.null(x.start.vals)){
			x.start <- c(runif(veclength)*tol)
		}else{
			if(reps >= 250){
				x.start <- c(runif(veclength)*tol)
			}else if(reps >= 150){
				x.start <- sapply(x.start.vals, function(x) abs(rnorm(mean = x, sd = x, n = 1)))
			}else if(reps >= 50){
				x.start <- sapply(x.start.vals, function(x) abs(rnorm(mean = x, sd = x/10, n = 1)))
			}else{
				x.start <- sapply(x.start.vals, function(x) abs(rnorm(mean = x, sd = x/100, n = 1)))
			}
		}
		res <- nleqslv(x.start, afun)
		if(res$message %in% c("Function criterion near zero", "x-values within tolerance 'xtol'") & 
			all(res$x > 0) & all(res$x < 1)) conv <- 1
		# if(reps == 60){conv <- 1; res <- rep(-1, veclength)}
	}
	res$x
}



