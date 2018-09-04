zero_finder_nleqslv <- function(afun, veclength, tol = 0.4, x.start.vals = NULL){
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



# I don't know what this is for....
zero_finder_nleqslv.2 <- function(afun, iter = 60)
{
	conv <- reps <- 0
	fres <- numeric(0)
	while(reps < iter){
		reps <- reps + 1
		x.start <- c(runif(1)*0.8)
		res <- tryCatch(nleqslv(x.start, afun), error = function(e) "error", warning = function(w) "warning")
		if(length(res) != 1){
			if(res$message == "Function criterion near zero" & all(res$x > 0) & all(res$x < 1)){
				if(res$x > 0 & res$x < 1){
					if(length(fres) == 0){
						fres <- c(fres, res$x)
					}else if(!any(sapply(fres, function(q) (all.equal(q, res$x, tolerance = 8e-7) == T)))){
						fres <- c(fres, res$x)
					}
				}
			}
		}
	}
	return(fres)
}
