ll_zero_finder <- function(af, ii = 4, lower = 0, upper = 1){
	# finds the zeros of a function between 0 and 1
	## next step is to reduce amount of calculations by only performing calculations near zero instead of the whole interval
	ii0 <- 0
	mdelta <- 1e-2
	xmin = lower
	xmax = upper
	non_consecutive_zeros <- F
	while(ii0 != ii){
		xl = seq(from = xmin - 10*mdelta, to = xmax + 10*mdelta, by = mdelta)
		xl <- xl[xl < upper & xl > lower]
		yl = sapply(xl, function(anx){
			warning = error = 0 #For error handling
			tryCatch(resy <- af(anx), 
				warning=function(warn){warning<<-1},
				error=function(err){error<<-1}
				)
			if(warning + error > 0) resy <- NA
			if(is.nan(resy)) resy <- NA
			return(resy)
		})
		if(any(diff(which(!is.na(yl))) != 1) & !(non_consecutive_zeros)) non_consecutive_zeros <- T
		xl0 <- xl[!is.na(yl)]
		yl0 <- yl[!is.na(yl)]
		xmin <- min(xl0)
		xmax <- max(xl0)
		ii0 <- ii0 + 1
		mdelta <- mdelta/10
		# first couple of iterations finds zeros, next few is just for increased accuracy, 
		#  so to reduce calculation time, we want to close in on the actual zero
		if(ii0 > 1){
			if(sum(diff(sign(yl0)) != 0) > 1) stop("after two iterations, still have multiple zeros")
			# xmin <- max(xl0[which.min(abs(yl0))] - 3*mdelta, xmin)
			# xmax <- min(xl0[which.min(abs(yl0))] + 3*mdelta, xmax)
			if(xmin >= xmax) stop("min/max calculation error")
		}
	}
	res <- list(range = c(xmin, xmax), zero = xl0[which.min(abs(yl0))], ncz = non_consecutive_zeros)
	return(res)
}


ll_zero_finder2 <- function(af, ii = 6, lower = 0, upper = 1, qdelta = 27){
	if(ii < 4) print("if iterations 'ii' set to less than 4, the function may not detect multiple zeros which are close together")
	# finds the zeros of a function between 0 and 1
	## next step is to reduce amount of calculations by only performing calculations near zero instead of the whole interval
	# want this one to be able to find multiple zeros
	# "qdelta" is the factor by which we multiply the y-value closest to zero and encompass all y-values in that range. 
	#  making it smaller should decrease calculation burden, but might obscure the solution
	ii0 <- 0
	mdelta <- 1e-2
	mxrange <- list(c(lower,upper))
	non_consecutive_zeros <- F
	t1 <- Sys.time()
	xl0<-yl0<-numeric(0)
	while(ii0 != ii){
		xl = sort(unique(unlist(lapply(mxrange, function(arange) seq(from = arange[1] - 10*mdelta, to = arange[2] + 10*mdelta, by = mdelta)))))
		xl <- xl[xl < upper & xl > lower]
		# print("0")
		yl = sapply(xl, function(anx){
			if(anx %in% xl0) return(yl0[xl0 == anx]) #reuse values we already calculated
			warning = error = 0 #For error handling
			tryCatch(resy <- af(anx), 
				warning=function(warn){warning<<-1},
				error=function(err){error<<-1}
				)
			if(warning + error > 0) resy <- NA
			if(any(is.nan(resy))) resy <- NA
			return(resy)
		})
		# print("1")
		# if(any(diff(which(!is.na(yl))) != 1) & !(non_consecutive_zeros)) non_consecutive_zeros <- T
		xl0 <- xl[!is.na(yl)]
		yl0 <- yl[!is.na(yl)]
		# go through the y-values with 2 values and pick the ones that are closest to their adjacent values
		# print("2")
		if(any(sapply(yl0, length) > 1)){
			# print("3")
			if(all(sapply(yl0, length) > 1)) stop("all x values have two y values")
			l1 <- min(which(sapply(yl0, length) == 1))
			for(zz in l1:length(yl0)){
				if(length(yl0[[zz]]) == 1) next
				yl0[[zz]] <- yl0[[zz]][which.min(abs(yl0[[zz]] - yl0[[zz - 1]]))]
			}
			l1 <- min(which(sapply(yl0, length) == 1))
			for(zz in l1:1){
				if(length(yl0[[zz]]) == 1) next
				yl0[[zz]] <- yl0[[zz]][which.min(abs(yl0[[zz]] - yl0[[zz + 1]]))]
			}
			if(any(sapply(yl0, length) > 1)){stop("double y value reduction to 1 did not work")}else  yl0 <- unlist(yl0)
		}
		ii0 <- ii0 + 1
		# print(ii0)
		mdelta <- mdelta/10
		# first couple of iterations finds zeros, next few is just for increased accuracy, 
		#  so to reduce calculation time, we want to close in on the actual zero
		if(ii0 == 1){
			if(length(xl0) > 0){
				mxrange <- list(c(min(xl0), max(xl0)))
			}else mxrange <- list(c(0,1))
		}
		if(ii0 > 1){
			# our "window" is defined by the x values for which the absolute value of the y-values is close the the minimum absolute y value
			iwin <- which(abs(yl0) < qdelta * mdelta)
			non_consecutive_zeros <- F
			if(length(table(sign(yl0))) == 1){
				iwin <- which(abs(yl0) < mdelta*100)
				if(length(iwin) == 0) iwin <- sort(order(abs(yl0))[1:min(100, length(yl0))])
				non_consecutive_zeros <- T
			}
			# looks for gaps in selected points. The length of "iwinv" should indicate how many windows exist in which the function is close to 0
			iwinv <- c(0, which(diff(iwin) != 1))
			# determine our x-ranges for which the y-value is close to 0
			# print(sprintf("iwinv: %s", iwinv))
			# print(sprintf("mdelta: %s", mdelta))
			# print(sprintf("iwin: %s", paste0(iwin,collapse=",")))
			mxrange <- lapply(data.frame(t(cbind(iwinv + 1, c(iwinv[-1], length(iwin))))), function(mvec) return(xl0[c(iwin[mvec[1]], iwin[mvec[2]])]))
			names(mxrange) <- NULL
			# print(sprintf("mxrange: %s", mxrange))
			# print(sprintf("ii0: %s", ii0))
			# cat("\n\nnext\n\n")
			# minima should be greater than maxima
			# iiA <<- ii0
			if(any(sapply(mxrange, function(xx) xx[1] > xx[2]))) stop("min/max calculation error")
		}
	}
	t2 <- Sys.time()
	if(min(abs(yl0) > 1e-3)) warning(sprintf(
		"minimum absolute y-value of log likelihood function is %s, which is not in the range [-0.001, 0.001], answer may not be correct", yl0[which.min(abs(yl0))]))
	res <- as.numeric(sapply(mxrange, function(arange) xl0[xl0 >= arange[1] & xl0 <= arange[2]][which.min(abs(yl0[xl0 >= arange[1] & xl0 <= arange[2]]))]))
	# res <- list(range = c(xmin, xmax), zero = xl0[which.min(abs(yl0))], ncz = non_consecutive_zeros)
	if(non_consecutive_zeros) warning("No sign change detected in log likelihood minimization function")
	return(res)
}

