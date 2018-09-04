#' Zero finding function
#'
#' Finds the zeros of a function af. Alternative to uniroot, designed specifically to work with the genpwr package.
#'
#' @param af The function for which to find the zero(s)
#' @param ii Number of iterations. The more iterations, the more accuracy. It is recommended that ii be at least 4.
#' @param lower Lower limit of region in which to find the zero
#' @param upper Upper limit of region in which to find the zero
#' @param qdelta Factor for finding intervals over which the function is close to zero
#' @return Points over the given interval at which the given function is approximately equal to zero
#'
#' @examples
#' ll_zero_finder2(function(x) (x-0.5)^2 - 0.1)
#' ll_zero_finder2(function(x) 8*x^3 - 11.2*x^2 + 4.56*x - 0.476)
#'
#' @export
#'


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
	lower0<-lower;upper0<-upper
	xl0<-yl0<-res<-llz<-numeric(0)
	while(ii0 != ii){
		xl <- sort(unique(unlist(lapply(mxrange, function(arange) seq(from = arange[1] - 10*mdelta, to = arange[2] + 10*mdelta, by = mdelta)))))
		xl <- sort(c(mdelta/100, mdelta/10, xl)) # in case there are ones near zero
		lower0 <- min(unlist(mxrange) - 10*mdelta)
		upper0 <- max(unlist(mxrange) + 10*mdelta)
		xl <- xl[xl <+ upper0 & xl >= lower0]
		yl <- lapply(xl, function(anx){
			if(anx %in% xl0) return(unlist(yl0[xl0 == anx])) #reuse values we already calculated
			warning = error = 0 #For error handling
			tryCatch(resy <- af(anx), 
				warning=function(warn){warning<<-1},
				error=function(err){error<<-1}
				)
			if(warning + error > 0) resy <- NA
			if(any(is.nan(resy))) resy <- NA
			# resy <- suppressWarnings(af(anx))
			return(resy)
		})
		# if(any(diff(which(!is.na(yl))) != 1) & !(non_consecutive_zeros)) non_consecutive_zeros <- T
		xl0 <- xl[sapply(yl, function(a_yl) !all(is.na(a_yl)))]
		yl0 <- yl[sapply(yl, function(a_yl) !all(is.na(a_yl)))]
		# go through the y-values with 2 values and pick the ones that are closest to their adjacent values
		if(any(sapply(yl0, length) > 1)){
			# if(any(sapply(yl0, length)) == 1 & any(sapply(yl0, length) == 2)) stop("a function where some points have 2 values and some have 1") #temporarily looking for this situation so I can debug it
			if(any(sapply(yl0, length) > 2)) stop("solving for a function with three values, don't know how to handle it")
			if(F){
				# going to try a different approach than this, but I don't want to delete this yet
				# print("3")
				if(all(sapply(yl0, length) > 1)){
					if(ii0 == ii - 1){stop("all x values have two y values")
					}else yl0<-numeric(0)
				}else{
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
			}
			xl00 <- list(xl0, xl0)
			yl00 <- list()
			yl00 <- c(yl00, list(sapply(yl0, function(x) x[1])))
			yl00 <- c(yl00, list(sapply(yl0, function(x) x[length(x)])))
			if(!all(sapply(yl0, length) == 2)){
				ll1 <- which(sapply(yl0, length) == 1)
				ll2 <- which(sapply(yl0, length) == 2)
				if(!all(sort(c(ll1,ll2)) == 1:length(yl0))) stop("every yl0 element should have 1 or 2 values")
				lld <- c(1, which(diff(ll1) != 1), length(ll1))
				ll0 <- list()
				for(lli in 1:(length(lld) - 1)) ll0 <- c(ll0, list(ll1[lld[lli]:lld[lli+1]]))
				# when you are in the situation where some but not all of the x coordinates map to more than 1, you have to match the singular y values
				# 	to one of the multiple y values. To do this, you find the number "llz" which indicates which of the two sets of y values does not match
				# 	with the y-values which are singular. Once this is achieved, it doesn't need to be done again, so llz is not recalculatd every loop.
				if(length(llz) == 0){
					for(lli in 1:length(ll0)){
						#want to determine which element in yl00 the region with only 1 solution belongs to
						#first we look at the left side
						llzl <- llzr <- NA
						if(ll0[[lli]][1] != 1){
							llzl <- which.max(c(abs(diff(yl00[[1]][(ll0[[lli]][1]-1):(ll0[[lli]][1])])), abs(diff(yl00[[2]][(ll0[[lli]][1]-1):(ll0[[lli]][1])]))))
						}
						if(ll0[[lli]][1] != length(yl00[[1]])){
							llzr <- which.max(c(abs(diff(yl00[[1]][(ll0[[lli]][length(ll0[[lli]])]):(ll0[[lli]][length(ll0[[lli]])]+1)])), 
								abs(diff(yl00[[2]][(ll0[[lli]][length(ll0[[lli]])]):(ll0[[lli]][length(ll0[[lli]])]+1)]))))
						}
						if(!any(is.na(c(llzl,llzr))) & llzl != llzr) stop("trying to match region with two outcomes and it matches to different lists on either side")
						llz <- unique(na.omit(c(llzl,llzr)))
						if(is.na(llz)) stop("matching region with two outcomes did not work")
					}
				}
				yl00[[llz]] <- yl00[[llz]][-ll0[[lli]]]
				xl00[[llz]] <- xl00[[llz]][-ll0[[lli]]]
			}
		}else {yl00 <- list(unlist(yl0)); xl00 <- list(xl0)}
		ii0 <- ii0 + 1
		mdelta <- mdelta/10
		for(yl0ii in length(yl00)){xl00[[yl0ii]] <- xl00[[yl0ii]][!is.na(yl00[[yl0ii]])]; yl00[[yl0ii]] <- yl00[[yl0ii]][!is.na(yl00[[yl0ii]])]}
		# first couple of iterations finds zeros, next few is just for increased accuracy, 
		#  so to reduce calculation time, we want to close in on the actual zero
		if(ii0 == 1){
			if(length(xl0) > 0){
				mxrange <- list(c(max(0, min(xl0) - mdelta), min(max(xl0) + mdelta, 1)))
			}else mxrange <- list(c(0,1))
		}else{
			if(length(xl0) > 0){
				# our "window" is defined by the x values for which the absolute value of the y-values is close the the minimum absolute y value
				# iwin <- which(abs(yl0) < qdelta * mdelta)
				if(!class(yl0) == "list") yl00 <- list(yl0)
				mxrange0 <- list()
				no_x_axis_cross <- rep(F, length(yl00))
				for(yl0ii in 1:length(yl00)){
					yl0a <- yl00[[yl0ii]]
					xl0a <- xl00[[yl0ii]]
					yl0a <- yl0a[!is.na(yl0a)]
					iwin <- unique(as.numeric(sapply(which(diff(sign(yl0a)) != 0), function(x) (x-9):(x+10))))# adjust this 9 and it might make it take less time
					iwin <- iwin[iwin > 1 & iwin <= length(yl0a)]
					if(length(table(sign(yl0a))) == 1){
						iwin <- which(abs(yl0a) < mdelta*100)
						if(length(iwin) == 0) iwin <- sort(order(abs(yl0a))[1:min(110, length(yl0a))])
						no_x_axis_cross[yl0ii] <- T
					}
					# looks for gaps in selected points. The length of "iwinv" should indicate how many windows exist in which the function is close to 0
					iwinv <- c(0, which(diff(iwin) != 1))
					# determine our x-ranges for which the y-value is close to 0
					mxrange <- lapply(data.frame(t(cbind(iwinv + 1, c(iwinv[-1], length(iwin))))), function(mvec) return(xl00[[yl0ii]][c(iwin[mvec[1]], iwin[mvec[2]])]))
					names(mxrange) <- NULL
					mxrange <- lapply(mxrange, function(mvec){
						if(mvec[1] == mvec[2]){mvec[1] <- mvec[1] - mdelta*10; mvec[2] <- mvec[2] + mdelta*10}
						return(mvec)
						})
					if(ii0 == ii & !no_x_axis_cross[yl0ii]){
						res <- c(res, as.numeric(sapply(mxrange, function(arange) 
							xl0a[xl0a >= arange[1] & xl0a <= arange[2]][which.min(abs(yl0a[xl0a >= arange[1] & xl0a <= arange[2]]))])))
					}
					if(!all(c(no_x_axis_cross[yl0ii], ii0 >= 3, min(abs(yl0a), na.rm = T) > 10*mdelta))) mxrange0 <- c(mxrange0, mxrange)
				}
				mxrange <- mxrange0
				if(any(sapply(mxrange, function(xx) xx[1] > xx[2]))) stop("min/max calculation error")
				if(length(mxrange) == 0){res <- NA; ii0 <- ii}
			}else if(ii0 < 3) {mxrange <- list(c(0,1))
			}else{res <- NA; ii0 <- ii; no_x_axis_cross <- NA}
		}
	}
	t2 <- Sys.time()
	# res <- res[!no_x_axis_cross]
	if(F){
		# thinking I want to refrain from returning an answer if the function does not cross the x axis, instead of returning a warning
		if(any(sapply(yl00, function(yl00xx) min(abs(yl00xx))) < 1e-3) & non_consecutive_zeros) 
			warning("minimum absolute y-value of log likelihood function is %s, which is not in the range [-0.001, 0.001], answer may not be correct")
		# res <- as.numeric(sapply(mxrange, function(arange) xl0[xl0 >= arange[1] & xl0 <= arange[2]][which.min(abs(yl0[xl0 >= arange[1] & xl0 <= arange[2]]))]))
		if(non_consecutive_zeros) warning("No sign change detected in log likelihood minimization function")
	}
	if(length(res) == 0) res <- NA
	res <- res[res >= lower & res <= upper]
	return(res)
}

