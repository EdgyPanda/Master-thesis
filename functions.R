# CONTAINS CONSTRUCTED FUNCTIONS THROUGHOUT THE THESIS. NOT COMMENTED YET. 
#
#
#
# DOES NOT CONTAIN CODE FROM LAST YEARS PROJECT. 
#
#
#
#############################################################################
ewma.filter <- function(x, half.life, correlation = F, realized = F, lambda = NULL){

	iT <- nrow(x)

	covar <- array(0L, dim = c(ncol(x),ncol(x), iT))

	#unconditional covariance for first month.
	if(realized){
		covar[,,1] <- realCov(x[1:100, ])
	}
	covar[,,1] <- cov(x[1:100, ])


	#calculation for half-life: lmd^t = 0.5 (half-life def) <=> t = ln(0.5)/ln(lmd) <=> lmd = exp(ln(0.5)/t)
	if(is.null(lambda)){

		lambda <- exp(log(0.5)/half.life)

	}

	print(sprintf("Lambda: %s", lambda))
	half.life <- log(0.5)/log(lambda)
	print(sprintf("half life: %s", half.life))

	for(i in 2:iT){


	covar[,,i] <-  lambda * covar[,,i-1] + (1-lambda) * t(x[i-1,])%*%x[i-1,]

	} 

	if(correlation){

		corr <- array(0L, dim = c(ncol(x),ncol(x), iT))
		#corr[,,1] <- cor(x[1:100, ])

		for(i in 1:iT){

		d <- sqrt(diag(covar[,,i]))
		d <- diag(d, ncol(covar), ncol(covar))

		corr[,,i] <- inv(d) %*% covar[,,i] %*% inv(d) 
		}

		return(corr)
	}


return(covar)

}

ewma.filter2006 <- function(data, correlation = F, realized = FALSE, tau1 = NULL, rho = NULL){

	#Reference values found in RiskMetrics 2006:
	#Rho and tau1 controls the half-life construction in the 1993 model. 
    #set Rho =1 and vary tau1 to vary the half-life. 
	if(is.null(tau1)){
    	tau1 <- 4
    }
    
    if(is.null(rho)){
    	rho <- sqrt(2)
	}

	tau0 <- 1560
    taumax <- 512
    if(is.null(tau1) && is.null(rho)){

    	kmax   <- round(log(taumax/tau1)/log(rho))

    }
    else{kmax <- 14 }
	#The one business day horizon at which
	#the data are available is denoted by dt.
	dt <- 1
	tau <- numeric(kmax)
	mu <- numeric(kmax)
	iT <- nrow(data)

	vol <- array(0L, dim = c(ncol(data), ncol(data), iT))
	volk <- list()
	weights <- numeric(kmax)

	if(realized){

		vol[,,1] <- realCov(data[1:100, ])
	}

	vol[,,1] <- cov(data[1:100, ])

	for(k in 1:kmax){
		#now it is initialized using sample covariance
		tau[k] <- tau1 * rho^(k-1)  
		mu[k] <- exp(-dt/tau[k])
		weights[k] <- 1-log(tau[k]/tau0)
		#Calculate recursion of kth EWMA.

		#need to initialize the vol formula. He initializes it using a backcast??

		for(t in 2:iT){

			vol[,,t] <- mu[k] * vol[,,t-1] + (1 - mu[k]) *  t(data[t-1,])%*%data[t-1, ]

		}

		#saves kth recursion in list. 
		volk[[k]] <- vol 
		

	}
	#we pre-calculated the weights and made them sum to 1 by normalizing by the sum of the weights:
	weigths <- weights / sum(weights)
	
	#efficient vol is now a weighted sum of the k EWMA:

	effvol <- array(0L, dim=c(ncol(data), ncol(data), iT))

	for(k in 1:kmax){

		effvol <- effvol + weights[k] * volk[[k]]

	}

	if(correlation){

		corr <- array(0L, dim = c(ncol(data),ncol(data), iT))
		#corr[,,1] <- cor(x[1:100, ])

		for(i in 1:iT){

		d <- sqrt(diag(effvol[,,i]))
		d <- diag(d, ncol(effvol), ncol(effvol))

		corr[,,i] <- inv(d) %*% effvol[,,i] %*% inv(d) 
		}

		return(corr)
	}

return(effvol)

}

ewma.filter2006.realized <- function(covariance, correlation = F, tau1 = NULL, rho = NULL){

	#Reference values found in RiskMetrics 2006:
	#Rho and tau1 controls the half-life construction in the 1993 model. 
    #set Rho =1 and vary tau1 to vary the half-life. 
	if(is.null(tau1)){
    	tau1 <- 4
    }
    
    if(is.null(rho)){
    	rho <- sqrt(2)
	}

	tau0 <- 1560
    taumax <- 512
    if(is.null(tau1) && is.null(rho)){

    	kmax   <- round(log(taumax/tau1)/log(rho))

    }
    else{kmax <- 14 }
	#The one business day horizon at which
	#the data are available is denoted by dt.
	dt <- 1
	tau <- numeric(kmax)
	mu <- numeric(kmax)
	iT <- length(covariance[1,1,])

	vol <- array(0L, dim = c(2, 2, iT))
	volk <- list()
	weights <- numeric(kmax)

	sampleR <- matrix(
	c(mean(covariance[1,1,1:100], na.rm = T),
	  mean(covariance[2,1,1:100], na.rm = T),
	  mean(covariance[2,1,1:100], na.rm = T),
	  mean(covariance[2,2,1:100], na.rm = T)),
	  ncol=2, nrow=2)

	#initializing by sample mean of realized covariance
	vol[,,1] <- sampleR
	
	for(k in 1:kmax){
		#now it is initialized using sample covariance
		tau[k] <- tau1 * rho^(k-1)  
		mu[k] <- exp(-dt/tau[k])
		weights[k] <- 1-log(tau[k]/tau0)
		#Calculate recursion of kth EWMA.

		#need to initialize the vol formula. He initializes it using a backcast??

		for(t in 2:iT){

			vol[,,t] <- mu[k] * vol[,,t-1] + (1 - mu[k]) *  covariance[,,t-1]

		}

		#saves kth recursion in list. 
		volk[[k]] <- vol 
		

	}
	#we pre-calculated the weights and made them sum to 1 by normalizing by the sum of the weights:
	weigths <- weights / sum(weights)
	
	#efficient vol is now a weighted sum of the k EWMA:

	effvol <- array(0L, dim=c(2, 2, iT))

	for(k in 1:kmax){

		effvol <- effvol + weights[k] * volk[[k]]

	}

	if(correlation){

		corr <- array(0L, dim = c(ncol(data),ncol(data), iT))
		#corr[,,1] <- cor(x[1:100, ])

		for(i in 1:iT){

		d <- sqrt(diag(effvol[,,i]))
		d <- diag(d, ncol(effvol), ncol(effvol))

		corr[,,i] <- inv(d) %*% effvol[,,i] %*% inv(d) 
		}

		return(corr)
	}

return(effvol)

}

ewma.filter.realized <- function(cov, half.life, correlation = F, lambda = NULL){

	iT <- length(cov[1,1,])

	covar <- array(0L, dim = c(2,2, iT))

	#unconditional covariance for first month.
	sampleR <- matrix(
	c(mean(cov[1,1,1:100], na.rm = T),
	  mean(cov[2,1,1:100], na.rm = T),
	  mean(cov[2,1,1:100], na.rm = T),
	  mean(cov[2,2,1:100], na.rm = T)),
	  ncol=2, nrow=2) 
	
	covar[,,1] <- sampleR


	#calculation for half-life: lmd^t = 0.5 (half-life def) <=> t = ln(0.5)/ln(lmd) <=> lmd = exp(ln(0.5)/t)
	if(is.null(lambda)){

		lambda <- exp(log(0.5)/half.life)

	}

	print(sprintf("Lambda: %s", lambda))
	half.life <- log(0.5)/log(lambda)
	print(sprintf("half life: %s", half.life))

	for(i in 2:iT){

	covar[,,i] <-  lambda * covar[,,i-1] + (1-lambda) * cov[,,i-1]

	} 

	if(correlation){

		corr <- array(0L, dim = c(2,2, iT))
		#corr[,,1] <- cor(x[1:100, ])

		for(i in 1:iT){

		d <- sqrt(diag(covar[,,i]))
		d <- diag(d, ncol(covar), ncol(covar))

		corr[,,i] <- inv(d) %*% covar[,,i] %*% inv(d) 
		}

		return(corr)
	}


return(covar)

}


realsemicov <- function(matrix,type, correlation = FALSE){
  # Realized SemiCovariance estimator for 1 day. 
  #
  # Args: 
  #
  #   matrix: A matrix of returns. Nxk 
  #
  #	  N: number of intraday returns, k: number of asset dimension, T: number of days. 
  realcovariance <- matrix()
  
  	positive <- matrix * (matrix >= 0)
  	negative <- matrix * (matrix <= 0)

  if(type =='P'){
  	realcovariance <- t(positive) %*% (positive)
  }
  if(type == 'N'){
  	realcovariance <- t(negative) %*% (negative)
  }
  if(type == 'M'){
  	realcovariance <- t(positive) %*% negative + t(negative) %*% positive
  }
  if(type == 'hM'){
  	realcovariance <- t(positive) %*% negative
  }

  if(correlation){

		d <- sqrt(diag(realcovariance))
		d <- diag(d, ncol(realcovariance), ncol(realcovariance))

		corr <- inv(d) %*% realcovariance %*% inv(d) 

		return(corr)
	}
  
 return( realcovariance ) 
}

bandwidthH <- function(list, sparsedata){

	n <- as.vector(sapply(list,length))

	c <- ((12)^2/0.269)^(1/5)

	w <- numeric()
	IV <- numeric()
	for(i in 1:length(list)){

		w[i] <- 1/(2*n[i]) * realCov(list[[i]])
		IV[i] <- realCov(sparsedata[[i]])

	}

	e <- w/(IV/n)

	H <- c * e^(4/5) * n^(3/5)

	return(H)

}

minvar <- function(Covar){

	ones <- matrix(rep(1, ncol(Covar)), ncol=1, nrow=ncol(Covar))

	t <- try(inv(Covar), silent = F)

	if("try-error" %in% class(t)){return(matrix(rep(NaN,ncol(Covar)), ncol=2, nrow =1))}
	else{
	w1 <- t %*% ones
	w2 <- t(ones) %*% t %*% (ones)

	w <- w1 %*% w2^-1

	return(w)
	}
}

#Kevin sheppard 2 sided hessian function:

hessian_2sided <- function(fun, theta, ...){

	f <- fun(params = theta, ...)
	h <- 1e-5 * abs(theta)
	thetah <- theta + h

	h <- thetah - theta 
	k <- length(theta)
	h <- diag(h)

	fp <- numeric(k)
	fm <- numeric(k)

	for (i in 1:k){

		fp[i] <- - fun(params =theta+h[i, ], ...)$fulldLLK #negating it since sum is
		fm[i] <- - fun(params =theta-h[i, ], ...)$fulldLLK
	}

	fpp <- matrix(0L, nrow = k, ncol = k)
	fmm <- matrix(0L, nrow = k, ncol = k)

	for(i in 1:k){
		for(j in i:k){
			fpp[i,j] <- fun(params = theta + h[i, ] + h[j, ], ...)$fulldLLK
			fpp[j,i] <- fpp[i,j]
			fmm[i,j] <- fun(params = theta - h[i, ] - h[j, ], ...)$fulldLLK
			fmm[j,i] <- fmm[i,j]

		}
	}


	hh <- diag(h)
	hh <- matrix(hh, nrow = k, ncol = 1)
	hh <- hh %*% t(hh)

	H <- matrix(0L, nrow = k, ncol = k)

	for(i in 1:k){
		for(j in i:k){
			H[i,j] <- (fpp[i,j] - fp[i] - fp[j] + f$fulldLLK  +  f$fulldLLK - fm[i] - fm[j] + fmm[i,j])/hh[i,j]/2
			H[j,i] <- H[i,j]
		}
	}

	return(H)

}

standardized.res <- function(returns, variances){
	#variances should be a matrix of ncol = number of assets, nrow = days. 

	lT <- ncol(returns)

	rT <- nrow(returns)

	res <- matrix(0L, nrow = rT, ncol = lT)

	for(i in 1:lT){
		for(j in 1:rT){

			res[j,i] <- returns[j,i] * (sqrt(variances[j,i]))^-1
		}
	}


	return(res)
}

standardized.res.intraday <- function(returns, variances){
	#variances should be a matrix of ncol = number of assets, nrow = days. 

	lT <- ncol(returns)

	rT <- nrow(returns)

	res <- matrix(0L, nrow = rT, ncol = lT)

	variances <- matrix(variances, nrow = 1, ncol = lT)

	variances <- matrix(rep(variances[1, ],rT), nrow =rT, ncol=lT, byrow = T)

	for(i in 1:lT){
		for(j in 1:rT){

			res[j,i] <- returns[j,i] * (sqrt(variances[j]))^-1
		}
	}


	return(res)
}


#UNIVARIATE INSANITY FILTER!
volatility.insanity.filter <- function(forecasts, lb, ub, adjusted_forecast){

	ut1 <- forecasts
	items <- sum(forecasts < lb) + sum( forecasts > ub)

	ut1[forecasts < lb] = adjusted_forecast
	ut1[forecasts > ub] = adjusted_forecast

	lOut <- list()
	lOut[["vol"]] <- ut1
	lOut[["items"]] <- items

	return(lOut)

}

RTQ <- function(rData) { 
  returns <- as.vector(as.numeric(rData))
  n <- length(returns)
  tq <- n * (n/(n-2)) *((2^(2/3) * gamma(7/6) * gamma(1/2)^(-1))^(-3)) *  sum(abs(returns[1:(n - 2)])^(4/3) * 
  	abs(returns[2:(n-1)])^(4/3) * abs(returns[3:n])^(4/3))
  return(tq)
} 


mlag <- function(x,n=1,init=0){

	xlag <- matrix(0L, ncol = ncol(x)*n, nrow = nrow(x))
	icnt <- 0 
	for(i in 1:ncol(x)){
		for(j in 1:n){
			xlag[(j+1):nrow(x), icnt+j] <- x[1:(nrow(x)-j),i]
		}
		icnt <- icnt+n
	}

	xlag[xlag == 0] <- init

	return(xlag)
}