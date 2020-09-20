# CONTAINS CONSTRUCTED FUNCTIONS THROUGHOUT THE THESIS. NOT COMMENTED YET. 
#
#
#
# DOES NOT CONTAIN CODE FROM LAST YEARS PROJECT. 
#
#
#
#############################################################################
ewma.filter <- function(x, half.life, correlation = F, lambda = NULL){

	iT <- nrow(x)

	covar <- array(0L, dim = c(ncol(x),ncol(x), iT))

	#unconditional covariance for first month.
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

ewma.filter2006 <- function(data, correlation = F, tau1 = NULL, rho = NULL){

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


	for(k in 1:kmax){
		#now it is initialized using sample covariance
		vol[,,1] <- cov(data[1:100, ])
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