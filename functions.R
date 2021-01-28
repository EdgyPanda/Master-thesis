# CONTAINS CONSTRUCTED FUNCTIONS THROUGHOUT THE THESIS. NOT COMMENTED YET. 
#
#
#
# DOES NOT CONTAIN CODE FROM LAST YEARS PROJECT. 
#
#
#
#############################################################################


########################################################################################################################
#
#						RISKMETRICS FILTERS
#
########################################################################################################################

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

ewma.filter.realized <- function(cov, half.life, correlation = F, lambda = NULL, trace){

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

	half.life <- log(0.5)/log(lambda)

	if(trace == 1){	
	print(sprintf("Lambda: %s", lambda))
	print(sprintf("half life: %s", half.life))
	}

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

min.qlike.riskmetrics <- function(cov, half.life, proxy, correlation = F, lambda){

	filter <- ewma.filter.realized(cov, half.life,F, lambda)

	qlikes.riskmetrics <- numeric()
	for(i in 1:length(filter[1,1,])){

		qlikes.riskmetrics[i] <- QLIKE(filter[,,i], proxy, 2)

	}

 return(mean(qlikes.riskmetrics))

}

########################################################################################################################
#
#						REALIZED AND INTERNAL FUNCTIONS FOR REALIZED MEASURES
#
########################################################################################################################

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


########################################################################################################################
#
#						GENERAL FUNCTIONS
#
########################################################################################################################

riskparity_2dim <- function(matrix, risktarget, rt = F){
  #
  #
  #
  #NOTE TO YOURSELF: Palomar uses non-sqrt portrisk, which gives 
  #reasonable risk for the unlevered risk-parity portfolio. Using sqrt
  #
  #REMEMBER TO USE COV(),SD()  etc. FOR CLOSE-TO-CLOSE RETURNS, SINCE REALCOV OVERESTIMATES
  #THE COVARIANCE FOR CLOSE-TO-CLOSE RETURNS
  
  bonds <- matrix[1,1]
  stocks <- matrix[2,2]
  
  w_1 <- sqrt(bonds)^-1 / (sqrt(stocks)^-1 + sqrt(bonds)^-1)
  w_2 <- sqrt(stocks)^-1 / (sqrt(stocks)^-1 + sqrt(bonds)^-1)
  
  w <- matrix(c(w_1, w_2), ncol=1, nrow=2) #dimnames = list(c(), c("Bond", "Stock"))
  
  #Palomar uses portfolio variance as portfolio risk, thus no sqrt. 
  portrisk <- as.numeric(sqrt(t(w) %*% (matrix) %*% w))
  
  riskcont <-  w * (matrix %*% w)/portrisk
  
  relativeriskcont <- (w * (matrix %*% w)) / portrisk^2
  
  if(rt){
    
    alpha <- risktarget / portrisk
    
    w_new <- w %*% alpha
    
    w_new <- matrix(c(w_new[1], w_new[2]), ncol=1, nrow=2)
    #here is sqrt, while palomar uses variance, no sqrt.
    portrisk <-  as.numeric(sqrt((t(w_new) %*% matrix %*% w_new)))  #gives marginal risk for each asset. 
    
    riskcont <- (w_new * (matrix %*% w_new)) / portrisk
    
    w_riskfree <- uniroot(function(x) colSums(w_new)+x-1, interval = c(-100,100), extendInt="yes")$root
    
    w_new <- matrix(c(w_new, w_riskfree), ncol=1, nrow=3)
    rownames(w_new) <- c("TLT", "SPY", "riskfree")
    
    lout <- list(w_new, portrisk, riskcont, alpha)
    
    names(lout) <- c("w", "portrisk", "riskcont", "alpha")
    
    return(lout)
  }
  w_riskfree <- uniroot(function(x) colSums(w)+x-1, interval = c(-100,100), extendInt="yes")$root

  w <- matrix(c(w_1, w_2, w_riskfree), ncol = 1, nrow = 3)

  lout <- list(w, portrisk, riskcont, relativeriskcont)
  names(lout) <- c("w", "portrisk", "riskcont", "relativeriskcont")
  
  return(lout)
  
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

#for univariate HAR models. 
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

#requires pracma
#takes entire array
invsquarerootmat <- function(Covar){

	iT <-  length(Covar[1,1,])
	d <- ncol(Covar[,,1]) 

	isrmat <- array(0L, dim = c(d,d,iT))

	for(i in 1:iT){

		isrmat[,,i] <- solve(sqrtm(t(chol(Covar[,,i])))$B %*% sqrtm(chol(Covar[,,i]))$B)

	}

	return(isrmat)

}

turnover <- function(weights, portret, assetret){

	assets <- ncol(weights)

	init.weight <- rep(0, assets)

	turnover <- numeric()
	turnover[1] <- abs(sum(weights[1, ])) #only when initial weights are zero. 

	for(i in 2:(nrow(portret)-1)){

		turnover[i] <- sum(abs( weights[i, ] -  weights[i-1, ] * (1 + assetret[i-1, ]) * (1+portret[i-1, ])^-1 ))


	}
	return(turnover)
}

portconcentration <- function(weights){

	concentration <- numeric()

	for(i in 1:(nrow(weights)-1)){

		concentration[i] <- sqrt(sum(weights[i, ]^2))

	}
	return(concentration)
}

########################################################################################################################
#
#						DRD-HAR MODEL AND CORRELATION MODELS. 
#
########################################################################################################################


minimizingfunc <- function(data, params){

	dA <- params[1]
	dB <- params[2]
	dC <- params[3]

	y <- as.numeric(data[[1]])
	x1 <- as.numeric(data[[2]])
	x2 <- as.numeric(data[[3]])
	x3 <- as.numeric(data[[4]])

	dint <- (1-dA-dB-dC) * mean(y)

	mini <- (y - dint - dA * x1 - dB * x2  - dC * x3)^2

	
	minis <- mini 
	mini <-  sum(mini)

	lOut <- list()
	lOut[["mini"]] <- mini 
	lOut[["minis"]] <- minis

	return(lOut)
}


ineqconstraint <- function(vPar, ...){

	dA <- vPar[1]  
	dB <- vPar[2]
	dC <- vPar[3]

	return(dA+dB+dC)

}

min.RSS <- function(vPar, data){

	rss <- minimizingfunc(data, vPar)$mini

	return(rss)

}


EstimatecorrHAR <- function(variances, correlation = NULL, proxy, trace=1, Forecasts = F, ineqfun = ineqconstraint, ineqLB = 0.00, ineqUB = 0.9999){
  #variances should be produced by the univariate models and not come from realized measures!
  #correlation should be found via eg. five min samples but adhere to the same frequency as variances from HAR models. 

  dA = 0.01
  dB  = 0.2185
  dC = 0.4
  ## vector of starting parameters
  vPar = c(dA, dB, dC)

  correlation <- as.matrix(correlation)
  corrday <- correlation
  corrweek <- rowMeans(cbind(correlation, mlag(correlation,4,mean(correlation))), na.rm = T)
  corrweek[is.na(corrweek)] <- mean(corrweek, na.rm = T)
  corrmonth <- rowMeans(cbind(correlation, mlag(correlation,21,mean(correlation))), na.rm = T) 
  corrmonth[is.na(corrmonth)] <- mean(corrmonth, na.rm = T)

  end <- length(correlation)

  data <- list(proxy[23:end], corrday[22:(end-1)], corrweek[22:(end-1)], corrmonth[22:(end-1)])
  optimizer = solnp(vPar, fun = min.RSS, data = data, 
                    ineqfun = ineqfun, #the inequality constraint
                    ineqLB  = ineqLB, ## the inequality lower bound
                    ineqUB = ineqUB, ## the inequality lower bound, i.e. 0.0 <= a + b + c < 0.9999
                    ## lower and upper bounds for all parameters
                    LB = c(0, 0, 0), UB = c(0.9999, 0.9999, 0.9999),
                    control = list(tol = 1e-22, outer.iter = 800, inner.iter = 1200, delta = 1e-8,
                    trace = trace))
                    

  params <- optimizer$pars

  min <- tail(optimizer$values, 1)

  hessian <- optimizer$hessian

  #calculating covariances: 

  hhatcorrHAR <- cbind(corrday[22:(end-1)], corrweek[22:(end-1)], corrmonth[22:(end-1)]) %*% matrix(params)

  if(Forecasts){ hhatcorrHAR <- cbind(corrday[end], corrweek[end], corrmonth[end]) %*% matrix(params) }

  covs <- hhatcorrHAR * sqrt(variances[,1]) * sqrt(variances[,2])

  vSigma2 <- array(0L, dim = c(2,2, (nrow(variances))))

  for(i in 1:(nrow(variances))){

	vSigma2[,,i] <- matrix(c(variances[i,1], covs[i], covs[i], variances[i,2]), ncol=2, nrow=2) 

  }

	  #R-squared

	Rsquared <- 1-var(correlation[23:end] - hhatcorrHAR)/var(correlation[23:end])



  lOut <- list()

  lOut[["vPar"]] <- params
  lOut[["vSigma2"]] <- vSigma2
  lOut[["R2"]] <- Rsquared
  lOut[["estcor"]] <- hhatcorrHAR
  lOut[["MSE"]] <- min/end
  lOut[["hessian"]] <- hessian

  return(lOut)

 }


EstimatecorrHAR2 <- function(returns, variances, proxy, trace=1, ineqfun = ineqconstraint, ineqLB = 0.00, ineqUB = 0.9999){
  #mEta is standardized residuals from univariate models!
  #variances should be produced by the univariate models and not come from realized measures!
  #correlation should be found via eg. five min samples but adhere to the same frequency as variances from HAR models. 
  mEta <- numeric()
  for(i in 1:nrow(variances)){
  										#returns should start at 22. 
  	mEta2 <- standardized.res.intraday(returns[[i]], variances[i, ])
  	mEta[i] <- realCov(mEta2, T)[2,1]

  }
  dA = 0.45 
  dB  = 0.09185
  dC = 0.26
  ## vector of starting parameters
  vPar = c(dA, dB, dC)

  mEta <- as.matrix(mEta)
  corrday <- mEta
  corrweek <- rowMeans(cbind(mEta, mlag(mEta,4,mean(mEta))))
  corrmonth <- rowMeans(cbind(mEta, mlag(mEta,21,mean(mEta)))) 

  data <- list(proxy, corrday, corrweek, corrmonth)

  optimizer = solnp(vPar, fun = min.RSS, data = data, 
                    ineqfun = ineqfun, #the inequality constraint
                    ineqLB  = ineqLB, ## the inequality lower bound
                    ineqUB = ineqUB, ## the inequality lower bound, i.e. 0.0 <= a + b + c < 0.9999
                    ## lower and upper bounds for all parameters
                    LB = c(0, 0, 0), UB = c(0.9999, 0.9999, 0.9999),
                    control = list(tol = 1e-18, outer.iter = 800, inner.iter = 1200, delta = 1e-7,
                    trace = trace))
                    

  params <- optimizer$pars

  min <- tail(optimizer$values, 1)

  hessian <- optimizer$hessian

  #calculating covariances: 
  
  hhatcorrHAR <- cbind(corrday, corrweek, corrmonth) %*% matrix(params)

  covs <- hhatcorrHAR * sqrt(variances[,1]) * sqrt(variances[,2])

  vSigma2 <- array(0L, dim = c(2,2, (nrow(variances))))

  for(i in 1:(nrow(variances))){

    vSigma2[,,i] <- matrix(c(variances[i,1], covs[i], covs[i], variances[i,2]), ncol=2, nrow=2) 

  }
  

  lOut <- list()

  lOut[["vPar"]] <- params
  lOut[["vSigma2"]] <- vSigma2
  lOut[["estcor"]] <- hhatcorrHAR
  lOut[["MSE"]] <- min/2516 
  lOut[["hessian"]] <- hessian
  lOut[["mEta"]] <- mEta

  return(lOut)

 }



########################################################################################################################
#
#					REALIZED SCALAR-BASED BIVARIATE GARCH MODEL (rBG)
#
########################################################################################################################


BivarGARCHFilter <- function(dailyret, params, covariance, inference = F){

	dAlpha <- params[1]
	dBeta <- params[2]

	samplecov <- cov(dailyret)


	if(inference){

		dAlpha <- params[1]
		dBeta <- params[2]
		samplecov <- params[3]

	}

	days <- length(covariance[1,1,])

	d <-  2 #being bivariate

	rcov <- array(covariance, dim = c(2,2,days))

	samplercov <- matrix(
		c(mean(rcov[1,1,]),
		  mean(rcov[2,1,]),
		  mean(rcov[2,1,]),
		  mean(rcov[2,2,])),
		  ncol=2, nrow=2)

	mSigma <- array(0L, dim = c(2, 2, days))

	#initializing  with unconditional mean. 
	id <- matrix(c(1,0,0,1), ncol = 2, nrow = 2)

	omega <-   samplecov - dAlpha * samplercov %*% id - dBeta * samplecov #  samplecov^(-1)

	mSigma[,,1] <- samplecov 

	#compute first observation of log-likelihood (we are minizing neg-log-likelihood):
	dLLK <-  log(det(mSigma[,,1])) + (dailyret[1, , drop = F]) %*% solve(mSigma[,,1]) %*% t(dailyret[1, , drop = F])


	dLLKs <- numeric()
	dLLKs[1] <-   dLLK
	for(i in 2:days){

		mSigma[,,i] <- omega + dBeta * mSigma[,,i-1] + dAlpha * rcov[,,i-1]

		#neglog collection for score calculation
		dLLKs[i] <-   (log(det(mSigma[,,i])) +  
		dailyret[i, , drop = F] %*% solve(mSigma[,,i]) %*% t(dailyret[i, , drop = F]))
		}

	fulldLLK <-  - 0.5 *(d* log(2*pi) * days +  sum(dLLKs)) #loglikelihood

	lOut <- list()

	lOut[["fulldLLK"]] <- fulldLLK
	lOut[["mSigma"]] <- mSigma
	lOut[["cov"]] <- omega
	lOut[["dLLKs"]] <- dLLKs
	lOut[["sampleH"]] <- samplecov

	return(lOut)
	
	
}

#objective function specific for scalar Bivariate GARCH
ObjFBivarGARCH <- function(vPar, dailyret, covariance) {
  
  dAlpha = vPar[1]
  dBeta  = vPar[2]
  dLLK = BivarGARCHFilter(dailyret, c(dAlpha, dBeta), covariance)$fulldLLK
  
  return(-as.numeric(dLLK))
}

ineqfun_GARCH_BIVAR <- function(vPar, ...) {
  dAlpha = vPar[1]
  dBeta  = vPar[2]
  
  return(dAlpha + dBeta)
}

EstimateBivarGARCH <- function(dailyret, covariance, bootstrap = FALSE, vPar=NULL, ineqfun_GARCH = ineqfun_GARCH_BIVAR, ineqLB = 0.00, ineqUB = 0.9999, tol = 1e-7){
  
  # We set starting value for alpha equal to 0.05, dBeta = 0.94, and chose omega to target
  # the empirical variance by targeting the unconditional variance of the 
  # GARCH model
  
  dAlpha = 0.20
  dBeta  = 0.50
  
  ## vector of starting parameters
  if(is.null(vPar)){
  vPar = c(dAlpha, dBeta)
  }
  # have a look at help(solnp)   ObjFBivarGARCH
  ##optimization step            
  optimizer = solnp(vPar, fun = ObjFBivarGARCH, dailyret = dailyret, covariance = covariance, 
                    ineqfun = ineqfun_GARCH, #the inequality constraint
                    ineqLB  = ineqLB, ## the inequality lower bound
                    ineqUB = ineqUB, ## the inequality lower bound, i.e. 0.0 < alpha + beta < 0.9999
                    ## lower and upper bounds for all parameters
                    LB = c(0.0001, 0.0001), UB = c(0.999, 0.999),
                    control = list(tol = tol, outer.iter = 800, inner.iter = 1200, delta = 1e-7)
                    ) 
  
  ## extract estimated parameters
  vPar = optimizer$pars
  
  ## extract the likelihood computed at its maximum
  dLLK = -tail(optimizer$values, 1)

  if(!bootstrap){
  #CALCULATING ROBUST STD. ERRORS WHEN NO COVARIANCE TARGETING:
  	  rse <- NULL
  	  se <- NULL
  	  scores <- NULL
	  if(is.null(vPar)){
		  scores <- matrix(0L, nrow=length(dailyret[,1]), ncol = length(vPar))

		  step <- 1e-5 * vPar
		  # This follows directly from Kevin Sheppard who uses percentage returns, Moreover scale open-to-close with 24/6.5.
		  for(i in 1:length(step)){

			h <- step[i]
		    delta <- rep(0, length(vPar))
		    delta[i] <- h
																		
			loglikeminus <- BivarGARCHFilter(dailyret, vPar-delta, covariance)$dLLKs
			loglikeplus <- BivarGARCHFilter(dailyret, vPar+delta, covariance,)$dLLKs

			scores[,i] <- (loglikeplus - loglikeminus)/(2*h)

		  }

		  J <- (t(scores) %*% scores)/length(dailyret[,1])

		  I <- optimizer$hessian/length(dailyret[,1])

		  I <- solve(I)[-1 ,-1]

		  vars <- (I * J * I)/length(dailyret[,1])
		  
		  rse <- sqrt(diag(vars))

		  ##Standard errors:
  		  se <- solve(optimizer$hessian)
  		  se <- matrix(sqrt(diag(se))[-1], ncol=length(vPar), nrow=1)/iT
  		  colnames(se) <- c("Alpha", "Beta")
	  }
  ## compute filtered volatility
  Filter <- BivarGARCHFilter(dailyret, c(vPar[1], vPar[2]), covariance)

  vSigma2 <- Filter$mSigma
  
  sampleH <- Filter$sampleH

  ## Compute the daily Average BIC
  iT = length(dailyret[,1])
  BIC = (-2 * dLLK + log(iT) * length(vPar))/iT

  ## return a list with estimated parameters, likelihood value and BIC
  lOut = list()
  lOut[["vPar"]] = vPar
  lOut[["dLLK"]] = dLLK
  lOut[["BIC"]] = BIC
  lOut[["vSigma2"]] = vSigma2
  lOut[["se"]] = se
  lOut[["Hessian"]] = optimizer$hessian
  lOut[["rse"]] = rse
  lOut[["sampleH"]] = sampleH
  lOut[["scores"]] = scores
  #lOut[["grad"]] = grad_Lk
    return(lOut)

  }
  
}





########################################################################################################################
#
#					CONTINUOUS ASSYMMETRICAL REALIZED SCALAR-BASED GARCH MODEL (crBG)
#
########################################################################################################################

BivarGARCHFilterContAsym <- function(dailyret, params, covariance){

	dAlphaP <- params[1]
	dAlphaN <- params[2]
	dAlphaM <- params[3]
	dBeta  <- params[4]

	days <- length(dailyret[,1])

	dailyret <-  matrix(dailyret, ncol=2, nrow=days)

	d <-  2 #being bivariate

	samplecov <- cov(dailyret)

	P <- covariance[[1]]
	N <- covariance[[2]]
	M <- covariance[[3]]

	sampleP <- matrix(
		c(mean(P[1,1,]),
		  mean(P[2,1,]),
		  mean(P[2,1,]),
		  mean(P[2,2,])),
		  ncol=2, nrow=2)

	sampleN <- matrix(
		c(mean(N[1,1,]),
		  mean(N[2,1,]),
		  mean(N[2,1,]),
		  mean(N[2,2,])),
		  ncol=2, nrow=2)

	sampleM <- matrix(
		c(mean(M[1,1,]),
		  mean(M[2,1,]),
		  mean(M[2,1,]),
		  mean(M[2,2,])),
		  ncol=2, nrow=2)

	mSigma <- array(0L, dim = c(2, 2, days))

	#initializing  with unconditional mean. 

	mSigma[,,1] <- samplecov

	#compute first observation of log-likelihood (we are minizing neg-log-likelihood):
	dLLK <-  log(det(mSigma[,,1])) + (dailyret[1, , drop = F]) %*% solve(mSigma[,,1]) %*% t(dailyret[1, , drop = F])

	id <- matrix(c(1,0,0,1), ncol = 2, nrow = 2)

	#astar <- (dAlphaP * sampleP + dAlphaN * sampleN + dAlphaM * sampleM) * samplecov^(-1)

	omega <- samplecov - (dAlphaP * sampleP) %*% id +  (dAlphaN * sampleN) %*% id + (dAlphaM * sampleM) %*% id - dBeta * samplecov

	dLLKs <- numeric()

	dLLKs[1] <-  dLLK

	for(i in 2:days){

		mSigma[,,i] <- omega  + dBeta * mSigma[,,i-1] + 
		dAlphaP * P[,,i-1] + dAlphaN * N[,,i-1] + dAlphaM * M[,,i-1]


		#if(!is.positive.definite(mSigma[,,i])){

		#	mSigma[,,i] <- matrix(nearPD(mSigma[,,i])$mat, ncol=2,nrow=2)

		#}

		#dLLK <- as.numeric(dLLK) + log(det(mSigma[,,i])) + 
		#dailyret[i, , drop = F] %*% solve(mSigma[,,i]) %*% t(dailyret[i, , drop = F]) 
		#full neg loglike
		
		dLLKs[i] <-  as.numeric(log(det(mSigma[,,i])) + 
		dailyret[i, , drop = F] %*% solve(mSigma[,,i]) %*% t(dailyret[i, , drop = F]))

		}

	#+ days * d * log(2*pi)
	fulldLLK <- - 0.5 * (d * log(2*pi) * days +  sum(dLLKs, na.rm=T))  

	lOut <- list()

	lOut[["fulldLLK"]] <- fulldLLK
	lOut[["mSigma"]] <- mSigma
	lOut[["dLLKs"]] <- dLLKs
	lOut[["cov"]] <- omega

	return(lOut)
	
	
}



ObjFBivarGARCHContAsym <- function(vPar, dailyret, covariance) {
  
  dAlphaP = vPar[1]
  dAlphaN = vPar[2]
  dAlphaM = vPar[3]
  dBeta  = vPar[4]
  dLLK = BivarGARCHFilterContAsym(dailyret, c(dAlphaP, dAlphaN, dAlphaM, dBeta), covariance)$fulldLLK
  
  return(-as.numeric(dLLK))
}


#imposing strong stationarity as implied from the paper. This is generalized and can be used for other imposed models.

ineqfun_GARCH_BIVARContAsym <- function(vPar, ...) {
  dAlpha = vPar[1] + vPar[2] + vPar[3]
  dBeta  = vPar[4]
  
  return(dAlpha + dBeta)
}




EstimateBivarGARCHContAsym <- function(dailyret, covariance, ineqfun_GARCH = ineqfun_GARCH_BIVARContAsym, ineqLB = 0.00, ineqUB = 0.9999, vPar = NULL){
  
  if(is.null(vPar)){
  dAlphaP = 0.09324267
  dAlphaN = 0.34553257
  dAlphaM = 0.04110211
  dBeta  = 0.52002253
  
  ## vector of starting parameters
  vPar = c(dAlphaP, dAlphaN, dAlphaM, dBeta)
  }
  # have a look at help(solnp)
  ##optimization step
  optimizer = solnp(vPar, fun = ObjFBivarGARCHContAsym, dailyret = dailyret, covariance = covariance, 
                    ineqfun = ineqfun_GARCH, #the inequality constraint
                    ineqLB  = ineqLB, ## the inequality lower bound
                    ineqUB = ineqUB, ## the inequality lower bound, i.e. 0.0 < alpha + beta < 0.9999
                    ## lower and upper bounds for all parameters
                    LB = c(0.0001, 0.0001, 0.0001, 0.0001), UB = c(0.999, 0.999, 0.999, 0.999)
                    ) #0.0001
  
  ## extract estimated parameters
  vPar = optimizer$pars
  
  ## extract the likelihood computed at its maximum
  dLLK = -tail(optimizer$values, 1)

  #CALCULATING ROBUST STD. ERRORS WHEN NO COVARIANCE TARGETING:
  se <- NULL
  rse <- NULL
  if(is.null(vPar)){
	  scores <- matrix(0L, nrow=length(dailyret[,1]), ncol = length(vPar))

	  step <- 1e-5 * vPar
	  # This follows directly from Kevin Sheppard who uses percentage returns, Moreover scale open-to-close with 24/6.5.
	  for(i in 1:length(step)){

		h <- step[i]
	    delta <- rep(0, length(vPar))
	    delta[i] <- h
																	
		loglikeminus <- BivarGARCHFilterContAsym(dailyret, vPar-delta, covariance)$dLLKs
		loglikeplus <- BivarGARCHFilterContAsym(dailyret, vPar+delta, covariance)$dLLKs

		scores[,i] <- (loglikeplus - loglikeminus)/(2*h)

	  }

	  J <- (t(scores) %*% scores)/2516

	  I <- optimizer$hessian/2516

	  I <- solve(I)[-1 ,-1]

	  vars <- (I * J * I)/2516
	  
	  rse <- sqrt(diag(vars))

	  ##Standard errors:
  	  se <- solve(optimizer$hessian)
  	  se <- matrix(sqrt(diag(se))[-1], ncol=length(vPar), nrow=1)
  	  colnames(se) <- c("AlphaP", "AlphaN", "AlphaM", "Beta")
  }
  
  ## compute filtered volatility
  vSigma2 = BivarGARCHFilterContAsym(dailyret, c(vPar[1], vPar[2], vPar[3], vPar[4]), covariance)$mSigma
  
  ## Compute the daily Average BIC
  iT = 2516
  BIC = (-2 * dLLK + log(iT) * length(vPar))/iT

  ## return a list with estimated parameters, likelihood value and BIC
  lOut = list()
  lOut[["vPar"]] = vPar
  lOut[["dLLK"]] = dLLK
  lOut[["BIC"]] = BIC
  lOut[["vSigma2"]] = vSigma2
  lOut[["se"]] = se
  lOut[["Hessian"]] = optimizer$hessian
  lOut[["rse"]] = rse
  #lOut[["grad"]] = grad_Lk
  
  return(lOut)
}









########################################################################################################################
#
#						REALIZED DCC MODEL (rDCC)
#
########################################################################################################################

rDCCFilter <- function(mEta, dA, dB, mQ, covariance) {
  
  iN <- ncol(mEta)
  iT <- nrow(mEta)
  
  # initialize the array for the correlations
  aCor <- array(0, dim = c(iN, iN, iT))
  # initialize the array for the Q matrices
  aQ <- array(0, dim = c(iN, iN, iT))
  
  #compute the realized correlations

  rcor <- array(0L, dim = c(iN, iN, iT))

  for(i in 1:iT){

  	rcor[,,i] <- diag(sqrt(1/diag(covariance[,,i]))) %*% covariance[,,i] %*% diag(sqrt(1/diag(covariance[,,i])))

  }
  ## initialization at the unconditional cor
  aCor[,, 1] <- mQ
  aQ[,,1] <- mQ
  
  #Compute the first likelihood contribution
  dLLK <- mEta[1, , drop = FALSE] %*% solve(aCor[,, 1]) %*% t(mEta[1, , drop = FALSE]) - 
    mEta[1, , drop = FALSE]%*% t(mEta[1, , drop = FALSE]) + log(det(aCor[,, 1]))
  
  #main loop
  for (t in 2:iT) {
    #update the Q matrix
    aQ[,, t] <- mQ * (1 - dA - dB) + dA * rcor[,,t - 1] + dB * aQ[,,t - 1]
    
    ## Compute the correlation as Q_tilde^{-1/2} Q Q_tilde^{-1/2} 
    aCor[,, t] <- diag(sqrt(1/diag(aQ[,, t]))) %*% aQ[,, t] %*% diag(sqrt(1/diag(aQ[,, t]))) 
    
    #augment the likelihood value
    dLLK <- dLLK + mEta[t, , drop = FALSE] %*% solve(aCor[,, t]) %*% t(mEta[t, , drop = FALSE]) - 
      mEta[t, , drop = FALSE] %*% t(mEta[t, , drop = FALSE]) + log(det(aCor[,, t]))
  }
  
  lOut <- list()
  #remember to include the -1/2 term in the likelihood evaluation
  #see the equations in the corresponding lecture
  lOut[["dLLK"]] <- -0.5 * dLLK
  lOut[["aCor"]] <- aCor
  lOut[["rcor"]] <- rcor
  
  return(lOut)
}


#testing with simple realGARCH(1,1) using rugarch package. 
#you need getDates for fit function in rugarch to work. 
#moreover the covariances you use in the DCC model is also used in the univariate models.
Estimate_rDCC <- function(mY, covariance, getDates, bootstrap = F, residuals = NULL, tol = 1e-7, vPar = NULL) {
  
  ## estimate the marginal models
  require(Rsolnp)
  require(rugarch) 
 
  #Marginal garch specification. THIS WORKS ONLY IN BIVARIATE SETUP. 
  
  #list where marginal models are stored
  #to get covariances from percentage log-returns you need to multiply by 10000 (DONE OUTSIDE ESTIMATOR)
  if(!bootstrap){
    cov1 <- (covariance[1,1,])
    cov2 <- (covariance[2,2,]) 

    spec <- ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), variance.model = 
  	list(model = 'realGARCH', garchOrder = c(2, 1)))

    setbounds(spec)<-list(alpha2=c(-1,1))


    spec1 <- ugarchfit(spec, mY[,1], solver = 'hybrid', realizedVol = 
    xts(cov1, order.by = as.Date(getDates)), 
    fit.control = list(stationarity = 1, fixed.se = 1))

    spec2 <- ugarchfit(spec, mY[,2], solver = 'hybrid', realizedVol = 
    xts(cov2, order.by = as.Date(getDates)), 
    fit.control = list(stationarity = 1, fixed.se = 1))

    #mspec <- multispec( replicate(spec, n=2) )
    						 #only used percentage log-returns, since it seems to converge better. 
    #lfit <- multifit(mspec, mY, solver = 'hybrid', realizedVol = 
    #xts(cbind(cov1, cov2), order.by = as.Date(getDates)), fit.control = list(stationarity = 1, fixed.se = 1))

    res1 <- residuals(spec1, standardize = T)
    res2 <- residuals(spec2, standardize = T)

    mEta <- cbind(res1, res2) #residuals(lfit, standardize = T)
  }else{mEta <- residuals}
  ####################################################
  
  ## maximization of the DCC likelihood
  
  #initial parameters
  if(is.null(vPar)){
  	vPar =  c(0.1974634, 0.8015366)
  }
  #unconditional correlation
  mQ = cor(mEta)
  
  #maximize the DCC likelihood
  optimizer = solnp(vPar, fun = function(vPar, mEta, mQ, covariance) {
    
    Filter = rDCCFilter(mEta, vPar[1], vPar[2], mQ, covariance)
    dNLLK = -as.numeric(Filter$dLLK)
    return(dNLLK)
    
  }, ineqfun = function(vPar, ...) {
    sum(vPar)
  }, ineqLB = 1e-4, ineqUB = 0.999, 
  LB = c(1e-4, 1e-4), UB = c(0.999, 0.999), 
  mEta = mEta, mQ = mQ, covariance = covariance, 
  control = list(tol = tol, outer.iter = 800, inner.iter = 1200, delta = 1e-7))
  
  #Extract the estimated parameters
  vPar = optimizer$pars
  
  #Extract the likelihood of the correlation part
  dLLK_C = -tail(optimizer$values, 1)
  
  #Filter the dynamic correlation using the estimated parameters
  Filter = rDCCFilter(mEta, vPar[1], vPar[2], mQ, covariance)

  if(bootstrap){

    return(vPar)
  
  }

  #extract univariate variances
  mSigma = cbind(sigma(spec1)^2, sigma(spec2)^2)
  


  #extract univariate estimated parameters
  mCoef = cbind(coef(spec1), coef(spec2))

  colnames(mCoef) <- colnames(mY)
  
  #compute the likelihood of the volatility  part
  dLLK_V = sum(-spec1@fit$partial.log.likelihoods) +sum(-spec2@fit$partial.log.likelihoods)
  
  dLLK_V2 = sum(likelihood(spec1)) + sum(likelihood(spec2))

  #compute the total likelihood
  dLLK = dLLK_V2 + dLLK_C
  
  ## Compute z_t
  aCor = Filter[["aCor"]]
  covs = sigma(spec1) * sigma(spec2) *  aCor[2,1,]
  

  iT = nrow(mY)
  mZ <- NULL
  se <- NULL
  if(is.null(vPar)){
  mZ = matrix(0, iT, ncol(mY))
  
  for (t in 1:iT) {
    mZ[t, ] = solve(chol(aCor[,,t])) %*% as.numeric(mEta[t, ])
  }
    
  #standard errors 
  se <- solve(optimizer$hessian)
  se <- matrix(sqrt(diag(se))[-1], ncol=length(vPar), nrow=1)
  }


  #compute estimated covariances:
  vSigma2 <- array(0L, dim = c(2,2, 2516))

  for(i in 1:length(mY[,1])){

    vSigma2[,,i] <- matrix(c(mSigma[i,1], covs[i], covs[i], mSigma[i,2]))

  }

  lOut = list()

   ## Compute the daily Average BIC
  iT = 2516
  BIC = (-2 * dLLK + log(iT) * (length(vPar) + length(coef(spec1)) + length(coef(spec2))))/iT
  
  #output the results
  lOut[["dLLK"]] = dLLK
  lOut[["dLLKV"]] = dLLK_V
  lOut[["BIC"]] = BIC
  lOut[["mCoef"]] = mCoef
  lOut[["vPar"]] = vPar
  lOut[["mSigma"]] = mSigma
  lOut[["aCor"]] = aCor
  lOut[["mEta"]] = mEta
  lOut[["mZ"]] = mZ
  lOut[["se"]] = se
  lOut[["dLLK_C"]] = dLLK_C
  lOut[["vSigma2"]] = vSigma2
  lOut[["covs"]] = covs

  return(lOut)
  
}




########################################################################################################################
#
#					CONTINUOUS ASSYMMETRICAL REALIZED DCC MODEL (crDCC)
#
########################################################################################################################

rcDCCFilter <- function(mEta, dAP, dAN, dAM , dB, mQ, covariance) {
  
  mEta <- as.matrix(mEta)

  iN <- ncol(mEta)
  iT <- nrow(mEta)
  
  # initialize the array for the correlations
  aCor <- array(0, dim = c(iN, iN, iT))
  # initialize the array for the Q matrices
  aQ <- array(0, dim = c(iN, iN, iT))
  
  #compute the realized correlations

  P <- covariance[[1]]
  N <- covariance[[2]]
  M <- covariance[[3]]

  Pcor <- array(0L, dim = c(2,2,iT))
  Ncor <- array(0L, dim = c(2,2,iT))
  Mcor <- array(0L, dim = c(2,2,iT))


  for(i in 1:iT){

  	Pcor[,,i] <- diag(sqrt(1/diag(P[,,i]))) %*% P[,,i] %*% diag(sqrt(1/diag(P[,,i])))
  	Ncor[,,i] <- diag(sqrt(1/diag(N[,,i]))) %*% N[,,i] %*% diag(sqrt(1/diag(N[,,i])))
  	Mcor[,,i] <- diag(sqrt(1/diag(N[,,i] + P[,,i] + M[,,i]))) %*% M[,,i] %*% diag(sqrt(1/diag(N[,,i] + P[,,i] + M[,,i]))) + diag(iN)


  	#ad-hog method for eliminating NaNs
  	if(is.nan(Pcor[2,1,i])){ Pcor[,,i] <- Pcor[,,i-1]}
  	if(is.nan(Ncor[2,1,i])){ Ncor[,,i] <- Ncor[,,i-1]}
  	if(is.nan(Mcor[2,1,i])){ Mcor[,,i] <- Mcor[,,i-1]}



  }

  sampleP <- matrix(
	c(mean(Pcor[1,1,], na.rm = T),
	  mean(Pcor[2,1,], na.rm = T),
	  mean(Pcor[2,1,], na.rm = T),
	  mean(Pcor[2,2,], na.rm = T)),
	  ncol=2, nrow=2)

  sampleN <- matrix(
	c(mean(Ncor[1,1,], na.rm = T),
	  mean(Ncor[2,1,], na.rm = T),
	  mean(Ncor[2,1,], na.rm = T),
	  mean(Ncor[2,2,], na.rm = T)),
	  ncol=2, nrow=2)

  sampleM <- matrix(
	c(mean(Mcor[1,1,], na.rm = T),
	  mean(Mcor[2,1,], na.rm = T),
	  mean(Mcor[2,1,], na.rm = T),
	  mean(Mcor[2,2,], na.rm = T)),
	  ncol=2, nrow=2)

  ## initialization at the unconditional cor
  aCor[,, 1] <- mQ
  aQ[,,1] <- mQ
  
  #Compute the first likelihood contribution
  dLLK <- mEta[1, , drop = FALSE] %*% solve(aCor[,, 1]) %*% t(mEta[1, , drop = FALSE]) - 
    mEta[1, , drop = FALSE]%*% t(mEta[1, , drop = FALSE]) + log(det(aCor[,, 1]))
  

  id <- matrix(c(1,0,0,1), ncol = 2, nrow = 2)


  omega <- mQ - ((dAP * sampleP)  + (dAN * sampleN)  + (dAM * sampleM)) %*% mQ %*% solve(mQ)  - dB * mQ

  #astar <- (dAP * sampleP + dAN * sampleN  + dAM * sampleM) %*% solve(mQ)
  #main loop
  for (t in 2:iT) {
    #update the Q matrix
    aQ[,, t] <- omega + dAP * Pcor[,,t - 1] + dAN * Ncor[,,t - 1]  + dAM * Mcor[,,t - 1]  + dB * aQ[,,t - 1]
    
    ## Compute the correlation as Q_tilde^{-1/2} Q Q_tilde^{-1/2} 
    aCor[,, t] <- diag(sqrt(1/diag(aQ[,, t]))) %*% aQ[,, t] %*% diag(sqrt(1/diag(aQ[,, t]))) 
    
    #augment the likelihood value
    dLLK <- dLLK + mEta[t, , drop = FALSE] %*% solve(aCor[,, t]) %*% t(mEta[t, , drop = FALSE]) - 
      mEta[t, , drop = FALSE] %*% t(mEta[t, , drop = FALSE]) + log(det(aCor[,, t]))
  }
  
  lOut <- list()
  #remember to include the -1/2 term in the likelihood evaluation
  #see the equations in the corresponding lecture
  lOut[["dLLK"]] <- -0.5 * dLLK
  lOut[["aCor"]] <- aCor
  lOut[["astar"]] <- Pcor
  lOut[["omega"]] <- omega
  
  return(lOut)
}



Estimate_rcDCC <- function(mY, covariance, cov2 = NULL, getDates, bootstrap = F, residuals = NULL, vPar = NULL, tol = 1e-7){
  
  ## estimate the marginal models
  require(Rsolnp)
  require(rugarch) 
 
  #Marginal garch specification. THIS WORKS ONLY IN BIVARIATE SETUP. 
  if(!bootstrap){
  	
    #theres some minor inconsistencies (e-15) and thus I used true covariance instead.
    #inconsistencies is enough to make the likelihood rough and thus end in a hole -133000. 
    if(is.null(cov2)){
      cov2 <- (covariance[[1]] + covariance[[2]] + covariance[[3]]) 
    }
  	#list where marginal models are stored
  	spec <- ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), variance.model = 
	list(model = 'realGARCH', garchOrder = c(2, 1)))
  	
    setbounds(spec)<-list(alpha2=c(-1,1))


    specforrobust1 <- ugarchfit(spec, mY[,1], solver = 'hybrid', realizedVol = 
  xts((cov2[1,1,]), order.by = as.Date(getDates)), fit.control = list(stationarity = 1, fixed.se = 0))

  	#specforrobust1 <- ugarchfit(spec, mY[ ,1], solver = 'solnp', realizedVol = 
	#xts(sqrt(cov[1,1,]), order.by = as.Date(getDates)), fit.control = list(stationarity = 1, 
    #fixed.se = 0))

  	specforrobust2 <- ugarchfit(spec, mY[ ,2], solver = 'hybrid', realizedVol = 
	xts((cov2[2,2,]), order.by = as.Date(getDates)), fit.control = list(stationarity = 1, fixed.se = 0))

  	#mspec <- multispec( replicate(spec, n=2) )

  	#theres some issues with the multifit procedure regarding parameter estimations
  	#lfit <- multifit(mspec, mY , solver = 'hybrid', realizedVol = 
  	#xts(cbind(sqrt(cov[1,1,]), sqrt(cov[2,2,])), order.by = as.Date(getDates), fit.control = list(stationarity = 1, fixed.se = 1)))
  
  	mEta <- cbind(residuals(specforrobust1, standardize = T), residuals(specforrobust2, standardize = T))   #residuals(lfit, standardize = T)
  }
  else{mEta <- residuals}
  ####################################################
  
  ## maximization of the DCC likelihood
  
  #initial parameters
  if(is.null(vPar)){
  	vPar = c(0.00000001384, 0.00000013727, 0.08197776333, 0.91702207742)
  }
  #unconditional correlation
  mQ = cor(mEta)
  
  #maximize the DCC likelihood
  optimizer = solnp(vPar, fun = function(vPar, mEta, mQ, covariance) {
    
    Filter = rcDCCFilter(mEta, vPar[1], vPar[2], vPar[3], vPar[4], mQ, covariance)
    dNLLK = -as.numeric(Filter$dLLK)
    return(dNLLK)
    
  }, ineqfun = function(vPar, ...) {
    sum(vPar)
  }, ineqLB = 1e-6, ineqUB = 0.999, 
  LB = c(0, 0, 0, 0), UB = c(0.999, 0.999, 0.999, 0.999), 
  mEta = mEta, mQ = mQ, covariance = covariance,
  control = list(tol = tol, outer.iter = 800, inner.iter = 1200, delta = 1e-7))
  
  #Extract the estimated parameters
  vPar = optimizer$pars
  
  if(bootstrap){return(vPar)}

  #Extract the likelihood of the correlation part
  dLLK_C = -tail(optimizer$values, 1)
  
  #Filter the dynamic correlation using the estimated parameters
  Filter = rcDCCFilter(mEta, vPar[1], vPar[2], vPar[3], vPar[4], mQ, covariance)

  #extract univariate volatilities
  mSigma = cbind(sigma(specforrobust1)^2, sigma(specforrobust2)^2)
  
  #extract univariate estimated parameters
  mCoef = cbind(coef(specforrobust1), coef(specforrobust2))

  colnames(mCoef) <- colnames(mY)
  
  #compute the likelihood of the volatility  part

  dLLK_V = sum(-specforrobust1@fit$partial.log.likelihoods) +sum(-specforrobust2@fit$partial.log.likelihoods)
  dLLK_V2 <- sum(likelihood(specforrobust1)) + sum(likelihood(specforrobust2))
  
  #compute the total likelihood
  dLLK = dLLK_V2 + dLLK_C
  
  ## Compute z_t
  aCor = Filter[["aCor"]]
  covs = sigma(specforrobust1) * sigma(specforrobust2) *  aCor[2,1,]
  

  iT = nrow(mY)
  mZ <- NULL
  if(is.null(vPar)){

	  mZ = matrix(0, iT, ncol(mY))
	  
	  for (t in 1:iT) {
	    mZ[t, ] = solve(chol(aCor[,,t])) %*% as.numeric(mEta[t, ])
	  }
    
  }
  #compute estimated covariances:
  vSigma2 <- array(0L, dim = c(2,2, 2516))

  for(i in 1:length(mY[,1])){

    vSigma2[,,i] <- matrix(c(mSigma[i,1], covs[i], covs[i], mSigma[i,2]))

  }

  lOut = list()
  allpars <- length(vPar) + length(coef(specforrobust1)) + (coef(specforrobust2))
   ## Compute the daily Average BIC
  iT = 2516
  BIC = (-2 * dLLK + log(iT) * (length(vPar) + length(coef(specforrobust1)) + length(coef(specforrobust2))))/iT
  lOut = list()

  #output the results
  lOut[["dLLK"]] = dLLK #can become unstable due to rugarch loglikes. should be around -4555 ish. 
  lOut[["mCoef"]] = mCoef
  lOut[["vPar"]] = vPar
  lOut[["mSigma"]] = mSigma
  lOut[["aCor"]] = aCor
  lOut[["mEta"]] = mEta
  lOut[["mZ"]] = mZ
  lOut[["seeall"]] = list(specforrobust1, specforrobust2) 
  lOut[["vSigma2"]] = vSigma2
  lOut[["BIC"]] = BIC
  lOut[["omega"]] = Filter[["omega"]]
  return(lOut)
  
}
