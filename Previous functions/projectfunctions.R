#                                                                                                  
#-------------------------------HIGH FREQUENCY PROJECT------------------------------------------            
# DESCRIPTION: This R script contains self-constructed functions used in the high frequency project.
# Some of the functions are inspired by other packages (highfrequency), but then tweaked to fit particular 
#data types or to do something else entirely. 
#
#
#
#
#------------------------------------REALIZED ESTIMATORS---------------------------------------

#' realCov
#'
#' Computes realized Covariance.
#' @param data of type xts.
#' @param correlation boolean value. 
#' @return realized estimates.
#' @import xts stats matlib
#' @export
realCov <- function(matrix, correlation = FALSE){
  # Realized Covariance estimator for 1 day. 
  #
  # Args: 
  #
  #   matrix: A matrix of returns. Nxk 
  #
  #	  N: number of intraday returns, k: number of asset dimension, T: number of days. 
  realcovariance <- matrix()
  	
  realcovariance <- t(matrix) %*% (matrix)

  if(correlation){

		d <- sqrt(diag(realcovariance))
		d <- diag(d, 2, 2)

		corr <- inv(d) %*% realcovariance %*% inv(d) 

		return(corr)
	}
  
 return( realcovariance ) 
}
#' preavCov
#'
#' Computes preaveraged realized Covariance.
#' @param data of type xts.
#' @param rescaled boolean value. Rescales the estimate as described in Christensen & Podolskij
#' @param biascorrected boolean value. Computes the bias-corrected estimate. 
#' @param correlation boolean value. 
#' @param theta parameter. 
#' @return realized estimates.
#' @import xts stats matlib
#' @export
preavCov <- function(matrix, rescaled = FALSE, biascorrected = FALSE, correlation = FALSE ,theta = 0.8){

	#require matlib for inverse function 

	psi2 <- 1 / 12

	n <- nrow(matrix)

	kn <- floor( theta * sqrt(n) ) 

	.weight <- function(x){ 

		ifelse( x > (1-x), (1-x), x)

	}

	psi1_kn <- kn*sum((.weight( (1:kn) / kn ) - .weight( ((1:kn) - 1) / kn))^2 )

	psi2_kn <- 1/kn * sum( (.weight( (1:kn) / kn))^2  )

	.preaveraging <- function(matrix, kn){

		kn <- as.numeric(kn)

		if( kn == 1 ){preavg = matrix} 
		else{

			weightedsum <- function(mat){sum( .weight( 1:(kn - 1) / kn ) * mat )}

			preavg <- rollapply(matrix, width = kn-1, weightedsum, align = "left")

			#refills the nan values in preavg with the values in  matrix.
			preavg[is.na(preavg[, 1]), ] = matrix[is.na(preavg[, 1])]
		}
		return(preavg)
	}


	#.univar_est <- function(matrix){

		#matrix <- as.matrix(as.numeric(matrix))

	#	r1 <- .preaveraging(matrix, kn=kn)

	#	est <- 1 / (sqrt(n)*theta*psi2_kn) * sum(r1^2, na.rm=TRUE) - psi1_kn * (1/n) / (2*theta^2*psi2_kn) * sum(matrix^2, na.rm=TRUE)

	 #return(est)
	#}

	.univar_est <- function(matrix){
		#Follows Christensen, Oomen (2014). 

		r1 <- .preaveraging(matrix, kn=kn)

		omega <- 1/(2*n) * realCov(matrix)

		constant1 <- n/(n-kn+2)

		constant2 <- 1/(kn*psi2_kn)

		bias <- psi1_kn/(theta^2 * psi2_kn) * omega^2

		covar <- constant1 * constant2 * realCov(r1) - bias

		return(covar)

	}

	if(ncol(matrix) == 1){

		MRC <- .univar_est(matrix)

		return(MRC)
	}
	else{


		.bivar_est <- function(matrix){

			#matrix[, i] <- as.matrix(as.numeric(matrix[, i]))
			#matrix[, j] <- as.matrix(as.numeric(matrix[, j]))

			r1 <- .preaveraging(matrix[, i], kn=kn)
			r2 <- .preaveraging(matrix[, j], kn=kn)

			if(biascorrected){

				bi_MRC <- n / (n-kn+2) * 1 / (psi2*kn) * (sum(t(r1)%*%r2, na.rm=TRUE)) - (psi1_kn / (theta^2 * psi2_kn)) * 1/(2*n)*(sum(t(matrix[, i])%*%matrix[, j], na.rm=TRUE))

				return(bi_MRC)
			}

			bi_MRC <- n / (n-kn+2)* 1 / (psi2*kn) * (sum(t(r1)%*%r2, na.rm=TRUE))

		 return(bi_MRC)
		}

		N <- ncol(matrix)

		covar <- matrix(rep(0, N^2), ncol = N)
		for (i in 1:N){

			diag(covar)[i] <- .univar_est(matrix[, i])

		}
		
		for (i in 2:N){
			for (j in 1:(i - 1)){

				covar[i, j] <- covar[j, i] <- .bivar_est(matrix)

			}
		}
		if(rescaled){

			covar <-  1 / (1 - (psi1_kn/(theta^2*psi2_kn)) * 1/(2*n) )  * covar

		}
	}
	if(correlation){

		d <- sqrt(diag(covar))
		d <- diag(d, 2, 2)

		corr <- inv(d) %*% covar %*% inv(d) 

		return(corr)
	}

	return(covar)
}

#' preavBPCov (Bias-corrected by default)
#'
#' Computes (preaveraged) bipower realized Covariance.
#' @param matrix of type xts.
#' @param preaveraging of type boolean. Tells you if you want to pre-average the returns. 
#' @param rescaled boolean value. Rescales the estimate as described in Christensen & Podolskij
#' @param correlation boolean value. 
#' @param theta parameter. 
#' @return realized estimates.
#' @import xts stats matlib
#' @export
preavBPCOV <- function(matrix, preaveraging = FALSE, correlation = FALSE, rescaled = FALSE, theta = 0.8){

	n <- as.numeric(nrow(matrix))

	N <- ncol(matrix)

	kn <- floor( theta * sqrt(n) )

	wk <- (1+2*kn^(-2)) / 12

	.weight <- function(x){ 

		ifelse( x > (1-x), (1-x), x)

	}

	psi2 <- 1 / 12

	
	psi1_kn <- kn*sum((.weight( (1:kn) / kn ) - .weight( ((1:kn) - 1) / kn))^2 )

	psi2_kn <- 1/kn * sum( (.weight( (1:kn) / kn))^2  )
	  
	
	.preaveraging <- function(matrix, kn){

		kn <- as.numeric(kn)
		if( kn == 1 ){preavg = matrix} 
		else{

			weightedsum <- function(mat){sum( .weight( 1:(kn - 1) / kn ) * mat )}

			preavg <- rollapply(matrix, width = kn-1, weightedsum, align = "left")

			#refills the nan values in preavg with the values in  matrix.
			preavg[is.na(preavg[, 1]), ] = matrix[is.na(preavg[, 1])]
		}
		return(preavg)
	}


	.uniBPCOV <- function(vector){

	
		if(preaveraging){

			r1 <- .preaveraging(vector, kn=kn)

			vector <- as.numeric(vector)
			
			ret <- matrix(vector, nrow = n, ncol = 1)

			omega <- 1/(2*n) * realCov(ret)
			
			#bpvar:

			mu <- sqrt(2/pi)

			constant1 <- n/(n-2*kn+2)

			constant2 <- 1/(kn*psi2_kn*mu^2)

			bias <- psi1_kn/(theta^2 * psi2_kn) * omega^2

			lead <- as.vector(r1[kn:n])

			lag <- as.vector(r1[1:(n-kn+1)])

			bpvar <- constant1 * constant2 * sum(abs(lead * lag)) - bias

		return(bpvar)

		}
		else{

			vector <- as.vector(as.numeric(vector))

			var <- (pi/2) * sum( abs( vector[1:(n-1)]) * abs(vector[2:n]) )
		}
	
	return(var)

	}

	if(ncol(matrix) == 1){


		covar <- .uniBPCOV(matrix)

		return(covar)

	}
	else{
		.bivarBPCOV <- function(matrix){

			if(preaveraging){

				r1 <- .preaveraging(matrix[, i], kn=kn)
				r2 <- .preaveraging(matrix[, j], kn=kn)


				add <- as.vector(as.numeric((r1 + r2)))

				subtract <- as.vector(as.numeric((r1 - r2)))

				omega <- 1/(2*n) * t(matrix[, i]) %*% matrix[, j]

				mu <- sqrt(8/pi)

				const1 <- n/(n-2*kn+2)

				const2 <- 1/(kn * psi2_kn * mu^2)

				biascorr <- psi1_kn / (theta^2 * psi2_kn) * omega^2

				covar <- const1 * const2 * sum(  (abs( add[1:(n-kn+1)] ) * abs( add[kn:n] )) - ( abs(subtract[1:(n-kn+1)]) * abs(subtract[kn:n]) ), na.rm = TRUE  ) - biascorr   

			}
			else{

				add <- as.numeric(matrix[, i] + matrix[, j])

				subtract <- as.numeric(matrix[, i] - matrix[, j])

				covar <- (pi/8) * sum(  (abs( add[1:(n-1)] * add[2:n] ) ) - ( abs( subtract[1:(n-1)] * subtract[2:n] ) )  )     
			}

		return(covar)

		}

		covar <- matrix(rep(0, N^2), ncol = N)

		for (i in 1:N){

			diag(covar)[i] <- .uniBPCOV(matrix[, i])

		}
		
		for (i in 2:N){
			for (j in 1:(i - 1)){

				covar[i, j] <- covar[j, i] <- .bivarBPCOV(matrix)

			}
		}
	}
	if(rescaled){

			covar <-  1 / (1 - (psi1_kn/(theta^2*psi2_kn)) * 1/(2*n) )  * covar

		}
	if(correlation){

		d <- sqrt(diag(covar))
		d <- diag(d, 2, 2)

		corr <- inv(d) %*% covar %*% inv(d) 

		return(corr)
	}
	return(covar)
}
#' preavthrCov (Bias-corrected by default)
#'
#' Computes (preaveraged) threshold realized Covariance.
#' @param matrix of type xts.
#' @param preaveraging of type boolean. Tells you if you want to pre-average the returns. 
#' @param getThreshold of boolean value. Prints the threshold values. 
#' @param theta parameter. 
#' @return realized estimates.
#' @import xts stats matlib
#' @export
preavthrCOV <- function(matrix, preaveraging = FALSE, theta = 0.8, getThreshold = FALSE){

	n <- nrow(matrix)

	kn <- floor( theta * sqrt(n) )

	.weight <- function(x){ 

		ifelse( x > (1-x), (1-x), x)

	}

	.preaveraging <- function(matrix, kn){

		kn <- as.numeric(kn)
		if( kn == 1 ){preavg = matrix} 
		else{

			weightedsum <- function(mat){sum( .weight( 1:(kn - 1) / kn ) * mat )}

			preavg <- rollapply(matrix, width = kn-1, weightedsum, align = "left")

			#refills the nan values in preavg with the values in  matrix.
			preavg[is.na(preavg[, 1]), ] = matrix[is.na(preavg[, 1])]
		}
		return(preavg)
	}

	.uniBPCOV <- function(vector){

	  n <- nrow(matrix)
	  
		vector <- as.vector(as.numeric(vector))

		var <- (pi/2) * sum( abs( vector[1:(n-1)]) * abs(vector[2:n]) )
		
	
	return(var)

	}

	delta <- 1/n

	threshold <- matrix(nrow = 1, ncol = ncol(matrix))

	for (i in 1:ncol(matrix)){

		threshold[, i] <- .uniBPCOV(matrix[, i])

	}

	#Following Jacod and Todorov (2009) p. 1814
	threshold <- 3*sqrt(threshold) * (delta^0.49)

	#Follows slides in Adv. fin econ:
	thresholdslides <- qnorm(0.999, mean = 0, sd = 1) * threshold * (delta^0.49)

	thresholdmat <- matrix(rep(threshold, n), ncol=ncol(threshold), nrow=n, byrow = TRUE)

	if(getThreshold){
	  
	  thresh <- matrix
	  thresh[ abs(thresh) > thresholdmat ] = 1e10
	  
	  thresh1 <-  table(thresh[,1])
	  thresh2 <- table(thresh[,2])
	  
	  thresh1 <- thresh1[names(thresh1) == 1e10]
	  thresh2 <- thresh2[names(thresh2) == 1e10]
	  
	  result <- matrix(c(thresh1,thresh2), nrow=1, ncol=2, byrow=T)
	  colnames(result) <- c("TLT", "SPY")
	  
	  return(result)
	  
	}
	
	
	
	#indicator funcion
	matrix[ abs(matrix) > thresholdmat ] = 0
	
	if(preaveraging){

		covar <- preavCov(matrix, FALSE, FALSE, theta = theta)

		return(covar)

	}
	else{

		covar <- realCov(matrix)

		return(covar)
	}
}
#' MLRV (Univariate Maximum Likelihood Realized Variance)
#'
#' Computes the Univariate Maximum Likelihood Realized Variance.
#' @param list of type list with list elements of type xts.
#' @param names string vector of names for each asset. 
#' @return realized estimates.
#' @import xts stats matlib
#' @export
MLRV <- function(list, names = NULL){
	#names: a vector of names of object string. 

	n <- sapply(list, length)

	.getMA <- function(list){
		mas <- list()

		for (i in 1:length(list)){

			mas[[i]] <- arima(list[[i]], order = c(0, 0, 1))

		}

		coef <- matrix(0L, nrow=length(mas), ncol = 2)

		for (j in 1:length(mas)){

			coef[j, 1] <- -mas[[j]]$model$theta
			coef[j, 2] <- mas[[j]]$sigma2

		}
		
		return(coef)
	}

	var <- matrix(0L, nrow = length(list), ncol = 1)

	args <- .getMA(list)

	for (i in 1:length(list)){

		var[i, 1] <- (n[i] * args[i, 2] * (1+args[i, 1])^2) 

	}

	rownames(var) <- names
	colnames(var) <- "MLRV"

	return(var)
}


#--------------------------------------------------------------------------------------------


#' noisetosignal (Oomen (2006))
#'
#' Computes the Univariate Maximum Likelihood Realized Variance.
#' @param list of type list with list elements of type xts.
#' @param names string vector of names for each asset. 
#' @return noise ratios.
#' @import xts stats matlib
#' @export
noisetosignal <- function(list, names = NULL){

	n <- as.vector(sapply(list, length))

	IVest <- MLRV(list, names)

	#a bar calculations
	condition <- as.numeric(vector())

	#sum only takes numeric vectors, therefore saving in list object.
	ret <- list() 
	for ( i in 1:length(list)){

		ret[[i]] <- as.vector(list[[i]])

	}

	for (i in 1:length(list)){

	condition[i] <- t(ret[[i]][2:n[i]]) %*% ret[[i]][1:(n[i]-1)]

	}

	abar <- matrix(0L, nrow = length(list), ncol=1)

	for(i in 1:length(list)){

		if(condition[i] < 0){

			abar[i, 1] <-  -1/(n[i]-1) * t(ret[[i]][1:(n[i]-1)]) %*% ret[[i]][2:n[i]]
		}
		else{

			abar[i, 1] <-  (1/(2*n[i]) * realCov(ret[[i]]))

		}

	}

	#eps <- abar %*% 1/(IVest %*% 1/n)
	eps <-  abar / (IVest / n)


	return(eps)

}


#' datesToPOSIXct
#'
#' Converges dates to POSIXct values used as index in multivariate timeseries, xts.
#' @param dates of type dates.
#' @return dates as POSIXct.
#' @import xts stats matlib
#' @export
datesToPOSIXct <- function(dates) {

    for (i in 1:length(dates)){

        dates[i] <- gsub('^([0-9]{4})([0-9]+)$', '\\1-\\2', dates[i])
        dates[i] <- gsub('([[:digit:]]{2,2})$', '-\\1\\2', dates[i])

    }

    dates <- as.POSIXct((dates))

    return(dates)
}


#' do.call.rbind
#'
#' Merges all of the days into one big xts.
#' @param list of type list with list elements being xts.
#' @return one big xts of dim (length(list) * nrow(list[[i]]))  x 1, for all i. 
#' @import xts stats matlib
#' @export 
do.call.rbind <- function(lst) {
  while(length(lst) > 1) {
    idxlst <- seq(from=1, to=length(lst), by=2)
    lst <- lapply(idxlst, function(i) {
      if(i==length(lst)) { return(lst[[i]]) }
      return(rbind(lst[[i]], lst[[i+1]]))
    })
  }
  lst[[1]]
}
#' riskparity_2dim
#' 
#' Calculates the riskparity portfolio utilizing a specific realized measure.
#' @param matrix of type xts.
#' @param risktarget double. Specify risk target. 
#' @param measure of type function. 
#' @return weights for riskparity portfolio in two dimensions.
#' @import xts stats matlib
#' @export
riskparity_2dim <- function(matrix, risktarget, measure){
  
  bonds <- matrix[,1]
  stocks <- matrix[,2]
  
  w_1 <- sqrt(measure(bonds))^-1 / (sqrt(measure(stocks))^-1 + sqrt(measure(bonds))^-1)
  w_2 <- sqrt(measure(stocks))^-1 / (sqrt(measure(stocks))^-1 + sqrt(measure(bonds))^-1)

  w <- matrix(c(w_1, w_2), ncol=1, nrow=2) #dimnames = list(c(), c("Bond", "Stock"))
  
  portrisk <- sqrt(t(w) %*% measure(matrix) %*% w) 
  
  alpha <- risktarget / portrisk
  
  w_new <- w %*% alpha
  
  return(w_new)
  
}

#' brownian
#' 
#' Simulates a correlated bivariate square-root process.
#' @param J integer. Specify paths. 
#' @param N integer. Specify mesh. 
#' @param sigma integer. Specify volatility.
#' @param alpha integer. specify drift.
#' @param beta integer.
#' @param gamma integer.
#' @return simulations of the bivariate square-root process. 
#' @import xts stats matlib
#' @export
brownian <- function(J, N, sigma, alpha, beta, gamma){ 
  dt = time/N
  
  w <- matrix(nrow = N+1, ncol=J)
  w2 <- matrix(nrow = N+1, ncol=J)
  
for (j in 1:J){
    phi <- matrix(rnorm(N, mean = 0, sd = 1),nrow=N, ncol=J)
    phi2 <- matrix(rnorm(N, mean = 0, sd = 1),nrow=N, ncol=J)
    noise <- matrix(rnorm(N, mean = 0, sd = 1),nrow=N, ncol=J)
    for (i in 1:N){
      w[1, j]<- 10
      w2[1,j] <- 100
      dw <- sqrt(dt)%*%phi[i,j]
      dw2 <- rho%*%dw + sqrt(1-rho^2)*sqrt(dt)%*%phi2[i,j]
      w[i+1,j] <- w[i,j] + (alpha*w[i,j]) * dt + sigma*w[i,j] * dw
      w2[i+1, j] <- w2[i,j] + (beta*w2[i,j])*dt+sigma*w2[i,j]^(gamma)*dw2 #+ noise[i,j]
      #w2[i,j]^(gamma)
  }
  print(sprintf("Current number of paths generated: %s, out of %s", j, J))
} 
sim <- cbind(w,w2)
return(sim)  
} 

#' owndmtest
#' 
#' Calculates the Diebold-Mariano test for equal predictive ability.
#' @param e1 vector of losses for benchmark. 
#' @param e2 vetor of losses for alternative forecasting model.
#' @param h integer. Number of forecasting steps. 
#' @return test statisic together with mean and variance for test statistic. 
#' @import xts stats matlib
#' @export 
owndmtest <- function(e1, e2, h=1){
  #d is mean loss differential
  #h is horizon 
  
  time <- length(e1)
  d <- e1 - e2
  
  dMean <- mean(d)
  dVar <- var(d)
  
  DM <- dMean/sqrt( (1/time) * dVar)
  
  temp <- list(DM, dVar, dMean)
  
  return(temp)
}

#' QLIKE
#' 
#' Calculates the bivariate QLIKE loss.
#' @param realized type function. Realized estimator function
#' @param proxy type function. Realized estimator function chosen as proxy.
#' @return loss
#' @import xts stats matlib
#' @export 
QLIKE <- function(realized, proxy){

	#k seems to be the asset dimension in the covariance structure
	#k <- ncol(realized)
  t <- try(inv(realized))
  if("try-error" %in% class(t)){return(0)}
  
	temp <- tr(inv(realized) %*% proxy) - log(Det(inv(realized) %*% proxy)) - 2

	return(temp[1])
  
}