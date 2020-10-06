#Volatility models 


source("Previous functions/projectfunctions.R")
source("functions.R")

library(xts)
library(highfrequency)
library(matlib)
library(MCS)
library(Rsolnp)
library(matrixcalc)
library(MASS)
library(Matrix)


calccov <- readRDS("calculatedcovariances.rds")

mergedfrequencies <- readRDS("mergedfrequencies.rds")


#-----------------------------------The scalar bivariate GARCH model estimation---------------------------

#lT is list of intraday returns as xts object with each list element being each day
#dailyret is daily returns either open-to-close or close-to-close, however, the model
#will target the unconditional realized covariance matrix accordingly, so it's a good idea
#to specify what you're using. Moreover, covariance needs to be an array with 3 dim being days.
#Start with open-to-close. And they should obviouslý have the same amount of days. 
#if covariance is not null, then realized covariance estimations will  be bypassed. 
#Currently you need to source projectfunctions to make it work, however you can set it up with
#highfrequency package. 
#
#
#IF STARTING VALS ARE CLOSE TO OPTIM PARAMETERS THEN THE ALGORITHM TAKES 10SEC! THEREFORE
#FOR A ROLLING FORECAST, SET STARTING VALS AS LAST OPTIM VALS!
#
#
#CODE BELOW IS IF DID NOT PRE-CALCULATE YOUR COVARIANCES. HOWEVER IT MAKES THE OPTIMIZER SLOW.
#IF YOU PRECALCULATE EXPECT 48SEC IF YOU DONT 98SEC. 
#
#if(is.null(covariance)){

#		rcov <- array(0L, dim = c(2, 2, days))

#		for(i in 1:days){

#			rcov[,,i] <- realCov(lT[[i]])
#		}
#	}

#	else{

#		rcov <- covariance 
#	}

BivarGARCHFilter <- function(lT, dailyret, dAlpha, dBeta, covariance){

	days <- length(lT)

	d <-  2 #being bivariate

	samplecov <- cov(dailyret)

	rcov <- array(covariance, dim = c(2,2,days))

	samplercov <- matrix(
		c(mean(rcov[1,1,]),
		  mean(rcov[2,1,]),
		  mean(rcov[2,1,]),
		  mean(rcov[2,2,])),
		  ncol=2, nrow=2)

	mSigma <- array(0L, dim = c(2, 2, days))

	#initializing  with unconditional mean. 

	mSigma[,,1] <- samplecov

	#compute first observation of log-likelihood (we are minizing neg-log-likelihood):
	dLLK <-   log(det(mSigma[,,1])) + (dailyret[1, , drop = F]) %*% solve(mSigma[,,1]) %*% t(dailyret[1, , drop = F])

	astar <- dAlpha * (samplercov * samplecov^(-1))

	for(i in 2:days){

		mSigma[,,i] <- samplecov * (1 - astar - dBeta)  + dBeta * mSigma[,,i-1] + dAlpha * rcov[,,i-1]

		dLLK <- as.numeric(dLLK) + log(det(mSigma[,,i])) + 
		dailyret[i, , drop = F] %*% solve(mSigma[,,i]) %*% t(dailyret[i, , drop = F]) 

		}

	#+ days * d * log(2*pi)
	fulldLLK <- -0.5 * days * d * log(2*pi) - 0.5 * dLLK  

	lOut <- list()

	lOut[["fulldLLK"]] <- fulldLLK
	lOut[["mSigma"]] <- mSigma
	lOut[["cov"]] <- astar

	return(lOut)
	
	
}

leltest <- BivarGARCHFilter(mergedfrequencies[[8]],dailyretotc, 0.33680, 0.66191, calccov[[1]][[8]])

is.null(calccov[[1]][[8]][,,1])


#covariance should be specified in a c() with three arrays.

#if(is.null(covariance)){
#
#		P <- array(0L, dim = c(2, 2, days))
#		N <- array(0L, dim = c(2, 2, days))
#		M <- array(0L, dim = c(2, 2, days))
#
#
#		for(i in 1:days){
#
#			P[,,i] <- realsemicov(lT[[i]], "P")
#			N[,,i] <- realsemicov(lT[[i]], "N")
#			M[,,i] <- realsemicov(lT[[i]], "M")
#		}
#	}
#	else{
#
#		P <- covariance[1]
#		N <- covariance[2]
#		M <- covariance[3] 
#	}


BivarGARCHFilterContAsym <- function(lT, dailyret, dAlphaP, dAlphaN, dAlphaM, dBeta, covariance){

	days <- length(lT)

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
	dLLK <-   log(det(mSigma[,,1])) + (dailyret[1, , drop = F]) %*% solve(mSigma[,,1]) %*% t(dailyret[1, , drop = F])

	astar <- (dAlphaP * sampleP + dAlphaN * sampleN + dAlphaM * sampleM) * samplecov^(-1)

	for(i in 2:days){

		mSigma[,,i] <- samplecov * (1 - astar - dBeta)  + dBeta * mSigma[,,i-1] + 
		dAlphaP * P[,,i-1] + dAlphaN * N[,,i-1] + dAlphaM * M[,,i-1]


		if(!is.positive.definite(mSigma[,,i])){

			mSigma[,,i] <- matrix(nearPD(mSigma[,,i])$mat, ncol=2,nrow=2)

		}

		dLLK <- as.numeric(dLLK) + log(det(mSigma[,,i])) + 
		dailyret[i, , drop = F] %*% solve(mSigma[,,i]) %*% t(dailyret[i, , drop = F]) 

		}

	#+ days * d * log(2*pi)
	fulldLLK <- -0.5 * days * d * log(2*pi) - 0.5 * dLLK  

	lOut <- list()

	lOut[["fulldLLK"]] <- fulldLLK
	lOut[["mSigma"]] <- mSigma
	lOut[["cov"]] <- astar

	return(lOut)
	
	
}

lel <- BivarGARCHFilterContAsym(mergedfrequencies[[8]], dailyretotc, 0.01, 0.005, 0.04, 0.94)
matrix(nearPD(lel$mSigma[,,87])$mat, ncol=2,nrow=2)

lel$mSigma[,,87]


test <- numeric()
for(i in 1:2516){

	test[i] <- !is.positive.definite(lel$mSigma[,,i])

}



#objective function specific for scalar Bivariate GARCH
ObjFBivarGARCH <- function(vPar, lT, dailyret, covariance) {
  
  dAlpha = vPar[1]
  dBeta  = vPar[2]
  dLLK = BivarGARCHFilter(lT, dailyret, dAlpha, dBeta, covariance)$fulldLLK
  
  return(-as.numeric(dLLK))
}

ObjFBivarGARCHContAsym <- function(vPar, lT, dailyret, covariance) {
  
  dAlphaP = vPar[1]
  dAlphaN = vPar[2]
  dAlphaM = vPar[3]
  dBeta  = vPar[4]
  dLLK = BivarGARCHFilterContAsym(lT, dailyret, dAlphaP, dAlphaN, dAlphaM, dBeta, covariance)$fulldLLK
  
  return(-as.numeric(dLLK))
}


#imposing strong stationarity as implied from the paper. This is generalized and can be used for other imposed models.

ineqfun_GARCH_BIVAR <- function(vPar, ...) {
  dAlpha = vPar[1]
  dBeta  = vPar[2]
  
  return(dAlpha + dBeta)
}

ineqfun_GARCH_BIVARContAsym <- function(vPar, ...) {
  dAlpha = vPar[1] + vPar[2] + vPar[3]
  dBeta  = vPar[2]
  
  return(dAlpha + dBeta)
}

#YOU NEED TO CALCULATE ROBUST STANDARD ERRORS. 
EstimateBivarGARCH <- function(lT, dailyret, covariance, ineqfun_GARCH = ineqfun_GARCH_BIVAR, ineqLB = 0.00, ineqUB = 0.9999){
  
  # We set starting value for alpha equal to 0.05, dBeta = 0.94, and chose omega to target
  # the empirical variance by targeting the unconditional variance of the 
  # GARCH model
  
  dAlpha = 0.33680 #
  dBeta  = 0.66190
  
  ## vector of starting parameters
  vPar = c(dAlpha, dBeta)
  
  # have a look at help(solnp)
  ##optimization step
  optimizer = solnp(vPar, fun = ObjFBivarGARCH, lT = lT, dailyret = dailyret, covariance = covariance, 
                    ineqfun = ineqfun_GARCH_BIVAR, #the inequality constraint
                    ineqLB  = ineqLB, ## the inequality lower bound
                    ineqUB = ineqUB, ## the inequality lower bound, i.e. 0.0 < alpha + beta < 0.9999
                    ## lower and upper bounds for all parameters
                    LB = c(0.0001, 0.0001), UB = c(0.999, 0.999)
                    ) 
  
  ## extract estimated parameters
  vPar = optimizer$pars
  
  ## extract the likelihood computed at its maximum
  dLLK = -tail(optimizer$values, 1)

  #gradient tryout,  finite difference:
  #h <- 0.001

  #grad_Lk <- (exp(ObjectiveFunction(c(vPar[1]+h,vPar[2],vPar[3]), vY)) -  
  #exp(ObjectiveFunction(c(vPar[1],vPar[2],vPar[3]), vY)))/h

  
  ## compute filtered volatility
  vSigma2 = BivarGARCHFilter(lT, dailyret, vPar[1], vPar[2], covariance)$mSigma
  
  ## Compute the daily Average BIC
  iT = length(lT)
  BIC = (-2 * dLLK + log(iT) * length(vPar))/iT
  
  ##Standard errors:
  se <- solve(optimizer$hessian)
  se <- matrix(sqrt(diag(se))[-1], ncol=length(vPar), nrow=1)
  colnames(se) <- c("Alpha", "Beta")

  ## return a list with estimated parameters, likelihood value and BIC
  lOut = list()
  lOut[["vPar"]] = vPar
  lOut[["dLLK"]] = dLLK
  lOut[["BIC"]] = BIC
  lOut[["vSigma2"]] = vSigma2
  lOut[["se"]] = se
  lOut[["Hessian"]] = optimizer$hessian
  #lOut[["grad"]] = grad_Lk
  
  return(lOut)
}


#YOU NEED TO CALCULATE ROBUST STANDARD ERRORS. 
EstimateBivarGARCHContAsym <- function(lT, dailyret, covariance, ineqfun_GARCH = ineqfun_GARCH_BIVARContAsym, ineqLB = 0.00, ineqUB = 0.9999){
  
  # We set starting value for alpha equal to 0.05, dBeta = 0.94, and chose omega to target
  # the empirical variance by targeting the unconditional variance of the 
  # GARCH model
  
  dAlphaP = 0.02
  dAlphaN = 0.01
  dAlphaM = 0.01
  dBeta  = 0.94
  
  ## vector of starting parameters
  vPar = c(dAlphaP, dAlphaN, dAlphaM, dBeta)
  
  # have a look at help(solnp)
  ##optimization step
  optimizer = solnp(vPar, fun = ObjFBivarGARCHContAsym, lT = lT, dailyret = dailyret, covariance = covariance, 
                    ineqfun = ineqfun_GARCH_BIVARContAsym, #the inequality constraint
                    ineqLB  = ineqLB, ## the inequality lower bound
                    ineqUB = ineqUB, ## the inequality lower bound, i.e. 0.0 < alpha + beta < 0.9999
                    ## lower and upper bounds for all parameters
                    LB = c(0.0001, 0.0001, -0.999, 0.0001), UB = c(0.999, 0.999, 0.999, 0.999)
                    ) 
  
  ## extract estimated parameters
  vPar = optimizer$pars
  
  ## extract the likelihood computed at its maximum
  dLLK = -tail(optimizer$values, 1)

  #gradient tryout,  finite difference:
  #h <- 0.001

  #grad_Lk <- (exp(ObjectiveFunction(c(vPar[1]+h,vPar[2],vPar[3]), vY)) -  
  #exp(ObjectiveFunction(c(vPar[1],vPar[2],vPar[3]), vY)))/h

  
  ## compute filtered volatility
  vSigma2 = BivarGARCHFilterContAsym(lT, dailyret, vPar[1], vPar[2], vPar[3], vPar[4], covariance)$mSigma
  
  ## Compute the daily Average BIC
  iT = length(lT)
  BIC = (-2 * dLLK + log(iT) * length(vPar))/iT
  
  ##Standard errors:
  se <- solve(optimizer$hessian)
  se <- matrix(sqrt(diag(se))[-1], ncol=length(vPar), nrow=1)
  colnames(se) <- c("AlphaP", "AlphaN", "AlphaM", "Beta")

  ## return a list with estimated parameters, likelihood value and BIC
  lOut = list()
  lOut[["vPar"]] = vPar
  lOut[["dLLK"]] = dLLK
  lOut[["BIC"]] = BIC
  lOut[["vSigma2"]] = vSigma2
  lOut[["se"]] = se
  lOut[["Hessian"]] = optimizer$hessian
  #lOut[["grad"]] = grad_Lk
  
  return(lOut)
}





library(numDeriv, MASS)
library(MASS)

tt <- jacobian(ObjFBivarGARCH, x = lel2$vPar, lT = mergedfrequencies[[8]], dailyret = dailyretotc)

t(tt) %*% ginv(lel2$Hessian)[2:3, 2:3] %*% (tt)



dataTLT <- readRDS("dataTLT.rds")

getDates <- unlist(lapply(dataTLT, function(x) as.character(index(x[1]))))

for(i in 1:length(getDates)){
	getDates[i] <- strsplit(getDates, " ")[[i]][1]
}






dailyretotc <- xts(t(sapply(mergedfrequencies[[10]], function(x) cbind(x[,1], x[,2]))), order.by = as.Date(getDates))

colnames(dailyretotc) <- c("TLT", "SPY")




library(tictoc)

tic()
lel2 <- EstimateBivarGARCH(mergedfrequencies[[8]], dailyretotc, calccov[[1]][[8]])
toc()


calccov[[1]][[8]][1,1,]


#now you have to precalculate the realized semi-covariances

P <- array(0L, c(2,2,2516))
N <- array(0L, c(2,2,2516))
M <- array(0L, c(2,2,2516))

for(i in 1:2516){

	P[,,i] <- realsemicov(mergedfrequencies[[9]][[i]], "P")
	N[,,i] <- realsemicov(mergedfrequencies[[9]][[i]], "N")
	M[,,i] <- realsemicov(mergedfrequencies[[9]][[i]], "M")

}


tic()
lel3 <- EstimateBivarGARCHContAsym(mergedfrequencies[[9]], dailyretotc, list(P, N, M))
toc()

lel2$BIC





test <- matrix(c(mean(calccov[[1]][[7]][1,1,]), mean(calccov[[1]][[7]][2,1,]), 
	mean(calccov[[1]][[7]][2,1,]), mean(calccov[[1]][[7]][2,2,])), ncol=2,nrow=2)

determinant(test, F)$modulus[1]
det(test)

test[1, , drop = F]
test

(mergedfrequencies[[7]][[1]][1, , drop = F]) %*% solve(test) %*% t(mergedfrequencies[[7]][[1]][1, , drop = F])


test2 <- array(rnorm(100), c(2,2,25))

test2[1, , , drop = F]