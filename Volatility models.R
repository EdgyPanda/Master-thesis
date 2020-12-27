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



acf(mergedfrequencies[[7]][[1]][,1], type = c("covariance"))$acf



dataTLT <- readRDS("dataTLT.rds")
dataSPY <- readRDS("dataSPY.rds")

getDates <- unlist(lapply(dataTLT, function(x) as.character(index(x[1]))))

for(i in 1:length(getDates)){
	getDates[i] <- strsplit(getDates, " ")[[i]][1]
}


dailyretotc <- xts(t(sapply(mergedfrequencies[[10]], function(x) cbind(x[,1], x[,2]))), order.by = as.Date(getDates))


acf(dailyretotc[,1]^2)

colnames(dailyretotc) <- c("TLT", "SPY")


P <- array(0L, dim =c(2,2,2516))
N <- array(0L, dim =c(2,2,2516))
M <- array(0L, dim =c(2,2,2516))

for(i in 1:2516){

	P[,,i] <- realsemicov(mergedfrequencies[[7]][[i]], "P")
	N[,,i] <- realsemicov(mergedfrequencies[[7]][[i]], "N")
	M[,,i] <- realsemicov(mergedfrequencies[[7]][[i]], "M")


}


#---------------------------------------------------------------------------------------

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
#
#
#
#
#NOTE TO SELF! SOLVERS CONVERGE BETTER WHEN LOG-RETURNS ARE IN PERCENTAGE!


################################################################################################################
#
#
#
#									Realized Bivarite Scalar GARCH (rBG)
#
#
#
################################################################################################################

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

	#omega under cov targeting
	omega <-   samplecov - dAlpha * samplercov %*% id - dBeta * samplecov #  samplecov^(-1)

	mSigma[,,1] <- samplecov 

	#compute first observation of log-likelihood (we are minizing neg-log-likelihood):
	dLLK <-  d* log(2*pi) + log(det(mSigma[,,1])) + (dailyret[1, , drop = F]) %*% solve(mSigma[,,1]) %*% t(dailyret[1, , drop = F])


	dLLKs <- numeric()
	dLLKs[1] <-   dLLK
	for(i in 2:days){

		mSigma[,,i] <- omega + dBeta * mSigma[,,i-1] + dAlpha * rcov[,,i-1]

		#neglog collection for score calculation
		dLLKs[i] <-   d* log(2*pi) +   (log(det(mSigma[,,i])) +  
		dailyret[i, , drop = F] %*% solve(mSigma[,,i]) %*% t(dailyret[i, , drop = F]))
		}

	fulldLLK <-  - 0.5 * sum(dLLKs) #loglikelihood

	lOut <- list()

	lOut[["fulldLLK"]] <- fulldLLK
	lOut[["mSigma"]] <- mSigma
	lOut[["cov"]] <- omega
	lOut[["dLLKs"]] <- dLLKs
	lOut[["sampleH"]] <- samplecov

	return(lOut)
	
	
}


tt <- BivarGARCHFilter(dailyretotc * 100, c(0.12,0.85), calccov[[1]][[7]] * 10000)

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

EstimateBivarGARCH <- function(dailyret, covariance, bootstrap = FALSE, vPar=NULL, ineqfun_GARCH = ineqfun_GARCH_BIVAR, ineqLB = 0.00, ineqUB = 0.9999){
  
  # We set starting value for alpha equal to 0.05, dBeta = 0.94, and chose omega to target
  # the empirical variance by targeting the unconditional variance of the 
  # GARCH model
  
  dAlpha = 0.39873  
  dBeta  = 0.57351
  
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
                    control = list(tol = 1e-7, outer.iter = 800, inner.iter = 1200, delta = 1e-7)
                    ) 
  
  ## extract estimated parameters
  vPar = optimizer$pars
  
  ## extract the likelihood computed at its maximum
  dLLK = -tail(optimizer$values, 1)

  if(!bootstrap){
  #CALCULATING ROBUST STD. ERRORS WHEN NO COVARIANCE TARGETING:

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

  ## compute filtered volatility
  Filter <- BivarGARCHFilter(dailyret, c(vPar[1], vPar[2]), covariance)

  vSigma2 <- Filter$mSigma
  
  sampleH <- Filter$sampleH

  ## Compute the daily Average BIC
  iT = length(dailyret[,1])
  BIC = (-2 * dLLK + log(iT) * length(vPar))/iT
  
  ##Standard errors:
  se <- solve(optimizer$hessian)
  se <- matrix(sqrt(diag(se))[-1], ncol=length(vPar), nrow=1)/iT
  colnames(se) <- c("Alpha", "Beta")

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

#it does not matter whether you use percentage returns with realcov calculated on percentage returns or 
#standard returns. but if you want the loglike you should not use percentage returns. 

tester <- EstimateBivarGARCH(dailyretotc * 100, calccov[[1]][[10]]*10000)

qlikes_rBG <- numeric()

for(i in 2:2516){
  qlikes_rBG[i] <-  QLIKE(tester$vSigma2[,,i-1],calccov[[1]][[7]][,,i]*10000, 2)
}

mean(qlikes_rBG, na.rm = T)

estacrossfreq_rBG <- list()

#should be done under reasonable sampling frequencies. 
for(i in 6:10){

  estacrossfreq_rBG[[i]] <- EstimateBivarGARCH(dailyretotc * 100, calccov[[1]][[i]] * 10000)

}

params <- matrix(unlist(lapply(estacrossfreq_rBG, function(x) x$vPar)), ncol = 2, byrow = T)

colMeans(params)
quantile(params[,1], 0.1)
quantile(params[,1], 0.9)
quantile(params[,2], 0.1)
quantile(params[,2], 0.9)

avgloglike <- mean(unlist(sapply(estacrossfreq_rBG, function(x) x$dLLK)))

avgavgbic <- mean(unlist(sapply(estacrossfreq_rBG, function(x) x$BIC)))

#calc QLIKES: 

sigmas_rBG <- lapply(estacrossfreq_rBG, function(x) x$vSigma2)

Qlikes <- matrix(0L, ncol = 10, nrow = 2516)

for(i in 6:10){
  for(j in 2:2516){
    Qlikes[j, i] <- QLIKE(sigmas_rBG[[i]][,,j-1],calccov[[1]][[7]][,,j]*10000, 2)
  }
}
mean(colMeans(Qlikes, na.rm = T)[6:10])

library(ggplot2)

ggplot() + geom_line(aes(as.Date(getDates), sigmas_rBG[[10]][2,2,])) + geom_line(aes(as.Date(getDates), calccov[[1]][[10]][2,2,]*10000))




#------------------------------------------bootstrapping----------------------------------------------
#
#
# bootstrapping is only done on a 5 minute realized covariance structure to preserve time. 

library(boot)

library(np)
#uses the row to construct the block length. Therefore you should transpose it. 
ll <- tsboot(dailyretotc[,1], mean ,R = 1000, l = 10, sim = "fixed", endcorr = TRUE)
indexation <- boot.array(ll)

#rearranging calculated covariances: correct. 

#calccov_rBG <- list(calccov[[1]], calccov[[4]], calccov[[5]], calccov[[6]], calccov[[7]], calccov[[8]]) 


estimates_5min_rBG <- matrix(0L, nrow = 1000, ncol = 2)

BS_dailyretotc <- matrix(dailyretotc, nrow=2516, ncol = 2)

BS_dailyretotc_l <- list()

for(i in 1:1000){

BS_dailyretotc_l[[i]] <- BS_dailyretotc[indexation[i, ], ]

}

#ONLY BOOTSTRAPPING FOR 5 MIN REALIZED COVARIANCE. 

for(j in 1:1000){
						#5min_covariances
	BS_calccov_rBG <- calccov_rBG[[1]][[7]][,,indexation[j, ]]
	BS_dailyretotc <- BS_dailyretotc_l[[j]]

	estimates_5min_rBG[j, ] <- EstimateBivarGARCH(BS_dailyretotc, BS_calccov_rBG, nobootstrap = F)
	print(sprintf("bootstrap %s",j))
}

#saveRDS(estimates_5min_rBG, "bootstrapestimates_5min_rBG.rds")

bootstraps_stimates_5min_rBG <- readRDS("bootstrapestimates_5min_rBG.rds") 

bootstrap_stderror_5min_rBG <- apply(bootstraps_stimates_5min_rBG, MARGIN = c(2), FUN = function(x) sd(x))



#------------------------------------------------------------------------------------------------------



#testing via simulation

leltest <- BivarGARCHFilter(mergedfrequencies[[8]],dailyretotc, c(0.2, 0.8), calccov[[1]][[8]]) 

leltest$mSigma[,,1]

sum(leltest$dLLKs)

set.seed(1234)
retsim <- matrix(0L, nrow=2516, ncol=2)
for(i in 1:2516){
	
	retsim[i,] <-  mvrnorm(1, c(0,0), leltest$mSigma[,,i])

}

#it is okay, but could be better. It very much depend on the choice of randomness. However, it converges. 
#it will converge better if you had more simulated data. 
#How to do it: construct 2 dim return data 100.000 points eg, from some covariance matrix.
#estimate 30 min covariance from 100.000/13 =7592 days. Split dataset and use endpoints as daily data. 
#now put into filter algoritm. Construct returns following mvrnorm. See estimation method. 
leltest2 <- EstimateBivarGARCH(mergedfrequencies[[8]],retsim, calccov[[1]][[8]])



################################################################################################################







#NOW IT NEEDS SEMI-COVARIANCES SPECIFIED AS ARRAY WITH 3RD DIM BEING DAYS. THEN PUT THE 3 SEMI-COVARIANCES
#INTO A LIST WITH REPSECTIVE ORDER P,N,M. 
BivarGARCHFilterContAsym <- function(dailyret, params, covariance){

	dAlphaP <- params[1]
	dAlphaN <- params[2]
	dAlphaM <- params[3]
	dBeta  <- params[4]

	days <- 2516

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
	dLLK <-   log(det(mSigma[,,1])) + (dailyret[1, , drop = F]) %*% solve(mSigma[,,1]) %*% t(dailyret[1, , drop = F])

	id <- matrix(c(1,0,0,1), ncol = 2, nrow = 2)

	#astar <- (dAlphaP * sampleP + dAlphaN * sampleN + dAlphaM * sampleM) * samplecov^(-1)
	#this is wrong see functions.R
	omega <- (dAlphaP * sampleP) %*% id +  (dAlphaN * sampleN) %*% id + (dAlphaM * sampleM) %*% id

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
		
		dLLKs[i] <- 0.5 * d * log(2*pi) + 0.5 * as.numeric(log(det(mSigma[,,i])) + 
		dailyret[i, , drop = F] %*% solve(mSigma[,,i]) %*% t(dailyret[i, , drop = F]))

		}

	#+ days * d * log(2*pi)
	fulldLLK <- - sum(dLLKs, na.rm=T)  

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




#YOU NEED TO CALCULATE ROBUST STANDARD ERRORS. 
EstimateBivarGARCHContAsym <- function(dailyret, covariance, ineqfun_GARCH = ineqfun_GARCH_BIVARContAsym, ineqLB = 0.00, ineqUB = 0.9999){
  
  # We set starting value for alpha equal to 0.05, dBeta = 0.94, and chose omega to target
  # the empirical variance by targeting the unconditional variance of the 
  # GARCH model
  
  dAlphaP = 0.09324267
  dAlphaN = 0.34553257
  dAlphaM = 0.04110211
  dBeta  = 0.52002253
  
  ## vector of starting parameters
  vPar = c(dAlphaP, dAlphaN, dAlphaM, dBeta)
  
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

  scores <- matrix(0L, nrow=2516, ncol = length(vPar))

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

  
  ## compute filtered volatility
  vSigma2 = BivarGARCHFilterContAsym(dailyret, c(vPar[1], vPar[2], vPar[3], vPar[4]), covariance)$mSigma
  
  ## Compute the daily Average BIC
  iT = 2516
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
  lOut[["rse"]] = rse
  #lOut[["grad"]] = grad_Lk
  
  return(lOut)
}


#5min frequency stemming from P, N, M.  I havent scaled in the estimator

leltest4 <- EstimateBivarGARCHContAsym(dailyretotc*100, list(P*10000,N*10000,M*10000))

qlikes_crBG <- numeric()

for(i in 2:2516){
  qlikes_crBG[i] <-  QLIKE(leltest4$vSigma2[,,i-1],calccov[[1]][[7]][,,i]*10000, 2)
}

mean(qlikes_crBG, na.rm = T)

#constructing semicovs from entire sampling frequency:
if(FALSE){
ptemp <- array(0L, dim =c(2,2,2516)) 
ntemp <- array(0L, dim =c(2,2,2516)) 
mtemp <- array(0L, dim =c(2,2,2516)) 

Pfreq <- list()
Nfreq <- list()
Mfreq <- list()

for(j in 1:10){
  for(i in 1:2516){

    ptemp[,,i] <- realsemicov(mergedfrequencies[[j]][[i]], "P")
    ntemp[,,i] <- realsemicov(mergedfrequencies[[j]][[i]], "N")
    mtemp[,,i] <- realsemicov(mergedfrequencies[[j]][[i]], "M")
    print(sprintf("%s", i))
  }

  Pfreq[[j]] <- ptemp 
  Nfreq[[j]] <- ntemp 
  Mfreq[[j]] <- mtemp 
  print(sprintf("%s", j))
}
}

#--------------------

semicov_acrossfreq <- list(Pfreq, Nfreq, Mfreq)

#saveRDS(semicov_acrossfreq, "semicov_acrossfreq.rds")

semicovs <- readRDS("semicov_acrossfreq.rds")

Pfreq <- semicovs[[1]]
Nfreq <- semicovs[[2]]
Mfreq <- semicovs[[3]]

estacrossfreq_crBG <- list()

for(i in 10:10){

  estacrossfreq_crBG[[i]] <- EstimateBivarGARCHContAsym(dailyretotc * 100, 
    list(Pfreq[[i]]*10000,Nfreq[[i]]*10000,Mfreq[[i]]*10000))

}

params <- matrix(unlist(lapply(estacrossfreq_crBG, function(x) x$vPar)), ncol = 4, byrow = T)

colMeans(params)
apply(params, MARGIN = c(2), FUN =function(x) quantile(x, 0.1))
apply(params, MARGIN = c(2), FUN =function(x) quantile(x, 0.9))


avgloglike <- mean(unlist(sapply(estacrossfreq_crBG, function(x) x$dLLK)))

avgavgbic <- mean(unlist(sapply(estacrossfreq_crBG, function(x) x$BIC)))

#calc QLIKES: 

sigmas_crBG <- lapply(estacrossfreq_crBG, function(x) x$vSigma2)

Qlikes <- matrix(0L, ncol = 10, nrow = 2516)

for(i in 6:10){
  for(j in 2:2516){
    Qlikes[j, i] <- QLIKE(sigmas_crBG[[i]][,,j-1],calccov[[1]][[7]][,,j]*10000, 2)
  }
}
mean(colMeans(Qlikes, na.rm = T)[6:10])




#------------------------------------------bootstrapping----------------------------------------------

library(boot)

library(np)
#uses the row to construct the block length. Therefore you should transpose it. 
ll <- tsboot(dailyretotc[,1], mean ,R = 1000, l = 10, sim = "fixed", endcorr = TRUE)
indexation <- boot.array(ll)

#rearranging calculated covariances: correct. 

#Due to the assymetric continuous garch model we omit the realized semicovariances from the estimation
#procedure of crBG:



estimates_5min_crBG <- matrix(0L, nrow = 1000, ncol = 4)


BS_dailyretotc <- matrix(dailyretotc, nrow=2516, ncol = 2)

BS_dailyretotc_l <- list()

for(i in 1:1000){

BS_dailyretotc_l[[i]] <- BS_dailyretotc[indexation[i, ], ]

}


for(j in 1:1000){
						#5min_covariances
	BS_P_crBG <- P[,,indexation[j, ]]
	BS_N_crBG <- N[,,indexation[j, ]]
	BS_M_crBG <- M[,,indexation[j, ]]

	BS_dailyretotc <- BS_dailyretotc_l[[j]]

	estimates_5min_crBG[j, ] <- EstimateBivarGARCHContAsym(BS_dailyretotc * 100, list(BS_P_crBG*10000, BS_N_crBG*10000, BS_M_crBG*10000))$vPar
	print(sprintf("Bootstrap %s", j))
	

}

#saveRDS(estimates_5min_crBG, "bootstapestimates_5min_crBG.rds")



estimates_5min_crBG <- readRDS("bootstapestimates_5min_crBG.rds")

#removing spurious results. 
rem1 <- which(estimates_5min_crBG[,2]<0.2)
rem2 <- which(estimates_5min_crBG[,4]<0.2)
rem3 <- which(estimates_5min_crBG[,3]>0.15)
rem4 <- which(estimates_5min_crBG[,1]>0.3)

removed <- unique(c(rem1, rem2, rem3, rem4))

estimates_5min_crBG <- estimates_5min_crBG[-c(removed), ]

bootstrap_stderror_5min_crBG <- apply(estimates_5min_crBG, MARGIN = c(2), FUN = function(x) sd(x))





#------------------------------------realized DCC model (rDCC) --------------------------------------

#mEta is the residuals from the univariate models. 

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
Estimate_rDCC <- function(mY, covariance, getDates, bootstrap = F, residuals = NULL) {
  
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
  vPar = c(0.1974634, 0.8015366)
  
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
  mEta = mEta, mQ = mQ, covariance = covariance)
  
  #Extract the estimated parameters
  vPar = optimizer$pars
  
  #Extract the likelihood of the correlation part
  dLLK_C = -tail(optimizer$values, 1)
  
  #Filter the dynamic correlation using the estimated parameters
  Filter = rDCCFilter(mEta, vPar[1], vPar[2], mQ, covariance)

  if(bootstrap){

    return(vPar)
  
  }
  #standard errors 
  se <- solve(optimizer$hessian)
  se <- matrix(sqrt(diag(se))[-1], ncol=length(vPar), nrow=1)

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
  
  mZ = matrix(0, iT, ncol(mY))
  
  for (t in 1:iT) {
    mZ[t, ] = solve(chol(aCor[,,t])) %*% as.numeric(mEta[t, ])
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


#solnp can (rarely) give an error, try and run all functions one more time and retry
rDCCests <- Estimate_rDCC(dailyretotc*100, calccov[[1]][[7]]*10000, getDates)


qlikes_rDCC <- numeric()

for(i in 2:2516){
  qlikes_rDCC[i] <-  QLIKE(rDCCests$vSigma2[,,i-1],calccov[[1]][[7]][,,i]*10000, 2)
}

mean(qlikes_rDCC, na.rm = T)


estacrossfreq_rDCC <- list()

for(i in 6:10){

  estacrossfreq_rDCC[[i]] <- Estimate_rDCC(dailyretotc * 100, calccov[[1]][[i]]*10000, getDates)

}


params <- matrix(unlist(lapply(estacrossfreq_rDCC, function(x) x$vPar)), ncol = 2, byrow = T)

colMeans(params)
apply(params, MARGIN = c(2), FUN =function(x) quantile(x, 0.1))
apply(params, MARGIN = c(2), FUN =function(x) quantile(x, 0.9))


avgloglike <- mean(unlist(sapply(estacrossfreq_rDCC, function(x) x$dLLK)))

avgavgbic <- mean(unlist(sapply(estacrossfreq_rDCC, function(x) x$BIC)))

#calc QLIKES: 

sigmas_rDCC <- lapply(estacrossfreq_rDCC, function(x) x$vSigma2)

Qlikes <- matrix(0L, ncol = 10, nrow = 2516)

for(i in 6:9){
  for(j in 2:2516){
    Qlikes[j, i] <- QLIKE(sigmas_rDCC[[i]][,,j-1],calccov[[1]][[7]][,,j]*10000, 2)
  }
}
mean(colMeans(Qlikes, na.rm = T)[6:9])









#------------------------testing if dcc likelihood is positive for my data following proper coded dcc model -----

# it is. 



library(rmgarch)

xspec <- ugarchspec()
xspec <- ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), variance.model = 
  list(model = 'realGARCH', garchOrder = c(2,1))) #he changed the order in the code see starch. 

setbounds(xspec)<-list(alpha2=c(-1,1))

uspec <- multispec(replicate(2, xspec))

spec1 = dccspec(uspec = uspec, dccOrder = c(1, 1), distribution = 'mvnorm')

fit <- dccfit(spec1, data = dailyretotc*100) #realizedVol = tt143)
cc
tt143 <- xts(cbind(sqrt(calccov[[1]][[7]][1,1,]*10000), sqrt(calccov[[1]][[7]][2,2,]*10000)), 
order.by = as.Date(1:2516))

tt12 <- ugarchfit(xspec, data = dailyretotc[,1] * 100)
tt13 <-  ugarchfit(xspec, data = dailyretotc[,2] * 100)

likelihood(tt12) + likelihood(tt13)

#implies that the optimizer values are negative in your case which is correct!
dcc.likelihood <-  likelihood(fit) - (likelihood(tt12) + likelihood(tt13))


#------------------------------------------bootstrapping----------------------------------------------

library(boot)

library(np)

#see patton "does anything beat rv5 min" says that the block length is driven by the persistence in the
# variable we are interested in testing. 


spec <- ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), variance.model = 
  list(model = 'realGARCH', garchOrder = c(2, 1)))

setbounds(spec)<-list(alpha2=c(-1,1))


spec1 <- ugarchfit(spec, dailyretotc[,1] * 100, solver = 'hybrid', realizedVol = 
xts(sqrt(calccov[[1]][[7]][1,1,] * 10000), order.by = as.Date(getDates)), 
fit.control = list(stationarity = 1, fixed.se = 1))

spec2 <- ugarchfit(spec, dailyretotc[,2] * 100, solver = 'hybrid', realizedVol = 
xts(sqrt(calccov[[1]][[7]][2,2,] * 10000), order.by = as.Date(getDates)), 
fit.control = list(stationarity = 1, fixed.se = 1))


mEta <- matrix(c(residuals(spec1, standardize = T), residuals(spec2, standardize = T)), ncol=2, byrow = T)

ts.plot(mEta[,1])

hist(mEta[,2], breaks = 60)

head(mEta)
#optimal blockk length:
b.star(calccov[[1]][[7]][2,1, ]/sqrt(calccov[[1]][[7]][1,1,  ] * calccov[[1]][[7]][2,2, ]))



		#correlation models are very persistent and therefore needs bigger block length. 
		#if a too small block length is chosen it can destroy the dependency in the data, also leading
		#to non-invertible hessians in the rugarch package. 
ll <- tsboot(dailyretotc[,1], mean ,R = 1000, l = 97, sim = "fixed", endcorr = TRUE)
indexation <- boot.array(ll)



estimates_5min_rDCC <- matrix(0L, nrow = 1000, ncol = 2)


BS_dailyretotc <- matrix(dailyretotc, nrow=2516, ncol = 2)

BS_dailyretotc_l <- list()

mEta_l <- list()

for(i in 1:1000){

BS_dailyretotc_l[[i]] <- xts(BS_dailyretotc[indexation[i, ], ], order.by = as.Date(getDates))
mEta_l[[i]] <- mEta[indexation[i, ], ]
}

head(mEta_l[[1]])
head(indexation[1, ])
head(mEta[1852: 1857, ])

#bootstrapping wrt residuals does not work at all. 

for(j in 1:1000){
	
	BS_calccov_rDCC <- calccov[[1]][[7]][,,indexation[j, ]]

	BS_dailyretotc <- BS_dailyretotc_l[[j]]

  mEta <- mEta_l[[j]]

	estimates_5min_rDCC[j, ] <- Estimate_rDCC(NULL, BS_calccov_rDCC*10000, getDates, TRUE, residuals = mEta)
	print(sprintf("Bootstrap %s", j))
	

}

#check
ts.plot(estimates_5min_rDCC[1:421,2])


ts.plot(estimates_5min_rDCC[-c(which(estimates_5min_rDCC[1:421,1]<0.01)),1])
estimates

#saveRDS(estimates_5min_rDCC, "bootstrapestimates_5min_rDCC_residuals.rds")

#saveRDS(bootstrap5min_est_complete, "bootstapestimates_5min_rDCC.rds")

bootstapestimates_5min_rDCC <- readRDS("bootstrapestimates_5min_rDCC_residuals.rds")


rem1 <- which(bootstapestimates_5min_rDCC[,1]<0.1 | bootstapestimates_5min_rDCC[,1]>0.4)
rem2 <- which(bootstapestimates_5min_rDCC[,2]<0.55)

rem <- unique(c(rem1,rem2))

bootstapestimates_5min_rDCC <- bootstapestimates_5min_rDCC[-c(rem), ]


bootstrap_stderror_5min_rDCC <- apply(bootstapestimates_5min_rDCC, MARGIN = c(2), FUN = function(x) sd(x))



#------------------------------------------------realGARCH 5min estimates and standard errors--------------------

library(rugarch)



spec <- ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), variance.model = 
	list(model = 'realGARCH', garchOrder = c(2,1))) #he changed the order in the code see starch. 

setbounds(spec)<-list(alpha2=c(-1,1)) #you need to reset the bounds for alpha2 otherwise it will be zero since it 
#originally is set between (0,1). 


#TLT
lfit2 = ugarchfit(spec, dailyretotc[,1]*100, solver = 'hybrid', realizedVol = 
	xts((calccov[[1]][[7]][1,1,]*10000), order.by = as.Date(getDates)), 
	fit.control = list(stationarity = 1, fixed.se = 0))


test2 <- sigma(lfit2)

rugarch.LL = c('logL' = sum(-lfit2@fit$log.likelihoods), 'pLogL' = sum(-lfit2@fit$partial.log.likelihoods))
sum(rugarch.LL)


#estimating across frequencies:
tt14 <- coef(lfit2)

estimatesacross_TLT <- matrix(0L, ncol = 9, nrow = 10)
estimatesacross_SPY <- matrix(0L, ncol = 9, nrow = 10)
tlt.loglike <- matrix(0L, ncol = 2, nrow = 10)
spy.loglike <- matrix(0L, ncol = 2, nrow = 10)



for(i in 1:10){
	#TLT

	tlt <-  ugarchfit(spec, dailyretotc[,1] * 100, solver = 'hybrid', realizedVol = 
	xts((calccov[[1]][[i]][1,1,]*10000), order.by = as.Date(getDates)), 
	fit.control = list(stationarity = 1, fixed.se = 1))

	spy <- ugarchfit(spec, dailyretotc[,2] * 100, solver = 'hybrid', realizedVol = 
	xts((calccov[[1]][[i]][2,2,]*10000), order.by = as.Date(getDates)), 
	fit.control = list(stationarity = 1, fixed.se = 1))

	estimatesacross_TLT[i, ] <- coef(tlt)

	tlt.loglike[i, ] <- c('logL' = sum(-(tlt@fit$log.likelihoods)), 'pLogL' = sum(-tlt@fit$partial.log.likelihoods))
	#SPY
	spy.loglike[i, ] <- c(sum(-(spy@fit$log.likelihoods)), sum(-spy@fit$partial.log.likelihoods))
	estimatesacross_SPY[i, ] <- coef(spy)

}


colnames(estimatesacross_TLT) <- names(coef(lfit2))
rownames(estimatesacross_TLT) <- c("1sec", "5sec", "15sec", "20sec", "30sec", "1min", "5min", "15min", "30min", "daily")
colnames(estimatesacross_SPY) <- colnames(estimatesacross_TLT)
rownames(estimatesacross_SPY) <- rownames(estimatesacross_TLT)
colnames(tlt.loglike) <- c("logL", "pLogL")
colnames(spy.loglike) <- colnames(tlt.loglike)

#needed to smooth out the daily scheme to ensure consistency
daily_SPY <- coef(ugarchfit(spec, dailyretotc[,2] * 100, solver = 'hybrid', realizedVol = 
	xts(sqrt(smooth(calccov[[1]][[10]][2,2,])) * 100, order.by = as.Date(getDates)), 
	fit.control = list(stationarity = 1, fixed.se = 1)))
daily_TLT <- coef(ugarchfit(spec, dailyretotc[,1] * 100, solver = 'hybrid', realizedVol = 
	xts(sqrt(smooth(calccov[[1]][[10]][1,1,])) * 100, order.by = as.Date(getDates)), 
	fit.control = list(stationarity = 1, fixed.se = 1)))

estimatesacross_TLT[10, ] <- daily_TLT
estimatesacross_SPY[10, ] <- daily_SPY

estimatesacross_TLT_mean <- colMeans(estimatesacross_TLT[-10, ])
estimatesacross_SPY_mean <- colMeans(estimatesacross_SPY[-10, ])

TLT_90percentquantile <- apply(estimatesacross_TLT, MARGIN = c(2), FUN = function(x) quantile(x, 0.9, na.rm=TRUE))
TLT_10percentquantile <- apply(estimatesacross_TLT, MARGIN = c(2), FUN = function(x) quantile(x, 0.1, na.rm=TRUE))

SPY_90percentquantile <- apply(estimatesacross_SPY[-10, ], MARGIN = c(2), FUN = function(x) quantile(x, 0.9, na.rm=TRUE))
SPY_10percentquantile <- apply(estimatesacross_SPY[-10, ], MARGIN = c(2), FUN = function(x) quantile(x, 0.1, na.rm=TRUE))



TLT_stats <- rbind(estimatesacross_TLT_mean, TLT_10percentquantile, TLT_90percentquantile)
SPY_stats <- rbind(estimatesacross_SPY_mean, SPY_10percentquantile, SPY_90percentquantile)




colMeans(tlt.loglike, na.rm = T)
quantile(tlt.loglike[,2], 0.9, na.rm = T)
quantile(tlt.loglike[,2], 0.1, na.rm = T)

mean(tlt.loglike[,1], na.rm = T)
quantile(tlt.loglike[,1], 0.9, na.rm = T)
quantile(tlt.loglike[,1], 0.1, na.rm = T)


colMeans(spy.loglike, na.rm = T)
quantile(spy.loglike[,2], 0.9, na.rm = T)
quantile(spy.loglike[,2], 0.1, na.rm = T)

mean(spy.loglike, na.rm = T)
quantile(spy.loglike[,1], 0.9, na.rm = T)
quantile(spy.loglike[,1], 0.1, na.rm = T)




#SPY
lfit = ugarchfit(spec, dailyretotc[,2] * 100, solver = 'hybrid', realizedVol = 
	(xts((calccov[[1]][[6]][2,2,]* 10000), order.by = as.Date(getDates))), 
	fit.control = list(stationarity = 0, fixed.se = 0))


rugarch.LL = c('logL' = sum(-lfit@fit$log.likelihoods), 'pLogL' = sum(-lfit@fit$partial.log.likelihoods))



lfit <- ugarchfit(spec, spyreal[,1] , solver = 'hybrid', realizedVol = 
	spyreal[,2], fit.control = list(stationarity = 1, fixed.se = 1))


spec = ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), variance.model = list(model = 'realGARCH', garchOrder = c(2, 1)))


fit = ugarchfit(spec, spyreal[, 1] * 100, solver = 'hybrid', realizedVol = spyreal[,2] * 100)

rugarch.LL = c('logL' = sum(-fit@fit$log.likelihoods[1:1492]), 'pLogL' = sum(-fit@fit$partial.log.likelihoods[1:1492]))

sum(rugarch.LL)


data(spyreal)
head(spyreal)

ts.plot(spyreal[,2])
ts.plot((calccov[[8]][[7]][2,2,]))
sqrt(head(calccov[[1]][[7]][2,2,]))

#---------------------------------realized continuous asymmetric DCC model (rcDCC)---------------------------

rcDCCFilter <- function(mEta, dAP, dAN, dAM , dB, mQ, covariance) {
  
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
  

  astar <- (dAP * sampleP + dAN * sampleN  + dAM * sampleM) * mQ^(-1)
  #main loop
  for (t in 2:iT) {
    #update the Q matrix
    aQ[,, t] <- mQ * (1 - astar - dB) + dAP * Pcor[,,t - 1] + dAN * Ncor[,,t - 1]  + dAM * Mcor[,,t - 1]  + dB * aQ[,,t - 1]
    
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
  
  return(lOut)
}



Estimate_rcDCC <- function(mY, covariance, cov2 = NULL, getDates, bootstrap = F, residuals = NULL){
  
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
  vPar = c(0.00000001384, 0.00000013727, 0.08197776333, 0.91702207742)
  
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
  mEta = mEta, mQ = mQ, covariance = covariance)
  
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
  
  mZ = matrix(0, iT, ncol(mY))
  
  for (t in 1:iT) {
    mZ[t, ] = solve(chol(aCor[,,t])) %*% as.numeric(mEta[t, ])
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
  return(lOut)
  
}

#I have not scaled inside the estimator. 
crDCCests <- Estimate_rcDCC(dailyretotc * 100, list(P * 10000, N * 10000, 
	M * 10000), NULL , getDates, F)


tt <- apply(crDCCests$vSigma, MARGIN = c(3), FUN = function(x) diag(sqrt(1/diag(x))))

tt[1,] <-  dailyretotc[,1]/tt[1, ] 
tt[4,] <-  dailyretotc[,2]/tt[4, ] 


ts.plot(tt[4,])
acf(tt[1,])

estacrossfreq_crDCC <- Estimate_rcDCC(dailyretotc[1:1000, ] * 100, 
    list(Pfreq[[1]][,,1:1000]*10000,Nfreq[[1]][,,1:1000]*10000, Mfreq[[1]][,,1:1000]*10000), NULL, getDates[1:1000])

qlikes_crDCC <- numeric()

for(i in 2:2516){
  qlikes_crDCC[i] <-  QLIKE(crDCCests$vSigma2[,,i-1],calccov[[1]][[7]][,,i]*10000, 2)
}

mean(qlikes_crDCC, na.rm = T)

semicovs <- readRDS("semicov_acrossfreq.rds")

Pfreq <- semicovs[[1]]
Nfreq <- semicovs[[2]]
Mfreq <- semicovs[[3]]

estacrossfreq_crDCC <- list()

#cannot do daily returns. 
for(i in 2:9){

  estacrossfreq_crDCC[[i]] <- Estimate_rcDCC(dailyretotc * 100, 
    list(Pfreq[[i]]*10000,Nfreq[[i]]*10000,Mfreq[[i]]*10000), NULL, getDates)

}

params <- matrix(unlist(lapply(estacrossfreq_crDCC, function(x) x$vPar)), ncol = 4, byrow = T)

options(digits = 6)
colMeans(params)
apply(params[-2, ], MARGIN = c(2), FUN =function(x) quantile(x, 0.1))
apply(params[-2, ], MARGIN = c(2), FUN =function(x) quantile(x, 0.9))


avgloglike <- mean(unlist(sapply(estacrossfreq_crDCC, function(x) x$dLLK))[-2])

avgavgbic <- mean(unlist(sapply(estacrossfreq_crDCC, function(x) x$BIC))[-2])

#calc QLIKES: 

sigmas_crDCC <- lapply(estacrossfreq_crDCC, function(x) x$vSigma2)

Qlikes <- matrix(0L, ncol = 10, nrow = 2516)

for(i in 6:9){
  for(j in 2:2516){
    Qlikes[j, i] <- QLIKE(sigmas_crDCC[[i]][,,j-1],calccov[[1]][[7]][,,j]*10000, 2)
  }
}

qlikes <- colMeans(Qlikes, na.rm = T)[6:9]

qlikes[2] <- mean(qlikes_crDCC, na.rm = T)
mean(qlikes, na.rm = T)








#------------------------------------------bootstrapping----------------------------------------------

library(boot)

library(np)
#uses the row to construct the block length. Therefore you should transpose it. 
ll <- tsboot(dailyretotc[,1], mean ,R = 1000, l = 97, sim = "fixed", endcorr = TRUE)
indexation <- boot.array(ll)

#rearranging calculated covariances: correct. 



#If we bootstrap the standardized residuals:
spec <- ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), variance.model = 
	list(model = 'realGARCH', garchOrder = c(2, 1)))
setbounds(spec)<-list(alpha2=c(-1,1))

lfittlt <- ugarchfit(spec, dailyretotc[,1]*100, solver = 'hybrid', realizedVol = xts(sqrt(calccov[[1]][[7]][1,1,]*10000), order.by = as.Date(getDates)))

lfitspy <- ugarchfit(spec, dailyretotc[,2]*100, solver = 'hybrid', realizedVol = 
  xts(sqrt(calccov[[1]][[7]][2,2,]*10000), order.by = as.Date(getDates)))


residuals <- matrix(cbind(residuals(lfittlt, standardize = T), residuals(lfitspy, standardize = T)), ncol=2)

hist(residuals[,2] ,breaks=50)

cor(residuals)

b.star(M[2,1,])

residualstest <- residuals[indexation[1, ], ]

head(residuals[indexation[1, ], ])


estimates_5min_crDCC <- matrix(0L, nrow = 1000, ncol = 4)


BS_dailyretotc <- matrix(dailyretotc, nrow=2516, ncol = 2)

BS_dailyretotc_l <- list()

resid <- list()

for(i in 1:1000){

BS_dailyretotc_l[[i]] <- BS_dailyretotc[indexation[i, ], ] 
resid[[i]] <- xts(residuals[indexation[i, ], ], order.by = as.Date(getDates))
}

head(resid[[1]])

#do it without bootstrapping standardized residuals first. 
for(j in 1:1000){
						#5min_covariances
	BS_P_crDCC <- P[,,indexation[j, ]] * 10000
	BS_N_crDCC <- N[,,indexation[j, ]] * 10000
	BS_M_crDCC <- M[,,indexation[j, ]] * 10000

  #ccc <- calccov[[1]][[7]][,,indexation[j, ]] * 10000
	BS_dailyretotc <- BS_dailyretotc_l[[j]] * 100

	estimates_5min_crDCC[j, ] <- Estimate_rcDCC(BS_dailyretotc, list(BS_P_crDCC, BS_N_crDCC, BS_M_crDCC),cov2 = NULL, 1:2516, T, resid[[j]])
	print(sprintf("Bootstrap %s", j))
	

}


#saveRDS(estimates_5min_crDCC, "bootstapestimates_5min_crDCC_stdresidualbootstrap.rds")





estimates_5min_crDCC_bootstraps <- readRDS("bootstapestimates_5min_crDCC_stdresidualbootstrap.rds")


rem1 <- which(estimates_5min_crDCC_bootstraps[,1] > 2.00067363e-08)
rem2 <- which(estimates_5min_crDCC_bootstraps[,2] > 0.00005)
rem3 <- which(estimates_5min_crDCC_bootstraps[,3] < 0.003)
rem4 <- which(estimates_5min_crDCC_bootstraps[,3]>0.65)
rem5 <- which(estimates_5min_crDCC_bootstraps[,4]<0.2)

inconsist <- unique(c(rem1, rem2, rem3, rem4, rem5))

estimates_5min_crDCC_bootstraps <- estimates_5min_crDCC_bootstraps[-inconsist, ]

demeanedbootstraps <- cbind(rep(crDCCests$vPar[1], 999), rep(crDCCests$vPar[2], 999),
  rep(crDCCests$vPar[3], 999), rep(crDCCests$vPar[4], 999))


demeanedbootstraps <- demeanedbootstraps - estimates_5min_crDCC_bootstraps

apply(demeanedbootstraps, MARGIN = c(2), FUN = function(x) sd(x))

bootstrap_stderror_5min_crDCC <- apply(estimates_5min_crDCC_bootstraps, MARGIN = c(2), FUN = function(x) sd(x))

bootstrap_stderror_5min_crDCC

mean(estimates_5min_crDCC[,3])
mean(estimates_5min_crDCC[,4])


tstat <- (crDCCests$vPar[1]/bootstrap_stderror_5min_crDCC[1])
lefttail <- 1-pnorm(abs(tstat))
righttail <- pnorm(-abs(tstat))
pobs <- lefttail + righttail
pobs #forkaster nul hypotesen om parameteren a^P er statistisk insignifikant ergo. H0: B = 0.

tstat <- (crDCCests$vPar[2]/bootstrap_stderror_5min_crDCC[2])
lefttail <- 1-pnorm(abs(tstat))
righttail <- pnorm(-abs(tstat))
pobs <- lefttail + righttail
pobs #accepterer hypotesen og a^N er statistisk insignificant

tstat <- (crDCCests$vPar[3]/bootstrap_stderror_5min_crDCC[3])
lefttail <- 1-pnorm(abs(tstat))
righttail <- pnorm(-abs(tstat))
pobs <- lefttail + righttail
pobs #forkaster hypotesen og a^M er statistisk significant

tstat <- (crDCCests$vPar[4]/bootstrap_stderror_5min_crDCC[4])
lefttail <- 1-pnorm(abs(tstat))
righttail <- pnorm(-abs(tstat))
pobs <- lefttail + righttail
pobs #forkaster hypotesen og b er statistisk significant



#---------------------------------------------------------------------------------------------------------------


#--------------------------SIMULATION TO SEE IF THE VOLATILITIES CONVERGES TOWARDS TRUE DATA------------------
#
#
#SIMULATION IS BIVARIATE WITH NO NOISE:



J <- 1

time <- 1000

#5min sampling
intradayticks <- (60 * 6.5)

N <- intradayticks * time

#vol
sigma <- 0.01 /252
sigma2 <- 0.98 /252
#drift
alpha <- 0 / 252
beta <- 0 / 252

#squareroot process
gamma <- 1
rho <- -0.4

simulation_prices <- brownian(J, N, sigma, sigma2, alpha, beta, gamma,rho, time)

dt <- length(simulation_prices[,1])/time


intodays <- list()

seq1 <- seq(0, length(simulation_prices[,1]), intradayticks)[-1]

intodays[[1]] <- cbind(simulation_prices[1:(intradayticks-1),1], simulation_prices[1:(intradayticks-1),2])


for(i in 2:(time-1)){


	intodays[[i]] <- cbind(simulation_prices[seq1[i]:seq1[i+1],1], simulation_prices[seq1[i]:seq1[i+1],2])

}

makesimret <- list()

for(i in 1:(time-1)){

	makesimret[[i]] <- cbind(diff(log(intodays[[i]][,1]))[-1] * 100, diff(log(intodays[[i]][,2]))[-1] * 100)

}

dailysimrets <- matrix(0L, nrow = time-1, ncol = 2)


for(i in 1:(time-1)){

	dailysimrets[i, ] <- cbind(log(intodays[[i]][1,1]) - log(intodays[[i]][length(intodays[[i]][,1]), 1]), log(intodays[[i]][1,2]) - log(intodays[[i]][length(intodays[[i]][,2]), 2])) * 100


}

rcovpos <- array(dim = c(2,2,time-1))
rcovneg <- array(dim = c(2,2,time-1))
rcovmix <- array(dim = c(2,2,time-1))
rcov <- array(dim = c(2,2,time-1))
prercov <- array(dim = c(2,2,time-1))


for(i in 2:(time-1)){

	rcovpos[,,i] <- realsemicov(makesimret[[i]], "P")
	rcovneg[,,i] <- realsemicov(makesimret[[i]], "N")
	rcovmix[,,i] <- realsemicov(makesimret[[i]], "M")
	rcov[,,i] <- realCov(makesimret[[i]])
	prercov[,,i] <- preavCov(matrix =xts(makesimret[[i]], order.by = as.Date(1:77)), T, T, F, 1)


}



simresults <- Estimate_rcDCC(mY = dailysimrets, covariance =list(rcovpos,rcovneg,rcovmix), getDates = 1:(time-1))

#getting correlations:

mean((simresults$aCor[2,1,] - rho)^2)

#cor
mean(simresults$aCor[2,1,])

#sigma
mean(sqrt(simresults$mSigma[,1])/252)
sigma
mean(sqrt(simresults$mSigma[,2])/252)
sigma2

#is in variances, now we need to transform it to vol/time
mean(sqrt(simresults$mSigma[,1])/time - sigma)^2 
mean(sqrt(simresults$mSigma[,2])/time - sigma2)^2 

#it converges towards the true parameters. 



mean(sqrt(prercov[1,1,-1]) / 252)
abs(mean(sqrt(prercov[2,2,-1]) / 252) - 0.98 /252)
abs(mean(sqrt(rcov[2,2,-1]) / 252) - 0.98 /252)










