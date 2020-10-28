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



dataTLT <- readRDS("dataTLT.rds")

getDates <- unlist(lapply(dataTLT, function(x) as.character(index(x[1]))))

for(i in 1:length(getDates)){
	getDates[i] <- strsplit(getDates, " ")[[i]][1]
}


dailyretotc <- xts(t(sapply(mergedfrequencies[[10]], function(x) cbind(x[,1], x[,2]))), order.by = as.Date(getDates))

colnames(dailyretotc) <- c("TLT", "SPY")


P <- array(0L, dim =c(2,2,2516))
N <- array(0L, dim =c(2,2,2516))
M <- array(0L, dim =c(2,2,2516))

for(i in 1:2516){

	P[,,i] <- realsemicov(mergedfrequencies[[7]][[i]], "P")
	N[,,i] <- realsemicov(mergedfrequencies[[7]][[i]], "N")
	M[,,i] <- realsemicov(mergedfrequencies[[7]][[i]], "M")


}

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

	astar <- dAlpha * (samplercov * samplecov^(-1))

	mSigma[,,1] <- samplecov 

	#compute first observation of log-likelihood (we are minizing neg-log-likelihood):
	dLLK <- log(det(mSigma[,,1])) + (dailyret[1, , drop = F]) %*% solve(mSigma[,,1]) %*% t(dailyret[1, , drop = F])


	dLLKs <- numeric()
	dLLKs[1] <-   dLLK
	for(i in 2:days){

		mSigma[,,i] <- samplecov * (1 - astar - dBeta)  + dBeta * mSigma[,,i-1] + dAlpha * rcov[,,i-1]

		#neglog collection for score calculation
		dLLKs[i] <- 0.5 * d* log(2*pi) + 0.5 * (log(det(mSigma[,,i])) + 
		dailyret[i, , drop = F] %*% solve(mSigma[,,i]) %*% t(dailyret[i, , drop = F]))
		}

	fulldLLK <- - sum(dLLKs) #loglikelihood

	lOut <- list()

	lOut[["fulldLLK"]] <- fulldLLK
	lOut[["mSigma"]] <- mSigma
	lOut[["cov"]] <- astar
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

EstimateBivarGARCH <- function(dailyret, covariance, nobootstrap = T, vPar=NULL, ineqfun_GARCH = ineqfun_GARCH_BIVAR, ineqLB = 0.00, ineqUB = 0.9999){
  
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

  if(nobootstrap){
  #CALCULATING ROBUST STD. ERRORS WHEN NO COVARIANCE TARGETING:

  scores <- matrix(0L, nrow=length(lT), ncol = length(vPar))

  step <- 1e-5 * vPar
  # This follows directly from Kevin Sheppard who uses percentage returns, Moreover scale open-to-close with 24/6.5.
  for(i in 1:length(step)){

	h <- step[i]
    delta <- rep(0, length(vPar))
    delta[i] <- h
																
	loglikeminus <- BivarGARCHFilter(lT,dailyret, vPar-delta, covariance)$dLLKs
	loglikeplus <- BivarGARCHFilter(lT,dailyret, vPar+delta, covariance,)$dLLKs

	scores[,i] <- (loglikeplus - loglikeminus)/(2*h)

  }

  J <- (t(scores) %*% scores)/length(lT)

  I <- optimizer$hessian/length(lT)

  I <- solve(I)[-1 ,-1]

  vars <- (I * J * I)/length(lT)
  
  rse <- sqrt(diag(vars))

  ## compute filtered volatility
  Filter <- BivarGARCHFilter(lT, dailyret, c(vPar[1], vPar[2]), covariance)

  vSigma2 <- Filter$mSigma
  
  sampleH <- Filter$sampleH

  ## Compute the daily Average BIC
  iT = length(lT)
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

  return(vPar)
  
}

#it does not matter whether you use percentage returns with realcov calculated on percentage returns or 
#standard returns. 

library(tictoc)
tic()
tester <- EstimateBivarGARCH(dailyretotc, calccov[[1]][[7]], nobootstrap = F)
toc()


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

saveRDS(estimates_5min_rBG, "bootstrapestimates_5min_rBG.rds")


bootstrap_stderror_5min_rBG <- apply(estimates_5min_rBG, MARGIN = c(2), FUN = function(x) sd(x))



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

	astar <- (dAlphaP * sampleP + dAlphaN * sampleN + dAlphaM * sampleM) * samplecov^(-1)

	dLLKs <- numeric()

	dLLKs[1] <-  dLLK

	for(i in 2:days){

		mSigma[,,i] <- samplecov * (1 - astar - dBeta)  + dBeta * mSigma[,,i-1] + 
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
	lOut[["cov"]] <- astar

	return(lOut)
	
	
}

leltest334 <- BivarGARCHFilterContAsym(mergedfrequencies[[8]],dailyretotc, c(0.02,0.01,0.01,0.94), list(P, N, M))




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
                    ) 
  
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


#5min frequency stemming from P, N, M. 
leltest4 <- EstimateBivarGARCHContAsym(dailyretotc, list(P,N,M))


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

	estimates_5min_crBG[j, ] <- EstimateBivarGARCHContAsym(BS_dailyretotc, list(BS_P_crBG, BS_N_crBG, BS_M_crBG))$vPar
	print(sprintf("Bootstrap %s", j))
	

}

saveRDS(estimates_5min_crBG, "bootstapestimates_5min_crBG.rds")



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






N[,,1292]
ts.plot(lel4$aCor[2,1,])
covariance <- list(P, N, M)

lel4 <- rcDCCFilter(dailyretotc,-0.06504, -0.05120,  0.25439,  0.86085,cor(dailyretotc), covariance)

lel4$aCor[2,1,]


mean(lel4$astar[2,1,])

is.nan(lel4$astar[2,1,])

lel15 <- rDCCFilter(dailyretotc, 0.1, 0.2, cor(dailyretotc), calccov[[1]][[8]])

ts.plot(lel15$aCor[2,1,])

ts.plot(lel4$aCor[2,1,])

0.1 * diag(sqrt(1/diag(calccov[[1]][[9]][,, 1]))) %*% calccov[[1]][[9]][,,1] %*% diag(sqrt(1/diag(calccov[[1]][[9]][,, 1])))

matrix(c(diag[1],0,0,diag[2]), nrow=2,ncol=2)^(-0.5) %*% calccov[[1]][[9]][,,1] %*% 
matrix(c(diag[1],0,0,diag[2]), nrow=2,ncol=2)^(-0.5)



#testing with simple realGARCH(1,1) using rugarch package. 
#you need getDates for fit function in rugarch to work. 
#moreover the covariances you use in the DCC model is also used in the univariate models.
Estimate_rDCC <- function(mY, covariance, getDates) {
  
  ## estimate the marginal models
  require(Rsolnp)
  require(rugarch) 
 
  #Marginal garch specification. THIS WORKS ONLY IN BIVARIATE SETUP. 
  
  #list where marginal models are stored

  cov1 <- sqrt(covariance[1,1,]) * 100
  cov2 <- sqrt(covariance[2,2,]) * 100

  spec <- ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), variance.model = 
	list(model = 'realGARCH', garchOrder = c(2, 1)))

  setbounds(spec)<-list(alpha2=c(-1,1))

  #spec1 <- ugarchfit(spec, mY[,1], solver = 'hybrid', realizedVol = 
   #xts(cov1, order.by = as.Date(getDates)), 
	#fit.control = list(stationarity = 1, fixed.se = 1))

  #spec2 <- ugarchfit(spec, mY[,2], solver = 'hybrid', realizedVol = 
	#xts(cov2, order.by = as.Date(getDates)), 
	#fit.control = list(stationarity = 1, fixed.se = 1))

  mspec <- multispec( replicate(spec, n=2) )

  lfit <- multifit(mspec, mY * 100, solver = 'hybrid', realizedVol = 
  xts(cbind(cov1, cov2), order.by = as.Date(getDates)))


  mEta <- residuals(lfit, standardize = T)
  ####################################################
  
  ## maximization of the DCC likelihood
  
  #initial parameters
  vPar = c(0.19, 0.79)
  
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

  #standard errors 
  se <- solve(optimizer$hessian)
  se <- matrix(sqrt(diag(se))[-1], ncol=length(vPar), nrow=1)

  #extract univariate volatilities
  mSigma = sigma(lfit)^2
  
  #extract univariate estimated parameters
  mCoef = coef(lfit)

  colnames(mCoef) <- colnames(mY)
  
  #compute the likelihood of the volatility  part
  dLLK_V = sum(likelihood(lfit)) 
  
  
  #compute the total likelihood
  dLLK = dLLK_V + dLLK_C
  
  ## Compute z_t
  aCor = Filter[["aCor"]]
  iT = nrow(mY)
  
  mZ = matrix(0, iT, ncol(mY))
  
  for (t in 1:iT) {
    mZ[t, ] = solve(chol(aCor[,,t])) %*% as.numeric(mEta[t, ])
  }
    
  lOut = list()

  #output the results
  lOut[["dLLK"]] = dLLK
  lOut[["mCoef"]] = mCoef
  lOut[["vPar"]] = vPar
  lOut[["mSigma"]] = mSigma
  lOut[["aCor"]] = aCor
  lOut[["mEta"]] = mEta
  lOut[["mZ"]] = mZ
  lOut[["se"]] = se
  lOut[["seeall"]] = coef(lfit)
  #you can get standard errors for GARCH model by fitting each univariately
  return(lOut)
  
}







lel5 <- Estimate_rDCC(dailyretotc, calccov[[1]][[7]], getDates)

#------------------------------------------------realGARCH 5min estimates and standard errors--------------------

library(rugarch)

spec <- ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), variance.model = 
	list(model = 'realGARCH', garchOrder = c(2,1))) #he changed the order in the code see starch. 

setbounds(spec)<-list(alpha2=c(-1,1)) #you need to reset the bounds for alpha2 otherwise it will be zero since it 
#originally is set between (0,1). 


mspec <- multispec(replicate(spec, n=2))

setbounds(mspec)<-replicate(list(alpha2=c(-1,1)), n=2)

lfit <- multifit(mspec, dailyretotc, solver = 'hybrid', realizedVol = 
  xts(cbind(sqrt(calccov[[1]][[7]][1,1,]), sqrt(calccov[[1]][[7]][2,2,])), order.by = as.Date(getDates)))

coef(lfit)


#TLT
lfit2 = ugarchfit(spec, dailyretotc[,1] * 100, solver = 'hybrid', realizedVol = 
	xts(sqrt(calccov[[1]][[7]][1,1,]) * 100, order.by = as.Date(getDates)), 
	fit.control = list(stationarity = 1, fixed.se = 1))


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
	xts(sqrt(calccov[[1]][[i]][1,1,]) * 100, order.by = as.Date(getDates)), 
	fit.control = list(stationarity = 1, fixed.se = 1))

	spy <- ugarchfit(spec, dailyretotc[,2] * 100, solver = 'hybrid', realizedVol = 
	xts(sqrt(calccov[[1]][[i]][2,2,]) * 100, order.by = as.Date(getDates)), 
	fit.control = list(stationarity = 1, fixed.se = 1))

	estimatesacross_TLT[i, ] <- coef(tlt)

	tlt.loglike[i, ] <- c('logL' = sum(-tlt@fit$log.likelihoods), 'pLogL' = sum(-tlt@fit$partial.log.likelihoods))
	#SPY
	spy.loglike <- c(sum(-spy@fit$log.likelihoods), sum(-spy@fit$partial.log.likelihoods))
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
estimatesacross_SPY_mean <- colMeans(estimatesacross_SPY)

TLT_90percentquantile <- apply(estimatesacross_TLT, MARGIN = c(2), FUN = function(x) quantile(x, 0.9, na.rm=TRUE))
TLT_10percentquantile <- apply(estimatesacross_TLT, MARGIN = c(2), FUN = function(x) quantile(x, 0.1, na.rm=TRUE))

TLT_stats <- rbind(estimatesacross_TLT_mean, TLT_10percentquantile, TLT_90percentquantile)

colMeans(tlt.loglike, na.rm = T)
quantile(tlt.loglike[,2], 0.9, na.rm = T)
quantile(tlt.loglike[,2], 0.1, na.rm = T)

tlt.completeloglike <- rowSums(tlt.loglike)
mean(tlt.completeloglike, na.rm = T)
quantile(tlt.completeloglike, 0.9, na.rm = T)
quantile(tlt.completeloglike, 0.1, na.rm = T)



#SPY
lfit = ugarchfit(spec, dailyretotc[,2] * 100, solver = 'hybrid', realizedVol = 
	(xts(sqrt(calccov[[1]][[7]][2,2,]), order.by = as.Date(getDates))) * 100, 
	fit.control = list(stationarity = 1, fixed.se = 1))


rugarch.LL = c('logL' = sum(-lfit@fit$log.likelihoods), 'pLogL' = sum(-lfit@fit$partial.log.likelihoods))



lfit <- ugarchfit(spec, spyreal[,1] , solver = 'hybrid', realizedVol = 
	spyreal[,2] , 
	fit.control = list(stationarity = 1, fixed.se = 1))


spec = ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), variance.model = list(model = 'realGARCH', garchOrder = c(2, 1)))


fit = ugarchfit(spec, spyreal[, 1] * 100, solver = 'hybrid', realizedVol = spyreal[,2] * 100)

rugarch.LL = c('logL' = sum(-fit@fit$log.likelihoods[1:1492]), 'pLogL' = sum(-fit@fit$partial.log.likelihoods[1:1492]))

sum(rugarch.LL)


data(spyreal)
head(spyreal)

sqrt(head(calccov[[1]][[7]][2,2,]))

#------------------------------------------------------------------------------------------------------------------


Estimate_rcDCC <- function(mY, covariance, getDates) {
  
  ## estimate the marginal models
  require(Rsolnp)
  require(rugarch) 
 
  #Marginal garch specification. THIS WORKS ONLY IN BIVARIATE SETUP. 
  
  cov <- covariance[[1]] + covariance[[2]] + covariance[[3]]

  #list where marginal models are stored
  spec <- ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), variance.model = 
	list(model = 'realGARCH', garchOrder = c(2, 1)))

  specforrobust1 <- ugarchfit(spec, dailyretotc[,1], solver = 'hybrid', realizedVol = 
	xts(sqrt(cov[1,1,]), order.by = as.Date(getDates)))

  specforrobust2 <- ugarchfit(spec, dailyretotc[,2], solver = 'hybrid', realizedVol = 
	xts(sqrt(cov[2,2,]), order.by = as.Date(getDates)))

  mspec <- multispec( replicate(spec, n=2) )

  lfit <- multifit(mspec, mY, solver = 'hybrid', realizedVol = 
  xts(cbind(sqrt(cov[1,1,]), sqrt(cov[2,2,])), order.by = as.Date(getDates)))
  
  mEta <- residuals(lfit, standardize = T)
  ####################################################
  
  ## maximization of the DCC likelihood
  
  #initial parameters
  vPar = c(0.001, 0.001, 0.20, 0.80)
  
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
  
  #Extract the likelihood of the correlation part
  dLLK_C = -tail(optimizer$values, 1)
  
  #Filter the dynamic correlation using the estimated parameters
  Filter = rcDCCFilter(mEta, vPar[1], vPar[2], vPar[3], vPar[4], mQ, covariance)

  #standard errors 
  se <- solve(optimizer$hessian)
  se <- matrix(sqrt(diag(se))[-1], ncol=length(vPar), nrow=1)

  #extract univariate volatilities
  mSigma = sigma(lfit)^2
  
  #extract univariate estimated parameters
  mCoef = coef(lfit)

  colnames(mCoef) <- colnames(mY)
  
  #compute the likelihood of the volatility  part
  dLLK_V = sum(likelihood(lfit)) 
  
  
  #compute the total likelihood
  dLLK = dLLK_V + dLLK_C
  
  ## Compute z_t
  aCor = Filter[["aCor"]]
  iT = nrow(mY)
  
  mZ = matrix(0, iT, ncol(mY))
  
  for (t in 1:iT) {
    mZ[t, ] = solve(chol(aCor[,,t])) %*% as.numeric(mEta[t, ])
  }
    
  lOut = list()

  #output the results
  lOut[["dLLK"]] = dLLK
  lOut[["mCoef"]] = mCoef
  lOut[["vPar"]] = vPar
  lOut[["mSigma"]] = mSigma
  lOut[["aCor"]] = aCor
  lOut[["mEta"]] = mEta
  lOut[["mZ"]] = mZ
  lOut[["se"]] = se
  lOut[["seeall"]] = list(specforrobust1, specforrobust2) 
  return(lOut)
  
}










library(numDeriv, MASS)
library(MASS)

tt <- jacobian(ObjFBivarGARCH, x = lel2$vPar, lT = mergedfrequencies[[8]], dailyret = dailyretotc)

t(tt) %*% ginv(lel2$Hessian)[2:3, 2:3] %*% (tt)





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