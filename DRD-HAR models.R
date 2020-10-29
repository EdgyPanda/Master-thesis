#DRD-HAR model


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


#They use percentage realized variances for HARQ even though they do not specify it. 

getDates <- readRDS("getDates.rds")

library(HARModel)


calccov <- readRDS("calculatedcovariances.rds")

mergedfrequencies <- readRDS("mergedfrequencies.rds")

#5min correlations

fivemincorr <- array(0L, dim = c(2,2,2516))

for(i in 1:2516){

	fivemincorr[,,i] <- realCov(mergedfrequencies[[7]][[i]]*100, T)

}



test <- vech(fivemincorr[,,1])
fivemincorr[,,1]

#manual "unvech":
matrix(c(1,test[2],test[2],1), nrow=2,ncol=2)

#no need for vech operator when working in a bivariate framework. 



#-----------------------------------------CORRELATION HAR USING VECH OPERATOR------------------------------------------


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

fivemincorr <- xts(fivemincorr[2,1,], order.by = as.Date(getDates))


#one-day lag is today since we are in F_t-1. 
corrday <- fivemincorr
corrweek <- rowMeans(cbind(fivemincorr, mlag(fivemincorr,4,mean(fivemincorr))))
corrmonth <- rowMeans(cbind(fivemincorr, mlag(fivemincorr,21,mean(fivemincorr))))

ones <- rep(1, length(corrday)-22)

#trying to get the intercept term:

intercept <-  mean(fivemincorr) * (1- corrday - corrweek - corrmonth)
#needs 23 data points to initialize. This is in agreement with HARestimate from HARmodel library when removing intercept part.
tester <- lm(fivemincorr[23:2516] ~  0 +  
	corrday[22:(2516-1)] + corrweek[22:(2516-1)] + corrmonth[22:(2516-1)])

summary(tester)

sum(coef(tester))
coef(tester)[coef(tester) == 0]

c(-1,2,3,4)[c(-1,2,3,4) < 0]

#rolling forecast to see if params violate restrictions when posing no restrictions:
window <- 1000

params <- matrix(NaN, ncol=3, nrow = 2516)
for(i in (window+1):2516){

	fivemincorr1 <- fivemincorr[(i-window):(i-1)]

	corrday1 <- fivemincorr1
	corrweek1 <- rowMeans(cbind(fivemincorr1, mlag(fivemincorr1,4,mean(fivemincorr1))))
	corrmonth1 <- rowMeans(cbind(fivemincorr1, mlag(fivemincorr1,21,mean(fivemincorr1))))

	fivemincorr1 <- fivemincorr1[22:length(fivemincorr1)]
	corrday1 <- corrday1[22:length(corrday1)]
	corrweek1 <- corrweek1[22:length(corrweek1)]
	corrmonth1 <- corrmonth1[22:length(corrmonth1)]


	#end <- length(corrday1)
	
	est <- EstimatecorrHAR(list(fivemincorr1[2:length(fivemincorr1)], corrday1[1:(length(corrday1)-1)], 
		corrweek1[1:(length(corrweek1)-1)], corrmonth1[1:(length(corrmonth1)-1)]), 0)

	params[i,] <- est$vPar

	print(sprintf("%s", i))

}

params2 <- params[!is.na(params[,1]), ]

#parameter changes for the correlation estimate under five minute sampling frequency. 
library(ggplot2)

ggplot() + geom_line(aes(as.Date(getDates[(window+1):2516]), params2[,1], col="daily")) + 
geom_line(aes(as.Date(getDates[(window+1):2516]), params2[,2], col="weekly")) + 
geom_line(aes(as.Date(getDates[(window+1):2516]), params2[,3], col="monthly")) 


#solving the above using quadratic programing, from solnp function:
#problems with y, it had date indexation and thus misaligned by other data. 
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
	mini <- sum(mini)

	lOut <- list()
	lOut[["mini"]] <- mini 
	lOut[["minis"]] <- minis

	return(lOut)
}


data <- list(fivemincorr[23:2516], corrday[22:(2516-1)], corrweek[22:(2516-1)], corrmonth[22:(2516-1)])

ineqconstraint <- function(vPar, ...){

	dA <- vPar[1]  
	dB <- vPar[2]
	dC <- vPar[3]

	return(dA+dB+dC)

}

min.RSS <- function(vPar, ...){

	rss <- minimizingfunc(data, vPar)$mini

	return(rss)

}

EstimatecorrHAR <- function(data, trace=1, ineqfun = ineqconstraint, ineqLB = 0.00, ineqUB = 0.9999){
  
  dA = 0.45 
  dB  = 0.03185
  dC = 0.03
  ## vector of starting parameters
  vPar = c(dA, dB, dC)

  optimizer = solnp(vPar, fun = min.RSS, data = data, 
                    ineqfun = ineqfun, #the inequality constraint
                    ineqLB  = ineqLB, ## the inequality lower bound
                    ineqUB = ineqUB, ## the inequality lower bound, i.e. 0.0 <= a + b + c < 0.9999
                    ## lower and upper bounds for all parameters
                    LB = c(0, 0, 0), UB = c(0.9999, 0.9999, 0.9999),
                    control = list(tol = 1e-18, outer.iter = 800, inner.iter = 1200, delta = 1e-7,
                    trace = trace))
                    

  vPar <- optimizer$pars

  min <- tail(optimizer$values, 1)

  hessian <- optimizer$hessian

  scores <- matrix(0L, nrow=length(data[[1]]), ncol = 3)

  step <- 1e-5 * vPar

  for(i in 1:length(step)){

	h <- step[i]
    delta <- rep(0, length(vPar))
    delta[i] <- h
																
	loglikeminus <- minimizingfunc(data, vPar-delta)$minis
	loglikeplus <- minimizingfunc(data, vPar+delta)$minis

	scores[,i] <- (loglikeplus - loglikeminus)/(2*h)

  }

  J <- (t(scores) %*% scores)/length(data[[1]])

  I <- optimizer$hessian/length(data[[1]])

  I <- solve(I)[-1 ,-1]

  vars <- (I * J * I)/length(data[[1]])
  
  rse <- sqrt(diag(vars))


  t.stat <- vPar/rse




  lOut <- list()

  lOut[["vPar"]] <- vPar
  lOut[["MSE"]] <- min/2516 
  lOut[["rse"]] <- rse
  lOut[["hessian"]] <- hessian
  lOut[["Tstats"]] <- t.stat

  return(lOut)

 }

dta <- list(fivemincorr[23:2516], corrday[22:(2516-1)], corrweek[22:(2516-1)], corrmonth[22:(2516-1)])

tt <- EstimatecorrHAR(dta, 0)




#when you get your intercept to be correct, then this is how you recover correlations. 
#A super great inspiration can be found at patton 2016 matlab code. 

corrHAR <- cbind(intercept[22:2515], corrday[22:2515], corrweek[22:2515], corrmonth[22:2515]) %*% matrix(coef(tester))

(1-coef(tester)[2]-coef(tester)[3]-coef(tester)[4]) 




#########################################################################################################
#
#
#							HAR models estimations 5min with standard errors
#
#
#########################################################################################################


#can be estimated using package for time efficiency:

library(HARModel)


HAREstimate(RM, BPV = NULL, RQ = NULL, periods = c(1,5,22),
periodsJ = NULL, periodsRQ = NULL, type = "HAR",
insanityFilter = TRUE, h = 1)


HAR_5min_TLT <- HAREstimate(sqrt(calccov[[1]][[7]][1,1,]))


fiveminvol <- matrix(sqrt(calccov[[1]][[7]][1,1,]))
volday <- fiveminvol
volweek <- rowMeans(cbind(fiveminvol, mlag(fiveminvol,4,mean(fiveminvol))))
volmonth <- rowMeans(cbind(fiveminvol, mlag(fiveminvol,21,mean(fiveminvol))))


#needs 23 data points to initialize. This is in agreement with HARestimate from HARmodel library when removing intercept part.
tester <- lm(fiveminvol[23:2516] ~    
	volday[22:(2516-1)] + volweek[22:(2516-1)] + volmonth[22:(2516-1)])

summary(tester)



t(mergedfrequencies[[5]][[8]] * 100) %*% mergedfrequencies[[5]][[8]] * 100

lll <- t(mergedfrequencies[[5]][[8]]) %*% mergedfrequencies[[5]][[8]]