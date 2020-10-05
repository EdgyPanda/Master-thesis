#Volatility models 


source("Previous functions/projectfunctions.R")
source("functions.R")

library(xts)
library(highfrequency)
library(matlib)
library(MCS)



calccov <- readRDS("calculatedcovariances.rds")

mergedfrequencies <- readRDS("mergedfrequencies.rds")


#The scalar bivariate GARCH model estimation


GARCHFilter <- function(vY, dOmega, dAlpha, dBeta) {
  
  ## number of observations
  iT = length(vY)
  
  ## initialize the variance
  vSigma2 = numeric(iT)
  
  ## set the variance at time t = 1 to the empirical variance of the first 0.1 * iT observations
  vSigma2[1] = var(vY[1:round(iT * 0.1)])
  
  # ## compute the log--likelihood of the first obs
  dLLK = dnorm(vY[1], 0, sqrt(vSigma2[1]), log = TRUE)
  # 
  ## loop over iT
  for (t in 2:iT) {
    # update the volatility
    vSigma2[t] = dOmega + dAlpha * vY[t-1]^2 + dBeta * vSigma2[t - 1]
    # add the log-likelihood contribution at time t
    dLLK = dLLK + dnorm(vY[t], 0, sqrt(vSigma2[t]), log = TRUE)
  }
  
  lOut = list()
  lOut[["dLLK"]] = dLLK
  lOut[["vSigma2"]] = vSigma2
  
  return(lOut)
  
}

#lT is list of intraday returns as xts object with each list element being each day
#dailyret is daily returns either open-to-close or close-to-close, however, the model
#will target the unconditional realized covariance matrix accordingly, so it's a good idea
#to specify what you're using. Start with open-to-close. And they should obviouslý have the same amount of days. 
#if covariance is not null, then realized covariance estimations will  be bypassed. 
#
#
#
#CURRENTLY IT IS MAXIMIZING LIKELIHOOD. YOU NEED TO NEGATE IT. 
BivarGARCHFilter <- function(lT, dailyret, dAlpha, dBeta, covariance = NULL){

	days <- length(lT)

	d <- ncol(lT[[1]])

	samplecov <- cov(dailyret)

	if(is.null(covariance)){

		rcov <- array(0L, dim = c(2, 2, days))

		for(i in 1:days){

			rcov[,,i] <- realCov(lT[[i]])

		}

		samplercov <- matrix(
			c(mean(rcov[1,1,]),
			  mean(rcov[2,1,]),
			  mean(rcov[2,1,]),
			  mean(rcov[2,2,])),
				ncol=2, nrow=2)

		mSigma <- array(0L, dim = c(2, 2, days))

		#initializing it with unconditional mean. 

		mSigma[,,1] <- samplercov

		#compute first observation of log-likelihood:
		dLLK <-  - log(det(mSigma[,,1])) - dailyret[1, , drop = F] %*% inv(mSigma[,,1]) %*% t(dailyret[1, , drop = F])

		for(i in 2:days){

			astar <- dAlpha * samplercov %*% inv(samplecov)

			mSigma[,,i] <- (1 - astar - dBeta) %*% samplecov + dBeta * mSigma[,,i-1] + dAlpha * rcov[,,i-1]

			dLLK <- dLLK - det(mSigma[,,i]) - dailyret[i, , drop = F] %*% inv(mSigma[,,i]) %*% t(dailyret[i, , drop = F])

		}

		fulldLLK <- -days * d/2 * log(2*pi) + dLLK  

		lOut <- list()

		lOut[["fulldLLK"]] <- fulldLLK
		lOut[["mSigma"]] <- mSigma

		return(lOut)

	}
}


dataTLT <- readRDS("dataTLT.rds")

getDates <- unlist(lapply(dataTLT, function(x) as.character(index(x[1]))))

for(i in 1:length(getDates)){
	getDates[i] <- strsplit(getDates, " ")[[i]][1]
}






dailyretotc <- xts(t(sapply(mergedfrequencies[[10]], function(x) cbind(x[,1], x[,2]))), order.by = as.Date(getDates))

colnames(dailyretotc) <- c("TLT", "SPY")

lel <- BivarGARCHFilter(mergedfrequencies[[8]], dailyretotc, 0.02, 0.98)


lel$fulldLLK


another <- numeric()
for(i in 1:2516){
another[i] <- det(lel$mSigma[,,i])
}


test <- matrix(c(mean(calccov[[1]][[7]][1,1,]), mean(calccov[[1]][[7]][2,1,]), 
	mean(calccov[[1]][[7]][2,1,]), mean(calccov[[1]][[7]][2,2,])), ncol=2,nrow=2)

test[1, , drop = F]
test

(mergedfrequencies[[7]][[1]][1, , drop = F]) %*% solve(test) %*% t(mergedfrequencies[[7]][[1]][1, , drop = F])


test2 <- array(rnorm(100), c(2,2,25))

test2[1, , , drop = F]