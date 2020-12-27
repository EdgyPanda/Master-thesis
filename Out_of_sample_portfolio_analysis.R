####################################################################
#
#				OUT OF SAMPLE PORTFOLIO ANALYSIS
#
####################################################################

library(PerformanceAnalytics)
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


turnover <- function(weights, portret, assetret){

	turnover <- numeric()
	turnover[1] <- abs(sum(weights[1, ])) #only when initial weights are zero. 

	#lagged <- weights * ((1 + assetret)/ (1+portret))

	for(i in 2:nrow(assetret)){

		turnover[i] <- sum(abs( weights[i, ] -  weights[i-1, ] * (1 + assetret[i-1, ])/(1+portret[i-1, ]) ))


	}
	return(turnover)
}


turnover.simple <- function(weights){

	turnover <- numeric()
	turnover[1] <- abs(sum(weights[1, ])) #only when initial weights are zero. 

	for(i in 2:nrow(weights)){

		turnover[i] <- sum(abs( weights[i, ] -  weights[i-1, ]))


	}
	return(turnover)
}



portconcentration <- function(weights){

	concentration <- numeric()

	for(i in 1:nrow(weights)){

		concentration[i] <- sqrt(sum(weights[i, ]^2))

	}
	return(concentration)
}

#####################################################################################################
#
#
#                     Calculating the risk-free rate from a 3-month T-bill:
#
#
#####################################################################################################


#Calculating the risk-free rate from a 3-month T-bill:
#3month T-bill. Has no coupons 
dtb3 <- read.csv("DTB3.csv", header =  T)

dtb3[,2] <- as.numeric(levels(dtb3[,2]))[dtb3[,2]]

dtb3[,2] <- na.approx(dtb3[,2])       #dtb3[!is.na(dtb3[,2]), ]



#this is how you should do it:
logriskfreerate <- log(1 + dtb3[,2]/(100*365)) 

logriskfreerate <- xts(logriskfreerate, order.by = as.Date(dtb3[,1]))

intersectMulti <- function(x=list()){
 for(i in 2:length(x)){
    if(i==2) foo <- x[[i-1]]
    foo <- intersect(foo,x[[i]]) #find intersection between ith and previous
 }
 return(x[[1]][match(foo, x[[1]])]) #get original to retain format
}

returns_TLT <- dailyretotc[,1]

indexes <- intersectMulti(list(index(logriskfreerate), index(returns_TLT)))

logriskfreerate <- logriskfreerate[indexes]


#you need to downscale the daily risk-free rate or upscale the percentage returns to daily close-to-close returns by multiplying by (24/6.5)

library(alphavantager)
library(ggplot2)
#source("functions.R")
source("APIKEY.R")

av_api_key(apikey)

TLT <- as.data.frame(av_get(symbol = "TLT", av_fun = "TIME_SERIES_DAILY", outputsize = "full"))
rownames(TLT) <- TLT$timestamp
SPY <- as.data.frame(av_get(symbol = "SPY", av_fun = "TIME_SERIES_DAILY", outputsize = "full"))
rownames(SPY) <- SPY$timestamp

returns_TLT <- as.xts(diff(log(TLT[,4])), order.by = as.Date(TLT[,1], format='%d/%m/%Y')[-1])
returns_SPY <- as.xts(diff(log(SPY[,4])), order.by = as.Date(SPY[,1])[-1])
returns_TLT <-  returns_TLT[seq(from= as.Date('2010-01-02'), to = as.Date('2019-12-31'), by=1), ]
returns_SPY <- returns_SPY[seq(from= as.Date('2010-01-02'), to = as.Date('2019-12-31'), by=1), ]
merged_ret <- cbind(returns_TLT, returns_SPY)

colMeans(merged_ret[1000:2516]) * 100 * 252

apply(merged_ret, MARGIN = c(2), FUN = function(x) sd(x*100)*sqrt(252))




excessret_5sec <- list()
excessret_5sec_2 <- list()
excessret_15sec_2 <- list()
excessret_5min_2 <- list()
excessret_daily <- list()

excessret_daily_2 <- list()
for(i in 1:2516){
  
  excessret_5sec[[i]] <- mergedfrequencies[[2]][[i]]  * (24/6.5) - as.numeric(logriskfreerate[i, ]) 

  excessret_5sec_2[[i]] <- mergedfrequencies[[2]][[i]]  - as.numeric(logriskfreerate[i, ] * (1/(24*12*60)))   
  excessret_15sec_2[[i]] <- mergedfrequencies[[3]][[i]]  - as.numeric(logriskfreerate[i, ] * (1/(24*4*60)))   

  excessret_5min_2[[i]] <- mergedfrequencies[[7]][[i]]  - as.numeric(logriskfreerate[i, ] * (1/(24*12)))  


  excessret_daily[[i]] <- dailyretotc[i, ]    - as.numeric(logriskfreerate[i, ])  * 6.5/24 #you have to scale it down to intraday hours
}



#merging intra-day returns and logfree-rate

mergedret_5sec <- list()
mergedret_15sec <- list()
mergedret_5min <- list()
mergedret_daily <- matrix(0L, ncol = 3, nrow = 2516)
mergedret_daily2 <- matrix(0L, ncol = 3, nrow = 2516)
#* 6.5/24) 
for(i in 1:2516){
  
  temp <- as.numeric(logriskfreerate[i, ] * ((6.5*12*60)/(24*12*60)))    
  temp2 <- as.numeric(logriskfreerate[i, ] * 1/(24*12)) 
  temp3 <-  as.numeric(logriskfreerate[i, ] * ((6.5*4*60)/(24*4*60))) 
  temp4 <- as.numeric(logriskfreerate[i, ]) * 6.5/24 


  mergedret_5sec[[i]] <- cbind(excessret_5sec_2[[i]], rep(temp, length(excessret_5sec_2[[i]][,1])))  
  mergedret_15sec[[i]] <- cbind(excessret_15sec_2[[i]], rep(temp3, length(excessret_15sec_2[[i]][,1])))  
  mergedret_5min[[i]] <- cbind(excessret_5min_2[[i]], rep(temp2, length(excessret_5min_2[[i]][,1])))  
  mergedret_daily[i, ] <- cbind(excessret_daily[[i]], temp4)  
  colnames(mergedret_5sec[[i]]) <- c("TLT", "SPY", "RF")
  colnames(mergedret_5min[[i]]) <- c("TLT", "SPY", "RF")
  colnames(mergedret_15sec[[i]]) <- c("TLT", "SPY", "RF")
  colnames(mergedret_daily) <- c("TLT", "SPY", "RF")

}



#####################################################################################################
#
#
#                     				Extracting covariance forecasts.				
#
#
#####################################################################################################



sigmafivemincovariances <- readRDS("covariances5minrBGcrBGrDCCcrDCC.rds")


sigmas_rBG <- sigmafivemincovariances$sigmas_rBG
sigmas_rDCC <- sigmafivemincovariances$sigmas_rDCC
sigmas_crBG <- sigmafivemincovariances$sigmas_crBG
sigmas_crDCC <- sigmafivemincovariances$sigmas_crDCC

sigmas_HARQ <- array(0L, dim = c(2,2,2516))
sigmas_HARQF <- array(0L, dim = c(2,2,2516))
sigmas_HAR <- array(0L, dim = c(2,2,2516))
sigmas_HARJ <- array(0L, dim = c(2,2,2516))
sigmas_HARQJ <- array(0L, dim = c(2,2,2516))
sigmas_HAR_daily <- array(0L, dim = c(2,2,2516))
sigmas_CHARQ <- array(0L, dim = c(2,2,2516))
sigmas_CHAR <- array(0L, dim = c(2,2,2516))
sigmas_RM_daily <- array(0L, dim = c(2,2,2516))

Correlations_measures <- readRDS("Correlations_measures.rds")
Quarticities <- readRDS("Quarticity_estimators.rds")
proxycorrelation <- array(0L, dim = c(2,2,2516))


for(i in 1:2516){
  
  proxycorrelation[,,i] <- realCov(mergedfrequencies[[7]][[i]]*100, T)
}

covDAILY <- rollapply(calccov[[1]][[10]][2,1, ], 2, mean, align = 'right', fill = mean(calccov[[1]][[10]][2,1, ])) 
TLTDAILY <- rollapply(calccov[[1]][[10]][1,1, ], 2, mean, align = 'right', fill = mean(calccov[[1]][[10]][1,1, ])) 
SPYDAILY <- rollapply(calccov[[1]][[10]][2,2, ], 2, mean, align = 'right', fill = mean(calccov[[1]][[10]][2,2, ])) 
correlationDAILY <- as.numeric(covDAILY / sqrt(TLTDAILY * SPYDAILY))


window <- 1000
for(i in (window + 1):2516){

	#############################################################################################
	#
	#
	#											PROXIES
	# 
	#
	#############################################################################################

	RV5min_TLT <- calccov[[1]][[7]][1,1,(i-window):(i-1)] * 10000
	RV5min_SPY <- calccov[[1]][[7]][2,2,(i-window):(i-1)] * 10000

	RV5min_TLT <- RV5min_TLT[22:length(RV5min_TLT)]
	RV5min_SPY <- RV5min_SPY[22:length(RV5min_SPY)]

	#proxy covariance for qlike losses
	RCov5min <- calccov[[1]][[7]][,,i] * 10000

	proxycor <- proxycorrelation[2,1, (i-window):(i-1)]
	#---------------------------------------------------------------------------------------------


	#############################################################################################
	#
	#
	#							PRELIMINARIES (LAGS AND QUARTICITIES)
	# 
	#
	#############################################################################################

	#THESE HAVE BEEN CHANGES ASWELL!!

	correlation <- Correlations_measures[[2]][(i-window):(i-1), 8]
	correlationjumprobust <- Correlations_measures[[3]][(i-window):(i-1), 6]
	end <- length(RV5min_TLT)
	

	#Lags

	#voltlt under RCOV 5min : 

	vol_TLT <- matrix((calccov[[1]][[7]][1,1,(i-window):(i-1)])) * 10000  

	#vol_TLT <- matrix((calccov[[8]][[2]][1,1,(i-window):(i-1)])) * 10000  
	volday_TLT <- vol_TLT
	volweek_TLT <- rowMeans(cbind(vol_TLT, mlag(vol_TLT,4,mean(vol_TLT))))
	volmonth_TLT <- rowMeans(cbind(vol_TLT, mlag(vol_TLT,21,mean(vol_TLT))))

	vol_TLT <- vol_TLT[22:length(vol_TLT)]
	volday_TLT <- volday_TLT[22:length(volday_TLT)]
	volweek_TLT <- volweek_TLT[22:length(volweek_TLT)]
	volmonth_TLT <- volmonth_TLT[22:length(volmonth_TLT)]



	#volspy under RCOV 5min : 

	vol_SPY <- matrix((calccov[[1]][[7]][2,2,(i-window):(i-1)])) * 10000  


	#vol_SPY <- matrix((calccov[[8]][[2]][2,2,(i-window):(i-1)])) * 10000  
	volday_SPY <- vol_SPY
	volweek_SPY <- rowMeans(cbind(vol_SPY, mlag(vol_SPY,4,mean(vol_SPY))))
	volmonth_SPY <- rowMeans(cbind(vol_SPY, mlag(vol_SPY,21,mean(vol_SPY))))

	vol_SPY <- vol_SPY[22:length(vol_SPY)]
	volday_SPY <- volday_SPY[22:length(volday_SPY)]
	volweek_SPY <- volweek_SPY[22:length(volweek_SPY)]
	volmonth_SPY <- volmonth_SPY[22:length(volmonth_SPY)]





	#DAILY LAGS FOR PROXY!!
	vol_TLTDAILY <- matrix((calccov[[1]][[10]][1,1,(i-window):(i-1)])) * 10000  
	volday_TLTDAILY <- vol_TLTDAILY
	volweek_TLTDAILY <- rowMeans(cbind(vol_TLTDAILY, mlag(vol_TLTDAILY,4,mean(vol_TLTDAILY))))
	volmonth_TLTDAILY <- rowMeans(cbind(vol_TLTDAILY, mlag(vol_TLTDAILY,21,mean(vol_TLTDAILY))))

	vol_TLTDAILY <- vol_TLTDAILY[22:length(vol_TLTDAILY)]
	volday_TLTDAILY <- volday_TLTDAILY[22:length(volday_TLTDAILY)]
	volweek_TLTDAILY <- volweek_TLTDAILY[22:length(volweek_TLTDAILY)]
	volmonth_TLTDAILY <- volmonth_TLTDAILY[22:length(volmonth_TLTDAILY)]



	vol_SPYDAILY <- matrix((calccov[[1]][[10]][2,2,(i-window):(i-1)])) * 10000  
	volday_SPYDAILY <- vol_SPYDAILY
	volweek_SPYDAILY <- rowMeans(cbind(vol_SPYDAILY, mlag(vol_SPYDAILY,4,mean(vol_SPYDAILY))))
	volmonth_SPYDAILY <- rowMeans(cbind(vol_SPYDAILY, mlag(vol_SPYDAILY,21,mean(vol_SPYDAILY))))

	vol_SPYDAILY <- vol_SPYDAILY[22:length(vol_SPYDAILY)]
	volday_SPYDAILY <- volday_SPYDAILY[22:length(volday_SPYDAILY)]
	volweek_SPYDAILY <- volweek_SPYDAILY[22:length(volweek_SPYDAILY)]
	volmonth_SPYDAILY <- volmonth_SPYDAILY[22:length(volmonth_SPYDAILY)]




	#jump robust lags:

	#bpcov 5min 

	ICvol_TLT <- matrix((calccov[[5]][[7]][1,1,(i-window):(i-1)])) * 10000  

	#ICvol_TLT <- matrix((calccov[[6]][[3]][1,1,(i-window):(i-1)])) * 10000  
	ICvolday_TLT <- ICvol_TLT
	ICvolweek_TLT <- rowMeans(cbind(ICvol_TLT, mlag(ICvol_TLT,4,mean(ICvol_TLT))))
	ICvolmonth_TLT <- rowMeans(cbind(ICvol_TLT, mlag(ICvol_TLT,21,mean(ICvol_TLT))))

	ICvol_TLT <- ICvol_TLT[22:length(ICvol_TLT)]
	ICvolday_TLT <- ICvolday_TLT[22:length(ICvolday_TLT)]
	ICvolweek_TLT <- ICvolweek_TLT[22:length(ICvolweek_TLT)]
	ICvolmonth_TLT <- ICvolmonth_TLT[22:length(ICvolmonth_TLT)]



	#bpcov 5min 

	ICvol_SPY <- matrix((calccov[[5]][[7]][2,2,(i-window):(i-1)])) * 10000  


	#ICvol_SPY <- matrix((calccov[[6]][[3]][2,2,(i-window):(i-1)])) * 10000  
	ICvolday_SPY <- ICvol_SPY
	ICvolweek_SPY <- rowMeans(cbind(ICvol_SPY, mlag(ICvol_SPY,4,mean(ICvol_SPY))))
	ICvolmonth_SPY <- rowMeans(cbind(ICvol_SPY, mlag(ICvol_SPY,21,mean(ICvol_SPY))))

	ICvol_SPY <- ICvol_SPY[22:length(ICvol_SPY)]
	ICvolday_SPY <- ICvolday_SPY[22:length(ICvolday_SPY)]
	ICvolweek_SPY <- ICvolweek_SPY[22:length(ICvolweek_SPY)]
	ICvolmonth_SPY <- ICvolmonth_SPY[22:length(ICvolmonth_SPY)]




	#RQ and TRQ from percentage returns
	rq_TLT <- Quarticities$rq_TLT[(i-window):(i-1),2]
	rq_SPY <- Quarticities$rq_SPY[(i-window):(i-1),2]

	rq_TLT <- matrix(rq_TLT)
	rq_SPY <- matrix(rq_SPY)

	trq_SPY <- Quarticities$trq_SPY[(i-window):(i-1),3]
	trq_TLT <- Quarticities$trq_TLT[(i-window):(i-1),3]

	trq_TLT <- trq_TLT[22:length(trq_TLT)]
	trq_SPY <- trq_SPY[22:length(trq_SPY)]


	sqrtrq_TLT <- sqrt(rq_TLT) - mean(sqrt(rq_TLT))
	sqrttrq_TLT <- sqrt(trq_TLT) - mean(sqrt(trq_TLT))

	sqrtrq_SPY <- sqrt(rq_SPY) - mean(sqrt(rq_SPY))
	sqrttrq_SPY<- sqrt(trq_SPY) - mean(sqrt(trq_SPY))


	#lagged for HARQF
	rqSPYweek <- rowMeans(cbind(rq_SPY, mlag(rq_SPY,4,mean(rq_SPY))))
	sqrtrq_SPYweek <- sqrt(rqSPYweek) - mean(sqrt(rqSPYweek))

	rqSPYmonth <- rowMeans(cbind(rq_SPY, mlag(rq_SPY,21,mean(rq_SPY))))
	sqrtrq_SPYmonth <- sqrt(rqSPYmonth) - mean(sqrt(rqSPYmonth))


	rqTLTweek <- rowMeans(cbind(rq_TLT, mlag(rq_TLT,4,mean(rq_TLT))))
	sqrtrq_TLTweek <- sqrt(rqTLTweek) - mean(sqrt(rqTLTweek))

	rqTLTmonth <- rowMeans(cbind(rq_TLT, mlag(rq_TLT,21,mean(rq_TLT))))
	sqrtrq_TLTmonth <- sqrt(rqTLTmonth) - mean(sqrt(rqTLTmonth))


	sqrtrq_TLT <- sqrtrq_TLT[22:length(sqrtrq_TLT)]
	sqrtrq_TLTweek <- sqrtrq_TLTweek[22:length(sqrtrq_TLTweek)]
	sqrtrq_TLTmonth <- sqrtrq_TLTmonth[22:length(sqrtrq_TLTmonth)]

	sqrtrq_SPY <- sqrtrq_SPY[22:length(sqrtrq_SPY)]
	sqrtrq_SPYweek <- sqrtrq_SPYweek[22:length(sqrtrq_SPYweek)]
	sqrtrq_SPYmonth <- sqrtrq_SPYmonth[22:length(sqrtrq_SPYmonth)]


	#jumpparams
	#jumpparam_TLT <- ifelse(calccov[[1]][[2]][1,1,(i-window):(i-1)]*10000 - calccov[[5]][[2]][1,1,(i-window):(i-1)]*10000>0, 
	#	calccov[[1]][[2]][1,1,(i-window):(i-1)]*10000 - calccov[[5]][[2]][1,1,(i-window):(i-1)]*10000, 0)

	#jumpparam_SPY <- ifelse(calccov[[1]][[2]][2,2,(i-window):(i-1)]*10000 - calccov[[5]][[2]][2,2,(i-window):(i-1)]*10000>0, 
	#	calccov[[1]][[2]][2,2,(i-window):(i-1)]*10000 - calccov[[5]][[2]][2,2,(i-window):(i-1)]*10000, 0)

	#under RCOV 5 min
	jumpparam_TLT <- ifelse(calccov[[1]][[7]][1,1,(i-window):(i-1)]*10000 - calccov[[5]][[7]][1,1,(i-window):(i-1)]*10000>0, 
		calccov[[1]][[7]][1,1,(i-window):(i-1)]*10000 - calccov[[5]][[7]][1,1,(i-window):(i-1)]*10000, 0)

	jumpparam_SPY <- ifelse(calccov[[1]][[7]][2,2,(i-window):(i-1)]*10000 - calccov[[5]][[7]][2,2,(i-window):(i-1)]*10000>0, 
		calccov[[1]][[7]][2,2,(i-window):(i-1)]*10000 - calccov[[5]][[7]][2,2,(i-window):(i-1)]*10000, 0)


	jumpparam_TLT <- jumpparam_TLT[22:length(jumpparam_TLT)]
	jumpparam_SPY <- jumpparam_SPY[22:length(jumpparam_SPY)]


	#############################################################################################
	#
	#
	#									DRD-HAR MODEL (MRK 5 sec)
	# 
	#
	#############################################################################################

	HAR_TLT <- lm(RV5min_TLT[2:end] ~ volday_TLT[1:(end-1)] + volweek_TLT[1:(end-1)] + volmonth_TLT[1:(end-1)])
	hhat_HAR_TLT <- cbind(1, volday_TLT[end], volweek_TLT[end], volmonth_TLT[end])  %*% matrix(coef(HAR_TLT))

	HAR_SPY <- lm(RV5min_SPY[2:end] ~ volday_SPY[1:(end-1)] + volweek_SPY[1:(end-1)] + volmonth_SPY[1:(end-1)])
	hhat_HAR_SPY <- cbind(1, volday_SPY[end], volweek_SPY[end], volmonth_SPY[end])  %*% matrix(coef(HAR_SPY))

	hhat_HAR_TLT <- volatility.insanity.filter(hhat_HAR_TLT, min(RV5min_TLT), max(RV5min_TLT), mean(RV5min_TLT))$vol
	hhat_HAR_SPY <- volatility.insanity.filter(hhat_HAR_SPY, min(RV5min_SPY), max(RV5min_SPY), mean(RV5min_SPY))$vol

	#warning is just that R does not like to add one 1 parameter array to a number eg. array(1, c(1,1,1)) + 2.  
	DRD_HAR <- suppressWarnings(EstimatecorrHAR(cbind(hhat_HAR_TLT, hhat_HAR_SPY), 
		correlation = correlation, proxy = proxycor, 0, T))

	sigmas_HAR[,,i] <-  DRD_HAR$vSigma2[,,1]


	#############################################################################################
	#
	#
	#									DRD-HAR MODEL (DAILY)
	# 
	#
	#############################################################################################
	HAR_TLTDAILY <- lm(RV5min_TLT[2:end] ~ volday_TLTDAILY[1:(end-1)] + volweek_TLTDAILY[1:(end-1)] + volmonth_TLTDAILY[1:(end-1)])
	hhat_HAR_TLTDAILY <- cbind(1, volday_TLTDAILY[end], volweek_TLTDAILY[end], volmonth_TLTDAILY[end])  %*% matrix(coef(HAR_TLTDAILY))

	HAR_SPYDAILY <- lm(RV5min_SPY[2:end] ~ volday_SPYDAILY[1:(end-1)] + volweek_SPYDAILY[1:(end-1)] + volmonth_SPYDAILY[1:(end-1)])
	hhat_HAR_SPYDAILY <- cbind(1, volday_SPYDAILY[end], volweek_SPYDAILY[end], volmonth_SPYDAILY[end]) %*% matrix(coef(HAR_SPYDAILY))
		
	hhat_HAR_TLTDAILY <- volatility.insanity.filter(hhat_HAR_TLTDAILY, min(RV5min_TLT), max(RV5min_TLT), mean(RV5min_TLT))$vol
	hhat_HAR_SPYDAILY <- volatility.insanity.filter(hhat_HAR_SPYDAILY, min(RV5min_SPY), max(RV5min_SPY), mean(RV5min_SPY))$vol

	DRD_HARDAILY <- suppressWarnings(EstimatecorrHAR(cbind(hhat_HAR_TLTDAILY, hhat_HAR_SPYDAILY), 
		correlation = correlationDAILY[(i-window):(i-1)], proxy = proxycor, 0, Forecasts = T))
	sigmas_HAR_daily[,,i] <- DRD_HARDAILY$vSigma2[,,1]


	#############################################################################################
	#
	#
	#									DRD-HARQ MODEL (RCOV, MRC, MRK)
	# 
	#
	#############################################################################################

	HARQ_TLT <- lm(RV5min_TLT[2:end] ~ volday_TLT[1:(end-1)] +  volweek_TLT[1:(end-1)] + volmonth_TLT[1:(end-1)] + I(volday_TLT[1:(end-1)] * sqrtrq_TLT[1:(end-1)]))
	hhat_HARQ_TLT <- cbind(1, volday_TLT[end], volweek_TLT[end], volmonth_TLT[end],volday_TLT[end] * sqrtrq_TLT[end])  %*% matrix(coef(HARQ_TLT))


	HARQ_SPY <- lm(RV5min_SPY[2:end] ~ volday_SPY[1:(end-1)] +  volweek_SPY[1:(end-1)] + volmonth_SPY[1:(end-1)] + I(volday_SPY[1:(end-1)] * sqrtrq_SPY[1:(end-1)]))
	hhat_HARQ_SPY <- cbind(1, volday_SPY[end], volweek_SPY[end], volmonth_SPY[end],volday_SPY[end] * sqrtrq_SPY[end])  %*% matrix(coef(HARQ_SPY))

	hhat_HARQ_TLT <- volatility.insanity.filter(hhat_HARQ_TLT, min(RV5min_TLT), max(RV5min_TLT), mean(RV5min_TLT))$vol
	hhat_HARQ_SPY <- volatility.insanity.filter(hhat_HARQ_SPY, min(RV5min_SPY), max(RV5min_SPY), mean(RV5min_SPY))$vol

	#SAME CORRELATION AS HAR MODELS. ONLY THING THAT CHANGED IS THE UNIVARIATE VOLS. 
	DRD_HARQ <- suppressWarnings(EstimatecorrHAR(cbind(hhat_HARQ_TLT, hhat_HARQ_SPY), correlation = correlation, proxy = proxycor, 0,T))

	sigmas_HARQ[,,i] <- DRD_HARQ$vSigma2[,,1]

	
	#############################################################################################
	#
	#
	#									DRD-HARQF MODEL (MRK 5 sec)
	# 
	#
	#############################################################################################
			
	HARQF_TLT <- lm(RV5min_TLT[2:end] ~ volday_TLT[1:(end-1)] +  volweek_TLT[1:(end-1)] + 
	volmonth_TLT[1:(end-1)] + I(volday_TLT[1:(end-1)] * sqrtrq_TLT[1:(end-1)]) + 
	I(volweek_TLT[1:(end-1)] * sqrtrq_TLTweek[1:(end-1)]) + I(volmonth_TLT[1:(end-1)] * sqrtrq_TLTmonth[1:(end-1)]))

	hhat_HARQF_TLT <- cbind(1, volday_TLT[end], volweek_TLT[end], volmonth_TLT[end],
	volday_TLT[end] * sqrtrq_TLT[end], volweek_TLT[end] * sqrtrq_TLTweek[end],
	volmonth_TLT[end] * sqrtrq_TLTmonth[end])  %*% matrix(coef(HARQF_TLT))


	HARQF_SPY <- lm(RV5min_SPY[2:end] ~ volday_SPY[1:(end-1)] +  volweek_SPY[1:(end-1)] + 
	volmonth_SPY[1:(end-1)] + I(volday_SPY[1:(end-1)] * sqrtrq_SPY[1:(end-1)]) + 
	I(volweek_SPY[1:(end-1)] * sqrtrq_SPYweek[1:(end-1)]) + I(volmonth_SPY[1:(end-1)] * sqrtrq_SPYmonth[1:(end-1)]))

	hhat_HARQF_SPY <- cbind(1, volday_SPY[end], volweek_SPY[end], volmonth_SPY[end],
	volday_SPY[end] * sqrtrq_SPY[end], volweek_SPY[end] * sqrtrq_SPYweek[end], 
	volmonth_SPY[end] * sqrtrq_SPYmonth[end])  %*% matrix(coef(HARQF_SPY))

	hhat_HARQF_SPY <- volatility.insanity.filter(hhat_HARQF_SPY, min(RV5min_SPY), max(RV5min_SPY), mean(RV5min_SPY))$vol
	hhat_HARQF_TLT <- volatility.insanity.filter(hhat_HARQF_TLT, min(RV5min_TLT), max(RV5min_TLT), mean(RV5min_TLT))$vol


	DRD_HARQF <- suppressWarnings(EstimatecorrHAR(cbind(hhat_HARQF_TLT, hhat_HARQF_SPY), correlation = correlation, 
		proxy = proxycor, 0, Forecasts = T))

	sigmas_HARQF[,,i] <- DRD_HARQF$vSigma2[,,1]


	#############################################################################################
	#
	#
	#									DRD-HARJ MODEL (RCOV, MRC, MRK)
	# 
	#
	#############################################################################################


	HARJ_TLT <- lm(RV5min_TLT[2:end] ~ volday_TLT[1:(end-1)] + volweek_TLT[1:(end-1)] + volmonth_TLT[1:(end-1)] + jumpparam_TLT[1:(end-1)])
	hhat_HARJ_TLT <- cbind(1, volday_TLT[end], volweek_TLT[end], volmonth_TLT[end],jumpparam_TLT[end])  %*% matrix(coef(HARJ_TLT))

	HARJ_SPY <- lm(RV5min_SPY[2:end] ~ volday_SPY[1:(end-1)] + volweek_SPY[1:(end-1)] + volmonth_SPY[1:(end-1)] + jumpparam_SPY[1:(end-1)])
	hhat_HARJ_SPY <- cbind(1, volday_SPY[end], volweek_SPY[end], volmonth_SPY[end],jumpparam_SPY[end])  %*% matrix(coef(HARJ_SPY))

	hhat_HARJ_SPY <- volatility.insanity.filter(hhat_HARJ_SPY, min(RV5min_SPY), max(RV5min_SPY), mean(RV5min_SPY))$vol
	hhat_HARJ_TLT <- volatility.insanity.filter(hhat_HARJ_TLT, min(RV5min_TLT), max(RV5min_TLT), mean(RV5min_TLT))$vol

	DRD_HARJ <- suppressWarnings(EstimatecorrHAR(cbind(hhat_HARJ_TLT, hhat_HARJ_SPY), correlation = correlation, proxy = proxycor, 0, T))
			
	sigmas_HARJ[,,i] <- DRD_HARJ$vSigma2[,,1]

	

	#############################################################################################
	#
	#
	#									DRD-HARQJ MODEL (RCOV, MRC, MRK)
	# 
	#
	#############################################################################################

	HARQJ_TLT <- lm(RV5min_TLT[2:end] ~ volday_TLT[1:(end-1)] +  volweek_TLT[1:(end-1)] + volmonth_TLT[1:(end-1)] + I(volday_TLT[1:(end-1)] * sqrtrq_TLT[1:(end-1)]) + jumpparam_TLT[1:(end-1)])
	hhat_HARQJ_TLT <- cbind(1, volday_TLT[end], volweek_TLT[end], volmonth_TLT[end],volday_TLT[end] * sqrtrq_TLT[end], jumpparam_TLT[end])  %*% matrix(coef(HARQJ_TLT))

	HARQJ_SPY <- lm(RV5min_SPY[2:end] ~ volday_SPY[1:(end-1)] +  volweek_SPY[1:(end-1)] + volmonth_SPY[1:(end-1)] + I(volday_SPY[1:(end-1)] * sqrtrq_SPY[1:(end-1)]) + jumpparam_SPY[1:(end-1)])
	hhat_HARQJ_SPY <- cbind(1, volday_SPY[end], volweek_SPY[end], volmonth_SPY[end],volday_SPY[end] * sqrtrq_SPY[end], jumpparam_SPY[end])  %*% matrix(coef(HARQJ_SPY))

	hhat_HARQJ_SPY <- volatility.insanity.filter(hhat_HARQJ_SPY, min(RV5min_SPY), max(RV5min_SPY), mean(RV5min_SPY))$vol
	hhat_HARQJ_TLT <- volatility.insanity.filter(hhat_HARQJ_TLT, min(RV5min_TLT), max(RV5min_TLT), mean(RV5min_TLT))$vol

	DRD_HARQJ <- suppressWarnings(EstimatecorrHAR(cbind(hhat_HARQJ_TLT, hhat_HARQJ_SPY), correlation = correlation, proxy = proxycor, 0, T))

	sigmas_HARQJ[,,i] <- DRD_HARQJ$vSigma2[,,1]


	#############################################################################################
	#
	#
	#									DRD-CHAR MODEL PBPCOV 15 SEC
	# 
	#
	#############################################################################################

	CHAR_TLT <- lm(RV5min_TLT[2:end] ~ ICvolday_TLT[1:(end-1)] + ICvolweek_TLT[1:(end-1)] + ICvolmonth_TLT[1:(end-1)])
	hhat_CHAR_TLT <- cbind(1, volday_TLT[end], ICvolweek_TLT[end], ICvolmonth_TLT[end])  %*% matrix(coef(CHAR_TLT))

	CHAR_SPY <- lm(RV5min_SPY[2:end] ~ ICvolday_SPY[1:(end-1)] + ICvolweek_SPY[1:(end-1)] + ICvolmonth_SPY[1:(end-1)])
	hhat_CHAR_SPY <- cbind(1, ICvolday_SPY[end], ICvolweek_SPY[end], ICvolmonth_SPY[end])  %*% matrix(coef(CHAR_SPY))

	hhat_CHAR_SPY <- volatility.insanity.filter(hhat_CHAR_SPY, min(RV5min_SPY), max(RV5min_SPY), mean(RV5min_SPY))$vol
	hhat_CHAR_TLT <- volatility.insanity.filter(hhat_CHAR_TLT, min(RV5min_TLT), max(RV5min_TLT), mean(RV5min_TLT))$vol

	DRD_CHAR <- suppressWarnings(EstimatecorrHAR(cbind(hhat_CHAR_TLT, hhat_CHAR_SPY), correlation = correlationjumprobust, proxy = proxycor, 0, T))

	sigmas_CHAR[,,i] <- DRD_CHAR$vSigma2[,,1]



	#############################################################################################
	#
	#
	#									DRD-CHARQ MODEL PBPCOV 15 SEC
	# 
	#
	#############################################################################################

	CHARQ_TLT <- lm(RV5min_TLT[2:end] ~ ICvolday_TLT[1:(end-1)] +  ICvolweek_TLT[1:(end-1)] + ICvolmonth_TLT[1:(end-1)] + I(ICvolday_TLT[1:(end-1)] * sqrttrq_TLT[1:(end-1)]))
	hhat_CHARQ_TLT <- cbind(1, ICvolday_TLT[end], ICvolweek_TLT[end], ICvolmonth_TLT[end], ICvolday_TLT[end] * sqrttrq_TLT[end])  %*% matrix(coef(CHARQ_TLT))

	CHARQ_SPY <- lm(RV5min_SPY[2:end] ~ ICvolday_SPY[1:(end-1)] +  ICvolweek_SPY[1:(end-1)] + ICvolmonth_SPY[1:(end-1)] + I(ICvolday_SPY[1:(end-1)] * sqrttrq_SPY[1:(end-1)]))
	hhat_CHARQ_SPY <- cbind(1, ICvolday_SPY[end], ICvolweek_SPY[end], ICvolmonth_SPY[end], ICvolday_SPY[end] * sqrttrq_SPY[end])  %*% matrix(coef(CHARQ_SPY))

	hhat_CHARQ_SPY <- volatility.insanity.filter(hhat_CHARQ_SPY, min(RV5min_SPY), max(RV5min_SPY), mean(RV5min_SPY))$vol
	hhat_CHARQ_TLT <- volatility.insanity.filter(hhat_CHARQ_TLT, min(RV5min_TLT), max(RV5min_TLT), mean(RV5min_TLT))$vol


	DRD_CHARQ <- suppressWarnings(EstimatecorrHAR(cbind(hhat_CHARQ_TLT, hhat_CHARQ_SPY), correlation = correlationjumprobust, 
		proxy = proxycor, 0, Forecasts = T))

	sigmas_CHARQ[,,i] <- DRD_CHARQ$vSigma2[,,1]

	
	#############################################################################################
	#
	#
	#									RM (RCov daily)
	# 
	#
	#############################################################################################
	
	Filter <- ewma.filter.realized(calccov[[1]][[10]][,,(i-window):(i-1)] * 10000, NULL, F, 0.94, 0)
	
	end2 <- length(Filter[1,1, ])
	sigmas_RM_daily[,,i] <- Filter[,,end2]
	
	print(sprintf("%s", i))

}

sigmas_HARQ <- sigmas_HARQ[,,1001:2516]
sigmas_HARQJ <- sigmas_HARQJ[,,1001:2516]
sigmas_HARQF <- sigmas_HARQF[,,1001:2516]
sigmas_HAR <- sigmas_HAR[,,1001:2516]
sigmas_HARJ <- sigmas_HARJ[,,1001:2516]
sigmas_HAR_daily <- sigmas_HAR_daily[,,1001:2516]
sigmas_CHAR <- sigmas_CHAR[,,1001:2516]
sigmas_CHARQ <- sigmas_CHARQ[,,1001:2516]
sigmas_RM_daily <- sigmas_RM_daily[,,1001:2516]

#check if they are correct via qlike losses:
qlikes <- numeric()
for(i in 1:1516){
  
  qlikes[i] <- QLIKE(sigmas_HARQF[,,i], calccov[[1]][[7]][,,i+1000]*10000, 2)
  
}
mean(qlikes)


#insample because out-of-sample resulted in bounded values majority of the time. 



sigmas_rBG_daily <- readRDS("sigmas_rBG_daily.rds")


#####################################################################################################
#
#
#                     					The Analysis 
#
#
#####################################################################################################


################################### WEIGHT CALCULATIONS #################################

#calculating weights of risky assets:
rp_weigths_HAR <- matrix(0L, nrow=1516, ncol = 3)
rp_weigths_HARJ <- matrix(0L, nrow=1516, ncol = 3)
rp_weigths_HARQJ <- matrix(0L, nrow=1516, ncol = 3)
rp_weigths_HAR_daily <- matrix(0L, nrow=1516, ncol = 3)
rp_weigths_CHAR <- matrix(0L, nrow=1516, ncol = 3)
rp_weigths_CHARQ <- matrix(0L, nrow=1516, ncol = 3)
rp_weigths_HARQ <- matrix(0L, nrow=1516, ncol = 3)
rp_weigths_HARQF <- matrix(0L, nrow=1516, ncol = 3)
rp_weigths_RM_daily <- matrix(0L, nrow=1516, ncol = 3)
rp_weigths_rBG <- matrix(0L, nrow=1516, ncol = 3)
rp_weigths_rBG_daily <- matrix(0L, nrow=1516, ncol = 3)
rp_weigths_crBG <- matrix(0L, nrow=1516, ncol = 3)
rp_weigths_rDCC <- matrix(0L, nrow=1516, ncol = 3)
rp_weigths_crDCC <- matrix(0L, nrow=1516, ncol = 3)



#alphas to divide turnover so its fixed between 0 and 1:
alpha_HAR <- numeric()
alpha_HARJ <- numeric()
alpha_HARQJ <- numeric()
alpha_HAR_daily <- numeric()
alpha_CHAR <- numeric()
alpha_CHARQ <- numeric()
alpha_HARQ <- numeric()
alpha_HARQF <- numeric()
alpha_RM_daily <- numeric()
alpha_rBG <- numeric()
alpha_rBG_daily <- numeric()
alpha_crBG <- numeric()
alpha_rDCC <- numeric()
alpha_crDCC <- numeric()

#risktarget should be annualized. therefore divide by sqrt(252). You could also scale it with 24/6.5. 

for(i in 1:1516){
  rp_weigths_HAR[i, ] <- t(riskparity_2dim(sigmas_HAR[,,i]*(1/10000), 0.1/(sqrt(252)), T)$w)
  rp_weigths_HARJ[i, ] <- t(riskparity_2dim(sigmas_HARJ[,,i]* (1/10000), 0.1/(sqrt(252)), T)$w)
  rp_weigths_HARQJ[i, ] <- t(riskparity_2dim(sigmas_HARQJ[,,i] * (1/10000), 0.1/(sqrt(252)), T)$w)
  rp_weigths_HAR_daily[i, ] <- t(riskparity_2dim(sigmas_HAR_daily[,,i]* (1/10000), 0.1/(sqrt(252)), T)$w)
  rp_weigths_CHAR[i, ] <- t(riskparity_2dim(sigmas_CHAR[,,i]* (1/10000), 0.1/(sqrt(252)), T)$w)
  rp_weigths_CHARQ[i, ] <- t(riskparity_2dim(sigmas_CHARQ[,,i]* (1/10000), 0.1/(sqrt(252)), T)$w)
  rp_weigths_HARQ[i, ] <- t(riskparity_2dim(sigmas_HARQ[,,i]* (1/10000), 0.1/(sqrt(252)), T)$w)
  rp_weigths_HARQF[i, ] <- t(riskparity_2dim(sigmas_HARQF[,,i]* (1/10000), 0.1/(sqrt(252)), T)$w)
  rp_weigths_RM_daily[i, ] <- t(riskparity_2dim(sigmas_RM_daily[,,i]* (1/10000), 0.1/sqrt(252), T)$w)
  rp_weigths_rBG[i, ]   <- t(riskparity_2dim(sigmas_rBG[,,i+1000]* (1/10000), 0.1/(sqrt(252)), T)$w)
  rp_weigths_crBG[i, ]  <- t(riskparity_2dim(sigmas_crBG[,,i+1000]* (1/10000), 0.1/(sqrt(252)), T)$w)
  rp_weigths_rDCC[i, ]  <- t(riskparity_2dim(sigmas_rDCC[,,i+1000]* (1/10000), 0.1/(sqrt(252)), T)$w)
  rp_weigths_crDCC[i, ]  <- t(riskparity_2dim(sigmas_crDCC[,,i+1000]* (1/10000), 0.1/(sqrt(252)), T)$w)
  rp_weigths_rBG_daily[i, ] <- t(riskparity_2dim(sigmas_rBG_daily[,,i+1000]* (1/10000), 0.1/(sqrt(252)), T)$w)


  alpha_HAR[i] <- t(riskparity_2dim(sigmas_HAR[,,i] * (1/10000), 0.1/(sqrt(252)), T)$alpha)
  alpha_HARJ[i] <- t(riskparity_2dim(sigmas_HARJ[,,i] * (1/10000), 0.1/(sqrt(252)), T)$alpha)
  alpha_HARQJ[i] <- t(riskparity_2dim(sigmas_HARQJ[,,i] * (1/10000), 0.1/(sqrt(252)), T)$alpha)
  alpha_HAR_daily[i] <- t(riskparity_2dim(sigmas_HAR_daily[,,i]* (1/10000), 0.1/(sqrt(252)), T)$alpha)
  alpha_CHAR[i] <- t(riskparity_2dim(sigmas_CHAR[,,i]* (1/10000), 0.1/(sqrt(252)), T)$alpha)
  alpha_CHARQ[i] <- t(riskparity_2dim(sigmas_CHARQ[,,i]* (1/10000), 0.1/(sqrt(252)), T)$alpha)
  alpha_HARQ[i] <- t(riskparity_2dim(sigmas_HARQ[,,i]* (1/10000), 0.1/(sqrt(252)), T)$alpha)
  alpha_HARQF[i] <- t(riskparity_2dim(sigmas_HARQF[,,i]* (1/10000), 0.1/(sqrt(252)), T)$alpha)
  alpha_RM_daily[i] <- t(riskparity_2dim(sigmas_RM_daily[,,i]* (1/10000), 0.1/(sqrt(252)), T)$alpha)
  alpha_rBG[i]   <- t(riskparity_2dim(sigmas_rBG[,,i+1000]* (1/10000), 0.1/(sqrt(252)), T)$alpha)
  alpha_crBG[i]  <- t(riskparity_2dim(sigmas_crBG[,,i+1000]* (1/10000), 0.1/(sqrt(252)), T)$alpha)
  alpha_rDCC[i]  <- t(riskparity_2dim(sigmas_rDCC[,,i+1000]* (1/10000), 0.1/(sqrt(252)), T)$alpha)
  alpha_crDCC[i]  <- t(riskparity_2dim(sigmas_crDCC[,,i+1000]* (1/10000), 0.1/(sqrt(252)), T)$alpha)
  alpha_rBG_daily[i] <- t(riskparity_2dim(sigmas_rBG_daily[,,i+1000]* (1/10000), 0.1/sqrt(252), T)$alpha)
}

############################## PORTFOLIO RETURNS AND REALIZED VOLATILITY ###########################


#calculaing portfolio excess returns each day:
portret_HAR <- numeric()
portvar_HAR <- numeric()

portret_HAR2 <- list()
portret_RM2 <- list()

portret_HARJ <- numeric()
portvar_HARJ <- numeric()

portret_HARQJ <- numeric()
portvar_HARQJ <- numeric()

portret_HARQ <- numeric()
portvar_HARQ <- numeric()

portret_HAR_daily <- numeric()
portvar_HAR_daily <- numeric()

portret_HARQF <- numeric()
portvar_HARQF <- numeric()

portret_CHAR <- numeric()
portvar_CHAR <- numeric()

portret_CHARQ <- numeric()
portvar_CHARQ <- numeric()

portret_rBG <- numeric()
portvar_rBG <- numeric()

portret_rBG_daily <- numeric()
portvar_rBG_daily <- numeric()

portret_crBG <- numeric()
portvar_crBG <- numeric()

portret_rDCC <- numeric()
portvar_rDCC <- numeric()

portret_crDCC <- numeric()
portvar_crDCC <- numeric()

portret_RM_daily <- numeric()
portvar_RM_daily <- numeric()


#I DONT KNOW HOW TO CALCULATE PORTRET FOR A COMPARABLE FRAMEWORK?? SHOULD THEY ALL BE FOR DAILY RETURNS?
for(i in 1:1516){
  
  #portret_HAR2[[i]] <-  (mergedret_5min[[i+1000]][,1:3]) %*% rp_weigths_HAR[i, 1:3] 
  #portret_RM2[[i]] <-  (mergedret_5min[[i+1000]][,1:3]) %*%  rp_weigths_RM_daily[i, ]


  portret_HAR[i] <- (mergedret_daily[i+1000, ]) %*% rp_weigths_HAR[i, ]
  portret_HARJ[i] <-mergedret_daily[i+1000, ] %*% rp_weigths_HARJ[i, ]
  portret_HARQJ[i] <-mergedret_daily[i+1000, ] %*% rp_weigths_HARQJ[i, ]
  portret_HAR_daily[i] <- mergedret_daily[i+1000, ] %*% rp_weigths_HAR_daily[i, ]
  portret_HARQF[i] <- mergedret_daily[i+1000, ] %*% rp_weigths_HARQF[i, ]  
  portret_CHAR[i] <- mergedret_daily[i+1000, ] %*% rp_weigths_CHAR[i, ]
  portret_CHARQ[i] <- mergedret_daily[i+1000, ] %*% rp_weigths_CHARQ[i, ]
  portret_HARQ[i] <- mergedret_daily[i+1000, ] %*% rp_weigths_HARQ[i, ]  
  portret_rBG[i] <- mergedret_daily[i+1000, ] %*% rp_weigths_rBG[i, ]
  portret_rBG_daily[i] <- mergedret_daily[i+1000, ] %*% rp_weigths_rBG_daily[i, ]
  portret_crBG[i] <- mergedret_daily[i+1000, ] %*% rp_weigths_crBG[i, ]
  portret_rDCC[i] <- mergedret_daily[i+1000, ] %*% rp_weigths_rDCC[i, ]
  portret_crDCC[i] <- mergedret_daily[i+1000, ] %*% rp_weigths_crDCC[i, ]
  portret_RM_daily[i] <- mergedret_daily[i+1000, ] %*% rp_weigths_RM_daily[i, ]
	

  #these are the realized portfolio variance. We see how the forecasted weights act under our proxy
  portvar_HAR[i] <-  sqrt(t(rp_weigths_HAR[i, 1:2]) %*% (calccov[[1]][[7]][,,(i+1000)]) %*% (rp_weigths_HAR[i, 1:2]))
  portvar_HARJ[i] <-  sqrt(t(rp_weigths_HARJ[i, 1:2]) %*% (calccov[[1]][[7]][,,(i+1000)]) %*% (rp_weigths_HARJ[i, 1:2]))
  portvar_HARQJ[i] <-  sqrt(t(rp_weigths_HARQJ[i, 1:2]) %*% (calccov[[1]][[7]][,,(i+1000)]) %*% (rp_weigths_HARQJ[i, 1:2]))
  portvar_HAR_daily[i] <-  sqrt(t(rp_weigths_HAR_daily[i, 1:2]) %*% (calccov[[1]][[7]][,,(i+1000)]) %*% (rp_weigths_HAR_daily[i, 1:2]))
  portvar_HARQF[i] <-  sqrt(t(rp_weigths_HARQF[i, 1:2]) %*% (calccov[[1]][[7]][,,(i+1000)]) %*% (rp_weigths_HARQF[i, 1:2]))
  portvar_HARQ[i] <-  sqrt(t(rp_weigths_HARQ[i, 1:2]) %*% (calccov[[1]][[7]][,,(i+1000)]) %*% (rp_weigths_HARQ[i, 1:2]))
  portvar_CHAR[i] <-  sqrt(t(rp_weigths_CHAR[i, 1:2]) %*% (calccov[[1]][[7]][,,(i+1000)]) %*% (rp_weigths_CHAR[i, 1:2]))
  portvar_CHARQ[i] <-  sqrt(t(rp_weigths_CHARQ[i, 1:2]) %*% (calccov[[1]][[7]][,,(i+1000)]) %*% (rp_weigths_CHARQ[i, 1:2]))
  portvar_rBG[i] <- sqrt(t(rp_weigths_rBG[i, 1:2]) %*% (calccov[[1]][[7]][,,(i+1000)]) %*% (rp_weigths_rBG[i, 1:2]))
  portvar_rBG_daily[i] <- sqrt(t(rp_weigths_rBG_daily[i, 1:2]) %*% (calccov[[1]][[7]][,,(i+1000)]) %*% (rp_weigths_rBG_daily[i, 1:2]))
  portvar_crBG[i] <- sqrt(t(rp_weigths_crBG[i, 1:2]) %*% (calccov[[1]][[7]][,,(i+1000)]) %*% (rp_weigths_crBG[i, 1:2]))
  portvar_rDCC[i] <- sqrt(t(rp_weigths_rDCC[i, 1:2]) %*% (calccov[[1]][[7]][,,(i+1000)]) %*% (rp_weigths_rDCC[i, 1:2]))
  portvar_crDCC[i] <- sqrt(t(rp_weigths_crDCC[i, 1:2]) %*% (calccov[[1]][[7]][,,(i+1000)]) %*% (rp_weigths_crDCC[i, 1:2]))

  portvar_RM_daily[i] <- sqrt(t(rp_weigths_RM_daily[i, 1:2]) %*% (calccov[[1]][[10]][,,(i+1000)]) %*% (rp_weigths_RM_daily[i, 1:2]))




}


library(ggplot2)

ggplot() + geom_line(aes(1:1516, sigmas_HARQF[2,2,], col = "HARQF")) + 
geom_line(aes(1:1516, sigmas_RM_daily[2,2,], col = "RM"))  
geom_line(aes(1:1516, calccov[[1]][[7]][1,1,1001:2516]*10000, col = "RCov"))

test <- do.call(rbind, portret_HAR2)
test2 <- do.call(rbind, portret_RM2)


#turntest <- turnover(rp_weigths_HAR, portret_HAR2[[i]][length(portret_HAR2[[i]])], mergedret_5min[[i]][length(mergedret_5min[[i]][,1]), ])

#sum(abs(rp_weigths_HAR[2, ]/alpha_HAR[2] - rp_weigths_HAR[1, ]/alpha_HAR[1] *  (1 + mergedret_5min[[i]][length(mergedret_5min[[i]][,1]),])/(1+ portret_HAR2[[1]][length(portret_HAR2[[1]])])))

#sum(abs(rp_weigths_RM_daily[2, ]/alpha_RM_daily[2] - rp_weigths_RM_daily[1, ]/alpha_RM_daily[1] *  (1 + mergedret_5min[[i]][length(mergedret_5min[[i]][,1]),])/(1+ portret_RM2[[1]][length(portret_RM2[[1]])])))



portret <- list(portret_HAR,
  portret_HAR_daily,
  portret_HARJ,
  portret_HARQJ,
  portret_HARQ,
  portret_HARQF,
  portret_CHAR,
   portret_CHARQ, 
  portret_rBG,
  portret_rBG_daily,
  portret_crBG,
  portret_rDCC,
  portret_crDCC,
  portret_RM_daily)



rp_weights <- list(rp_weigths_HAR = rp_weigths_HAR,
 	rp_weigths_HAR_daily = rp_weigths_HAR_daily,
 	rp_weigths_HARJ,
  	rp_weigths_HARQJ, 
	rp_weigths_HARQ = rp_weigths_HARQ, 
	rp_weigths_HARQF = rp_weigths_HARQF,
	rp_weigths_CHAR = rp_weigths_CHAR, 
	rp_weigths_CHARQ = rp_weigths_CHARQ, 
	rp_weights_rBG = rp_weigths_rBG, 
	rp_weigths_rBG_daily = rp_weigths_rBG_daily, 
	rp_weights_crBG = rp_weigths_crBG,
	rp_weights_rDCC = rp_weigths_rDCC, 
	rp_weights_crDCC = rp_weigths_crDCC, 
	rp_weights_RM_daily = rp_weigths_RM_daily)

head(rp_weights[[14]]) - head( rp_weigths_RM_daily)



alphas <- list(alpha_HAR = alpha_HAR, 
	alpha_HAR_daily = alpha_HAR_daily, 
	alpha_HARJ,
  	alpha_HARQJ,
	alpha_HARQ = alpha_HARQ, 
	alpha_HARQF = alpha_HARQF,
	alpha_CHAR = alpha_CHAR, 
	alpha_CHARQ = alpha_CHARQ, 
	rp_weights_rBG = alpha_rBG, 
	alpha_rBG_daily = alpha_rBG_daily, 
	alpha_crBG = alpha_crBG,
	alpha_rDCC = alpha_rDCC, 
	alpha_crDCC = alpha_crDCC, 
	alpha_RM_daily = alpha_RM_daily)

#######NOTE: MAYBE YOU SHOULD NOT DO PERCENTAGE LOG-RETURNS IN THE PORTFOLIO RETURNS AND JUST MULTIPLY YOUR 
#PERFORMANCE NUMBERS AFTERWARDS. 
ts.plot(cumsum(portret_RM_daily), col ="red") + lines(cumsum(portret_HAR))

#now it should be average daily return. 
mean((portret_rBG)) *100 * 252
mean((portret_rBG_daily)) * 252 *100  
mean((portret_HAR))   *100 * 252 
mean((portret_HAR_daily))* 252  *100
mean((portret_HARJ))* 252  *100
mean((portret_HARQJ))* 252  *100
mean((portret_HARQ))* 252  *100
mean((portret_HARQF))* 252  *100
mean((portret_CHAR))* 252  *100
mean((portret_CHARQ))* 252  *100
mean((portret_crBG)) * 252  *100
mean((portret_rDCC)) * 252  *100
mean((portret_crDCC)) * 252  *100
mean((portret_RM_daily)) * 252  *100





#expected standard deviation (dont use it in table it does not matter)
mean(portvar_rBG_daily) * sqrt(252)  *100
mean(portvar_rBG)* sqrt(252)  *100
mean(portvar_HAR)* sqrt(252)  *100
mean(portvar_HAR_daily)* sqrt(252)  *100
mean(portvar_HARQF)* sqrt(252)  *100
mean(portvar_CHAR)* sqrt(252)* 100
mean(portvar_CHARQ)* sqrt(252)  *100
mean(portvar_crBG)* sqrt(252)  *100
mean(portvar_rDCC)* sqrt(252)  *100
mean(portvar_crDCC)* sqrt(252)  *100
mean(portvar_RM_daily)* sqrt(252)  *100




################vol of vol############
sd(portvar_rBG_daily) * sqrt(252)  * 100 
sd(portvar_rBG) * sqrt(252)  * 100
sd(portvar_HAR) * sqrt(252)  * 100
sd(portvar_HARJ) * sqrt(252)  * 100
sd(portvar_HARQJ) * sqrt(252)  * 100
sd(portvar_HAR_daily) * sqrt(252)  * 100
sd(portvar_HARQ) * sqrt(252)  * 100
sd(portvar_HARQF) * sqrt(252)  * 100
sd(portvar_CHAR) * sqrt(252)  * 100
sd(portvar_CHARQ) * sqrt(252)  * 100
sd(portvar_crBG) * sqrt(252)  * 100
sd(portvar_rDCC) * sqrt(252)  * 100
sd(portvar_crDCC) * sqrt(252)  * 100
sd(portvar_RM_daily) * sqrt(252)  * 100

################ ES #####################
ES <- numeric()
for(i in 1:14){

	ES[i] <- ES(portret[[i]], p=.95, method="historical")

}

#under the assumption of normality you can scale with sqrt(252)
ES * 100


########## avg alpha ############
avgalpha <- numeric()

for(i in 1:14){


	avgalpha[i] <- mean(alphas[[i]])

}

avgalpha

#############abs distance#############
mean(abs(portvar_rBG_daily*sqrt(252) - 0.1)) *100
mean(abs(portvar_rBG*sqrt(252) - 0.1))*100
mean(abs(portvar_HAR*sqrt(252) - 0.1))*100
mean(abs(portvar_HARJ*sqrt(252) - 0.1))*100
mean(abs(portvar_HARQJ*sqrt(252) - 0.1))*100
mean(abs(portvar_HAR_daily*sqrt(252) - 0.1))*100
mean(abs(portvar_HARQ*sqrt(252) - 0.1))*100
mean(abs(portvar_HARQF*sqrt(252) - 0.1))*100
mean(abs(portvar_CHAR*sqrt(252) - 0.1))*100
mean(abs(portvar_CHARQ*sqrt(252) - 0.1))*100
mean(abs(portvar_crBG*sqrt(252) - 0.1))*100
mean(abs(portvar_rDCC*sqrt(252) - 0.1))*100
mean(abs(portvar_crDCC*sqrt(252) - 0.1))*100
mean(abs(portvar_RM_daily*sqrt(252) - 0.1))*100




########### TURN OVER ###############


#ncol is the amount of models
Turnovers <- matrix(0L, ncol = 14, nrow=1516)
Turnovers.simple <- matrix(0L, ncol = 14, nrow=1516)


oosmergedret_daily <- matrix(mergedret_daily[1001:2516, ], nrow = 1516, ncol = 3)


Turnovers[, 1] <- turnover(rp_weigths_HAR[,1:2], matrix(portret_HAR), oosmergedret_daily[,1:2])
Turnovers[, 2] <- turnover(rp_weigths_HAR_daily[,1:2], matrix(portret_HAR_daily), oosmergedret_daily[,1:2])
Turnovers[, 3] <- turnover(rp_weigths_HARJ[,1:2], matrix(portret_HARJ), oosmergedret_daily[,1:2])
Turnovers[, 4] <- turnover(rp_weigths_HARQJ[,1:2], matrix(portret_HARQJ), oosmergedret_daily[,1:2])
Turnovers[, 5] <- turnover(rp_weigths_HARQ[,1:2], matrix(portret_HARQ), oosmergedret_daily[,1:2]) 
Turnovers[, 6] <- turnover(rp_weigths_HARQF[,1:2], matrix(portret_HARQF), oosmergedret_daily[,1:2]) 
Turnovers[, 7] <- turnover(rp_weigths_CHAR[,1:2], matrix(portret_CHAR), oosmergedret_daily[,1:2]) 
Turnovers[, 8] <- turnover(rp_weigths_CHARQ[,1:2], matrix(portret_CHARQ), oosmergedret_daily[,1:2]) 
Turnovers[, 9] <- turnover(rp_weigths_rBG[,1:2], matrix(portret_rBG), oosmergedret_daily[,1:2]) 
Turnovers[, 10] <- turnover(rp_weigths_rBG_daily[,1:2], matrix(portret_rBG_daily), oosmergedret_daily[,1:2]) 
Turnovers[, 11] <- turnover(rp_weigths_crBG[,1:2], matrix(portret_crBG), oosmergedret_daily[,1:2]) 
Turnovers[, 12] <- turnover(rp_weigths_rDCC[,1:2], matrix(portret_rDCC), oosmergedret_daily[,1:2]) 
Turnovers[, 13] <- turnover(rp_weigths_crDCC[,1:2], matrix(portret_crDCC), oosmergedret_daily[,1:2])
Turnovers[, 14] <- turnover(rp_weigths_RM_daily[,1:2], matrix(portret_RM_daily), oosmergedret_daily[,1:2]) 

colMeans(Turnovers) 


tt <- Turnovers[,1] * (252/(2 * 2 * rowSums(rp_weigths_HAR[,1:2])))

notionalcost <- rowSums(rp_weigths_HAR[,1:2]) * Turnovers[,1] * 0.00005

turnover(rp_weigths_HAR[,1:2]/alpha_HAR, matrix(portret_HAR), oosmergedret_daily[,1:2])


pturnoverDN  = function(weight,rets,port_ret){ 
  weight[is.na(weight)] <- 0  # NAs = 0
  weighteop = weight*(1+rets)/(1+port_ret)
  dweight = abs(weight-lag(weighteop,1))
  out = (rowSums( dweight) )
  return(out)
}





Turnovers.simple[, 1] <- turnover.simple(rp_weigths_HAR[,1:2])
Turnovers.simple[, 2] <- turnover.simple(rp_weigths_HAR_daily[,1:2])
Turnovers.simple[, 3] <- turnover.simple(rp_weigths_HARJ[,1:2])
Turnovers.simple[, 4] <- turnover.simple(rp_weigths_HARQJ[,1:2])
Turnovers.simple[, 5] <- turnover.simple(rp_weigths_HARQ[,1:2]) 
Turnovers.simple[, 6] <- turnover.simple(rp_weigths_HARQF[,1:2]) 
Turnovers.simple[, 7] <- turnover.simple(rp_weigths_CHAR[,1:2]) 
Turnovers.simple[, 8] <- turnover.simple(rp_weigths_CHARQ[,1:2]) 
Turnovers.simple[, 9] <- turnover.simple(rp_weigths_rBG[,1:2]) 
Turnovers.simple[, 10] <- turnover.simple(rp_weigths_rBG_daily[,1:2]) 
Turnovers.simple[, 11] <- turnover.simple(rp_weigths_crBG[,1:2]) 
Turnovers.simple[, 12] <- turnover.simple(rp_weigths_rDCC[,1:2]) 
Turnovers.simple[, 13] <- turnover.simple(rp_weigths_crDCC[,1:2])
Turnovers.simple[, 14] <- turnover.simple(rp_weigths_RM_daily[,1:2]) 


colMeans(Turnovers.simple)


mean(rowSums(rp_weigths_RM_daily[,1:2]))
mean(rowSums(rp_weigths_HARQF[,1:2]))

0.0001 * 2.21 * ((0.25554792*252)/(2*2.21))


#cost
c <- c(0, 0.01, 0.02)

portret_transcost <- list()
portret_transcost.simple <- list()
temp <- matrix(0L, ncol = 3, nrow = 1516)
temp2 <- matrix(0L, ncol = 3, nrow = 1516)
sharpes <- matrix(0L, ncol = 14, nrow = 3)
sharpes.simple <- matrix(0L, ncol = 14, nrow = 3)
for(j in 1:14){
	for(i in 1:3){

		temp[, i] <- portret[[j]]  - c[i]/252 * Turnovers[, j]
		temp2[, i] <- portret[[j]] - c[i]/252 * Turnovers.simple[, j]
	}

	portret_transcost[[j]] <- temp 

	portret_transcost.simple[[j]] <- temp2
									#*252*(24/6.5)															#*sqrt(252)*sqrt(24/6.5)
	sharpes[, j] <- t(( ((colMeans(portret_transcost[[j]]))*252) / (apply(portret_transcost[[j]],MARGIN = c(2), FUN = function(x) sd(x)*sqrt(252))) ))
	sharpes.simple[, j] <- t((colMeans(portret_transcost.simple[[j]]*252*(24/6.5)) / apply(portret_transcost.simple[[j]],MARGIN = c(2), FUN = function(x) sd(x)*sqrt(252)*sqrt(24/6.5))))

	colnames(sharpes) <- c("HAR", "HAR_daily", "HARJ", "HARQJ", "HARQ", "HARQF", "CHAR", "CHARQ", "rBG", "rBG_daily", "crBG", "rDCC", "crDCC", "RM_daily")
	colnames(sharpes.simple) <- colnames(sharpes)
}

#sharpes are in annualized (percentage) terms. Scaling expected ret and vol with 100, 
#eliminates eachother in denominator and numerator  
sharpes 

names(sort(sharpes[1, ], decreasing =  T))
library(PerformanceAnalytics)


sharpes.simple
(mean(portret_transcost[[1]][,1])*252*24/6.5) / (sd(portret_transcost[[1]][,1])*sqrt(252)*24/6.5)


tt <- portret[[14]]*(24/6.5) - c[3] * Turnovers[, 14] * (24/6.5)


mean(tt) *(252 * 24/6.5)/(sd(tt) * sqrt(252 * 24/6.5))

################## CONCENTRATION ##############

concentration <- matrix(0L, ncol = 14, nrow=1516)

concentration[, 1] <- portconcentration(rp_weigths_HAR)
concentration[, 2] <- portconcentration(rp_weigths_HAR_daily)
concentration[, 3] <- portconcentration(rp_weigths_HARJ)
concentration[, 4] <- portconcentration(rp_weigths_HARQJ)
concentration[, 5] <- portconcentration(rp_weigths_HARQ) 
concentration[, 6] <- portconcentration(rp_weigths_HARQF) 
concentration[, 7] <- portconcentration(rp_weigths_CHAR) 
concentration[, 8] <- portconcentration(rp_weigths_CHARQ) 
concentration[, 9] <- portconcentration(rp_weigths_rBG) 
concentration[, 10] <- portconcentration(rp_weigths_rBG_daily) 
concentration[, 11] <- portconcentration(rp_weigths_crBG) 
concentration[, 12] <- portconcentration(rp_weigths_rDCC) 
concentration[, 13] <- portconcentration(rp_weigths_crDCC)
concentration[, 14] <- portconcentration(rp_weigths_RM_daily) 


colMeans(concentration)




######################## TEST ON SMOOTHED harqf sigmas #######################

#alihn "right" gives the prior indexes

temp1 <- rollapply(sigmas_HARQF[1,1,], 100, mean, align = "right")
temp2 <- rollapply(sigmas_HARQF[2,1,], 100, mean, align = "right")
temp3 <- rollapply(sigmas_HARQF[2,2,], 100, mean, align = "right")


tempnew2 <- ewma.filter.realized(sigmas_HARQF, 20/252, F, NULL, 0)

days <- 1516 #- 100 + 1

tempnew <- array(0L, dim = c(2,2,days))
for(i in 1:days){

tempnew[,,i] <- matrix(c(temp1[i], temp2[i], temp2[i], temp3[i]), ncol =2, nrow = 2)

}

rp_weights_harsmooth <- matrix(0L, ncol = 3, nrow = days)
alpha_harsmooth <- numeric()


for(i in 1:days){

	rp_weights_harsmooth[i, ] <- t(riskparity_2dim(tempnew2[,,i] * (1/10000), 0.1/sqrt(252), T)$w)
	alpha_harsmooth[i] <- t(riskparity_2dim(tempnew2[,,i] * (1/10000), 0.1, T)$alpha)


}


portret_harqfsmooth <- numeric()
portvar_harqfsmooth <- numeric()
for(i in 1:days){


  portret_harqfsmooth[i] <- mergedret_daily[i+1000, ] %*% rp_weights_harsmooth[i, ]

  portvar_harqfsmooth[i] <-  sqrt(t(rp_weights_harsmooth[i, 1:2]) %*% (calccov[[1]][[7]][,,(i+1000)]) %*% (rp_weights_harsmooth[i, 1:2]))

}

mean(portret_harqfsmooth) * 100 * 252

mean((portret_HARQF))* (252) * 100 

sd(portvar_harqfsmooth)* sqrt(252) * 100

sd(portvar_HARQF) * sqrt(252) * 100


smoothed <- turnover(rp_weights_harsmooth[,1:2], matrix(portret_harqfsmooth), oosmergedret_daily[20:1516,1:2])

mean(smoothed)

############################### REALIZED PORTFOLIO VOL PLOT ######################
library(ggplot2)


#TO INFER SOME UNDERSTANDING IN THE DIFFERENCE OF THE PORTFOLIO VOLATILITY SMOOTH YOUR ESTIMATES.
#HOWEVER, IN APPENDIX PUT THE NON-SMOOTHED GRAPHS ALSO. 

#realized portfolio volatility
ggplot() + geom_line(aes(as.Date(getDates[1001:2516]), portvar_RM_daily*sqrt(252)*100, col = "RM")) + 
geom_line(aes(as.Date(getDates[1001:2516]), portvar_HARQF*sqrt(252)*100, col = "HARQF")) + ylab("Realized portfolio volatility") 
geom_line(aes(as.Date(getDates[1001:2516]), portvar_rBG_daily*sqrt(252)*100, col = "rBG daily")) 

ggplot() + geom_line(aes(as.Date(getDates[1001:2516]), portvar_harqfsmooth*sqrt(252)*100, col ="smoothed")) +
geom_line(aes(as.Date(getDates[1001:2516]), portvar_HARQF*sqrt(252)*100, col = "HARQF"))



tt <- rollapply(portvar_RM_daily*sqrt(252)  *100, 5, mean, align = "right")
tt2 <- rollapply(portvar_HARQF*sqrt(252)  *100, 5, mean, align = "right")
tt3 <- rollapply(portvar_CHAR*sqrt(252)  *100, 5, mean, align = "right")


library(ggplot2)

p1 <- ggplot() + geom_line(aes(as.Date(getDates[1005:2516]), tt, col ="RM")) + 
geom_line(aes(as.Date(getDates[1005:2516]), tt2, col ="HARQF")) + geom_hline(yintercept = 10, linetype="dashed", 
	col="steelblue", lwd = 1) + ylab("Portfolio Volatility") + xlab("") + 
theme(legend.title = element_blank(),legend.position = c(0.90, 0.83), legend.background = element_rect(fill="lightblue", size=0.5,
 linetype="solid"), 
	plot.title = element_text(hjust = 0.5, face = "bold"),  axis.title=element_text(size=12))



#PLOTTING CORRELATIONS:

cor_HARQF <- sigmas_HARQF[2,1, ] / sqrt(sigmas_HARQF[1,1, ] * sigmas_HARQF[2,2, ])
cor_RM <- sigmas_RM_daily[2,1, ] / sqrt(sigmas_RM_daily[1,1, ] * sigmas_RM_daily[2,2, ])

smooth_cor_HARQF <- rollapply(cor_HARQF, 5, mean, align = "right")
smooth_cor_RM <- rollapply(cor_RM, 5, mean, align = "right")


p2 <- ggplot() + geom_line(aes(as.Date(getDates[1005:2516]), smooth_cor_RM,col ="RM")) + 
geom_line(aes(as.Date(getDates[1005:2516]), smooth_cor_HARQF, col ="HARQF")) + geom_hline(yintercept = 0, linetype="dashed", 
	col="steelblue", lwd = 1) + ylab("Correlation") + xlab("Dates") + theme(legend.position = "none", text = element_text(size=8))



library(gridExtra)

grid.arrange(p1, p2, ncol = 1)

ggsave("correlationandvolplot.eps", device = "eps")




###########################non smoothed######################################
tt <- rollapply(portvar_RM_daily*sqrt(252)  *100, 1, mean, align = "right")
tt2 <- rollapply(portvar_HARQF*sqrt(252)  *100, 1, mean, align = "right")
tt3 <- rollapply(portvar_CHAR*sqrt(252)  *100, 1, mean, align = "right")


library(ggplot2)

p1 <- ggplot() + geom_line(aes(as.Date(getDates[1001:2516]), tt, col ="RM")) + 
geom_line(aes(as.Date(getDates[1001:2516]), tt2, col ="HARQF")) + geom_hline(yintercept = 10, linetype="dashed", 
	col="steelblue", lwd = 1) + ylab("Portfolio Volatility") + xlab("") + 
theme(legend.title = element_blank(),legend.position = c(0.90, 0.83), legend.background = element_rect(fill="lightblue", size=0.5,
 linetype="solid"), 
	plot.title = element_text(hjust = 0.5, face = "bold"),  axis.title=element_text(size=12))



#PLOTTING CORRELATIONS:

cor_HARQF <- sigmas_HARQF[2,1, ] / sqrt(sigmas_HARQF[1,1, ] * sigmas_HARQF[2,2, ])
cor_RM <- sigmas_RM_daily[2,1, ] / sqrt(sigmas_RM_daily[1,1, ] * sigmas_RM_daily[2,2, ])

smooth_cor_HARQF <- rollapply(cor_HARQF, 1, mean, align = "right")
smooth_cor_RM <- rollapply(cor_RM, 1, mean, align = "right")


p2 <- ggplot() + geom_line(aes(as.Date(getDates[1001:2516]), smooth_cor_RM,col ="RM")) + 
geom_line(aes(as.Date(getDates[1001:2516]), smooth_cor_HARQF, col ="HARQF")) + geom_hline(yintercept = 0, linetype="dashed", 
	col="steelblue", lwd = 1) + ylab("Correlation") + xlab("Dates") + theme(legend.position = "none", text = element_text(size=8))


library(gridExtra)

grid.arrange(p1, p2, ncol = 1)

ggsave("correlationandvolplot.eps", device = "eps")


###########################################################################
#
#					correlation sensitivity  
#
###########################################################################

#getting average vol for tlt and spy.

#both are already in percentage. 

meanvolTLT_HARQF <- mean(sqrt(sigmas_HARQF[1,1,]*252)) 

meanvolSPY_HARQF <-  mean(sqrt(sigmas_HARQF[2,2,]*252))

meanvolTLT_RM <- mean(sqrt(sigmas_RM_daily[1,1,]*252))

meanvolSPY_RM <-  mean(sqrt(sigmas_RM_daily[2,2,]*252))

correlations <- seq(-0.9,0.9,0.01)

sensitivity_HARQF <- numeric()

sensitivity_RM <- numeric()

mean5rcov <- matrix(c(mean(calccov[[1]][[7]][1,1,1001:2516]), mean(calccov[[1]][[7]][2,1,1001:2516]), mean(calccov[[1]][[7]][2,1,1001:2516]), 
	mean(calccov[[1]][[7]][2,2,1001:2516])), nrow = 2, ncol = 2)


for(i in 1:length(correlations)){

	newcovariance <- matrix(c(meanvolTLT_HARQF^2, meanvolTLT_HARQF*meanvolSPY_HARQF*correlations[i],
		meanvolTLT_HARQF*meanvolSPY_HARQF*correlations[i], meanvolSPY_HARQF^2), ncol = 2, nrow = 2)

	newcovariance2 <- matrix(c(meanvolTLT_RM^2, meanvolTLT_RM*meanvolSPY_RM*correlations[i],
		meanvolTLT_RM*meanvolSPY_RM*correlations[i], meanvolSPY_RM^2), ncol = 2, nrow = 2)

	w_HARQF <- t(riskparity_2dim(newcovariance,10,T)$w)

	w_RM <- t(riskparity_2dim(newcovariance2,10,T)$w)

	sensitivity_HARQF[i] <-  (w_HARQF[1]*w_HARQF[2]*meanvolTLT_HARQF*meanvolSPY_HARQF)  / (sqrt( t(w_HARQF[, 1:2]) %*% (mean5rcov) %*% (w_HARQF[, 1:2])) * sqrt(252)*100)
	sensitivity_RM[i] <-  (w_RM[1]*w_RM[2]*meanvolTLT_RM*meanvolSPY_RM)/(sqrt( t(w_RM[, 1:2]) %*% (mean5rcov) %*% (w_RM[, 1:2])) * sqrt(252)*100)


}



#Shows the volatility's sensitivity to correlation for the unlevered risk-parity portfolio. In essence, 
#upscaling and downscaling the portfolios using a leverage parameter only shifts the graph.

ggplot() + geom_line(aes(correlations, sensitivity_HARQF, col = "HARQF"), lwd = 1) + 
geom_line(aes(correlations, sensitivity_RM, col = "RM"), lwd = 1) + 
scale_x_continuous(breaks = round(seq(-0.9,0.9, by = 0.1),1)) + ylab("portfolio volatility (%)") + xlab("Correlation")



################FOR EACH DAY###################

correlations <- seq(-0.9,0.9,0.01)

sensitivity_HARQF <- matrix(0L, nrow = length(correlations), ncol = 1516)

sensitivity_RM <- matrix(0L, nrow = length(correlations), ncol = 1516)

volTLT_HARQF <- (sqrt(sigmas_HARQF[1,1,]*252)) 

volSPY_HARQF <-  (sqrt(sigmas_HARQF[2,2,]*252))

volTLT_RM <- (sqrt(sigmas_RM_daily[1,1,]*252))

volSPY_RM <-  (sqrt(sigmas_RM_daily[2,2,]*252))

for(j in 1:1516){
	for(i in 1:length(correlations)){

		newcovariance <- matrix(c(volTLT_HARQF[j]^2, volTLT_HARQF[j]*volSPY_HARQF[j]*correlations[i],
			volTLT_HARQF[j]*volSPY_HARQF[j]*correlations[i], volSPY_HARQF[j]^2), ncol = 2, nrow = 2)

		newcovariance2 <- matrix(c(volTLT_RM[j]^2, volTLT_RM[j]*volSPY_RM[j]*correlations[i],
			volTLT_RM[j]*volSPY_RM[j]*correlations[i], volSPY_RM[j]^2), ncol = 2, nrow = 2) 

		w_HARQF <- t(riskparity_2dim(newcovariance,10,T)$w)

		w_RM <- t(riskparity_2dim(newcovariance2,10,T)$w)

			#*volTLT_HARQF[j]*volSPY_HARQF[j]
			#*volTLT_RM[j]*volSPY_RM[j]
		realized_vol <- sqrt(calccov[[1]][[7]][1,1,j+1000]*10000*252)*sqrt(calccov[[1]][[7]][2,2,j+1000]*10000*252)
		sensitivity_HARQF[i,j] <-  (w_HARQF[1]*w_HARQF[2]*realized_vol)  / (sqrt( t(w_HARQF[, 1:2]) %*% calccov[[1]][[7]][,,1000+j] %*% (w_HARQF[, 1:2])) * sqrt(252)*100)
		sensitivity_RM[i,j] <-  (w_RM[1]*w_RM[2]*realized_vol)/(sqrt( t(w_RM[, 1:2]) %*% calccov[[1]][[7]][,,1000+j] %*% (w_RM[, 1:2])) * sqrt(252)*100)
	}
	print(sprintf("%s", j))
}

sensitivity_HARQF2 <- rowMeans(sensitivity_HARQF)
sensitivity_RM2 <- rowMeans(sensitivity_RM)

ggplot() + geom_line(aes(correlations, sensitivity_HARQF2, col = "HARQF"), lwd = 1) + 
geom_line(aes(correlations, sensitivity_RM2, col = "RM"), lwd = 1) + 
scale_x_continuous(breaks = round(seq(-0.9,0.9, by = 0.1),1)) + ylab("portfolio volatility (%)") + xlab("Correlation")








yt <- arima.sim(list(order=c(2,0,0), ar=c(.5, 0.49)), n=500)
ts.plot(yt)


yt <- merged_ret[2:2516,2]^2 + 0.4 * merged_ret[5:2511,2]^2 + 1.43 * merged_ret[22:2494,2]^2

ts.plot(yt)

tseries::adf.test(yt)













#testing whether you can multiply by -1 to flip the tail calculation 

VaR(dailyretotc[,1], p=.95, method="historical")

VaR(-1*dailyretotc[,1], p=.95, method="historical")


#constructing synthetic data..
#original lower tail has var around -1 and -2. Upper tail between 0.001 and 0.02. 
synth <- c(rep(-2,15), rep(-1,15), rep(0,30), rep(0.001,50), rep(0.02,10))

hist(synth, breaks = 80)
VaR(synth, p=.95, method="historical")

VaR(-1*synth, p=.95, method="historical")

#it seems that you can flip the histogram by multiplying with -1 



###################testing under minvar portfolio ###################
rp_weigths_HAR <- matrix(0L, nrow=1516, ncol = 2)
rp_weigths_HAR_daily <- matrix(0L, nrow=1516, ncol = 2)
rp_weigths_CHAR <- matrix(0L, nrow=1516, ncol = 2)
rp_weigths_CHARQ <- matrix(0L, nrow=1516, ncol = 2)
rp_weigths_HARQF <- matrix(0L, nrow=1516, ncol = 2)
rp_weigths_RM_daily <- matrix(0L, nrow=1516, ncol = 2)
rp_weigths_rBG <- matrix(0L, nrow=1516, ncol = 2)
rp_weigths_rBG_daily <- matrix(0L, nrow=1516, ncol = 2)
rp_weigths_crBG <- matrix(0L, nrow=1516, ncol = 2)
rp_weigths_rDCC <- matrix(0L, nrow=1516, ncol = 2)
rp_weigths_crDCC <- matrix(0L, nrow=1516, ncol = 2)

#run above code 
for(i in 1:1516){

  rp_weigths_HAR[i, ] <- t(minvar(sigmas_HAR[,,i]* (1/10000))) #* (1/10000)
  rp_weigths_HAR_daily[i, ] <- t(minvar(sigmas_HAR_daily[,,i]* (1/10000)))
  rp_weigths_CHAR[i, ] <- t(minvar(sigmas_CHAR[,,i]* (1/10000)))
  rp_weigths_CHARQ[i, ] <- t(minvar(sigmas_CHARQ[,,i]* (1/10000)))
  rp_weigths_HARQF[i, ] <- t(minvar(sigmas_HARQF[,,i]* (1/10000)))
  rp_weigths_RM_daily[i, ] <- t(minvar(sigmas_RM_daily[,,i]* (1/10000)))
  rp_weigths_rBG[i, ]   <- t(minvar(sigmas_rBG[,,i+1000]* (1/10000)))
  rp_weigths_crBG[i, ]  <- t(minvar(sigmas_crBG[,,i+1000]* (1/10000)))
  rp_weigths_rDCC[i, ]  <- t(minvar(sigmas_rDCC[,,i+1000]* (1/10000)))
  rp_weigths_crDCC[i, ]  <- t(minvar(sigmas_crDCC[,,i+1000]* (1/10000)))
  rp_weigths_rBG_daily[i, ] <- t(minvar(sigmas_rBG_daily[,,i+1000]* (1/10000)))
}

for(i in 1:1516){
  
  #portret_HAR[i] <- colMeans(mergedret_5min[[i+1000]]) %*% rp_weigths_HAR[i, ]
  #portret_HAR_daily[i] <- colMeans(mergedret_5min[[i+1000]]) %*% rp_weigths_HAR_daily[i, ]
  #portret_HARQF[i] <- colMeans(mergedret_5min[[i+1000]]) %*% rp_weigths_HARQF[i, ]  
  #portret_CHAR[i] <- colMeans(mergedret_5min[[i+1000]]) %*% rp_weigths_CHAR[i, ]
  #portret_CHARQ[i] <- colMeans(mergedret_5min[[i+1000]]) %*% rp_weigths_CHARQ[i, ]
  #portret_rBG[i] <- colMeans(mergedret_5min[[i+1000]]) %*% rp_weigths_rBG[i, ]
  #portret_rBG_daily[i] <- colMeans(mergedret_5min[[i+1000]]) %*% rp_weigths_rBG_daily[i, ]
  #portret_crBG[i] <- colMeans(mergedret_5min[[i+1000]]) %*% rp_weigths_crBG[i, ]
  #portret_rDCC[i] <- colMeans(mergedret_5min[[i+1000]]) %*% rp_weigths_rDCC[i, ]
  #portret_crDCC[i] <- colMeans(mergedret_5min[[i+1000]]) %*% rp_weigths_crDCC[i, ]
  #portret_RM_daily[i] <- colMeans(mergedret_5min[[i+1000]]) %*% rp_weigths_RM_daily[i, ]

  portret_HAR[i] <-mergedret_daily[i+1000, 1:2] %*% rp_weigths_HAR[i, ]
  portret_HAR_daily[i] <- mergedret_daily[i+1000, 1:2] %*% rp_weigths_HAR_daily[i, ]
  portret_HARQF[i] <- mergedret_daily[i+1000, 1:2] %*% rp_weigths_HARQF[i, ]  
  portret_CHAR[i] <- mergedret_daily[i+1000, 1:2] %*% rp_weigths_CHAR[i, ]
  portret_CHARQ[i] <- mergedret_daily[i+1000, 1:2] %*% rp_weigths_CHARQ[i, ]
  portret_rBG[i] <- mergedret_daily[i+1000, 1:2] %*% rp_weigths_rBG[i, ]
  portret_rBG_daily[i] <- mergedret_daily[i+1000, 1:2] %*% rp_weigths_rBG_daily[i, ]
  portret_crBG[i] <- mergedret_daily[i+1000, 1:2] %*% rp_weigths_crBG[i, ]
  portret_rDCC[i] <- mergedret_daily[i+1000, 1:2] %*% rp_weigths_rDCC[i, ]
  portret_crDCC[i] <- mergedret_daily[i+1000, 1:2] %*% rp_weigths_crDCC[i, ]
  portret_RM_daily[i] <- mergedret_daily[i+1000, 1:2] %*% rp_weigths_RM_daily[i, ]
	

  #these are the realized portfolio variance. We see how the forecasted weights act under our proxy
  portvar_HAR[i] <-  sqrt(t(rp_weigths_HAR[i, 1:2]) %*% (calccov[[1]][[7]][,,(i+1000)]*10000) %*% (rp_weigths_HAR[i, 1:2]))
  portvar_HAR_daily[i] <-  sqrt(t(rp_weigths_HAR_daily[i, 1:2]) %*% (calccov[[1]][[7]][,,(i+1000)]*10000) %*% (rp_weigths_HAR_daily[i, 1:2]))
  portvar_HARQF[i] <-  sqrt(t(rp_weigths_HARQF[i, 1:2]) %*% (calccov[[1]][[7]][,,(i+1000)]*10000) %*% (rp_weigths_HARQF[i, 1:2]))
  portvar_CHAR[i] <-  sqrt(t(rp_weigths_CHAR[i, 1:2]) %*% (calccov[[1]][[7]][,,(i+1000)]*10000) %*% (rp_weigths_CHAR[i, 1:2]))
  portvar_CHARQ[i] <-  sqrt(t(rp_weigths_CHARQ[i, 1:2]) %*% (calccov[[1]][[7]][,,(i+1000)]*10000) %*% (rp_weigths_CHARQ[i, 1:2]))
  portvar_rBG[i] <- sqrt(t(rp_weigths_rBG[i, 1:2]) %*% (calccov[[1]][[7]][,,(i+1000)]*10000) %*% (rp_weigths_rBG[i, 1:2]))
  portvar_rBG_daily[i] <- sqrt(t(rp_weigths_rBG_daily[i, 1:2]) %*% (calccov[[1]][[7]][,,(i+1000)]*10000) %*% (rp_weigths_rBG_daily[i, 1:2]))
  portvar_crBG[i] <- sqrt(t(rp_weigths_crBG[i, 1:2]) %*% (calccov[[1]][[7]][,,(i+1000)]*10000) %*% (rp_weigths_crBG[i, 1:2]))
  portvar_rDCC[i] <- sqrt(t(rp_weigths_rDCC[i, 1:2]) %*% (calccov[[1]][[7]][,,(i+1000)]*10000) %*% (rp_weigths_rDCC[i, 1:2]))
  portvar_crDCC[i] <- sqrt(t(rp_weigths_crDCC[i, 1:2]) %*% (calccov[[1]][[7]][,,(i+1000)]*10000) %*% (rp_weigths_crDCC[i, 1:2]))

  portvar_RM_daily[i] <- sqrt(t(rp_weigths_RM_daily[i, 1:2]) %*% (calccov[[1]][[7]][,,(i+1000)]*10000) %*% (rp_weigths_RM_daily[i, 1:2]))

}



tun_HAR <- turnover(rp_weigths_HAR, matrix(portret_HAR), mergedret_daily[1001:2516, 1:2]) 
tun_HARQF <- turnover(rp_weigths_HARQF, matrix(portret_HARQF), mergedret_daily[1001:2516, 1:2]) 
tun_RM <- turnover(rp_weigths_RM_daily, matrix(portret_RM_daily), mergedret_daily[1001:2516, 1:2]) 


mean(tun_HAR)
mean(tun_HARQF)
mean(tun_RM)


ttt <- portret_HAR*252 - tun_HAR*(252/2) * 0.01/252










##########################################################################################
#
#
# 			RANKING BASED ON SHARPE RATIO
#
#
##########################################################################################


sigmas_HARQ <- array(0L, dim = c(2,2,2516))
sigmas_HARQF <- array(0L, dim = c(2,2,2516))
sigmas_HAR <- array(0L, dim = c(2,2,2516))
sigmas_HARJ <- array(0L, dim = c(2,2,2516))
sigmas_HARQJ <- array(0L, dim = c(2,2,2516))
sigmas_HAR_daily <- array(0L, dim = c(2,2,2516))
sigmas_CHARQ <- array(0L, dim = c(2,2,2516))
sigmas_CHAR <- array(0L, dim = c(2,2,2516))

rp_weigths_HAR <- matrix(0L, nrow=1516, ncol = 3)
rp_weigths_HARQ <- matrix(0L, nrow=1516, ncol = 3)
rp_weigths_HARJ <- matrix(0L, nrow=1516, ncol = 3)
rp_weigths_HARQJ <- matrix(0L, nrow=1516, ncol = 3)
rp_weigths_CHAR <- matrix(0L, nrow=1516, ncol = 3)
rp_weigths_CHARQ <- matrix(0L, nrow=1516, ncol = 3)
rp_weigths_HARQF <- matrix(0L, nrow=1516, ncol = 3)

portret_HAR <- numeric()
portret_HARJ <- numeric()
portret_HARQJ <- numeric()
portret_HARQ <- numeric()
portret_HAR_daily <- numeric()
portret_HARQF <- numeric()
portret_CHAR <- numeric()
portret_CHARQ <- numeric()
portret_rBG <- numeric()

Turnovers <- matrix(0L, ncol = 14, nrow=1516)

#REMEMBER DAILIES!
rank1_names <- matrix(0L, ncol = 3, nrow = 9)
rank1_vals <- matrix(0L, ncol = 3, nrow = 9)

rank2_names <- matrix(0L, ncol = 3, nrow = 9)
rank2_vals <- matrix(0L, ncol = 3, nrow = 9)

rank3_names <- matrix(0L, ncol = 3, nrow = 9)
rank3_vals <- matrix(0L, ncol = 3, nrow = 9)

rank1_names_list <- list()
rank1_vals_list <- list()

rank2_names_list <- list()
rank2_vals_list <- list()

rank3_names_list <- list()
rank3_vals_list <- list()


c <- c(0, 0.01, 0.02)

portret_transcost <- list()
temp <- matrix(0L, ncol = 3, nrow = 1516)
sharpes <- matrix(0L, ncol = 7, nrow = 3)

Correlations_measures <- readRDS("Correlations_measures.rds")
Quarticities <- readRDS("Quarticity_estimators.rds")
window <- 1000
oosmergedret_daily <- matrix(mergedret_daily[1001:2516, ], nrow = 1516, ncol = 3)

sharpes_raw <- list()
sharpes_complete <- list()

for(k in 2:3){
	for(j in 1:9){ 
		for(i in (window+1):2516){

			#measures
			QCchoice <- c(1,7,8)
			ICchoice <- c(4,5,6)



			#############################################################################################
			#
			#
			#											PROXIES
			# 
			#
			#############################################################################################

			RV5min_TLT <- calccov[[1]][[7]][1,1,(i-window):(i-1)] * 10000
			RV5min_SPY <- calccov[[1]][[7]][2,2,(i-window):(i-1)] * 10000

			RV5min_TLT <- RV5min_TLT[22:length(RV5min_TLT)]
			RV5min_SPY <- RV5min_SPY[22:length(RV5min_SPY)]

			#proxy covariance for qlike losses
			RCov5min <- calccov[[1]][[7]][,,i] * 10000

			proxycor <- proxycorrelation[2,1, (i-window):(i-1)]
			#---------------------------------------------------------------------------------------------


			#############################################################################################
			#
			#
			#							PRELIMINARIES (LAGS AND QUARTICITIES)
			# 
			#
			#############################################################################################




			#Lags
			vol_TLT <- matrix((calccov[[QCchoice[k]]][[j]][1,1,(i-window):(i-1)])) * 10000  
			volday_TLT <- vol_TLT
			volweek_TLT <- rowMeans(cbind(vol_TLT, mlag(vol_TLT,4,mean(vol_TLT))))
			volmonth_TLT <- rowMeans(cbind(vol_TLT, mlag(vol_TLT,21,mean(vol_TLT))))

			vol_TLT <- vol_TLT[22:length(vol_TLT)]
			volday_TLT <- volday_TLT[22:length(volday_TLT)]
			volweek_TLT <- volweek_TLT[22:length(volweek_TLT)]
			volmonth_TLT <- volmonth_TLT[22:length(volmonth_TLT)]





			vol_SPY <- matrix((calccov[[QCchoice[k]]][[j]][2,2,(i-window):(i-1)])) * 10000  
			volday_SPY <- vol_SPY
			volweek_SPY <- rowMeans(cbind(vol_SPY, mlag(vol_SPY,4,mean(vol_SPY))))
			volmonth_SPY <- rowMeans(cbind(vol_SPY, mlag(vol_SPY,21,mean(vol_SPY))))

			vol_SPY <- vol_SPY[22:length(vol_SPY)]
			volday_SPY <- volday_SPY[22:length(volday_SPY)]
			volweek_SPY <- volweek_SPY[22:length(volweek_SPY)]
			volmonth_SPY <- volmonth_SPY[22:length(volmonth_SPY)]



			#jump robust lags:

			ICvol_TLT <- matrix((calccov[[ICchoice[k]]][[j]][1,1,(i-window):(i-1)])) * 10000  
			ICvolday_TLT <- ICvol_TLT
			ICvolweek_TLT <- rowMeans(cbind(ICvol_TLT, mlag(ICvol_TLT,4,mean(ICvol_TLT))))
			ICvolmonth_TLT <- rowMeans(cbind(ICvol_TLT, mlag(ICvol_TLT,21,mean(ICvol_TLT))))

			ICvol_TLT <- ICvol_TLT[22:length(ICvol_TLT)]
			ICvolday_TLT <- ICvolday_TLT[22:length(ICvolday_TLT)]
			ICvolweek_TLT <- ICvolweek_TLT[22:length(ICvolweek_TLT)]
			ICvolmonth_TLT <- ICvolmonth_TLT[22:length(ICvolmonth_TLT)]





			ICvol_SPY <- matrix((calccov[[ICchoice[k]]][[j]][2,2,(i-window):(i-1)])) * 10000  
			ICvolday_SPY <- ICvol_SPY
			ICvolweek_SPY <- rowMeans(cbind(ICvol_SPY, mlag(ICvol_SPY,4,mean(ICvol_SPY))))
			ICvolmonth_SPY <- rowMeans(cbind(ICvol_SPY, mlag(ICvol_SPY,21,mean(ICvol_SPY))))

			ICvol_SPY <- ICvol_SPY[22:length(ICvol_SPY)]
			ICvolday_SPY <- ICvolday_SPY[22:length(ICvolday_SPY)]
			ICvolweek_SPY <- ICvolweek_SPY[22:length(ICvolweek_SPY)]
			ICvolmonth_SPY <- ICvolmonth_SPY[22:length(ICvolmonth_SPY)]



			#RQ and TRQ from percentage returns
			rq_TLT <- Quarticities$rq_TLT[(i-window):(i-1),j]
			rq_SPY <- Quarticities$rq_SPY[(i-window):(i-1),j]

			rq_TLT <- matrix(rq_TLT)
			rq_SPY <- matrix(rq_SPY)

			trq_SPY <- Quarticities$trq_SPY[(i-window):(i-1),j]
			trq_TLT <- Quarticities$trq_TLT[(i-window):(i-1),j]

			trq_TLT <- trq_TLT[22:length(trq_TLT)]
			trq_SPY <- trq_SPY[22:length(trq_SPY)]


			sqrtrq_TLT <- sqrt(rq_TLT) - mean(sqrt(rq_TLT))
			sqrttrq_TLT <- sqrt(trq_TLT) - mean(sqrt(trq_TLT))

			sqrtrq_SPY <- sqrt(rq_SPY) - mean(sqrt(rq_SPY))
			sqrttrq_SPY<- sqrt(trq_SPY) - mean(sqrt(trq_SPY))


			#lagged for HARQF
			rqSPYweek <- rowMeans(cbind(rq_SPY, mlag(rq_SPY,4,mean(rq_SPY))))
			sqrtrq_SPYweek <- sqrt(rqSPYweek) - mean(sqrt(rqSPYweek))

			rqSPYmonth <- rowMeans(cbind(rq_SPY, mlag(rq_SPY,21,mean(rq_SPY))))
			sqrtrq_SPYmonth <- sqrt(rqSPYmonth) - mean(sqrt(rqSPYmonth))


			rqTLTweek <- rowMeans(cbind(rq_TLT, mlag(rq_TLT,4,mean(rq_TLT))))
			sqrtrq_TLTweek <- sqrt(rqTLTweek) - mean(sqrt(rqTLTweek))

			rqTLTmonth <- rowMeans(cbind(rq_TLT, mlag(rq_TLT,21,mean(rq_TLT))))
			sqrtrq_TLTmonth <- sqrt(rqTLTmonth) - mean(sqrt(rqTLTmonth))


			sqrtrq_TLT <- sqrtrq_TLT[22:length(sqrtrq_TLT)]
			sqrtrq_TLTweek <- sqrtrq_TLTweek[22:length(sqrtrq_TLTweek)]
			sqrtrq_TLTmonth <- sqrtrq_TLTmonth[22:length(sqrtrq_TLTmonth)]

			sqrtrq_SPY <- sqrtrq_SPY[22:length(sqrtrq_SPY)]
			sqrtrq_SPYweek <- sqrtrq_SPYweek[22:length(sqrtrq_SPYweek)]
			sqrtrq_SPYmonth <- sqrtrq_SPYmonth[22:length(sqrtrq_SPYmonth)]


			#jumpparams
			jumpparam_TLT <- ifelse(calccov[[1]][[j]][1,1,(i-window):(i-1)]*10000 - calccov[[5]][[j]][1,1,(i-window):(i-1)]*10000>0, 
				calccov[[1]][[j]][1,1,(i-window):(i-1)]*10000 - calccov[[5]][[j]][1,1,(i-window):(i-1)]*10000, 0)

			jumpparam_SPY <- ifelse(calccov[[1]][[j]][2,2,(i-window):(i-1)]*10000 - calccov[[5]][[j]][2,2,(i-window):(i-1)]*10000>0, 
				calccov[[1]][[j]][2,2,(i-window):(i-1)]*10000 - calccov[[5]][[j]][2,2,(i-window):(i-1)]*10000, 0)


			jumpparam_TLT <- jumpparam_TLT[22:length(jumpparam_TLT)]
			jumpparam_SPY <- jumpparam_SPY[22:length(jumpparam_SPY)]


			#############################################################################################
			#
			#
			#									DRD-HAR MODEL (RCOV, MRC, MRK)
			# 
			#
			#############################################################################################



			correlation <- Correlations_measures[[j]][(i-window):(i-1), QCchoice[k]]
			correlationjumprobust <- Correlations_measures[[j]][(i-window):(i-1), ICchoice[k]]


			#They should all have the same length:
			end <- length(RV5min_TLT)

			#hhat is your one period ahead forecast, therefore we use "end" since these are the true "t" values, while
			#coef(HAR_...) are estimated until time t-1

			HAR_TLT <- lm(RV5min_TLT[2:end] ~ volday_TLT[1:(end-1)] + volweek_TLT[1:(end-1)] + volmonth_TLT[1:(end-1)])
			hhat_HAR_TLT <- cbind(1, volday_TLT[end], volweek_TLT[end], volmonth_TLT[end])  %*% matrix(coef(HAR_TLT))

			HAR_SPY <- lm(RV5min_SPY[2:end] ~ volday_SPY[1:(end-1)] + volweek_SPY[1:(end-1)] + volmonth_SPY[1:(end-1)])
			hhat_HAR_SPY <- cbind(1, volday_SPY[end], volweek_SPY[end], volmonth_SPY[end])  %*% matrix(coef(HAR_SPY))

			hhat_HAR_TLT <- volatility.insanity.filter(hhat_HAR_TLT, min(RV5min_TLT), max(RV5min_TLT), mean(RV5min_TLT))$vol
			hhat_HAR_SPY <- volatility.insanity.filter(hhat_HAR_SPY, min(RV5min_SPY), max(RV5min_SPY), mean(RV5min_SPY))$vol

			#warning is just that R does not like to add one 1 parameter array to a number eg. array(1, c(1,1,1)) + 2.  
			DRD_HAR <- suppressWarnings(EstimatecorrHAR(cbind(hhat_HAR_TLT, hhat_HAR_SPY), 
				correlation = correlation, proxy = proxycor, 0, T))

			sigmas_HAR[,,i] <- DRD_HAR$vSigma2

			#############################################################################################
			#
			#
			#									DRD-HARQ MODEL (RCOV, MRC, MRK)
			# 
			#
			#############################################################################################

			HARQ_TLT <- lm(RV5min_TLT[2:end] ~ volday_TLT[1:(end-1)] +  volweek_TLT[1:(end-1)] + volmonth_TLT[1:(end-1)] + I(volday_TLT[1:(end-1)] * sqrtrq_TLT[1:(end-1)]))
			hhat_HARQ_TLT <- cbind(1, volday_TLT[end], volweek_TLT[end], volmonth_TLT[end],volday_TLT[end] * sqrtrq_TLT[end])  %*% matrix(coef(HARQ_TLT))


			HARQ_SPY <- lm(RV5min_SPY[2:end] ~ volday_SPY[1:(end-1)] +  volweek_SPY[1:(end-1)] + volmonth_SPY[1:(end-1)] + I(volday_SPY[1:(end-1)] * sqrtrq_SPY[1:(end-1)]))
			hhat_HARQ_SPY <- cbind(1, volday_SPY[end], volweek_SPY[end], volmonth_SPY[end],volday_SPY[end] * sqrtrq_SPY[end])  %*% matrix(coef(HARQ_SPY))

			hhat_HARQ_TLT <- volatility.insanity.filter(hhat_HARQ_TLT, min(RV5min_TLT), max(RV5min_TLT), mean(RV5min_TLT))$vol
			hhat_HARQ_SPY <- volatility.insanity.filter(hhat_HARQ_SPY, min(RV5min_SPY), max(RV5min_SPY), mean(RV5min_SPY))$vol

			#SAME CORRELATION AS HAR MODELS. ONLY THING THAT CHANGED IS THE UNIVARIATE VOLS. 
			DRD_HARQ <- suppressWarnings(EstimatecorrHAR(cbind(hhat_HARQ_TLT, hhat_HARQ_SPY), correlation = correlation, proxy = proxycor, 0,T))

			sigmas_HARQ[,,i] <- DRD_HARQ$vSigma2

			#############################################################################################
			#
			#
			#									DRD-HARQF MODEL (RCOV, MRC, MRK)
			# 
			#
			#############################################################################################
			
			HARQF_TLT <- lm(RV5min_TLT[2:end] ~ volday_TLT[1:(end-1)] +  volweek_TLT[1:(end-1)] + 
			volmonth_TLT[1:(end-1)] + I(volday_TLT[1:(end-1)] * sqrtrq_TLT[1:(end-1)]) + 
			I(volweek_TLT[1:(end-1)] * sqrtrq_TLTweek[1:(end-1)]) + I(volmonth_TLT[1:(end-1)] * sqrtrq_TLTmonth[1:(end-1)]))

			hhat_HARQF_TLT <- cbind(1, volday_TLT[end], volweek_TLT[end], volmonth_TLT[end],
			volday_TLT[end] * sqrtrq_TLT[end], volweek_TLT[end] * sqrtrq_TLTweek[end],
			 volmonth_TLT[end] * sqrtrq_TLTmonth[end])  %*% matrix(coef(HARQF_TLT))


			HARQF_SPY <- lm(RV5min_SPY[2:end] ~ volday_SPY[1:(end-1)] +  volweek_SPY[1:(end-1)] + 
			volmonth_SPY[1:(end-1)] + I(volday_SPY[1:(end-1)] * sqrtrq_SPY[1:(end-1)]) + 
			I(volweek_SPY[1:(end-1)] * sqrtrq_SPYweek[1:(end-1)]) + I(volmonth_SPY[1:(end-1)] * sqrtrq_SPYmonth[1:(end-1)]))

			hhat_HARQF_SPY <- cbind(1, volday_SPY[end], volweek_SPY[end], volmonth_SPY[end],
			volday_SPY[end] * sqrtrq_SPY[end], volweek_SPY[end] * sqrtrq_SPYweek[end], 
			volmonth_SPY[end] * sqrtrq_SPYmonth[end])  %*% matrix(coef(HARQF_SPY))

			hhat_HARQF_SPY <- volatility.insanity.filter(hhat_HARQF_SPY, min(RV5min_SPY), max(RV5min_SPY), mean(RV5min_SPY))$vol
			hhat_HARQF_TLT <- volatility.insanity.filter(hhat_HARQF_TLT, min(RV5min_TLT), max(RV5min_TLT), mean(RV5min_TLT))$vol


			DRD_HARQF <- suppressWarnings(EstimatecorrHAR(cbind(hhat_HARQF_TLT, hhat_HARQF_SPY), correlation = correlation, proxy = proxycor, 0, T))

			sigmas_HARQF[,,i] <- DRD_HARQF$vSigma2

			#############################################################################################
			#
			#
			#									DRD-HARJ MODEL (RCOV, MRC, MRK)
			# 
			#
			#############################################################################################


			HARJ_TLT <- lm(RV5min_TLT[2:end] ~ volday_TLT[1:(end-1)] + volweek_TLT[1:(end-1)] + volmonth_TLT[1:(end-1)] + jumpparam_TLT[1:(end-1)])
			hhat_HARJ_TLT <- cbind(1, volday_TLT[end], volweek_TLT[end], volmonth_TLT[end],jumpparam_TLT[end])  %*% matrix(coef(HARJ_TLT))

			HARJ_SPY <- lm(RV5min_SPY[2:end] ~ volday_SPY[1:(end-1)] + volweek_SPY[1:(end-1)] + volmonth_SPY[1:(end-1)] + jumpparam_SPY[1:(end-1)])
			hhat_HARJ_SPY <- cbind(1, volday_SPY[end], volweek_SPY[end], volmonth_SPY[end],jumpparam_SPY[end])  %*% matrix(coef(HARJ_SPY))

			hhat_HARJ_SPY <- volatility.insanity.filter(hhat_HARJ_SPY, min(RV5min_SPY), max(RV5min_SPY), mean(RV5min_SPY))$vol
			hhat_HARJ_TLT <- volatility.insanity.filter(hhat_HARJ_TLT, min(RV5min_TLT), max(RV5min_TLT), mean(RV5min_TLT))$vol

			DRD_HARJ <- suppressWarnings(EstimatecorrHAR(cbind(hhat_HARJ_TLT, hhat_HARJ_SPY), correlation = correlation, proxy = proxycor, 0, T))
			
			sigmas_HARJ[,,i] <- DRD_HARJ$vSigma2

			#############################################################################################
			#
			#
			#									DRD-HARQJ MODEL (RCOV, MRC, MRK)
			# 
			#
			#############################################################################################

			HARQJ_TLT <- lm(RV5min_TLT[2:end] ~ volday_TLT[1:(end-1)] +  volweek_TLT[1:(end-1)] + volmonth_TLT[1:(end-1)] + I(volday_TLT[1:(end-1)] * sqrtrq_TLT[1:(end-1)]) + jumpparam_TLT[1:(end-1)])
			hhat_HARQJ_TLT <- cbind(1, volday_TLT[end], volweek_TLT[end], volmonth_TLT[end],volday_TLT[end] * sqrtrq_TLT[end], jumpparam_TLT[end])  %*% matrix(coef(HARQJ_TLT))

			HARQJ_SPY <- lm(RV5min_SPY[2:end] ~ volday_SPY[1:(end-1)] +  volweek_SPY[1:(end-1)] + volmonth_SPY[1:(end-1)] + I(volday_SPY[1:(end-1)] * sqrtrq_SPY[1:(end-1)]) + jumpparam_SPY[1:(end-1)])
			hhat_HARQJ_SPY <- cbind(1, volday_SPY[end], volweek_SPY[end], volmonth_SPY[end],volday_SPY[end] * sqrtrq_SPY[end], jumpparam_SPY[end])  %*% matrix(coef(HARQJ_SPY))

			hhat_HARQJ_SPY <- volatility.insanity.filter(hhat_HARQJ_SPY, min(RV5min_SPY), max(RV5min_SPY), mean(RV5min_SPY))$vol
			hhat_HARQJ_TLT <- volatility.insanity.filter(hhat_HARQJ_TLT, min(RV5min_TLT), max(RV5min_TLT), mean(RV5min_TLT))$vol

			DRD_HARQJ <- suppressWarnings(EstimatecorrHAR(cbind(hhat_HARQJ_TLT, hhat_HARQJ_SPY), correlation = correlation, proxy = proxycor, 0, T))

			sigmas_HARQJ[,,i] <- DRD_HARQJ$vSigma2

			#############################################################################################
			#
			#
			#									DRD-CHAR MODEL (TCOV, BPCOV, PBPCOV)
			# 
			#
			#############################################################################################

			CHAR_TLT <- lm(RV5min_TLT[2:end] ~ ICvolday_TLT[1:(end-1)] + ICvolweek_TLT[1:(end-1)] + ICvolmonth_TLT[1:(end-1)])
			hhat_CHAR_TLT <- cbind(1, volday_TLT[end], ICvolweek_TLT[end], ICvolmonth_TLT[end])  %*% matrix(coef(CHAR_TLT))

			CHAR_SPY <- lm(RV5min_SPY[2:end] ~ ICvolday_SPY[1:(end-1)] + ICvolweek_SPY[1:(end-1)] + ICvolmonth_SPY[1:(end-1)])
			hhat_CHAR_SPY <- cbind(1, ICvolday_SPY[end], ICvolweek_SPY[end], ICvolmonth_SPY[end])  %*% matrix(coef(CHAR_SPY))

			hhat_CHAR_SPY <- volatility.insanity.filter(hhat_CHAR_SPY, min(RV5min_SPY), max(RV5min_SPY), mean(RV5min_SPY))$vol
			hhat_CHAR_TLT <- volatility.insanity.filter(hhat_CHAR_TLT, min(RV5min_TLT), max(RV5min_TLT), mean(RV5min_TLT))$vol

			DRD_CHAR <- suppressWarnings(EstimatecorrHAR(cbind(hhat_CHAR_TLT, hhat_CHAR_SPY), correlation = correlationjumprobust, proxy = proxycor, 0, T))

			sigmas_CHAR[,,i] <- DRD_CHAR$vSigma2

			#############################################################################################
			#
			#
			#									DRD-CHARQ MODEL (TCOV, BPCOV, PBPCOV)
			# 
			#
			#############################################################################################

			CHARQ_TLT <- lm(RV5min_TLT[2:end] ~ ICvolday_TLT[1:(end-1)] +  ICvolweek_TLT[1:(end-1)] + ICvolmonth_TLT[1:(end-1)] + I(ICvolday_TLT[1:(end-1)] * sqrttrq_TLT[1:(end-1)]))
			hhat_CHARQ_TLT <- cbind(1, ICvolday_TLT[end], ICvolweek_TLT[end], ICvolmonth_TLT[end], ICvolday_TLT[end] * sqrttrq_TLT[end])  %*% matrix(coef(CHARQ_TLT))

			CHARQ_SPY <- lm(RV5min_SPY[2:end] ~ ICvolday_SPY[1:(end-1)] +  ICvolweek_SPY[1:(end-1)] + ICvolmonth_SPY[1:(end-1)] + I(ICvolday_SPY[1:(end-1)] * sqrttrq_SPY[1:(end-1)]))
			hhat_CHARQ_SPY <- cbind(1, ICvolday_SPY[end], ICvolweek_SPY[end], ICvolmonth_SPY[end], ICvolday_SPY[end] * sqrttrq_SPY[end])  %*% matrix(coef(CHARQ_SPY))

			hhat_CHARQ_SPY <- volatility.insanity.filter(hhat_CHARQ_SPY, min(RV5min_SPY), max(RV5min_SPY), mean(RV5min_SPY))$vol
			hhat_CHARQ_TLT <- volatility.insanity.filter(hhat_CHARQ_TLT, min(RV5min_TLT), max(RV5min_TLT), mean(RV5min_TLT))$vol


			DRD_CHARQ <- suppressWarnings(EstimatecorrHAR(cbind(hhat_CHARQ_TLT, hhat_CHARQ_SPY), correlation = correlationjumprobust, proxy = proxycor, 0, T))

			sigmas_CHARQ[,,i] <- DRD_CHARQ$vSigma2


			print(sprintf("Forecast: %s, Frequency: %s, Measure: %s", i,j,k))

		}


		sigmas_HARQ_temp <- sigmas_HARQ[,,1001:2516]
		sigmas_HARQJ_temp <- sigmas_HARQJ[,,1001:2516]
		sigmas_HARQF_temp <- sigmas_HARQF[,,1001:2516]
		sigmas_HAR_temp <- sigmas_HAR[,,1001:2516]
		sigmas_HARJ_temp <- sigmas_HARJ[,,1001:2516]
		sigmas_CHAR_temp <- sigmas_CHAR[,,1001:2516]
		sigmas_CHARQ_temp <- sigmas_CHARQ[,,1001:2516]


	for(i in 1:1516){
  		rp_weigths_HAR[i, ] <- t(riskparity_2dim(sigmas_HAR_temp[,,i]*(1/10000), 0.1/(sqrt(252)), T)$w)
 		rp_weigths_HARJ[i, ] <- t(riskparity_2dim(sigmas_HARJ_temp[,,i]* (1/10000), 0.1/(sqrt(252)), T)$w)
  		rp_weigths_HARQJ[i, ] <- t(riskparity_2dim(sigmas_HARQJ_temp[,,i] * (1/10000), 0.1/(sqrt(252)), T)$w)
  		rp_weigths_CHAR[i, ] <- t(riskparity_2dim(sigmas_CHAR_temp[,,i]* (1/10000), 0.1/(sqrt(252)), T)$w)
  		rp_weigths_CHARQ[i, ] <- t(riskparity_2dim(sigmas_CHARQ_temp[,,i]* (1/10000), 0.1/(sqrt(252)), T)$w)
  		rp_weigths_HARQ[i, ] <- t(riskparity_2dim(sigmas_HARQ_temp[,,i]* (1/10000), 0.1/(sqrt(252)), T)$w)
  		rp_weigths_HARQF[i, ] <- t(riskparity_2dim(sigmas_HARQF_temp[,,i]* (1/10000), 0.1/(sqrt(252)), T)$w)
	}

	for(i in 1:1516){
  
  		portret_HAR[i] <- (mergedret_daily[i+1000, ]) %*% rp_weigths_HAR[i, ]
  		portret_HARJ[i] <-mergedret_daily[i+1000, ] %*% rp_weigths_HARJ[i, ]
  		portret_HARQJ[i] <-mergedret_daily[i+1000, ] %*% rp_weigths_HARQJ[i, ]
 		portret_HARQF[i] <- mergedret_daily[i+1000, ] %*% rp_weigths_HARQF[i, ]  
  		portret_CHAR[i] <- mergedret_daily[i+1000, ] %*% rp_weigths_CHAR[i, ]
  		portret_CHARQ[i] <- mergedret_daily[i+1000, ] %*% rp_weigths_CHARQ[i, ]
  		portret_HARQ[i] <- mergedret_daily[i+1000, ] %*% rp_weigths_HARQ[i, ]  
	}

	portret <- list(portret_HAR,
  				portret_HARJ,
  				portret_HARQJ,
  				portret_HARQ,
  				portret_HARQF,
  				portret_CHAR,
   				portret_CHARQ)


	Turnovers[, 1] <- turnover(rp_weigths_HAR[,1:2], matrix(portret_HAR), oosmergedret_daily[,1:2])
	Turnovers[, 2] <- turnover(rp_weigths_HARJ[,1:2], matrix(portret_HARJ), oosmergedret_daily[,1:2])
	Turnovers[, 3] <- turnover(rp_weigths_HARQJ[,1:2], matrix(portret_HARQJ), oosmergedret_daily[,1:2])
	Turnovers[, 4] <- turnover(rp_weigths_HARQ[,1:2], matrix(portret_HARQ), oosmergedret_daily[,1:2]) 
	Turnovers[, 5] <- turnover(rp_weigths_HARQF[,1:2], matrix(portret_HARQF), oosmergedret_daily[,1:2]) 
	Turnovers[, 6] <- turnover(rp_weigths_CHAR[,1:2], matrix(portret_CHAR), oosmergedret_daily[,1:2]) 
	Turnovers[, 7] <- turnover(rp_weigths_CHARQ[,1:2], matrix(portret_CHARQ), oosmergedret_daily[,1:2]) 




	for(f in 1:7){
		for(i in 1:3){

			temp[, i] <- portret[[f]]  - c[i]/252 * Turnovers[, f]
		}

		portret_transcost[[f]] <- temp 

		sharpes[, f] <- t(( (colMeans(portret_transcost[[f]])*252*(24/6.5)) / 
			(apply(portret_transcost[[f]],MARGIN = c(2), FUN = function(x) sd(x))*sqrt(252)*sqrt(24/6.5)) ))

		colnames(sharpes) <- c("HAR", "HARJ", "HARQJ", "HARQ", "HARQF", "CHAR", "CHARQ")
		}	

	rank1_names[j, ] <- names(sort(sharpes[1, ], decreasing =  T)[1:3])
	rank1_vals[j, ] <- sort(sharpes[1, ], decreasing =  T)[1:3]

	rank2_names[j, ] <- names(sort(sharpes[2, ], decreasing =  T)[1:3])
	rank2_vals[j, ] <- sort(sharpes[2, ], decreasing =  T)[1:3]

	rank3_names[j, ] <- names(sort(sharpes[3, ], decreasing =  T)[1:3])
	rank3_vals[j, ] <- sort(sharpes[3, ], decreasing =  T)[1:3]
	
	sharpes_raw[[j]] <- sharpes

	}

	sharpes_complete[[k]] <- sharpes_raw

	rank1_names_list[[k]] <- rank1_names
	rank1_vals_list[[k]] <- rank1_vals

	rank2_names_list[[k]] <- rank2_names
	rank2_vals_list[[k]] <- rank2_vals

	rank3_names_list[[k]] <- rank3_names
	rank3_vals_list[[k]] <- rank3_vals

	
}


notransactioncost <- do.call(rbind, rank1_names_list)

rank1instances <-matrix(c(
sum(colSums(notransactioncost == "HAR")),
sum(colSums(notransactioncost == "HARQ")),
sum(colSums(notransactioncost == "HARQF")),
sum(colSums(notransactioncost == "HARJ")),
sum(colSums(notransactioncost == "HARQJ")),
sum(colSums(notransactioncost == "CHAR")),
sum(colSums(notransactioncost == "CHARQ"))), nrow = 1, ncol = 7)

colnames(rank1instances) <- c("HAR", "HARQ", "HARQF", "HARJ", "HARQJ", "CHAR", "CHARQ")
rank1instances

onepercenttransactioncost <- do.call(rbind, rank2_names_list)

rank2instances <-matrix(c(
sum(colSums(onepercenttransactioncost == "HAR")),
sum(colSums(onepercenttransactioncost == "HARQ")),
sum(colSums(onepercenttransactioncost == "HARQF")),
sum(colSums(onepercenttransactioncost == "HARJ")),
sum(colSums(onepercenttransactioncost == "HARQJ")),
sum(colSums(onepercenttransactioncost == "CHAR")),
sum(colSums(onepercenttransactioncost == "CHARQ"))), nrow = 1, ncol = 7)

colnames(rank2instances) <- c("HAR", "HARQ", "HARQF", "HARJ", "HARQJ", "CHAR", "CHARQ")
rank2instances


twopercenttransactioncost <- do.call(rbind, rank3_names_list)

rank3instances <-matrix(c(
sum(colSums(twopercenttransactioncost == "HAR")),
sum(colSums(twopercenttransactioncost == "HARQ")),
sum(colSums(twopercenttransactioncost == "HARQF")),
sum(colSums(twopercenttransactioncost == "HARJ")),
sum(colSums(twopercenttransactioncost == "HARQJ")),
sum(colSums(twopercenttransactioncost == "CHAR")),
sum(colSums(twopercenttransactioncost == "CHARQ"))), nrow = 1, ncol = 7)

colnames(rank3instances) <- c("HAR", "HARQ", "HARQF", "HARJ", "HARQJ", "CHAR", "CHARQ")
rank3instances


colSums(matrix(c(rank1instances, rank2instances, rank3instances), ncol=7,nrow=3))


vals1 <- do.call(cbind, rank1_vals_list)

sharpes_complete_l1 <-  do.call(cbind, sharpes_complete[[1]])
sharpes_complete_l2 <-  do.call(cbind, sharpes_complete[[2]])
sharpes_complete_l3 <-  do.call(cbind, sharpes_complete[[3]])

all <- cbind(sharpes_complete_l1, sharpes_complete_l2, sharpes_complete_l3)

sort(all[1, ], decreasing = T)[1:10]
sort(all[2, ], decreasing = T)[1:10]
sort(all[3, ], decreasing = T)[1:10]
