

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

calccov <- readRDS("calculatedcovariances.rds")

mergedfrequencies <- readRDS("mergedfrequencies.rds")

#5min proxy correlations

proxycorrelation <- array(0L, dim = c(2,2,2516))

for(i in 1:2516){

	proxycorrelation[,,i] <- realCov(mergedfrequencies[[7]][[i]]*100, T)

}




########################################################################################################################
#
#										OUT OF SAMPLE FORECAST ANALYSIS (ROLLING WINDOW)
#
########################################################################################################################



#[[1]],.. contains sample frequencies from 1sec to 30 min. and each column in [[i]] contains the measure following same 
#style as calccov. 
Correlations_measures <- readRDS("Correlations_measures.rds")
Quarticities <- readRDS("Quarticity_estimators.rds")
#works. But for one measure. Thus you need to wrap it into another for-loop looping over frequencies. And then another looping over
#measures.

#DUE TO INSTABILITY ISSUES, CALCULATING OUT-OF-SAMPLE LOSSES FOR DAILY FREQUENCIES WAS DONE SEPARATELY. 

#in here you can define HAR, HARQ, HARQF, HARJ, HARQJ, CHAR, CHARQ

#HAR uses Rcov, MRC, and MRK
#HARQ: Rcov, MRC, MRK
#HARQF: Rcov, MRC, MRK
#HARJ: RCov, MRC, MRK, 
#HARQJ: RCov, MRC, MRK
#CHAR: TCov, BPCov, PBPCov
#CHARQ: TCov, BPCov, PBPCov

#Do SHAR separately. 

window <- 1000

Qlikes <- list()

temp <- matrix(0L, ncol = 9, nrow = 2516)
temp_HARQ <- matrix(0L, ncol = 9, nrow = 2516)
temp_HARQF <- matrix(0L, ncol = 9, nrow = 2516)
temp_HARJ <- matrix(0L, ncol = 9, nrow = 2516)
temp_HARQJ <- matrix(0L, ncol = 9, nrow = 2516)
temp_CHAR <- matrix(0L, ncol = 9, nrow = 2516)
temp_CHARQ <- matrix(0L, ncol = 9, nrow = 2516)

#minvar 

mvpvariance_DRD <- list()

temp_mvp <- matrix(0L, ncol = 9, nrow = 2516)
temp_HARQ_mvp <- matrix(0L, ncol = 9, nrow = 2516)
temp_HARQF_mvp <- matrix(0L, ncol = 9, nrow = 2516)
temp_HARJ_mvp <- matrix(0L, ncol = 9, nrow = 2516)
temp_HARQJ_mvp <- matrix(0L, ncol = 9, nrow = 2516)
temp_CHAR_mvp <- matrix(0L, ncol = 9, nrow = 2516)
temp_CHARQ_mvp <- matrix(0L, ncol = 9, nrow = 2516)



for(k in 1:3){
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

			#temp[i,j] <- QLIKE(DRD_HAR$vSigma2[,,1], RCov5min, 2)

			weight <- minvar(DRD_HAR$vSigma2[,,1])

			temp_mvp[i,j] <- t(weight) %*% RCov5min %*% weight

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

			#temp_HARQ[i,j] <- QLIKE(DRD_HARQ$vSigma2[,,1], RCov5min, 2)
			weight <- minvar(DRD_HARQ$vSigma2[,,1])

			temp_HARQ_mvp[i,j] <- t(weight) %*% RCov5min %*% weight


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

			#temp_HARQF[i,j] <- QLIKE(DRD_HARQF$vSigma2[,,1], RCov5min, 2)
			weight <- minvar(DRD_HARQF$vSigma2[,,1])

			temp_HARQF_mvp[i,j] <- t(weight) %*% RCov5min %*% weight

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
			
			#temp_HARJ[i,j] <- QLIKE(DRD_HARJ$vSigma2[,,1], RCov5min, 2)

			weight <- minvar(DRD_HARJ$vSigma2[,,1])

			temp_HARJ_mvp[i,j] <- t(weight) %*% RCov5min %*% weight

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

			#temp_HARQJ[i,j] <- QLIKE(DRD_HARQJ$vSigma2[,,1], RCov5min, 2)

			weight <- minvar(DRD_HARQJ$vSigma2[,,1])

			temp_HARQJ_mvp[i,j] <- t(weight) %*% RCov5min %*% weight

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

			temp_CHAR[i,j] <- QLIKE(DRD_CHAR$vSigma2[,,1], RCov5min, 2)

			weight <- minvar(DRD_CHAR$vSigma2[,,1])

			temp_CHAR_mvp[i,j] <- t(weight) %*% RCov5min %*% weight

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

			temp_CHARQ[i,j] <- QLIKE(DRD_CHARQ$vSigma2[,,1], RCov5min, 2)
			weight <- minvar(DRD_CHARQ$vSigma2[,,1])

			temp_CHARQ_mvp[i,j] <- t(weight) %*% RCov5min %*% weight

			print(sprintf("Forecast: %s, Frequency: %s, Measure: %s", i,j,k))

		}
	}
	#Qlikes[[k]] <- temp[1001:nrow(temp), ]
	#Qlikes[[k+3]] <- temp_HARQ[1001:nrow(temp_HARQ), ]
	#Qlikes[[k+6]] <- temp_HARQF[1001:nrow(temp_HARQF), ]
	#Qlikes[[k+9]] <- temp_HARJ[1001:nrow(temp_HARJ), ]
	#Qlikes[[k+12]] <- temp_HARQJ[1001:nrow(temp_HARQJ), ]
	#putting SHAR model at 16th list item 
	#Qlikes[[k+16]] <- temp_CHAR[1001:nrow(temp_CHAR), ]
	#Qlikes[[k+19]] <- temp_CHARQ[1001:nrow(temp_CHARQ), ]

	mvpvariance_DRD[[k]] <- temp_mvp[1001:nrow(temp_mvp), ]
	mvpvariance_DRD[[k+3]] <- temp_HARQ_mvp[1001:nrow(temp_HARQ_mvp), ]
	mvpvariance_DRD[[k+6]] <- temp_HARQF_mvp[1001:nrow(temp_HARQF_mvp), ]
	mvpvariance_DRD[[k+9]] <- temp_HARJ_mvp[1001:nrow(temp_HARJ_mvp), ]
	mvpvariance_DRD[[k+12]] <- temp_HARQJ_mvp[1001:nrow(temp_HARQJ_mvp), ]
	mvpvariance_DRD[[k+16]] <- temp_CHAR_mvp[1001:nrow(temp_CHAR_mvp), ]
	mvpvariance_DRD[[k+19]] <- temp_CHARQ_mvp[1001:nrow(temp_CHARQ_mvp), ]


}


temp_SHAR <- matrix(0L, ncol = 9, nrow = 2516)
temp_SHAR_mvp <- matrix(0L, ncol = 9, nrow = 2516)

for(j in 1:9){
	for(i in (window+1):2516){

		voldayposvar_TLT <- matrix((calccov[[2]][[j]][1,1,(i-window):(i-1)])) * 10000 #sqrt
		voldaynegvar_TLT <- matrix((calccov[[3]][[j]][1,1,(i-window):(i-1)])) * 10000 #sqrt
		voldayposvar_SPY <- matrix((calccov[[2]][[j]][2,2,(i-window):(i-1)])) * 10000 #sqrt
		voldaynegvar_SPY <- matrix((calccov[[3]][[j]][2,2,(i-window):(i-1)])) * 10000 #sqrt

		voldayposvar_TLT <- voldayposvar_TLT[22:length(voldayposvar_TLT)]
		voldaynegvar_TLT <- voldaynegvar_TLT[22:length(voldaynegvar_TLT)]
		voldayposvar_SPY <- voldayposvar_SPY[22:length(voldayposvar_SPY)]
		voldaynegvar_SPY <- voldaynegvar_SPY[22:length(voldaynegvar_SPY)]

		correlation <- Correlations_measures[[j]][(i-window):(i-1), 1]

		proxycor <- proxycorrelation[2,1, (i-window):(i-1)]

		RV5min_TLT <- calccov[[1]][[7]][1,1,(i-window):(i-1)] * 10000
		RV5min_SPY <- calccov[[1]][[7]][2,2,(i-window):(i-1)] * 10000

		RV5min_TLT <- RV5min_TLT[22:length(RV5min_TLT)]
		RV5min_SPY <- RV5min_SPY[22:length(RV5min_SPY)]

		#proxy covariance for qlike losses
		RCov5min <- calccov[[1]][[7]][,,i] * 10000

		#############################################################################################
		#
		#
		#											DRD-SHAR MODEL 
		# 
		#
		#############################################################################################


		SHAR_TLT  <- lm(RV5min_TLT[2:end] ~ voldayposvar_TLT[1:(end-1)] + voldaynegvar_TLT[1:(end-1)] + volweek_TLT[1:(end-1)] + volmonth_TLT[1:(end-1)])
		hhat_SHAR_TLT <- cbind(1, voldayposvar_TLT[end], voldaynegvar_TLT[end], volweek_TLT[end], volmonth_TLT[end])  %*% matrix(coef(SHAR_TLT))

		SHAR_SPY  <- lm(RV5min_SPY[2:end] ~ voldayposvar_SPY[1:(end-1)] + voldaynegvar_SPY[1:(end-1)] + volweek_SPY[1:(end-1)] + volmonth_SPY[1:(end-1)])
		hhat_SHAR_SPY <- cbind(1, voldayposvar_SPY[end], voldaynegvar_SPY[end], volweek_SPY[end], volmonth_SPY[end])  %*% matrix(coef(SHAR_SPY))

		hhat_SHAR_SPY <- volatility.insanity.filter(hhat_SHAR_SPY, min(RV5min_SPY), max(RV5min_SPY), mean(RV5min_SPY))$vol
		hhat_SHAR_TLT <- volatility.insanity.filter(hhat_SHAR_TLT, min(RV5min_TLT), max(RV5min_TLT), mean(RV5min_TLT))$vol

		DRD_SHAR <- suppressWarnings(EstimatecorrHAR(cbind(hhat_SHAR_TLT, hhat_SHAR_SPY), correlation = correlation, proxy = proxycor, 0, T))

		#temp_SHAR[i,j] <- QLIKE(DRD_SHAR$vSigma2[,,1], RCov5min, 2)
		weights <- minvar(DRD_SHAR$vSigma2[,,1])

		temp_SHAR_mvp[i,j] <- t(weights) %*% RCov5min %*% weights

	}
	print(sprintf("%s", j))

}

mvpvariance_DRD[[16]] <- temp_SHAR_mvp[1001:nrow(temp_SHAR_mvp), ]
#Qlikes[[16]] <- temp_SHAR[1001:nrow(temp_SHAR), ]

#saveRDS(mvpvariance_DRD, "mvpvariance_DRD.rds")




#saveRDS(Qlikes, "Qlikes_Forecast_DRDmodels_intraday.rds")





means <- numeric()

for(i in 1:22){

means[i] <- mean(rowMeans(mvpvariance_DRD[[i]]))

}







#############################################################################################
#
#
#										RISKMETRICS 
# 
#
#############################################################################################
Qlikes_riskmetrics  <- list()
mvpvariance_RM <- list()
temp_variance_RM <- matrix(0L, ncol = 9, nrow = 2516)

temp_riskmetrics <- matrix(0L, ncol = 9, nrow = 2516)
for(k in 1:6){
	for(j in 1:9){
		for(i in (window+1):2516){

			choice <- c(1,4,5,6,7,8)

			Filter <- ewma.filter.realized(calccov[[choice[k]]][[j]][,,(i-window):(i-1)], NULL, F, 0.94, 0)

			end <- length(Filter[1,1, ])
			#temp_riskmetrics[i,j] <- QLIKE(Filter[,,end], calccov[[1]][[7]][,,i], 2)

			weight <- minvar(Filter[,,end])

			temp_variance_RM[i,j] <- t(weight) %*% (calccov[[1]][[7]][,,i]*10000) %*% weight
		}
		print(sprintf("%s", j))
	}
	#Qlikes_riskmetrics[[k]] <- temp_riskmetrics[1001:nrow(temp_riskmetrics), ]
	mvpvariance_RM[[k]] <- temp_variance_RM[1001:nrow(temp_variance_RM), ]
}

#saveRDS(mvpvariance, "mvpvariance_RM.rds")

ttt <- readRDS("mvpvariance_RM.rds")

means_riskmetrics <- numeric()

for(i in 1:6){

means_riskmetrics[i] <- mean(rowMeans(Qlikes_riskmetrics[[i]]))

}

#saveRDS(Qlikes_riskmetrics, "Qlikes_riskmetrics.rds")

#1second frequency for Tcov and Bpcov are huge!, thats why you get 7.48 and 1.10  in 2 and 3. 


#############################################################################################
#
#
#										rBG and rDCC
# 
#
#############################################################################################


#---------------------------------------PRELIMINARIES-----------------------------------------
dataTLT <- readRDS("dataTLT.rds")

getDates <- unlist(lapply(dataTLT, function(x) as.character(index(x[1]))))

for(i in 1:length(getDates)){
	getDates[i] <- strsplit(getDates, " ")[[i]][1]
}


dailyretotc <- xts(t(sapply(mergedfrequencies[[10]], function(x) cbind(x[,1], x[,2]))), order.by = as.Date(getDates))

colnames(dailyretotc) <- c("TLT", "SPY")

if(FALSE){
temp_mixed <- array(0L, dim = c(2,2,2516))
Mixed <- list()

for(j in 1:9){
	for(i in 1:2516){
		temp_mixed[,,i] <- realsemicov(mergedfrequencies[[j]][[i]], "M", F)
	}
	Mixed[[j]] <- temp_mixed
}
#saveRDS(Mixed, "Mixed_semicovariances.rds")
}


#--------------------------------------------------------------------------------------------


#instability issues often arises when you have chosen bad starting params, try and choose others. 



Qlikes_rBG<- list()

#you need mixed for crBG and crDCC. 



temp_rBG <- matrix(0L, ncol = 9, nrow = 2516)

window <- 1000
#rBG
rBG <- list(vPar = c(0.308691, 0.573546))
for(k in 1:6){
	for(j in 1:9){
		for(i in (window+1):2516){

			choice <- c(1,4,5,6,7,8)
			
			pars <- rBG$vPar 

			rBG <- EstimateBivarGARCH(dailyretotc[(i-window):(i-1), ] * 100, 
				calccov[[choice[k]]][[j]][,,(i-window):(i-1)]*10000, vPar = pars, tol = 1e-5)
	
			end <- 1000

			temp_rBG[i,j] <- QLIKE(rBG$vSigma2[,,end], calccov[[1]][[7]][,,i]*10000, 2)

			print(sprintf("Forecast %s, in frequency %s with measure %s",i,j,k))
		}
	}
	Qlikes_rBG[[k]] <- temp_rBG[1001:nrow(temp_rBG), ]
}


#saveRDS(Qlikes_rBG, "Qlikes_rBG.rds")


Qlikes_rDCC <- list()
temp_rDCC <- matrix(0L, ncol = 9, nrow = 2516)


#rDCC
#Tcov and pbpcov has numerically instable covariances under 1 second, so we exclude 1 sec and calculate it later
for(k in 1:6){
	for(j in 2:9){
		for(i in (window+1):2516){

			end <- 1000

			choice <- c(1,4,5,6,7,8)

			rDCC <- Estimate_rDCC(dailyretotc[(i-window):(i-1), ] * 100, calccov[[choice[k]]][[j]][,,(i-window):(i-1)]*10000, 
				getDates[(i-window):(i-1)], tol = 1e-5)

			temp_rDCC[i,j] <- QLIKE(rDCC$vSigma2[,,end], calccov[[1]][[7]][,,i]*10000, 2)
			
			print(sprintf("Forecast %s, in frequency %s with measure %s",i,j,k))
		}
	}
	Qlikes_rDCC[[k]] <- temp_rDCC
}


#one-second calculations are done in a separate setup together with daily frequencies


#saveRDS(Qlikes_rDCC, "Qlikes_rDCC.rds")



#############################################################################################
#
#
#										crBG and crDCC
# 
#
#############################################################################################
Qlikes_crBG <- matrix(0L, ncol = 9, nrow = 2516)
Qlikes_crDCC <- matrix(0L, ncol = 9, nrow = 2516)

#initiating
crBGpars <- c(0.09, 0.34, 0.04, 0.52)
crDCCpars <- c(0.0001, 0.0001, 0.08, 0.91)
#Mixed <- readRDS("Mixed_semicovariances.rds")

semicovs <- readRDS("semicov_acrossfreq.rds")

Pfreq <- semicovs[[1]]
Nfreq <- semicovs[[2]]
Mfreq <- semicovs[[3]]

window <- 1000
for(j in 1:9){
	for(i in (window+1):2516){

		decompcov <- list(Pfreq[[j]][,,(i-window):(i-1)] * 10000, 
		                  Nfreq[[j]][,,(i-window):(i-1)] * 10000, Mfreq[[j]][,,(i-window):(i-1)] * 10000)

		crBG <- EstimateBivarGARCHContAsym(dailyretotc[(i-window):(i-1), ] * 100, decompcov, vPar = crBGpars)

		crBGpars <-  crBG$vPar 
		
		end <- 1000
		
		Qlikes_crBG[i,j] <- QLIKE(crBG$vSigma2[,,end], calccov[[1]][[7]][,,i]*10000, 2)
		
		#cov2 <- calccov[[1]][[j]][,,(i-window):(i-1)] * 10000

		crDCC <- Estimate_rcDCC(dailyretotc[(i-window):(i-1), ] * 100, decompcov, NULL, getDates = getDates[(i-window):(i-1)], tol = 1e-7, vPar = crDCCpars)
		
		crDCCpars <- crDCC$vPar
		
		Qlikes_crDCC[i,j] <- QLIKE(crDCC$vSigma2[,,end], calccov[[1]][[7]][,,i]*10000, 2)
		print(sprintf("Forecast %s, in frequency %s",i,j))
		
	}
}


#saveRDS(Qlikes_crBG, "Qlikes_crBGcorrect.rds")
#saveRDS(Qlikes_crDCC, "Qlikes_crDCCcorrect.rds")



#############################################################################################
#
#
#										DAILY FREQUENCIES
# 
#
#############################################################################################


RV5min_SPY <- calccov[[1]][[7]][2,2,]*10000
RV5min_TLT <- calccov[[1]][[7]][1,1,]*10000


spy_1sec <- rollapply(calccov[[1]][[1]][2,2,]*10000, 25, mean, align = "left", fill = mean(calccov[[1]][[1]][2,2,]*10000))
tlt_1sec <- rollapply(calccov[[1]][[1]][1,1,]*10000, 25, mean, align = "left", fill = mean(calccov[[1]][[1]][1,1,]*10000))
cor_1sec <- rollapply(calccov[[1]][[1]][2,1,]*10000, 25, mean, align = "left", fill = mean(calccov[[1]][[1]][2,1,]*10000))

tempcov <- calccov[[1]][[1]]*10000

tempcov[1,1,] <- tlt_1sec
tempcov[2,2,] <- spy_1sec
tempcov[2,1,] <- cor_1sec


temptlt <- rollapply(dailyretotc[, 1]*100, 7, mean, align = "left", fill = mean(dailyretotc[, 1]*100))
tempspy <-   rollapply(dailyretotc[, 2]*100, 7, mean, align = "left", fill = mean(dailyretotc[, 2]*100))
tempretotc <- cbind(temptlt, tempspy)

