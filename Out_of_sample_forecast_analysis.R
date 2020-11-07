

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


#you need correlations for all 434 models...

window <- 1000

hhat_HAR_SPY_RV_5min <- numeric()
hhat_HAR_TLT_RV_5min <- numeric()
Qlikes <- numeric()


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

hhat_HAR_SPY_RV_5min <- numeric()
hhat_HAR_TLT_RV_5min <- numeric()
Qlikes <- list()

temp <- matrix(0L, ncol = 9, nrow = 2516)
temp_HARQ <- matrix(0L, ncol = 9, nrow = 2516)
temp_HARQF <- matrix(0L, ncol = 9, nrow = 2516)
temp_HARJ <- matrix(0L, ncol = 9, nrow = 2516)
temp_HARQJ <- matrix(0L, ncol = 9, nrow = 2516)
temp_CHAR <- matrix(0L, ncol = 9, nrow = 2516)
temp_CHARQ <- matrix(0L, ncol = 9, nrow = 2516)


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



			#RQ and TRQ
			rq_TLT <- Quarticities$rq_TLT[,j]
			rq_SPY <- Quarticities$rq_SPY[,j]

			rq_TLT <- matrix(rq_TLT)
			rq_SPY <- matrix(rq_SPY)

			trq_SPY <- Quarticities$trq_SPY[,j]
			trq_TLT <- Quarticities$trq_TLT[,j]



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




			#jumpparams
			jumpparam_TLT <- ifelse(calccov[[1]][[j]][1,1,] - calccov[[5]][[j]][1,1,]>0, 
				calccov[[1]][[j]][1,1,] - calccov[[5]][[j]][1,1,], 0)

			jumpparam_SPY <- ifelse(calccov[[1]][[j]][2,2,] - calccov[[5]][[j]][2,2,]>0, 
				calccov[[1]][[j]][2,2,] - calccov[[5]][[j]][2,2,], 0)


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
			hhat_HAR_TLT[i] <- cbind(1, volday_TLT[end], volweek_TLT[end], volmonth_TLT[end])  %*% matrix(coef(HAR_TLT))

			HAR_SPY <- lm(RV5min_SPY[2:end] ~ volday_SPY[1:(end-1)] + volweek_SPY[1:(end-1)] + volmonth_SPY[1:(end-1)])
			hhat_HAR_SPY[i] <- cbind(1, volday_SPY[end], volweek_SPY[end], volmonth_SPY[end])  %*% matrix(coef(HAR_SPY))

			#warning is just that R does not like to add one 1 parameter array to a number eg. array(1, c(1,1,1)) + 2.  
			DRD_HAR <- suppressWarnings(EstimatecorrHAR(cbind(hhat_HAR_TLT[i], hhat_HAR_SPY[i]), 
				correlation = correlation, proxy = proxycor, 0, T))

			temp[i,j] <- QLIKE(DRD_HAR$vSigma2[,,1], RCov5min, 2)



			#############################################################################################
			#
			#
			#									DRD-HARQ MODEL (RCOV, MRC, MRK)
			# 
			#
			#############################################################################################

			HARQ_TLT <- lm(RV5min_TLT[2:end] ~ volday_TLT[1:(end-1)] +  volweek_TLT[1:(end-1)] + volmonth_TLT[1:(end-1)] + I(volday_TLT[1:(end-1)] * sqrtrq_TLT[1:(end-1)]))
			hhat_HARQ_TLT <- cbind(1, volday_TLT[end], volweek_TLT[end], volmonth_TLT[end],volday_TLT[end] * sqrtrq_TLT[end])  %*% matrix(coef(HARQ_TLT))

			hhat_HARQ_TLT <- volatility.insanity.filter(hhat_HARQ_TLT, min(RV5min_TLT), max(RV5min_TLT), mean(RV5min_TLT))$vol

			HARQ_SPY <- lm(RV5min_SPY[2:end] ~ volday_SPY[1:(end-1)] +  volweek_SPY[1:(end-1)] + volmonth_SPY[1:(end-1)] + I(volday_SPY[1:(end-1)] * sqrtrq_SPY[1:(end-1)]))
			hhat_HARQ_SPY <- cbind(1, volday_SPY[end], volweek_SPY[end], volmonth_SPY[end],volday_SPY[end] * sqrtrq_SPY[end])  %*% matrix(coef(HARQ_SPY))

			hhat_HARQ_SPY <- volatility.insanity.filter(hhat_HARQ_SPY, min(RV5min_SPY), max(RV5min_SPY), mean(RV5min_SPY))$vol

			#SAME CORRELATION AS HAR MODELS. ONLY THING THAT CHANGED IS THE UNIVARIATE VOLS. 
			DRD_HARQ <- suppressWarnings(EstimatecorrHAR(cbind(hhat_HARQ_TLT, hhat_HARQ_SPY), correlation = correlation, proxy = proxycor, 0,T))

			temp_HARQ[i,j] <- QLIKE(DRD_HARQ$vSigma2[,,1], RCov5min, 2)


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

			temp_HARQF[i,j] <- QLIKE(DRD_HARQF$vSigma2[,,1], RCov5min, 2)

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

			DRD_HARJ <- suppressWarnings(EstimatecorrHAR(cbind(hhat_HARJ_TLT, hhat_HARJ_SPY), correlation = correlation, proxy = proxycor, 0, T))
			
			temp_HARJ[i,j] <- QLIKE(DRD_HARJ$vSigma2[,,1], RCov5min, 2)

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

			temp_HARQJ[i,j] <- QLIKE(DRD_HARQJ$vSigma2[,,1], RCov5min, 2)

			#############################################################################################
			#
			#
			#									DRD-CHAR MODEL (TCOV, BPCOV, PBPCOV)
			# 
			#
			#############################################################################################

			CHAR_TLT <- lm(RV5min_TLT[2:end] ~ volday_TLT[1:(end-1)] + volweek_TLT[1:(end-1)] + volmonth_TLT[1:(end-1)])
			hhat_CHAR_TLT <- cbind(1, volday_TLT[end], volweek_TLT[end], volmonth_TLT[end])  %*% matrix(coef(CHAR_TLT))

			CHAR_SPY <- lm(RV5min_SPY[2:end] ~ volday_SPY[1:(end-1)] + volweek_SPY[1:(end-1)] + volmonth_SPY[1:(end-1)])
			hhat_CHAR_SPY <- cbind(1, volday_SPY[end], volweek_SPY[end], volmonth_SPY[end])  %*% matrix(coef(CHAR_SPY))

			DRD_CHAR <- suppressWarnings(EstimatecorrHAR(cbind(hhat_CHAR_TLT, hhat_CHAR_SPY), correlation = correlationjumprobust, proxy = proxycor, 0, T))

			temp_CHAR[i,j] <- QLIKE(DRD_CHAR$vSigma2[,,1], RCov5min, 2)


			#############################################################################################
			#
			#
			#									DRD-CHARQ MODEL (TCOV, BPCOV, PBPCOV)
			# 
			#
			#############################################################################################

			CHARQ_TLT <- lm(RV5min_TLT[2:end] ~ volday_TLT[1:(end-1)] +  volweek_TLT[1:(end-1)] + volmonth_TLT[1:(end-1)] + I(volday_TLT[1:(end-1)] * sqrttrq_TLT[1:(end-1)]))
			hhat_CHARQ_TLT <- cbind(1, volday_TLT[end], volweek_TLT[end], volmonth_TLT[end], volday_TLT[end] * sqrttrq_TLT[end])  %*% matrix(coef(CHARQ_TLT))

			CHARQ_SPY <- lm(RV5min_SPY[2:end] ~ volday_SPY[1:(end-1)] +  volweek_SPY[1:(end-1)] + volmonth_SPY[1:(end-1)] + I(volday_SPY[1:(end-1)] * sqrttrq_SPY[1:(end-1)]))
			hhat_CHARQ_SPY <- cbind(1, volday_SPY[end], volweek_SPY[end], volmonth_SPY[end], volday_SPY[end] * sqrttrq_SPY[end])  %*% matrix(coef(CHARQ_SPY))

			hhat_CHARQ_SPY <- volatility.insanity.filter(hhat_CHARQ_SPY, min(RV5min_SPY), max(RV5min_SPY), mean(RV5min_SPY))$vol
			hhat_CHARQ_TLT <- volatility.insanity.filter(hhat_CHARQ_TLT, min(RV5min_TLT), max(RV5min_TLT), mean(RV5min_TLT))$vol


			DRD_CHARQ <- suppressWarnings(EstimatecorrHAR(cbind(hhat_CHARQ_TLT, hhat_CHARQ_SPY), correlation = correlationjumprobust, proxy = proxycor, 0, T))

			temp_CHARQ[i,j] <- QLIKE(DRD_CHARQ$vSigma2[,,1], RCov5min, 2)

			#############################################################################################
			#
			#
			#											DRD-SHAR MODEL 
			# 
			#
			#############################################################################################
			

		}
		print(sprintf("%s", j))
	}
	Qlikes[[k]] <- temp[1001:nrow(temp), ]
	Qlikes[[k+3]] <- temp_HARQ[1001:nrow(temp_HARQ), ]
	Qlikes[[k+6]] <- temp_HARQF[1001:nrow(temp_HARQF), ]
	Qlikes[[k+9]] <- temp_HARJ[1001:nrow(temp_HARJ), ]
	Qlikes[[k+12]] <- temp_HARQJ[1001:nrow(temp_HARQJ), ]
	#putting SHAR model at 16th list item 
	Qlikes[[k+16]] <- temp_CHAR[1001:nrow(temp_CHAR), ]
	Qlikes[[k+19]] <- temp_CHARQ[1001:nrow(temp_CHAR), ]


}


for(j in 1:9){
		for(i in (window+1):1050){



			#############################################################################################
			#
			#
			#											DRD-SHAR MODEL 
			# 
			#
			#############################################################################################


		}

	}