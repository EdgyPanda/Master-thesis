##########out of sample volatility models###########



###################################################################################
#
#
#                             OUT OF SAMPLE QLIKE LOSSES
#
#
###################################################################################

QLikes_DRD <- readRDS("Qlikes_Forecast_DRDmodels_intraday.rds")

QLikes_RM <- readRDS("Qlikes_riskmetrics.rds")

QLikes_rBG <- readRDS("Qlikes_rBG.rds")

QLikes_crBG <- readRDS("Qlikes_crBG.rds")

QLikes_crDCC <- readRDS("Qlikes_crDCC.rds")

QLikes_rDCC <- readRDS("Qlikes_rDCC.rds")

colmeanQLikes_DRD <- lapply(QLikes_DRD, FUN = function(x) colMeans(x)) 
colmeanQLikes_RM <- lapply(QLikes_RM, FUN = function(x) colMeans(x)) 

colmeanQLikes_rBG <- lapply(QLikes_rBG, FUN = function(x) colMeans(x)) 

colmeanQLikes_crBG <- colMeans(QLikes_crBG, na.rm = T)
colmeanQLikes_crDCC <- colMeans(QLikes_crDCC, na.rm = T)

colmeanQLikes_rDCC <- lapply(QLikes_rDCC, FUN = function(x) if(!is.null(x)){colMeans(x, na.rm =T)}) 

#Minute frequency:
DRD_HAR_Rcov_min <- colmeanQLikes_DRD[[1]][7]
               
DRD_all_minute <- lapply(colmeanQLikes_DRD, FUN = function(x) mean(x[6:9])/DRD_HAR_Rcov_min)

RM_all_minute <- lapply(colmeanQLikes_RM, FUN = function(x) mean(x[6:9])/DRD_HAR_Rcov_min)

rBG_all_minute <- lapply(colmeanQLikes_rBG, FUN = function(x) mean(x[6:9])/DRD_HAR_Rcov_min)

rBG_all_minute_2 <- lapply(colmeanQLikes_rBG, FUN = function(x) mean(x[6:7])/DRD_HAR_Rcov_min)


rBG_all_minute[[4]] <- rBG_all_minute_2[[4]]
rBG_all_minute[[5]] <- rBG_all_minute_2[[5]]
rBG_all_minute[[6]] <- rBG_all_minute_2[[6]]

rBG_all_minute

crBG_all_minute <- mean(colmeanQLikes_crBG[6:9]/DRD_HAR_Rcov_min)

crDCC_all_minute <- mean(colmeanQLikes_crDCC[6:9]/DRD_HAR_Rcov_min)

rDCC_all_minute <- lapply(colmeanQLikes_rDCC, FUN = function(x) if(!is.null(x)){mean(x[6:9], na.rm =T)/DRD_HAR_Rcov_min}) 




#second frequency
DRD_all_second <- lapply(colmeanQLikes_DRD, FUN = function(x) mean(x[1:5])/DRD_HAR_Rcov_min)

RM_all_second <- lapply(colmeanQLikes_RM, FUN = function(x) mean(x[1:5])/DRD_HAR_Rcov_min)

rBG_all_second <- lapply(colmeanQLikes_rBG, FUN = function(x) mean(x[1:5])/DRD_HAR_Rcov_min)

crBG_all_second <- mean(colmeanQLikes_crBG[1:5]/DRD_HAR_Rcov_min)

crDCC_all_second <- mean(colmeanQLikes_crDCC[1:5]/DRD_HAR_Rcov_min)

#have zeros in first column
rDCC_all_second <- lapply(colmeanQLikes_rDCC, FUN = function(x) if(!is.null(x)){mean(x[2:5], na.rm =T)/DRD_HAR_Rcov_min}) 

###################################################################################
#
#
#                            OUT-OF-SAMPLE STRATIFIED LOSSES 
#
#
###################################################################################


######## GETTING THE ERROR VARIANCES FOR THE HARQ MODEL.
########## JUST DO AS IN THE ORIGINAL ARTICLE AND USE RQ INSTEAD OR ERROR VARIANCES
EV_TLT <- list()
errorvariance_TLT <- matrix(0L, ncol = 9, nrow = 2516)

EV_SPY <- list()
errorvariance_SPY <- matrix(0L, ncol = 9, nrow = 2516)
for(k in 1:3){
	for(j in 1:9){
		for(i in (window+1):2516){


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




			rq_TLT <- Quarticities$rq_TLT[(i-window):(i-1),j]
			rq_SPY <- Quarticities$rq_SPY[(i-window):(i-1),j]
			rq_SPY <- rq_SPY[22:length(rq_SPY)]
			rq_TLT <- rq_TLT[22:length(rq_TLT)]
			rq_TLT[rq_TLT>10] <- mean(rq_TLT)
			rq_SPY[rq_SPY>10] <- mean(rq_SPY)


			rq_TLT <- matrix(rq_TLT)

			rq_SPY <- matrix(rq_SPY)



			sqrtrq_TLT <- sqrt(rq_TLT) - mean(sqrt(rq_TLT))

			sqrtrq_SPY <- sqrt(rq_SPY) - mean(sqrt(rq_SPY))


			correlation <- Correlations_measures[[j]][(i-window):(i-1), QCchoice[k]]

			
			
			#############################################################################################
			#
			#
			#									DRD-HARQ MODEL (RCOV, MRC, MRK)
			# 
			#
			#############################################################################################

			HARQ_TLT <- lm(RV5min_TLT[2:end] ~ volday_TLT[1:(end-1)] +  volweek_TLT[1:(end-1)] + volmonth_TLT[1:(end-1)] + I(volday_TLT[1:(end-1)] * sqrtrq_TLT[1:(end-1)]))
			hhat_HARQ_TLT <- cbind(1, volday_TLT[end], volweek_TLT[end], volmonth_TLT[end],volday_TLT[end] * sqrtrq_TLT[end])  %*% matrix(coef(HARQ_TLT))


			HARQ_SPY <- lm(RV5min_SPY[2:end] ~ volday_SPY[1:(end-1)] +  volweek_SPY[1:(end-1)] + volmonth_SPY[1:(end-1)] 
				+ I(volday_SPY[1:(end-1)] * sqrtrq_SPY[1:(end-1)]))
			hhat_HARQ_SPY <- cbind(1, volday_SPY[end], volweek_SPY[end], volmonth_SPY[end],volday_SPY[end] * sqrtrq_SPY[end])  %*% matrix(coef(HARQ_SPY))

			hhat_HARQ_TLT <- volatility.insanity.filter(hhat_HARQ_TLT, min(RV5min_TLT), max(RV5min_TLT), mean(RV5min_TLT))$vol
			hhat_HARQ_SPY <- volatility.insanity.filter(hhat_HARQ_SPY, min(RV5min_SPY), max(RV5min_SPY), mean(RV5min_SPY))$vol

			errorvariance_TLT[i,j] <- coef(HARQ_TLT)[2]  + coef(HARQ_TLT)[5] * sqrtrq_TLT[end]
			errorvariance_SPY[i,j] <- coef(HARQ_SPY)[2]  + coef(HARQ_SPY)[5] * sqrtrq_SPY[end]

			
		}
		print(sprintf("%s", j))
	}
	EV_TLT[[k]] <- errorvariance_TLT[1001:2516, ]
	EV_SPY[[k]] <- errorvariance_SPY[1001:2516, ]

}

frob.norm <- function(Covar){

	frob <- sqrt(sum(diag(t(Covar)%*%Covar)))

	return(frob)

}

Quarticities <- readRDS("Quarticity_estimators.rds")

#Remember its error variances and not only quarticities. So get the params of your model. 
qmat <- array(0L, dim = c(2,2,2516))
trqmat <- array(0L, dim = c(2,2,2516))
Qmats <- list()
TRQmats <- list()

for(j in 1:9){
	for(i in 1:2516){
				tltq <- Quarticities$rq_TLT[i,j]
				tltq[tltq>30] <- mean(Quarticities$rq_TLT[,j])
				spyq <- Quarticities$rq_SPY[i,j]
				spyq[spyq>30] <- mean(Quarticities$rq_SPY[,j])
				tlttrq <- Quarticities$trq_TLT[i,j]
				tlttrq[tlttrq>30] <- mean(Quarticities$trq_TLT[,j])
				spytrq <- Quarticities$trq_SPY[i,j]
				spytrq[spytrq>30] <- mean(Quarticities$trq_SPY[,j])
		tempq <- matrix(c(tltq, 0,0, spyq), ncol = 2, nrow = 2)
		temptrq <- matrix(c(tlttrq, 0, 0, spytrq), nrow = 2, ncol = 2)
		qmat[,,i] <- tempq
		trqmat[,,i] <- temptrq
	}
	Qmats[[j]] <- qmat
	TRQmats[[j]] <- trqmat
}


frobloss.Qmats <- list()
frobs <- numeric()
frobloss.TRQmats <- list()
frobs2 <- numeric()

for(j in 1:9){
    for(i in 1:2516){
    	frobs[i] <- frob.norm(Qmats[[j]][,,i])
    	frobs2[i] <- frob.norm(TRQmats[[j]][,,i])
    }
    frobloss.Qmats[[j]] <- frobs
    frobloss.TRQmats[[j]] <- frobs2
}

frobloss.Qmats <- lapply(frobloss.Qmats, function(x) x[1001:2516])
frobloss.TRQmats <- lapply(frobloss.TRQmats, function(x) x[1001:2516])



findHighQs <- lapply(frobloss.Qmats, function(x) which(x > quantile(x, 0.95)))
findLowQs <- lapply(frobloss.Qmats, function(x) which(x <= quantile(x, 0.95)))
findHighTRQs <- lapply(frobloss.TRQmats, function(x) which(x > quantile(x, 0.95)))
findLowTRQs <- lapply(frobloss.TRQmats, function(x) which(x <= quantile(x, 0.95)))



#Finding Qlikes for the days of high and low quarticities.

 QLikes_DRD_highQs <- list()
 tempq <- matrix(0L, ncol = 9, nrow = 1516)
  QLikes_DRD_lowQs <- list()
 tempq2 <- matrix(0L, ncol = 9, nrow = 1516)
 QLikes_DRD_highTRQs <- list()
 temptrq <- matrix(0L, ncol = 9, nrow = 1516)
 QLikes_DRD_lowTRQs <- list()
 temptrq2 <- matrix(0L, ncol = 9, nrow = 1516)

 #RM 
 QLikes_RM_highQs <- list()
 tempqRM <- matrix(0L, ncol = 9, nrow = 1516)
  QLikes_RM_lowQs <- list()
 tempq2RM <- matrix(0L, ncol = 9, nrow = 1516)
 QLikes_RM_highTRQs <- list()
 temptrqRM <- matrix(0L, ncol = 9, nrow = 1516)
 QLikes_RM_lowTRQs <- list()
 temptrq2RM <- matrix(0L, ncol = 9, nrow = 1516)

 #rBG 
  QLikes_rBG_highQs <- list()
 tempqrBG <- matrix(0L, ncol = 9, nrow = 1516)
  QLikes_rBG_lowQs <- list()
 tempq2rBG <- matrix(0L, ncol = 9, nrow = 1516)
 QLikes_rBG_highTRQs <- list()
 temptrqrBG <- matrix(0L, ncol = 9, nrow = 1516)
 QLikes_rBG_lowTRQs <- list()
 temptrq2rBG <- matrix(0L, ncol = 9, nrow = 1516)

  #rDCC 
  QLikes_rDCC_highQs <- list()
 tempqrDCC <- matrix(0L, ncol = 9, nrow = 1516)
  QLikes_rDCC_lowQs <- list()
 tempq2rDCC <- matrix(0L, ncol = 9, nrow = 1516)
 QLikes_rDCC_highTRQs <- list()
 temptrqrDCC <- matrix(0L, ncol = 9, nrow = 1516)
 QLikes_rDCC_lowTRQs <- list()
 temptrq2rDCC <- matrix(0L, ncol = 9, nrow = 1516)

   #crDCC 
 tempqcrDCC <- matrix(0L, ncol = 9, nrow = 76)
 tempq2crDCC <- matrix(0L, ncol = 9, nrow = 1440)
 temptrqcrDCC <- matrix(0L, ncol = 9, nrow = 76)
 temptrq2crDCC <- matrix(0L, ncol = 9, nrow = 1440)

 #crBG 
 tempqcrBG <- matrix(0L, ncol = 9, nrow = 76)
 tempq2crBG <- matrix(0L, ncol = 9, nrow = 1440)
 temptrqcrBG <- matrix(0L, ncol = 9, nrow = 76)
 temptrq2crBG <- matrix(0L, ncol = 9, nrow = 1440)


for(k in 1:22){
	for(j in 1:9){
		tempq[1:76,j] <- QLikes_DRD[[k]][findHighQs[[j]],j]
		tempq2[1:1440,j] <- QLikes_DRD[[k]][findLowQs[[j]],j]
		temptrq[1:76,j] <- QLikes_DRD[[k]][findHighTRQs[[j]],j]
		temptrq2[1:1440,j] <- QLikes_DRD[[k]][findLowTRQs[[j]],j]
	}
	QLikes_DRD_highQs[[k]] <- tempq[1:76, ]
	QLikes_DRD_lowQs[[k]] <- tempq2[1:1440, ]
	QLikes_DRD_highTRQs[[k]] <- temptrq[1:76, ]
	QLikes_DRD_lowTRQs[[k]] <- temptrq2[1:1440, ]

}
#RM
for(k in 1:6){
	for(j in 1:9){
		tempqRM[1:76,j] <- QLikes_RM[[k]][findHighQs[[j]],j]
		tempq2RM[1:1440,j] <- QLikes_RM[[k]][findLowQs[[j]],j]
		temptrqRM[1:76,j] <- QLikes_RM[[k]][findHighTRQs[[j]],j]
		temptrq2RM[1:1440,j] <- QLikes_RM[[k]][findLowTRQs[[j]],j]
	}
	QLikes_RM_highQs[[k]] <- tempqRM[1:76, ]
	QLikes_RM_lowQs[[k]] <- tempq2RM[1:1440, ]
	QLikes_RM_highTRQs[[k]] <- temptrqRM[1:76, ]
	QLikes_RM_lowTRQs[[k]] <- temptrq2RM[1:1440, ]

	
}
#rBG
#for k > 3 j >7 is invalid. 
for(k in 1:6){
	if(k<=3){
		for(j in 1:9){
			tempqrBG[1:76,j] <- QLikes_rBG[[k]][findHighQs[[j]],j]
			tempq2rBG[1:1440,j] <- QLikes_rBG[[k]][findLowQs[[j]],j]
			temptrqrBG[1:76,j] <- QLikes_rBG[[k]][findHighTRQs[[j]],j]
			temptrq2rBG[1:1440,j] <- QLikes_rBG[[k]][findLowTRQs[[j]],j]
		}
	}
	else{
		for(j in 1:7){
			tempqrBG[1:76,j] <- QLikes_rBG[[k]][findHighQs[[j]],j]
			tempq2rBG[1:1440,j] <- QLikes_rBG[[k]][findLowQs[[j]],j]
			temptrqrBG[1:76,j] <- QLikes_rBG[[k]][findHighTRQs[[j]],j]
			temptrq2rBG[1:1440,j] <- QLikes_rBG[[k]][findLowTRQs[[j]],j]
		}
	}
	QLikes_rBG_highQs[[k]] <- tempqrBG[1:76, ]
	QLikes_rBG_lowQs[[k]] <- tempq2rBG[1:1440, ]
	QLikes_rBG_highTRQs[[k]] <- temptrqrBG[1:76, ]
	QLikes_rBG_lowTRQs[[k]] <- temptrq2rBG[1:1440, ]
}

#rDCC
for(k in 1:2){
	for(j in 1:9){
		tempqrDCC[1:76,j] <- QLikes_rDCC[[k]][findHighQs[[j]],j]
		tempq2rDCC[1:1440,j] <- QLikes_rDCC[[k]][findLowQs[[j]],j]
		temptrqrDCC[1:76,j] <- QLikes_rDCC[[k]][findHighTRQs[[j]],j]
		temptrq2rDCC[1:1440,j] <- QLikes_rDCC[[k]][findLowTRQs[[j]],j]
	}
	QLikes_rDCC_highQs[[k]] <- tempqrDCC[1:76, ]
	QLikes_rDCC_lowQs[[k]] <- tempq2rDCC[1:1440, ]
	QLikes_rDCC_highTRQs[[k]] <- temptrqrDCC[1:76, ]
	QLikes_rDCC_lowTRQs[[k]] <- temptrq2rDCC[1:1440, ]
}
#only goes to 1001 ......
for(k in 5:6){
	for(j in 1:7){
		tempqrDCC[1:length(which(findHighQs[[j]]<1001)),j] <- QLikes_rDCC[[k]][     findHighQs[[j]][which(findHighQs[[j]]<1001)]  ,j]
		tempq2rDCC[1:length(which(findLowQs[[j]]<1001)),j] <- QLikes_rDCC[[k]][findLowQs[[j]][which(findLowQs[[j]]<1001)],j]
		temptrqrDCC[1:length(which(findHighTRQs[[j]]<1001)),j] <- QLikes_rDCC[[k]][findHighTRQs[[j]][which(findHighTRQs[[j]]<1001)],j]
		temptrq2rDCC[1:length(which(findLowTRQs[[j]]<1001)),j] <- QLikes_rDCC[[k]][findLowTRQs[[j]][which(findLowTRQs[[j]]<1001)],j]

	}
	QLikes_rDCC_highQs[[k]] <- tempqrDCC
	QLikes_rDCC_lowQs[[k]] <- tempq2rDCC
	QLikes_rDCC_highTRQs[[k]] <- temptrqrDCC
	QLikes_rDCC_lowTRQs[[k]] <- temptrq2rDCC
}

for(j in 1:9){
	tempqcrDCC[1:76,j] <- QLikes_crDCC[findHighQs[[j]],j]
	tempq2crDCC[1:1440,j] <- QLikes_crDCC[findLowQs[[j]],j]
	temptrqcrDCC[1:76,j] <- QLikes_crDCC[findHighTRQs[[j]],j]
	temptrq2crDCC[1:1440,j] <- QLikes_crDCC[findLowTRQs[[j]],j]
	tempqcrBG[1:76,j] <- QLikes_crBG[findHighQs[[j]],j]
	tempq2crBG[1:1440,j] <- QLikes_crBG[findLowQs[[j]],j]
	temptrqcrBG[1:76,j] <- QLikes_crBG[findHighTRQs[[j]],j]
	temptrq2crBG[1:1440,j] <- QLikes_crBG[findLowTRQs[[j]],j]



}

#################### HIGH RQ ###############
library(rugarch)

DRD_highQs <- do.call(cbind, QLikes_DRD_highQs)
RM_highQs <- do.call(cbind, QLikes_RM_highQs)
rBG_highQs <- do.call(cbind, QLikes_rBG_highQs[1:3])
rBG_highQs <- cbind(rBG_highQs, QLikes_rBG_highQs[[4]][,1:7], QLikes_rBG_highQs[[5]][,1:7], QLikes_rBG_highQs[[6]][,1:7])

#rDCC_highQs dont care, just write qlikes down but no mcs
crBG_highQs <- tempqcrBG[,1:8]
crDCC_highQs <- tempqcrDCC[,1:8]

highQs_all <- cbind(DRD_highQs, RM_highQs, rBG_highQs, crBG_highQs, crDCC_highQs)

mcs_highQs <- mcsTest(highQs_all, 0.05, nboot = 1000, nblock = 10, boot = c("stationary"))

sort(mcs_highQs$includedR)


duplicated(t(highQs_all))

#saveRDS(mcs_highQs, "mcs_highQs.rds")

ncol(DRD_highQs) + ncol(RM_highQs)


#minutes

DRD_highQs_minutes <- lapply(QLikes_DRD_highQs, function(x) mean(colMeans(x[,6:9])))

DRD_highs_minutes_HAR <- mean(unlist(DRD_highQs_minutes[1:3]))
DRD_highs_minutes_SHAR <- mean(unlist(DRD_highQs_minutes[16]))

DRD_highs_minutes_HARQ <- mean(unlist(DRD_highQs_minutes[4:6]))
DRD_highs_minutes_CHAR <- mean(unlist(DRD_highQs_minutes[17:19]))
DRD_highs_minutes_CHARQ <- mean(unlist(DRD_highQs_minutes[20:22]))

#seconds 
DRD_highQs_seconds <- lapply(QLikes_DRD_highQs, function(x) mean(colMeans(x[,1:5])))
DRD_highs_seconds_HAR <- mean(unlist(DRD_highQs_seconds[1:3]))
DRD_highs_seconds_SHAR <- mean(unlist(DRD_highQs_seconds[16]))

DRD_highs_seconds_HARQ <- mean(unlist(DRD_highQs_seconds[4:6]))
DRD_highs_seconds_CHAR <- mean(unlist(DRD_highQs_seconds[17:19]))
DRD_highs_seconds_CHARQ <- mean(unlist(DRD_highQs_seconds[20:22]))


#RM

RM_highQs_minutes <- lapply(QLikes_RM_highQs, function(x) mean(colMeans(x[,6:9])))
RM_highQs_minutes <- mean(unlist(RM_highQs_minutes))


RM_highQs_seconds <- lapply(QLikes_RM_highQs, function(x) mean(colMeans(x[,1:5])))
RM_highQs_seconds <- mean(unlist(RM_highQs_seconds)[c(1,4,5,6)])


#rBG 
rBG_minutes  <- lapply(QLikes_rBG_highQs, function(x) colMeans(x[, 6:9]))

rBG_minutes <- mean(unlist(rBG_minutes))

rBG_seconds  <- lapply(QLikes_rBG_highQs, function(x) colMeans(x[, 1:5]))

rBG_seconds <- mean(unlist(rBG_seconds))

#crBG
mean(colMeans(crBG_highQs)[6:8])
mean(colMeans(crBG_highQs)[1:5])


#crDCC 
mean(colMeans(crDCC_highQs)[6:8])
mean(colMeans(crDCC_highQs)[1:5])

#rDCC
mean(c(colMeans(QLikes_rDCC_highQs[[1]])[6:9], colMeans(QLikes_rDCC_highQs[[2]])[6:9],
colMeans(QLikes_rDCC_highQs[[5]][1:76, 6:7]),colMeans(QLikes_rDCC_highQs[[6]][1:76, 6:7])))


mean(c(colMeans(QLikes_rDCC_highQs[[1]])[2:5],colMeans(QLikes_rDCC_highQs[[2]])[2:5], 
mean(QLikes_rDCC_highQs[[5]][1:45, 1]),  mean(QLikes_rDCC_highQs[[6]][1:45, 1]), colMeans(QLikes_rDCC_highQs[[5]][1:76, 2:5]), 
colMeans(QLikes_rDCC_highQs[[6]][1:76, 1:5])))

colMeans(DRD_highQs)
colMeans(RM_highQs)
colMeans(rBG_highQs)



colMeans(QLikes_rDCC_highQs[[1]])
colMeans(QLikes_rDCC_highQs[[2]])
mean(QLikes_rDCC_highQs[[5]][1:45, 1])
colMeans(QLikes_rDCC_highQs[[5]][1:76, ])
mean(QLikes_rDCC_highQs[[6]][1:45, 1])
colMeans(QLikes_rDCC_highQs[[6]][1:76, ])

colMeans(crBG_highQs)
colMeans(crDCC_highQs)

################# HIGH TRQ ##################
DRD_highTRQs <- do.call(cbind, QLikes_DRD_highTRQs)
RM_highTRQs <- do.call(cbind, QLikes_RM_highTRQs)
rBG_highTRQs <- do.call(cbind, QLikes_rBG_highTRQs[1:3])
rBG_highTRQs <- cbind(rBG_highTRQs, QLikes_rBG_highTRQs[[4]][,1:7], QLikes_rBG_highTRQs[[5]][,1:7], QLikes_rBG_highTRQs[[6]][,1:7])

crBG_highTRQs <- temptrqcrBG
crDCC_highTRQs <- temptrqcrDCC
highTRQs_all <- cbind(DRD_highTRQs, RM_highTRQs, rBG_highTRQs, crBG_highTRQs, crDCC_highTRQs)

mcs_highTRQs <- mcsTest(highTRQs_all, 0.05, nboot = 1000, nblock = 10, boot = c("stationary"))

colMeans(DRD_highTRQs)
colMeans(RM_highTRQs)
colMeans(rBG_highTRQs)

colMeans(QLikes_rDCC_highTRQs[[1]])
colMeans(QLikes_rDCC_highTRQs[[2]])
mean(QLikes_rDCC_highTRQs[[5]][1:45, 1])
colMeans(QLikes_rDCC_highTRQs[[5]][1:76, ])
mean(QLikes_rDCC_highTRQs[[6]][1:45, 1])
colMeans(QLikes_rDCC_highTRQs[[6]][1:76, ])

colMeans(crBG_highTRQs)
colMeans(crDCC_highTRQs)


#################### LOW RQ ###############
DRD_lowQs <- do.call(cbind, QLikes_DRD_lowQs)
RM_lowQs <- do.call(cbind, QLikes_RM_lowQs)
rBG_lowQs <- do.call(cbind, QLikes_rBG_lowQs[1:3])
rBG_lowQs <- cbind(rBG_lowQs, QLikes_rBG_lowQs[[4]][,1:7], QLikes_rBG_lowQs[[5]][,1:7], QLikes_rBG_lowQs[[6]][,1:7])

#rDCC_lowQs dont care, just write qlikes down but no mcs
crBG_lowQs <- tempq2crBG[,1:8]
crDCC_lowQs <- tempq2crDCC[,1:8]

lowQs_all <- cbind(DRD_lowQs, RM_lowQs, rBG_lowQs, crBG_lowQs, crDCC_lowQs)

duplicated(t(lowQs_all))

mcs_lowQs <- mcsTest(lowQs_all, 0.05, nboot = 2000, nblock = 10, boot = c("stationary"))

sort(mcs_lowQs$includedR)

#saveRDS(mcs_lowQs, "mcs_lowQs.rds")



#minutes

DRD_lowQs_minutes <- lapply(QLikes_DRD_lowQs, function(x) mean(colMeans(x[,6:9])))

DRD_lows_minutes_HAR <- mean(unlist(DRD_lowQs_minutes[1:3]))
DRD_lows_minutes_SHAR <- mean(unlist(DRD_lowQs_minutes[16]))
DRD_lows_minutes_HARQ <- mean(unlist(DRD_lowQs_minutes[4:6]))
DRD_lows_minutes_CHAR <- mean(unlist(DRD_lowQs_minutes[17:19]))
DRD_lows_minutes_CHARQ <- mean(unlist(DRD_lowQs_minutes[20:22]))

#seconds 
DRD_lowQs_seconds <- lapply(QLikes_DRD_lowQs, function(x) mean(colMeans(x[,1:5])))
DRD_lows_seconds_HAR <- mean(unlist(DRD_lowQs_seconds[1:3]))
DRD_lows_seconds_SHAR <- mean(unlist(DRD_lowQs_seconds[16]))
DRD_lows_seconds_HARQ <- mean(unlist(DRD_lowQs_seconds[4:6]))
DRD_lows_seconds_CHAR <- mean(unlist(DRD_lowQs_seconds[17:19]))
DRD_lows_seconds_CHARQ <- mean(unlist(DRD_lowQs_seconds[20:22]))

#RM

RM_lowQs_minutes <- lapply(QLikes_RM_lowQs, function(x) mean(colMeans(x[,6:9])))
RM_lowQs_minutes <- mean(unlist(RM_lowQs_minutes)[c(1,4,5,6)])


RM_lowQs_seconds <- lapply(QLikes_RM_lowQs, function(x) mean(colMeans(x[,1:5])))
RM_lowQs_seconds <- mean(unlist(RM_lowQs_seconds)[c(1,4,5,6)])


#rBG 
rBG_minutes  <- lapply(QLikes_rBG_lowQs, function(x) colMeans(x[, 6:9]))

rBG_minutes <- mean(unlist(rBG_minutes))

rBG_seconds  <- lapply(QLikes_rBG_lowQs, function(x) colMeans(x[, 1:5]))

rBG_seconds <- mean(unlist(rBG_seconds))

#crBG
mean(colMeans(crBG_lowQs)[6:8])
mean(colMeans(crBG_lowQs)[1:5])


#crDCC 
mean(colMeans(crDCC_lowQs)[6:8])
mean(colMeans(crDCC_lowQs)[1:5])

#rDCC
mean(c(colMeans(QLikes_rDCC_lowQs[[1]])[6:9], colMeans(QLikes_rDCC_lowQs[[2]])[6:9],
colMeans(QLikes_rDCC_lowQs[[5]][1:76, 6:7]),colMeans(QLikes_rDCC_lowQs[[6]][1:76, 6:7])))


mean(c(colMeans(QLikes_rDCC_lowQs[[1]])[2:5],colMeans(QLikes_rDCC_lowQs[[2]])[2:5], 
mean(QLikes_rDCC_lowQs[[5]][1:45, 1]),  mean(QLikes_rDCC_lowQs[[6]][1:45, 1]), colMeans(QLikes_rDCC_lowQs[[5]][1:76, 2:5]), 
colMeans(QLikes_rDCC_lowQs[[6]][1:76, 1:5])))







################# LOW TRQ ##################
DRD_lowTRQs <- do.call(cbind, QLikes_DRD_lowTRQs)
RM_lowTRQs <- do.call(cbind, QLikes_RM_lowTRQs)
rBG_lowTRQs <- do.call(cbind, QLikes_rBG_lowTRQs)

crBG_lowTRQs <- temptrq2crBG
crDCC_lowTRQs <- temptrq2crDCC
lowTRQs_all <- cbind(DRD_lowTRQs, RM_lowTRQs, rBG_lowTRQs, crBG_lowTRQs, crDCC_lowTRQs)

mcs_lowTRQs <- mcsTest(lowTRQs_all, 0.05, nboot = 1000, nblock = 10, boot = c("stationary"))






colMeans(DRD_lowTRQs)
colMeans(RM_lowTRQs)
colMeans(rBG_lowTRQs)

colMeans(QLikes_rDCC_lowTRQs[[1]])
colMeans(QLikes_rDCC_lowTRQs[[2]])
mean(QLikes_rDCC_lowTRQs[[5]][1:45, 1])
colMeans(QLikes_rDCC_lowTRQs[[5]][1:76, ])
mean(QLikes_rDCC_lowTRQs[[6]][1:45, 1])
colMeans(QLikes_rDCC_lowTRQs[[6]][1:76, ])

colMeans(crBG_lowQs)
colMeans(crDCC_lowQs)



###################################################################################
#
#
#                            		MCS PROCEDURE 
#
#
###################################################################################

#replacing NA with unconditional mean
rdcc6matrix <- matrix(colMeans(QLikes_rDCC[[6]]), ncol = 7, nrow = 515, byrow = T)
rdcc5matrix <- matrix(colMeans(QLikes_rDCC[[5]]), ncol = 7, nrow = 515, byrow = T)


temprdcc <- list(QLikes_rDCC[[1]], QLikes_rDCC[[2]][,-c(1:2)], rbind(QLikes_rDCC[[5]],rdcc5matrix), rbind(QLikes_rDCC[[6]],rdcc6matrix))


library(rugarch)

crBG_models <- QLikes_crBG

crBG_models[836:1516, 9] <- mean(crBG_models[,9], na.rm = T)

crDCC_models <- QLikes_crDCC

crDCC_models[835:1516, 9] <- mean(crDCC_models[,9], na.rm = T)


#constructing loss matrix. 

DRD_models <- do.call(cbind, QLikes_DRD)
RM_models <- do.call(cbind, QLikes_RM)
rBG_models <- do.call(cbind, QLikes_rBG)
rDCC_models <- do.call(cbind, temprdcc)[, -c(1,10)]


loss_matrix <- cbind(DRD_models, RM_models, rBG_models, rDCC_models, crBG_models, crDCC_models)

#TO BE USED WITHIN MCS PROCEDURE OF SHEPPARD!
#write.table(loss_matrix, file="lossmatrix.csv") 

ncol(loss_matrix)

#skipped: DRD = daily (8), RM = daily (2), rBG = daily (2), rBG = misses last two rows of PBPCov, MRC and MRK (6)
#crBG = daily (1), crDCC = daily (1), rDCC = daily (2), misses BPCov and PBPCov (18) and 2 rows in TCov,MRC,MRK (6) and
#missing 1 and 10 (2). 

#so the total on the models: 346 + 8 + 2 + 2 + 6 + 1 + 1+ 2 + 18 + 6 + 2 = 394 (CORRECT)! 


ncol(DRD_models) #+ ncol(RM_models) + ncol(rBG_models)

library(tictoc)



tic()
mcs_realized <- mcsTest(loss_matrix, 0.05, nboot = 1000, nblock = 10, boot = c("stationary"))
toc()

mcs_realized <- readRDS("MCS_results_volmodels.rds")
#saveRDS(mcs_realized, "MCS_results_volmodels.rds")

#The last integer in $includedSQ or $includedR gives you the one with highest probability in MCS procedure. 

superiorlosses <- colMeans(loss_matrix[, mcs_realized$includedSQ])
superiorlosses

length(mcs_realized$includedR)

mcs_realized$includedSQ
mcs_realized$includedR


#HARQF-MRK-5SEC
mean(QLikes_DRD[[9]][,2])


###################################################################################
#
#
#                            		MINVAR MCS PROCEDURE 
#
#
###################################################################################

test <- t(matrix(c(0.1, 0.9))) %*% (calccov[[1]][[7]][,,1] *10000)%*% matrix(c(0.1, 0.9))


minvar_RM <- readRDS("mvpvariance_RM.rds")
minvar_DRD <- readRDS("mvpvariance_DRD.rds")
minvar_rBG <- readRDS("insample_rBGminvar.rds")
minvar_rDCC <- readRDS("insample_rDCCminvar.rds")

losses_DRD <- do.call(cbind, minvar_DRD)
losses_RM <- do.call(cbind, minvar_RM)
losses_rBG <- do.call(cbind, minvar_rBG)
losses_rDCC <- do.call(cbind, minvar_rDCC)

temprBG <- losses_rBG[1000:2515,-c(10, 19, 28, 42)]
temprDCC <-  losses_rDCC[1000:2515,-c(10, 19, 28, 42) ]

loss_matrix_minvar <- cbind(losses_DRD, losses_RM, temprBG, temprDCC) #1.15 

ncol(losses_DRD) + ncol(losses_RM) 

loss_matrix_minvar <- loss_matrix_minvar[,-c(284,328)] * 1e-4 #downscaled by 1e-4 so we dont work with percentage variance.

#write.table(loss_matrix_minvar, "loss_matrix_minvar.csv")

ncol(losses_DRD) + ncol(losses_RM)

duplicated(t(loss_matrix_minvar))
ncol(loss_matrix_minvar)

mcs_realized_minvar <- mcsTest(loss_matrix_minvar, 0.05, nboot = 1000, nblock = 10, boot = c("stationary"))

#saveRDS(mcs_realized_minvar, "MCS_results_minvar_volmodels.rds")


mcs_realized_minvar <- readRDS("MCS_results_minvar_volmodels.rds")

#if you want to construct tables you can just do ascend since p-vals does not matter. Only which model are contained within.

sortedinclusions <- sort(mcs_realized_minvar$includedSQ)
sortedinclusions
length(sortedinclusions)

colMeans(loss_matrix_minvar[, mcs_realized_minvar$includedSQ]) * 10000


#insamples since I lost out of sample covariances.....
#However, insample still performs much worse than har models. So this only helps
#narrow the superior set down. 

rBG_minvar <- list()
rDCC_minvar <- list()
temp_rBG_minvar <- matrix(0L, nrow = 2516, ncol = 7)
temp_rDCC_minvar <- matrix(0L, nrow = 2516, ncol = 7)


for(k in 5:6){
	for(j in 2:7){
			choice <- c(1,4,5,6,7,8)

		rBGpars <- rBG$vPar
		rDCCpars <- rDCC$vPar
		rBG <- EstimateBivarGARCH(dailyretotc * 100, calccov[[choice[k]]][[j]] * 10000, tol = 1e-5, vPar = rBGpars)
		rDCC <- Estimate_rDCC(dailyretotc * 100, calccov[[choice[k]]][[j]] * 10000, getDates, tol = 1e-5, vPar = rDCCpars)

		for(i in 1:2515){

			weights <- minvar(rBG$vSigma2[,,i])
			weightsrDCC <- minvar(rDCC$vSigma2[,,i])

			temp_rBG_minvar[i, j] <- t(weights) %*% (calccov[[1]][[7]][,,i+1] * 10000) %*% weights
			temp_rDCC_minvar[i, j] <- t(weightsrDCC) %*% (calccov[[1]][[7]][,,i+1] * 10000) %*% weightsrDCC

		}
		print(sprintf("%s", j))
	}
	rBG_minvar[[k]] <- temp_rBG_minvar
	rDCC_minvar[[k]] <- temp_rDCC_minvar
}


#saveRDS(rBG_minvar, "insample_rBGminvar.rds")
#saveRDS(rDCC_minvar, "insample_rDCCminvar.rds")


###################################################################################
#
#
#                       DRD-HAR(Q)-RCOV PLOTS OVER A TWO WEEK PERIOD
#
#
###################################################################################

library(xts)
library(matlib)
source("Previous functions/projectfunctions.R")
source("functions.R")

mergedfrequencies <- readRDS("mergedfrequencies.rds")
calccov <- readRDS("calculatedcovariances.rds")
dataTLT <- readRDS("dataTLT.rds")

getDates <- unlist(lapply(dataTLT, function(x) as.character(index(x[1]))))

for(i in 1:length(getDates)){
	getDates[i] <- strsplit(getDates, " ")[[i]][1]
}


getDates[1415:1429]

ts.plot(calccov[[1]][[7]][2,2,1415:1429] * 10000)

proxycorrelation <- array(0L, dim = c(2,2,2516))

for(i in 1:2516){

	proxycorrelation[,,i] <- realCov(mergedfrequencies[[7]][[i]]*100, T)

}



library(xts)
library(highfrequency)
library(matlib)
library(MCS)
library(Rsolnp)
library(matrixcalc)
library(MASS)
library(Matrix)


Correlations_measures <- readRDS("Correlations_measures.rds")
Quarticities <- readRDS("Quarticity_estimators.rds")

#SKETCH

sigmas_HAR <- array(0L, dim = c(2,2,2516))
sigmas_HARQ <- array(0L, dim = c(2,2,2516))

harq_params <- matrix(0L, ncol=2, nrow = 2516)
window <- 1000

for(i in (window+1):2516){




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
	vol_TLT <- matrix((calccov[[7]][[2]][1,1,(i-window):(i-1)])) * 10000  
	volday_TLT <- vol_TLT
	volweek_TLT <- rowMeans(cbind(vol_TLT, mlag(vol_TLT,4,mean(vol_TLT))))
	volmonth_TLT <- rowMeans(cbind(vol_TLT, mlag(vol_TLT,21,mean(vol_TLT))))

	vol_TLT <- vol_TLT[22:length(vol_TLT)]
	volday_TLT <- volday_TLT[22:length(volday_TLT)]
	volweek_TLT <- volweek_TLT[22:length(volweek_TLT)]
	volmonth_TLT <- volmonth_TLT[22:length(volmonth_TLT)]





	vol_SPY <- matrix((calccov[[7]][[2]][2,2,(i-window):(i-1)])) * 10000  
	volday_SPY <- vol_SPY
	volweek_SPY <- rowMeans(cbind(vol_SPY, mlag(vol_SPY,4,mean(vol_SPY))))
	volmonth_SPY <- rowMeans(cbind(vol_SPY, mlag(vol_SPY,21,mean(vol_SPY))))

	vol_SPY <- vol_SPY[22:length(vol_SPY)]
	volday_SPY <- volday_SPY[22:length(volday_SPY)]
	volweek_SPY <- volweek_SPY[22:length(volweek_SPY)]
	volmonth_SPY <- volmonth_SPY[22:length(volmonth_SPY)]




	rq_TLT <- Quarticities$rq_TLT[(i-window):(i-1),2]
	rq_SPY <- Quarticities$rq_SPY[(i-window):(i-1),2]
	rq_SPY <- rq_SPY[22:length(rq_SPY)]
	rq_TLT <- rq_TLT[22:length(rq_TLT)]
	rq_TLT[rq_TLT>10] <- mean(rq_TLT)
	rq_SPY[rq_SPY>10] <- mean(rq_SPY)


	rq_TLT <- matrix(rq_TLT)

	rq_SPY <- matrix(rq_SPY)



	sqrtrq_TLT <- sqrt(rq_TLT) - mean(sqrt(rq_TLT))

	sqrtrq_SPY <- sqrt(rq_SPY) - mean(sqrt(rq_SPY))



	#############################################################################################
	#
	#
	#									DRD-HAR MODEL (RCOV, MRC, MRK)
	# 
	#
	#############################################################################################
	correlation <- Correlations_measures[[2]][(i-window):(i-1), 7]

	end <- length(RV5min_TLT)

	HAR_TLT <- lm(RV5min_TLT[2:end] ~ volday_TLT[1:(end-1)] + volweek_TLT[1:(end-1)] + volmonth_TLT[1:(end-1)])
	hhat_HAR_TLT <- cbind(1, volday_TLT[end], volweek_TLT[end], volmonth_TLT[end])  %*% matrix(coef(HAR_TLT))

	HAR_SPY <- lm(RV5min_SPY[2:end] ~ volday_SPY[1:(end-1)] + volweek_SPY[1:(end-1)] + volmonth_SPY[1:(end-1)])
	hhat_HAR_SPY <- cbind(1, volday_SPY[end], volweek_SPY[end], volmonth_SPY[end])  %*% matrix(coef(HAR_SPY))

	hhat_HAR_TLT <- volatility.insanity.filter(hhat_HAR_TLT, min(RV5min_TLT), max(RV5min_TLT), mean(RV5min_TLT))$vol
	hhat_HAR_SPY <- volatility.insanity.filter(hhat_HAR_SPY, min(RV5min_SPY), max(RV5min_SPY), mean(RV5min_SPY))$vol

	DRD_HAR <- suppressWarnings(EstimatecorrHAR(cbind(hhat_HAR_TLT, hhat_HAR_SPY), 
		correlation = correlation, proxy = proxycor, 0, T))


	sigmas_HAR[,,i] <- DRD_HAR$vSigma2[,,1]
	
	#############################################################################################
	#
	#
	#									DRD-HARQ MODEL (RCOV, MRC, MRK)
	# 
	#
	#############################################################################################

	HARQ_TLT <- lm(RV5min_TLT[2:end] ~ volday_TLT[1:(end-1)] +  volweek_TLT[1:(end-1)] + volmonth_TLT[1:(end-1)] + I(volday_TLT[1:(end-1)] * sqrtrq_TLT[1:(end-1)]))
	hhat_HARQ_TLT <- cbind(1, volday_TLT[end], volweek_TLT[end], volmonth_TLT[end],volday_TLT[end] * sqrtrq_TLT[end])  %*% matrix(coef(HARQ_TLT))


	HARQ_SPY <- lm(RV5min_SPY[2:end] ~ volday_SPY[1:(end-1)] +  volweek_SPY[1:(end-1)] + volmonth_SPY[1:(end-1)] 
		+ I(volday_SPY[1:(end-1)] * sqrtrq_SPY[1:(end-1)]))
	hhat_HARQ_SPY <- cbind(1, volday_SPY[end], volweek_SPY[end], volmonth_SPY[end],volday_SPY[end] * sqrtrq_SPY[end])  %*% matrix(coef(HARQ_SPY))

	hhat_HARQ_TLT <- volatility.insanity.filter(hhat_HARQ_TLT, min(RV5min_TLT), max(RV5min_TLT), mean(RV5min_TLT))$vol
	hhat_HARQ_SPY <- volatility.insanity.filter(hhat_HARQ_SPY, min(RV5min_SPY), max(RV5min_SPY), mean(RV5min_SPY))$vol

	DRD_HARQ <- suppressWarnings(EstimatecorrHAR(cbind(hhat_HARQ_TLT, hhat_HARQ_SPY), 
		correlation = correlation, proxy = proxycor, 0,T))

	sigmas_HARQ[,,i] <- DRD_HARQ$vSigma2[,,1]

	harq_params[i, 1] <- matrix(coef(HARQ_TLT))[5,]
	harq_params[i, 2] <- matrix(coef(HARQ_SPY))[5,]

	print(sprintf("%s", i))
}




#---------------------insample estimations:
RC5min <- calccov[[1]][[7]][2,1,]/sqrt(calccov[[1]][[7]][1,1,] * calccov[[1]][[7]][2,2,])

RV5min_TLT <-  matrix(calccov[[1]][[7]][1,1,] * 10000)
RV5min_SPY <-  matrix(calccov[[1]][[7]][2,2,] * 10000)

#RV
oneminvol_TLT <- matrix((calccov[[1]][[6]][1,1,])) * 10000  #sqrt
volday_TLT <- oneminvol_TLT
volweek_TLT <- rowMeans(cbind(oneminvol_TLT, mlag(oneminvol_TLT,4,mean(oneminvol_TLT))))
volmonth_TLT <- rowMeans(cbind(oneminvol_TLT, mlag(oneminvol_TLT,21,mean(oneminvol_TLT))))

oneminvol_SPY <- matrix((calccov[[1]][[6]][2,2,])) * 10000  #sqrt
volday_SPY <- oneminvol_SPY
volweek_SPY <- rowMeans(cbind(oneminvol_SPY, mlag(oneminvol_SPY,4,mean(oneminvol_SPY))))
volmonth_SPY <- rowMeans(cbind(oneminvol_SPY, mlag(oneminvol_SPY,21,mean(oneminvol_SPY))))

rq_TLT <- Quarticities$rq_TLT[,6]
rq_TLT[rq_TLT>10] <- mean(rq_TLT)
rq_SPY <- Quarticities$rq_SPY[,6]
rq_SPY[rq_SPY>10] <- mean(rq_SPY)
sqrtrq_TLT_1min <- sqrt(rq_TLT) - mean(sqrt(rq_TLT))
sqrtrq_SPY_1min <- sqrt(rq_SPY) - mean(sqrt(rq_SPY))

ones <- matrix(rep(1, 2494))
#DRD-HAR-RV
HAR_TLT_RV <- lm(RV5min_TLT[23:2516] ~ volday_TLT[22:(2516-1)] + volweek_TLT[22:(2516-1)] + volmonth_TLT[22:(2516-1)])
hhat_HAR_TLT_RV_1min <- cbind(ones, volday_TLT[22:2515], volweek_TLT[22:2515], volmonth_TLT[22:2515])  %*% matrix(coef(HAR_TLT_RV))

HAR_SPY_RV <- lm(RV5min_SPY[23:2516] ~ volday_SPY[22:(2516-1)] + volweek_SPY[22:(2516-1)] + volmonth_SPY[22:(2516-1)])
hhat_HAR_SPY_RV_1min <- cbind(ones, volday_SPY[22:2515], volweek_SPY[22:2515], volmonth_SPY[22:2515])  %*% matrix(coef(HAR_SPY_RV))

DRD_HAR_RV_1min <- EstimatecorrHAR(cbind(hhat_HAR_TLT_RV_1min, hhat_HAR_SPY_RV_1min), correlation = Correlations_measures[[1]][, 7], proxy = RC5min, 1)

#DRD-HARQ-RV
HARQ_TLT_RV <- lm(RV5min_TLT[23:2516] ~ volday_TLT[22:(2516-1)] +  volweek_TLT[22:(2516-1)] + volmonth_TLT[22:(2516-1)] + I(volday_TLT[22:(2516-1)] * sqrtrq_TLT_1min[22:(2516-1)]))
hhat_HARQ_TLT_RV <- cbind(ones, volday_TLT[22:2515], volweek_TLT[22:2515], volmonth_TLT[22:2515],volday_TLT[22:(2516-1)] * sqrtrq_TLT_1min[22:(2516-1)])  %*% matrix(coef(HARQ_TLT_RV))

HARQ_SPY_RV <- lm(RV5min_SPY[23:2516] ~ volday_SPY[22:(2516-1)] +  volweek_SPY[22:(2516-1)] + volmonth_SPY[22:(2516-1)] + I(volday_SPY[22:(2516-1)] * sqrtrq_SPY_1min[22:(2516-1)]))
hhat_HARQ_SPY_RV <- cbind(ones, volday_SPY[22:2515], volweek_SPY[22:2515], volmonth_SPY[22:2515],volday_SPY[22:(2516-1)] * sqrtrq_SPY_1min[22:(2516-1)])  %*% matrix(coef(HARQ_SPY_RV))

#SAME CORRELATION AS HAR MODELS. ONLY THING THAT CHANGED IS THE UNIVARIATE VOLS. 
DRD_HARQ_RV_1min <- EstimatecorrHAR(cbind(hhat_HARQ_TLT_RV, hhat_HARQ_SPY_RV), correlation = Correlations_measures[[1]][, 7], proxy = RC5min, 1)

sigmas_HARQ <- DRD_HARQ_RV_1min$vSigma2
sigmas_HAR <- DRD_HAR_RV_1min$vSigma2

library(ggplot2)

#one day ahead forecasts so, 1000 forecasts 1001. 

#RV CI:

days <- 1840:1865 #500:520  #2260:2270 #1440:1450

CI_TLT <- 1.96*sqrt(2*390*Quarticities$rq_TLT[,6]/3)/sqrt(390)

CI_SPY <- 1.96*sqrt(2*390*Quarticities$rq_SPY[,6]/3)/sqrt(390) 

CI_COV <- 1.96*sqrt(2*390*(10000*calccov[[1]][[6]][2,1,])^4/3)/sqrt(390) 


spy_5min <- calccov[[1]][[6]][2,2,days+1] *10000 

tlt_5min <- calccov[[1]][[6]][1,1,days+1] *10000

cov_5min <- calccov[[1]][[6]][2,1,days+1] *10000



#spy

ggplot() + geom_errorbar(aes(as.Date(getDates[days]), ymin=(spy_5min)-(CI_SPY[days+1]), 
ymax=(spy_5min)+(CI_SPY[days+1])), colour="blue", width=.4)  + 
geom_point(aes(as.Date(getDates[days]), (spy_5min), col ="RV")) + 
geom_line(aes(as.Date(getDates[days]), (sigmas_HARQ[2,2,days]), col = "HARQ"), lwd =1) + 
geom_line(aes(as.Date(getDates[days]), (sigmas_HAR[2,2,days]), col = "HAR"), lwd =1)


#tlt
ggplot() + geom_errorbar(aes(as.Date(getDates[days]), ymin=tlt_5min-CI_TLT[days-1], 
ymax=tlt_5min +CI_TLT[days-1]), colour="blue", width=.4)  + 
geom_point(aes(as.Date(getDates[days]), tlt_5min , col ="RV")) + 
geom_line(aes(as.Date(getDates[days]), (sigmas_HARQ[1,1,days]), col = "HARQ"), lwd =1) + 
geom_line(aes(as.Date(getDates[days]), (sigmas_HAR[1,1,days]), col = "HAR"), lwd =1)


#covariance
ggplot() + geom_errorbar(aes(as.Date(getDates[days]), ymin=cov_5min-CI_COV[days-1], 
ymax=cov_5min +CI_COV[days-1]), colour="blue", width=.4)  + 
geom_point(aes(as.Date(getDates[days]), cov_5min , col ="RV")) + 
geom_line(aes(as.Date(getDates[days]), (sigmas_HARQ[2,1,days]), col = "HARQ"), lwd =1) + 
geom_line(aes(as.Date(getDates[days]), (sigmas_HAR[2,1,days]), col = "HAR"), lwd =1)


ts.plot(sigmas_HAR[2,2,1001:2270] - sigmas_HARQ[2,2, 1001:2270])

max(abs(sigmas_HAR[2,2,1001:2516] - sigmas_HARQ[2,2, 1001:2516]))

ts.plot(Quarticities$rq_SPY[1000:1200,7]) 

ts.plot(sigmas_HARQ[2,1,1001:2516])



#Doing correlations. 

corr_HARQ <- sigmas_HARQ[2,1,1001:2516]/sqrt(sigmas_HARQ[1,1,1001:2516] * sigmas_HARQ[2,2,1001:2516])

corr_HAR <- sigmas_HAR[2,1,1001:2516]/sqrt(sigmas_HAR[1,1,1001:2516] * sigmas_HAR[2,2,1001:2516])



dailyretotc <- xts(t(sapply(mergedfrequencies[[10]], function(x) cbind(x[,1], x[,2]))), order.by = as.Date(getDates))

colnames(dailyretotc) <- c("TLT", "SPY")



rBG_insample <- EstimateBivarGARCH(dailyretotc * 100, calccov[[7]][[2]]*10000) 
rDCC_insample <-

rBG_insample_corr <- rBG_insample$vSigma2[2,1,]/sqrt(rBG_insample$vSigma2[1,1,] * rBG_insample$vSigma2[2,2,])

ggplot() + geom_line(aes(as.Date(getDates)[2014:2265], corr_HARQ[1015:1266], col = "HARQ"), lwd = 1) +
geom_line(aes(as.Date(getDates)[2014:2265], rBG_insample_corr[2014:2265] , col = "rBG"), lwd = 1)







###
library(highfrequency)



Qlikes_quarticity_spy <- numeric()
Qlikes_quarticity_tlt <- numeric()

for(i in 1:2515){
	Qlikes_quarticity_spy[i] <- QLIKE(Quarticities$rq_SPY[i,6], Quarticities$rq_SPY[i+1,7], 1)
	Qlikes_quarticity_tlt[i] <- QLIKE(Quarticities$rq_TLT[i,6], Quarticities$rq_TLT[i+1,7], 1)

}

ts.plot(Qlikes_quarticity[1000:2516])


which(Qlikes_quarticity_spy[1000:2516] < quantile(Qlikes_quarticity_spy[1000:2516], 0.05, na.rm = T))

which(Qlikes_quarticity_spy[1000:2516] > quantile(Qlikes_quarticity_spy[1000:2516], 0.95, na.rm = T))

which(Qlikes_quarticity_tlt[1000:2516] < quantile(Qlikes_quarticity_tlt[1000:2516], 0.05, na.rm = T))

which(Qlikes_quarticity_tlt[1000:2516] > quantile(Qlikes_quarticity_tlt[1000:2516], 0.95, na.rm = T))



###################################################################################
#
#
#            TRYING TO FIGURE OUT WHY GARCH AND DCC MODELS PERFORM POORLY
#
#
###################################################################################








#---------------------prelims----------------------------------------------

library(xts)
library(matlib)
source("Previous functions/projectfunctions.R")
source("functions.R")
library(highfrequency)
library(MCS)
library(Rsolnp)
library(matrixcalc)
library(MASS)
library(Matrix)

mergedfrequencies <- readRDS("mergedfrequencies.rds")
calccov <- readRDS("calculatedcovariances.rds")
dataTLT <- readRDS("dataTLT.rds")
Mixed <- readRDS("Mixed_semicovariances.rds")

getDates <- unlist(lapply(dataTLT, function(x) as.character(index(x[1]))))

for(i in 1:length(getDates)){
	getDates[i] <- strsplit(getDates, " ")[[i]][1]
}

proxycorrelation <- array(0L, dim = c(2,2,2516))

for(i in 1:2516){

	proxycorrelation[,,i] <- realCov(mergedfrequencies[[7]][[i]]*100, T)

}


dailyretotc <- xts(t(sapply(mergedfrequencies[[10]], function(x) cbind(x[,1], x[,2]))), order.by = as.Date(getDates))



Correlations_measures <- readRDS("Correlations_measures.rds")
Quarticities <- readRDS("Quarticity_estimators.rds")


semicovs <- readRDS("semicov_acrossfreq.rds")

Pfreq <- semicovs[[1]]
Nfreq <- semicovs[[2]]
Mfreq <- semicovs[[3]]
#SKETCH

sigmas_HAR <- array(0L, dim = c(2,2,2516))
sigmas_HARQ <- array(0L, dim = c(2,2,2516))
sigmas_rBG <- array(0L, dim = c(2,2,2516))
sigmas_rDCC <- array(0L, dim = c(2,2, 2516))
sigmas_crBG <- array(0L, dim = c(2,2,2516))
sigmas_crDCC <- array(0L, dim = c(2,2, 2516))


harq_params <- matrix(0L, ncol=2, nrow = 2516)
window <- 1000


#KEEPING IT SIMPLE AND CLEAR THEREFORE USE A RCOV1MIN. 

#using pars from in-sample to initialize models.

rBG <- EstimateBivarGARCH(dailyretotc * 100, 
				calccov[[1]][[6]]*10000, tol = 1e-5)

rDCC <- Estimate_rDCC(dailyretotc * 100, calccov[[1]][[7]]*10000, 
				getDates, tol = 1e-5)

decompcov <- list(calccov[[2]][[6]] * 10000, 
		calccov[[3]][[6]] * 10000, Mixed[[6]] * 10000)

crBG <- EstimateBivarGARCHContAsym(dailyretotc * 100, decompcov)

crDCC <- Estimate_rcDCC(dailyretotc * 100, decompcov, getDates = getDates, tol = 1e-5)

crBGpars <- c(0.09, 0.34, 0.04, 0.50)
crDCCpars <- c(0.01, 0.01, 0.303243942, 0.695750518)

for(i in (window+1+1292):2516){




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
	vol_TLT <- matrix((calccov[[1]][[6]][1,1,(i-window):(i-1)])) * 10000  
	volday_TLT <- vol_TLT
	volweek_TLT <- rowMeans(cbind(vol_TLT, mlag(vol_TLT,4,mean(vol_TLT))))
	volmonth_TLT <- rowMeans(cbind(vol_TLT, mlag(vol_TLT,21,mean(vol_TLT))))

	vol_TLT <- vol_TLT[22:length(vol_TLT)]
	volday_TLT <- volday_TLT[22:length(volday_TLT)]
	volweek_TLT <- volweek_TLT[22:length(volweek_TLT)]
	volmonth_TLT <- volmonth_TLT[22:length(volmonth_TLT)]





	vol_SPY <- matrix((calccov[[1]][[6]][2,2,(i-window):(i-1)])) * 10000  
	volday_SPY <- vol_SPY
	volweek_SPY <- rowMeans(cbind(vol_SPY, mlag(vol_SPY,4,mean(vol_SPY))))
	volmonth_SPY <- rowMeans(cbind(vol_SPY, mlag(vol_SPY,21,mean(vol_SPY))))

	vol_SPY <- vol_SPY[22:length(vol_SPY)]
	volday_SPY <- volday_SPY[22:length(volday_SPY)]
	volweek_SPY <- volweek_SPY[22:length(volweek_SPY)]
	volmonth_SPY <- volmonth_SPY[22:length(volmonth_SPY)]




	rq_TLT <- Quarticities$rq_TLT[(i-window):(i-1),6]
	rq_SPY <- Quarticities$rq_SPY[(i-window):(i-1),6]
	rq_SPY <- rq_SPY[22:length(rq_SPY)]
	rq_TLT <- rq_TLT[22:length(rq_TLT)]
	rq_TLT[rq_TLT>10] <- mean(rq_TLT)
	rq_SPY[rq_SPY>10] <- mean(rq_SPY)


	rq_TLT <- matrix(rq_TLT)

	rq_SPY <- matrix(rq_SPY)



	sqrtrq_TLT <- sqrt(rq_TLT) - mean(sqrt(rq_TLT))

	sqrtrq_SPY <- sqrt(rq_SPY) - mean(sqrt(rq_SPY))



	#############################################################################################
	#
	#
	#									DRD-HAR MODEL (RCOV, MRC, MRK)
	# 
	#
	#############################################################################################
	correlation <- Correlations_measures[[2]][(i-window):(i-1), 7]

	end <- length(RV5min_TLT)

	HAR_TLT <- lm(RV5min_TLT[2:end] ~ volday_TLT[1:(end-1)] + volweek_TLT[1:(end-1)] + volmonth_TLT[1:(end-1)])
	hhat_HAR_TLT <- cbind(1, volday_TLT[end], volweek_TLT[end], volmonth_TLT[end])  %*% matrix(coef(HAR_TLT))

	HAR_SPY <- lm(RV5min_SPY[2:end] ~ volday_SPY[1:(end-1)] + volweek_SPY[1:(end-1)] + volmonth_SPY[1:(end-1)])
	hhat_HAR_SPY <- cbind(1, volday_SPY[end], volweek_SPY[end], volmonth_SPY[end])  %*% matrix(coef(HAR_SPY))

	hhat_HAR_TLT <- volatility.insanity.filter(hhat_HAR_TLT, min(RV5min_TLT), max(RV5min_TLT), mean(RV5min_TLT))$vol
	hhat_HAR_SPY <- volatility.insanity.filter(hhat_HAR_SPY, min(RV5min_SPY), max(RV5min_SPY), mean(RV5min_SPY))$vol

	DRD_HAR <- suppressWarnings(EstimatecorrHAR(cbind(hhat_HAR_TLT, hhat_HAR_SPY), 
		correlation = correlation, proxy = proxycor, 0, T))


	sigmas_HAR[,,i] <- DRD_HAR$vSigma2[,,1]
	
	#############################################################################################
	#
	#
	#									DRD-HARQ MODEL (RCOV, MRC, MRK)
	# 
	#
	#############################################################################################

	HARQ_TLT <- lm(RV5min_TLT[2:end] ~ volday_TLT[1:(end-1)] +  volweek_TLT[1:(end-1)] + volmonth_TLT[1:(end-1)] + I(volday_TLT[1:(end-1)] * sqrtrq_TLT[1:(end-1)]))
	hhat_HARQ_TLT <- cbind(1, volday_TLT[end], volweek_TLT[end], volmonth_TLT[end],volday_TLT[end] * sqrtrq_TLT[end])  %*% matrix(coef(HARQ_TLT))


	HARQ_SPY <- lm(RV5min_SPY[2:end] ~ volday_SPY[1:(end-1)] +  volweek_SPY[1:(end-1)] + volmonth_SPY[1:(end-1)] 
		+ I(volday_SPY[1:(end-1)] * sqrtrq_SPY[1:(end-1)]))
	hhat_HARQ_SPY <- cbind(1, volday_SPY[end], volweek_SPY[end], volmonth_SPY[end],volday_SPY[end] * sqrtrq_SPY[end])  %*% matrix(coef(HARQ_SPY))

	hhat_HARQ_TLT <- volatility.insanity.filter(hhat_HARQ_TLT, min(RV5min_TLT), max(RV5min_TLT), mean(RV5min_TLT))$vol
	hhat_HARQ_SPY <- volatility.insanity.filter(hhat_HARQ_SPY, min(RV5min_SPY), max(RV5min_SPY), mean(RV5min_SPY))$vol

	DRD_HARQ <- suppressWarnings(EstimatecorrHAR(cbind(hhat_HARQ_TLT, hhat_HARQ_SPY), 
		correlation = correlation, proxy = proxycor, 0,T))

	sigmas_HARQ[,,i] <- DRD_HARQ$vSigma2[,,1]

	harq_params[i, 1] <- matrix(coef(HARQ_TLT))[5,]
	harq_params[i, 2] <- matrix(coef(HARQ_SPY))[5,]

	#############################################################################################
	#
	#
	#									rBG (Rcov 1min)
	# 
	#
	#############################################################################################
	#end is 979 long.........

	

	pars <- rBG$vPar 

	rBG <- EstimateBivarGARCH(dailyretotc[(i-window):(i-1), ] * 100, 
				calccov[[1]][[6]][,,(i-window):(i-1)]*10000, vPar = pars, tol = 1e-5)

	sigmas_rBG[,,i] <- rBG$vSigma2[,,1000]

	#############################################################################################
	#
	#
	#									rDCC (Rcov 1min)
	# 
	#
	#############################################################################################

	if(parsDCC[1] > 0.5){

		parsDCC <- c(0.20551, 0.79348)
	}


	rDCC <- Estimate_rDCC(dailyretotc[(i-window):(i-1), ] * 100, calccov[[1]][[7]][,,(i-window):(i-1)]*10000, 
				getDates[(i-window):(i-1)], tol = 1e-5, vPar = parsDCC)


	parsDCC <- rDCC$vPar 

	sigmas_rDCC[,,i] <- rDCC$vSigma2[,,1000]

	#############################################################################################
	#
	#
	#									crBG (Rcov 1min)
	# 
	#
	#############################################################################################

		decompcov <- list(Pfreq[[7]][,,(i-window):(i-1)] * 10000, 
		Nfreq[[7]][,,(i-window):(i-1)] * 10000, Mfreq[[7]][,,(i-window):(i-1)] * 10000)

		crBG <- EstimateBivarGARCHContAsym(dailyretotc[(i-window):(i-1), ] * 100, decompcov, vPar = crBGpars)

		crBGpars <-  crBG$vPar 
		
		sigmas_crBG[,,i] <- crBG$vSigma2[,,1000]

	#############################################################################################
	#
	#
	#									crDCC (Rcov 1min)
	# 
	#
	#############################################################################################
		if(crDCCpars[1]>0.5){

			crDCCpars <- c(0.01, 0.01, 0.303243942, 0.695750518)

		}

		cov2 <- calccov[[1]][[7]][,,(i-window):(i-1)] * 10000

		crDCC <- Estimate_rcDCC(dailyretotc[(i-window):(i-1), ] * 100, decompcov, cov2, 
			getDates = getDates[(i-window):(i-1)], tol = 1e-5, vPar = crDCCpars)
		
		crDCCpars <- crDCC$vPar
		
	
		sigmas_crDCC[,,i] <- crDCC$vSigma2[,,1000]



	print(sprintf("%s", i))
}


sigmasonemincovariances <- list(sigmas_rBG = sigmas_rBG, sigmas_rDCC = sigmas_rDCC, sigmas_crBG = sigmas_crBG, sigmas_crDCC = sigmas_crDCC)

#saveRDS(sigmasonemincovariances, "covariances5minrBGcrBGrDCCcrDCC.rds")

#getDates(2014) till getDates(2264) is 2018 period

ts.plot(sigmas_rBG[2,2,2014:2264])

days <-  1001:2516 #2014:2264

library(ggplot2)

correlation_rBG <- sigmas_rBG[2,1,]/sqrt(sigmas_rBG[1,1,]*sigmas_rBG[2,2,])
correlation_rDCC <- sigmas_rDCC[2,1,]/sqrt(sigmas_rDCC[1,1,]*sigmas_rDCC[2,2,])
correlation_crBG <- sigmas_crBG[2,1,]/sqrt(sigmas_crBG[1,1,]*sigmas_crBG[2,2,])
correlation_crDCC <- sigmas_crDCC[2,1,]/sqrt(sigmas_crDCC[1,1,]*sigmas_crDCC[2,2,])
correlation_HAR <- sigmas_HAR[2,1,]/sqrt(sigmas_HAR[1,1,]*sigmas_HAR[2,2,])
correlation_HARQ <- sigmas_HARQ[2,1,]/sqrt(sigmas_HARQ[1,1,]*sigmas_HARQ[2,2,])


ggplot() + geom_line(aes(as.Date(getDates[days]), sigmas_rBG[1,1,days]), col = "red") + 
geom_line(aes(as.Date(getDates[days]), sigmas_rBG[2,2,days]))


ggplot() + geom_line(aes(as.Date(getDates[days]), correlation_rBG[days], col = "rBG")) + 
geom_line(aes(as.Date(getDates[days]), correlation_rDCC[days], col = "rDCC")) +
geom_line(aes(as.Date(getDates[days]), correlation_HAR[days], col = "DRD-HAR")) 


ggplot() + geom_line(aes(as.Date(getDates[1001:2516]), correlation_crBG[1001:2516], col = "crBG")) + 
geom_line(aes(as.Date(getDates[1001:2516]), correlation_crDCC[1001:2516], col = "crDCC")) +
geom_line(aes(as.Date(getDates[days]), correlation_HAR[days], col = "DRD-HAR")) 



head(days-22)



p1 <- ggplot() + geom_line(aes(as.Date(getDates[1001:2516]), sigmas_crDCC[1,1,1001:2516], col = "crDCC")) + 
geom_line(aes(as.Date(getDates[1001:2516]), sigmas_crBG[1,1,1001:2516], col = "crBG")) + 
geom_line(aes(as.Date(getDates[1001:2516]), sigmas_HAR[1,1,1001:2516], col = "DRD-HAR")) + xlab("") + ylab("Variance")


p2 <- ggplot() + geom_line(aes(as.Date(getDates[1001:2516]), sigmas_crDCC[2,2,1001:2516], col = "crDCC")) + 
geom_line(aes(as.Date(getDates[1001:2516]), sigmas_crBG[2,2,1001:2516], col = "crBG")) + 
geom_line(aes(as.Date(getDates[1001:2516]), sigmas_HAR[2,2,1001:2516], col = "DRD-HAR")) +xlab("Dates") +  ylab("Variance")
#geom_line(aes(as.Date(getDates[days]), calccov[[1]][[7]][2,2,days]*10000, col = "Rcov")) 




library(gridExtra)

grid.arrange(p1, p2, ncol = 1)


ggplot() + geom_line(aes(as.Date(getDates[days]), sigmas_crDCC[2,1,days], col = "crDCC")) + 
geom_line(aes(as.Date(getDates[days]), sigmas_crBG[2,1,days], col = "crBG")) + 
geom_line(aes(as.Date(getDates[days]), sigmas_HAR[2,1,days], col = "DRD-HAR")) 


ggplot() + geom_line(aes(as.Date(getDates[days]), sigmas_HAR[2,2,days], col = "DRD-HAR")) +
geom_line(aes(as.Date(getDates[days]), calccov[[1]][[7]][2,2,days]*10000, col = "Rcov")) 

ggplot() + geom_line(aes(as.Date(getDates[1001:2516]), calccov[[1]][[7]][2,2,1001:2516]*10000, col = "Rcov"))  +
geom_line(aes(as.Date(getDates[days-22]), sigmas_crDCC[2,2,days], col = "crBG")) 



