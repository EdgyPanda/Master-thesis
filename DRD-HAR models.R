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

calccov <- readRDS("calculatedcovariances.rds")

mergedfrequencies <- readRDS("mergedfrequencies.rds")

#5min correlations

fivemincorr <- array(0L, dim = c(2,2,2516))

for(i in 1:2516){

	fivemincorr[,,i] <- realCov(mergedfrequencies[[7]][[i]]*100, T)

}



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

#-----------------------------rolling forecast to see if params violate restrictions when posing no restrictions:
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
#-------------------------------------------------------------------------------------------------------------------

params2 <- params[!is.na(params[,1]), ]

#parameter changes for the correlation estimate under five minute sampling frequency. 
library(ggplot2)

ggplot() + geom_line(aes(as.Date(getDates[(window+1):2516]), params2[,1], col="daily")) + 
geom_line(aes(as.Date(getDates[(window+1):2516]), params2[,2], col="weekly")) + 
geom_line(aes(as.Date(getDates[(window+1):2516]), params2[,3], col="monthly")) 


#solving the above using quadratic programing, from solnp function:
#problems with y, it had date indexation and thus misaligned by other data. 






#################################################################################################################
#
#
#								HAR models estimations 5min with standard errors
#
#
#################################################################################################################


#can be estimated using package for time efficiency:

library(highfrequency)

#calculating realized quarticity and tripower quarticity for 5 min interval:

#from highfrequency. Could not find it in the package, so found in source.
RTQ <- function(rData) { 
  returns <- as.vector(as.numeric(rData))
  n <- length(returns)
  tq <- n * (n/(n-2)) *((2^(2/3) * gamma(7/6) * gamma(1/2)^(-1))^(-3)) *  sum(abs(returns[1:(n - 2)])^(4/3) * 
  	abs(returns[2:(n-1)])^(4/3) * abs(returns[3:n])^(4/3))
  return(tq)
} 

rq_TLT <- numeric() 
rq_SPY <- numeric() 

trq_TLT <- numeric()
trq_SPY <- numeric()

for(i in 1:length(mergedfrequencies[[7]])){
				 #realized quarticity
	rq_TLT[i] <- rQuar(mergedfrequencies[[7]][[i]][,1] * 100)
	rq_SPY[i] <- rQuar(mergedfrequencies[[7]][[i]][,2] * 100)
				#realized tripower quarticity
	trq_TLT[i] <- RTQ(mergedfrequencies[[7]][[i]][,1] * 100)
	trq_SPY[i] <- RTQ(mergedfrequencies[[7]][[i]][,2] * 100)

}


rq_TLT <- matrix(rq_TLT)
rq_SPY <- matrix(rq_SPY)
#in HARestimation he both squareroots the realized quarticity and demeans it. 




#-------------------------------------------------TLT---------------------------------------------------------
#
#
#
#
#
#TRYING WITHOUT SQRT BECAUSE IT FOLLOWS PATTONS DATA. THEREFORE TO GET VARIANCES ON PERCENTAGE LOG-RETURN
#YOU HAVE TO SCALE THEM BY A FACTOR 10000!



#RV
fiveminvol_TLT <- matrix((calccov[[1]][[7]][1,1,])) * 10000  #sqrt
volday_TLT <- fiveminvol_TLT
voldayposvar_TLT <- matrix((calccov[[2]][[7]][1,1,])) * 10000 #sqrt
voldaynegvar_TLT <- matrix((calccov[[3]][[7]][1,1,])) * 10000 #sqrt
volweek_TLT <- rowMeans(cbind(fiveminvol_TLT, mlag(fiveminvol_TLT,4,mean(fiveminvol_TLT))))
volmonth_TLT <- rowMeans(cbind(fiveminvol_TLT, mlag(fiveminvol_TLT,21,mean(fiveminvol_TLT))))

#BPV
fiveminbpvol_TLT <- matrix((calccov[[5]][[7]][1,1,])*10000) #sqrt
bpvolday_TLT <- fiveminbpvol_TLT
bpvolweek_TLT <- rowMeans(cbind(fiveminbpvol_TLT, mlag(fiveminbpvol_TLT,4,mean(fiveminbpvol_TLT))))
bpvolmonth_TLT <- rowMeans(cbind(fiveminbpvol_TLT, mlag(fiveminbpvol_TLT,21,mean(fiveminbpvol_TLT))))

#RQ
sqrtrq_TLT <- sqrt(rq_TLT) - mean(sqrt(rq_TLT))
sqrttrq_TLT <- sqrt(trq_TLT) - mean(sqrt(trq_TLT))

rqTLTweek <- rowMeans(cbind(rq_TLT, mlag(rq_TLT,4,mean(rq_TLT))))
sqrtrq_TLTweek <- sqrt(rqTLTweek) - mean(sqrt(rqTLTweek))

rqTLTmonth <- rowMeans(cbind(rq_TLT, mlag(rq_TLT,21,mean(rq_TLT))))
sqrtrq_TLTmonth <- sqrt(rqTLTmonth) - mean(sqrt(rqTLTmonth))


jumpparam_TLT <- ifelse(fiveminvol_TLT - fiveminbpvol_TLT>0, fiveminvol_TLT - fiveminbpvol_TLT, 0)

ones <- matrix(rep(1, 2516-22))



#HAR
HAR_TLT <- lm(fiveminvol_TLT[23:2516] ~ volday_TLT[22:(2516-1)] + volweek_TLT[22:(2516-1)] + volmonth_TLT[22:(2516-1)])
hhat_HAR_TLT <- cbind(ones, volday_TLT[22:2515], volweek_TLT[22:2515], volmonth_TLT[22:2515])  %*% matrix(coef(HAR_TLT))
QLIKE_HAR_TLT <- mean(QLIKE(hhat_HAR_TLT, fiveminvol_TLT[23:2516], 1))
#R-squared
Rsquared_HAR_TLT <-  1- var(fiveminvol_TLT[23:2516] - hhat_HAR_TLT)/var(fiveminvol_TLT[23:2516])


#HARQ 
HARQ_TLT <- lm(fiveminvol_TLT[23:2516] ~ volday_TLT[22:(2516-1)] +  volweek_TLT[22:(2516-1)] + volmonth_TLT[22:(2516-1)] + I(volday_TLT[22:(2516-1)] * sqrtrq_TLT[22:(2516-1)]))
hhat_HARQ_TLT <- cbind(ones, volday_TLT[22:2515], volweek_TLT[22:2515], volmonth_TLT[22:2515],volday_TLT[22:(2516-1)] * sqrtrq_TLT[22:(2516-1)])  %*% matrix(coef(HARQ_TLT))
QLIKE_HARQ_TLT <- mean(QLIKE(hhat_HARQ_TLT, fiveminvol_TLT[23:2516], 1))

#HARQF_TLT
HARQF_TLT <- lm(fiveminvol_TLT[23:2516] ~ volday_TLT[22:(2516-1)] +  volweek_TLT[22:(2516-1)] + 
	volmonth_TLT[22:(2516-1)] + I(volday_TLT[22:(2516-1)] * sqrtrq_TLT[22:(2516-1)]) + 
	I(volweek_TLT[22:(2516-1)] * sqrtrq_TLTweek[22:(2516-1)]) + I(volmonth_TLT[22:(2516-1)] * sqrtrq_TLTmonth[22:(2516-1)]))
hhat_HARQF_TLT <- cbind(ones, volday_TLT[22:2515], volweek_TLT[22:2515], volmonth_TLT[22:2515],
	volday_TLT[22:(2516-1)] * sqrtrq_TLT[22:(2516-1)], volweek_TLT[22:(2516-1)] * sqrtrq_TLTweek[22:(2516-1)], volmonth_TLT[22:(2516-1)] * sqrtrq_TLTmonth[22:(2516-1)])  %*% matrix(coef(HARQF_TLT))
QLIKE_HARQF_TLT <- mean(QLIKE(hhat_HARQF_TLT, fiveminvol_TLT[23:2516], 1))



#HARJ TLT
HARJ_TLT <- lm(fiveminvol_TLT[23:2516] ~ volday_TLT[22:(2516-1)] + volweek_TLT[22:(2516-1)] + volmonth_TLT[22:(2516-1)] + jumpparam_TLT[22:(2516-1)])
hhat_HARJ_TLT <- cbind(ones, volday_TLT[22:2515], volweek_TLT[22:2515], volmonth_TLT[22:2515],jumpparam_TLT[22:(2516-1)])  %*% matrix(coef(HARJ_TLT))
QLIKE_HARJ_TLT <- mean(QLIKE(hhat_HARJ_TLT, fiveminvol_TLT[23:2516], 1))

#HARQJ TLT
HARQJ_TLT <- lm(fiveminvol_TLT[23:2516] ~ volday_TLT[22:(2516-1)] +  volweek_TLT[22:(2516-1)] + volmonth_TLT[22:(2516-1)] + I(volday_TLT[22:(2516-1)] * sqrtrq_TLT[22:(2516-1)]) + jumpparam_TLT[22:(2516-1)])
hhat_HARQJ_TLT <- cbind(ones, volday_TLT[22:2515], volweek_TLT[22:2515], volmonth_TLT[22:2515],volday_TLT[22:(2516-1)] * sqrtrq_TLT[22:(2516-1)], jumpparam_TLT[22:(2516-1)])  %*% matrix(coef(HARQJ_TLT))
QLIKE_HARQJ_TLT <- mean(QLIKE(hhat_HARQJ_TLT, fiveminvol_TLT[23:2516], 1))

#CHAR TLT
CHAR_TLT <- lm(fiveminvol_TLT[23:2516] ~ bpvolday_TLT[22:(2516-1)] + bpvolweek_TLT[22:(2516-1)] + bpvolmonth_TLT[22:(2516-1)])
hhat_CHAR_TLT <- cbind(ones, bpvolday_TLT[22:2515], bpvolweek_TLT[22:2515], bpvolmonth_TLT[22:2515])  %*% matrix(coef(CHAR_TLT))
QLIKE_CHAR_TLT <- mean(QLIKE(hhat_CHAR_TLT, fiveminvol_TLT[23:2516], 1))

#CHARQ TLT
CHARQ_TLT <- lm(fiveminvol_TLT[23:2516] ~ bpvolday_TLT[22:(2516-1)] +  bpvolweek_TLT[22:(2516-1)] + bpvolmonth_TLT[22:(2516-1)] + I(bpvolday_TLT[22:(2516-1)] * sqrttrq_TLT[22:(2516-1)]))
hhat_CHARQ_TLT <- cbind(ones, bpvolday_TLT[22:2515], bpvolweek_TLT[22:2515], bpvolmonth_TLT[22:2515], bpvolday_TLT[22:(2516-1)] * sqrttrq_TLT[22:(2516-1)])  %*% matrix(coef(CHARQ_TLT))
QLIKE_CHARQ_TLT <- mean(QLIKE(hhat_CHARQ_TLT, fiveminvol_TLT[23:2516], 1))
#SHAR TLT
SHAR_TLT  <- lm(fiveminvol_TLT[23:2516] ~ voldayposvar_TLT[22:(2516-1)] + voldaynegvar_TLT[22:(2516-1)] + volweek_TLT[22:(2516-1)] + volmonth_TLT[22:(2516-1)])
hhat_SHAR_TLT <- cbind(ones, voldayposvar_TLT[22:(2516-1)], voldaynegvar_TLT[22:(2516-1)], volweek_TLT[22:2515], volmonth_TLT[22:2515])  %*% matrix(coef(SHAR_TLT))
QLIKE_SHAR_TLT <- mean(QLIKE(hhat_SHAR_TLT, fiveminvol_TLT[23:2516], 1))


QLIKES_TLT <- c(QLIKE_HAR_TLT, QLIKE_HARQ_TLT, QLIKE_HARQF_TLT, QLIKE_HARJ_TLT, QLIKE_HARQJ_TLT, QLIKE_CHAR_TLT, QLIKE_CHARQ_TLT, QLIKE_SHAR_TLT)

library(sandwich)
#ROBUST STANDARD ERRORS:
vcv1 <- sqrt(diag(vcovHC(HAR_TLT)))
vcv2 <- sqrt(diag(vcovHC(HARQ_TLT)))
vcv3 <- sqrt(diag(vcovHC(HARQF_TLT)))
vcv4 <- sqrt(diag(vcovHC(HARJ_TLT)))
vcv5 <- sqrt(diag(vcovHC(HARQJ_TLT)))
vcv6 <- sqrt(diag(vcovHC(CHAR_TLT)))
vcv7 <- sqrt(diag(vcovHC(CHARQ_TLT)))
vcv8 <- sqrt(diag(vcovHC(SHAR_TLT)))


HAR_TLT.show <- rbind(coef(HAR_TLT), vcv1)
HARQ_TLT.show <- rbind(coef(HARQ_TLT), vcv2)
HARQF_TLT.show <- rbind(coef(HARQF_TLT), vcv3)
HARJ_TLT.show <- rbind(coef(HARJ_TLT), vcv4)
HARQJ_TLT.show <- rbind(coef(HARQJ_TLT), vcv5)
CHAR_TLT.show <- rbind(coef(CHAR_TLT), vcv6)
CHARQ_TLT.show <- rbind(coef(CHARQ_TLT), vcv7)
SHAR_TLT.show <- rbind(coef(SHAR_TLT), vcv8)

(cbind(HAR_TLT.show, SHAR_TLT.show, HARQ_TLT.show, HARQF_TLT.show, HARJ_TLT.show, HARQJ_TLT.show, CHAR_TLT.show, CHARQ_TLT.show))

#--------------------------------------------SPY-----------------------------------------------------

#RV
fiveminvol_SPY <- matrix((calccov[[1]][[7]][2,2,]*10000)) #sqrt
volday_SPY <- fiveminvol_SPY
voldayposvar_SPY <- matrix((calccov[[2]][[7]][2,2,]* 10000)) 
voldaynegvar_SPY <- matrix((calccov[[3]][[7]][2,2,]* 10000)) 
volweek_SPY <- rowMeans(cbind(fiveminvol_SPY, mlag(fiveminvol_SPY,4,mean(fiveminvol_SPY)))) #sqrt
volmonth_SPY <- rowMeans(cbind(fiveminvol_SPY, mlag(fiveminvol_SPY,21,mean(fiveminvol_SPY)))) #sqrt

#BPV
fiveminbpvol_SPY <- matrix((calccov[[5]][[7]][2,2,]*10000)) #sqrt
bpvolday_SPY <- fiveminbpvol_SPY
bpvolweek_SPY <- rowMeans(cbind(fiveminbpvol_SPY, mlag(fiveminbpvol_SPY,4,mean(fiveminbpvol_SPY))))
bpvolmonth_SPY <- rowMeans(cbind(fiveminbpvol_SPY, mlag(fiveminbpvol_SPY,21,mean(fiveminbpvol_SPY))))

#RQ
sqrtrq_SPY <- sqrt(rq_SPY) - mean(sqrt(rq_SPY))
sqrttrq_SPY <- sqrt(trq_SPY) - mean(sqrt(trq_SPY))

rqSPYweek <- rowMeans(cbind(rq_SPY, mlag(rq_SPY,4,mean(rq_SPY))))
sqrtrq_SPYweek <- sqrt(rqSPYweek) - mean(sqrt(rqSPYweek))

rqSPYmonth <- rowMeans(cbind(rq_SPY, mlag(rq_SPY,21,mean(rq_SPY))))
sqrtrq_SPYmonth <- sqrt(rqSPYmonth) - mean(sqrt(rqSPYmonth))


jumpparam_SPY <- ifelse(fiveminvol_SPY - fiveminbpvol_SPY>0, fiveminvol_SPY - fiveminbpvol_SPY, 0)

#HAR
HAR_SPY <- lm(fiveminvol_SPY[23:2516] ~ volday_SPY[22:(2516-1)] + volweek_SPY[22:(2516-1)] + volmonth_SPY[22:(2516-1)])
hhat_HAR_SPY <- cbind(ones, volday_SPY[22:2515], volweek_SPY[22:2515], volmonth_SPY[22:2515])  %*% matrix(coef(HAR_SPY))
QLIKE_HAR_SPY <- colMeans(QLIKE(hhat_HAR_SPY, fiveminvol_SPY[23:2516], 1))
#HARQ 
HARQ_SPY <- lm(fiveminvol_SPY[23:2516] ~ volday_SPY[22:(2516-1)] +  volweek_SPY[22:(2516-1)] + volmonth_SPY[22:(2516-1)] + I(volday_SPY[22:(2516-1)] * sqrtrq_SPY[22:(2516-1)]))
hhat_HARQ_SPY <- cbind(ones, volday_SPY[22:2515], volweek_SPY[22:2515], volmonth_SPY[22:2515],volday_SPY[22:(2516-1)] * sqrtrq_SPY[22:(2516-1)])  %*% matrix(coef(HARQ_SPY))
QLIKE_HARQ_SPY <- colMeans(QLIKE(hhat_HARQ_SPY, fiveminvol_SPY[23:2516], 1))


#HARQF_SPY
HARQF_SPY <- lm(fiveminvol_SPY[23:2516] ~  volday_SPY[22:(2516-1)] +  volweek_SPY[22:(2516-1)] + 
	volmonth_SPY[22:(2516-1)] + I(volday_SPY[22:(2516-1)] * sqrtrq_SPY[22:(2516-1)]) + 
	I(volweek_SPY[22:(2516-1)] * sqrtrq_SPYweek[22:(2516-1)]) + I(volmonth_SPY[22:(2516-1)] * sqrtrq_SPYmonth[22:(2516-1)]))
hhat_HARQF_SPY <- cbind(ones, volday_SPY[22:2515], volweek_SPY[22:2515], volmonth_SPY[22:2515],
	volday_SPY[22:(2516-1)] * sqrtrq_SPY[22:(2516-1)], volweek_SPY[22:(2516-1)] * sqrtrq_SPYweek[22:(2516-1)], volmonth_SPY[22:(2516-1)] * sqrtrq_SPYmonth[22:(2516-1)])  %*% matrix(coef(HARQF_SPY))
hhat_HARQF_SPY[hhat_HARQF_SPY<0] <- mean(hhat_HARQF_SPY[80:110])
QLIKE_HARQF_SPY <- mean(QLIKE(hhat_HARQF_SPY, fiveminvol_SPY[23:2516], 1))


#HARJ SPY
HARJ_SPY <- lm(fiveminvol_SPY[23:2516] ~ volday_SPY[22:(2516-1)] + volweek_SPY[22:(2516-1)] + volmonth_SPY[22:(2516-1)] + jumpparam_SPY[22:(2516-1)])
hhat_HARJ_SPY <- cbind(ones, volday_SPY[22:2515], volweek_SPY[22:2515], volmonth_SPY[22:2515],jumpparam_SPY[22:(2516-1)])  %*% matrix(coef(HARJ_SPY))
QLIKE_HARJ_SPY <- mean(QLIKE(hhat_HARJ_SPY, fiveminvol_SPY[23:2516], 1))

#HARQJ SPY
HARQJ_SPY <- lm(fiveminvol_SPY[23:2516] ~ volday_SPY[22:(2516-1)] +  volweek_SPY[22:(2516-1)] + volmonth_SPY[22:(2516-1)] + I(volday_SPY[22:(2516-1)] * sqrtrq_SPY[22:(2516-1)]) + jumpparam_SPY[22:(2516-1)])
hhat_HARQJ_SPY <- cbind(ones, volday_SPY[22:2515], volweek_SPY[22:2515], volmonth_SPY[22:2515],volday_SPY[22:(2516-1)] * sqrtrq_SPY[22:(2516-1)], jumpparam_SPY[22:(2516-1)])  %*% matrix(coef(HARQJ_SPY))
QLIKE_HARQJ_SPY <- mean(QLIKE(hhat_HARQJ_SPY, fiveminvol_SPY[23:2516], 1))

#CHAR SPY
CHAR_SPY <- lm(fiveminvol_SPY[23:2516] ~ bpvolday_SPY[22:(2516-1)] + bpvolweek_SPY[22:(2516-1)] + bpvolmonth_SPY[22:(2516-1)])
hhat_CHAR_SPY <- cbind(ones, bpvolday_SPY[22:2515], bpvolweek_SPY[22:2515], bpvolmonth_SPY[22:2515])  %*% matrix(coef(CHAR_SPY))
QLIKE_CHAR_SPY <- mean(QLIKE(hhat_CHAR_SPY, fiveminvol_SPY[23:2516], 1))


#CHARQ SPY
CHARQ_SPY <- lm(fiveminvol_SPY[23:2516] ~ bpvolday_SPY[22:(2516-1)] +  bpvolweek_SPY[22:(2516-1)] + bpvolmonth_SPY[22:(2516-1)] + I(bpvolday_SPY[22:(2516-1)] * sqrttrq_SPY[22:(2516-1)]))
hhat_CHARQ_SPY <- cbind(ones, bpvolday_SPY[22:2515], bpvolweek_SPY[22:2515], bpvolmonth_SPY[22:2515], bpvolday_SPY[22:(2516-1)] * sqrttrq_SPY[22:(2516-1)])  %*% matrix(coef(CHARQ_SPY))
QLIKE_CHARQ_SPY <- mean(QLIKE(hhat_CHARQ_SPY, fiveminvol_SPY[23:2516], 1))
#SHAR SPY
SHAR_SPY  <- lm(fiveminvol_SPY[23:2516] ~ voldayposvar_SPY[22:(2516-1)] + voldaynegvar_SPY[22:(2516-1)] + volweek_SPY[22:(2516-1)] + volmonth_SPY[22:(2516-1)])
hhat_SHAR_SPY <- cbind(ones, voldayposvar_SPY[22:(2516-1)], voldaynegvar_SPY[22:(2516-1)], volweek_SPY[22:2515], volmonth_SPY[22:2515])  %*% matrix(coef(SHAR_SPY))
QLIKE_SHAR_SPY <- mean(QLIKE(hhat_SHAR_SPY, fiveminvol_SPY[23:2516], 1))

QLIKES_SPY <- c(QLIKE_HAR_SPY, QLIKE_HARQ_SPY, QLIKE_HARQF_SPY, QLIKE_HARJ_SPY, QLIKE_HARQJ_SPY, QLIKE_CHAR_SPY, QLIKE_CHARQ_SPY, QLIKE_SHAR_SPY)


library(sandwich)

#ROBUST STANDARD ERRORS:
vcv1s <- sqrt(diag(vcovHC(HAR_SPY)))
vcv2s <- sqrt(diag(vcovHC(HARQ_SPY)))
vcv3s <- sqrt(diag(vcovHC(HARQF_SPY)))
vcv4s <- sqrt(diag(vcovHC(HARJ_SPY)))
vcv5s <- sqrt(diag(vcovHC(HARQJ_SPY)))
vcv6s <- sqrt(diag(vcovHC(CHAR_SPY)))
vcv7s <- sqrt(diag(vcovHC(CHARQ_SPY)))
vcv8s <- sqrt(diag(vcovHC(SHAR_SPY)))


HAR_SPY.show <- rbind(coef(HAR_SPY), vcv1s)
HARQ_SPY.show <- rbind(coef(HARQ_SPY), vcv2s)
HARQF_SPY.show <- rbind(coef(HARQF_SPY), vcv3s)
HARJ_SPY.show <- rbind(coef(HARJ_SPY), vcv4s)
HARQJ_SPY.show <- rbind(coef(HARQJ_SPY), vcv5s)
CHAR_SPY.show <- rbind(coef(CHAR_SPY), vcv6s)
CHARQ_SPY.show <- rbind(coef(CHARQ_SPY), vcv7s)
SHAR_SPY.show <- rbind(coef(SHAR_SPY), vcv8s)

#########################################################################################################
#
#
#					PERSISTENCE IN SEMIVARIANCES FOR TLT AS SEEN BY AUTOCORRELATIONS
#
#
#########################################################################################################


acf_TLT <-  matrix(0L, ncol = 7, nrow = 101)
acf_posvol_TLT <- matrix(0L, ncol = 7, nrow = 101)
acf_negvol_TLT <- matrix(0L, ncol = 7, nrow = 101)

for(i in 1:7){
	temp <- matrix((calccov[[1]][[i]][1,1,])) * 10000
	temp2 <- matrix((calccov[[2]][[i]][1,1,])) * 10000 #sqrt
	temp3 <- matrix((calccov[[3]][[i]][1,1,])) * 10000 #sqrt

	acf_TLT[,i] <- acf(temp, lag.max  = 100, plot =  F)$acf
	acf_posvol_TLT[,i] <- acf(temp2, lag.max  = 100, plot =  F)$acf
	acf_negvol_TLT[,i] <- acf(temp3, lag.max  = 100, plot =  F)$acf

}

acfpos_TLT <- rowMeans(acf_posvol_TLT)
acfneg_TLT <- rowMeans(acf_negvol_TLT)
acf_TLT <- rowMeans(acf_TLT)
p1 <- ggplot() + geom_line(aes(1:100, acfpos_TLT[-1], col ="RVPos"), lwd = 1) + geom_line(aes(1:100, acfneg_TLT[-1], col ="RVneg"), lwd = 1) +
geom_point(aes(1:100, acf_TLT[-1], col="RV"), size = 1.4, shape = 16) + ylim(c(0,0.6)) + xlab("Lag") + ylab("Autocorrelation") +
theme(legend.title = element_blank(), legend.position = c(0.8, 0.8), legend.background = element_rect(fill="lightblue", 
size=0.5, linetype="solid", colour = 'darkblue'), 
legend.text = element_text(colour="black", size=9, face="bold"), legend.text.align = 0) + 
scale_fill_manual(name ='', values = c("blue", "red")) 

#ggsave("TLT_Autocorrelation.eps", device = "eps")

#See below for combined graph with QLIKE losses under different DRD-type models. 

#########################################################################################################
#
#
#					SIGNATURE VOLATILITY PLOT OF REALIZED QUARTICITY ESTIMATORS
#
#
#########################################################################################################

library(highfrequency)
dataTLT <- readRDS("dataTLT.rds")
dataSPY <- readRDS("dataSPY.rds")


secs <-c(1,5,15,20,30,60,300,900,1800)


frequenciesTLT <- list()
frequenciesSPY <- list()

tempTLT <- list()
tempSPY <- list()

days <- 2516

for(j in 1:length(secs)){
	for(i in 1:days){

		tempTLT[[i]] <- aggregatets(dataTLT[[i]], on="seconds", k=secs[j])
		tempSPY[[i]] <- aggregatets(dataSPY[[i]], on="seconds", k=secs[j])
	}

	frequenciesTLT[[j]] <- tempTLT
	frequenciesSPY[[j]] <- tempSPY
	print(sprintf("%s",j))
}

for(j in 1:length(secs)){
	for(i in 1:days){

		frequenciesTLT[[j]][[i]] <- diff(log(frequenciesTLT[[j]][[i]]))[-1] * 100
		frequenciesSPY[[j]][[i]] <- diff(log(frequenciesSPY[[j]][[i]]))[-1] * 100

	}
}


RQ_TLT <- matrix(0L, ncol = length(secs), nrow= days)
RQ_SPY <- matrix(0L, ncol = length(secs), nrow= days)

TRQ_TLT <- matrix(0L, ncol = length(secs), nrow= days)
TRQ_SPY <- matrix(0L, ncol = length(secs), nrow= days)

#contains outliers
for(j in 1:length(secs)){
	for(i in 1:days){

		RQ_TLT[i,j] <- rQuar(frequenciesTLT[[j]][[i]])
		RQ_SPY[i,j] <- rQuar(frequenciesSPY[[j]][[i]])
		TRQ_TLT[i,j] <- RTQ(frequenciesTLT[[j]][[i]])
		TRQ_SPY[i,j] <- RTQ(frequenciesSPY[[j]][[i]])

	}
}

#removing outliers

rem <- which(TRQ_TLT[,1]>11)
rem2 <- which(TRQ_TLT[,2]>11)
rem3 <- which(TRQ_TLT[,3]>11)
rem4 <- which(TRQ_TLT[,4]>11)
rem5 <- which(TRQ_TLT[,5]>11)
rem6 <- which(TRQ_SPY[,1]>9)
rem7 <- which(TRQ_SPY[,2]>9)
rem8 <- which(RQ_SPY[,1]>11)
rem9 <- which(RQ_SPY[,2]>11)
remove <- unique(c(rem, rem2, rem3, rem4, rem5, rem6, rem7,rem8, rem9))

RQ_TLT <- RQ_TLT[-remove, ]
RQ_SPY <- RQ_SPY[-remove, ]

TRQ_SPY <- TRQ_SPY[-remove, ]
TRQ_TLT <- TRQ_TLT[-remove, ]
RQ_TLT <- colMeans(RQ_TLT)
RQ_SPY <- colMeans(RQ_SPY)
RQ_SPY[1] <- RQ_SPY[1] + 2.5 #removing too many values
TRQ_TLT <- colMeans(TRQ_TLT)
TRQ_SPY <- colMeans(TRQ_SPY)


library(ggplot2)

ggplot() + geom_line(aes(secs, RQ_TLT, col ="RQ_TLT"), lwd =1) + geom_line(aes(secs, RQ_SPY, col ="RQ_SPY"), lwd=1) + 
geom_line(aes(secs, TRQ_TLT, col ="TRQ_TLT"), lwd=1) + geom_line(aes(secs, TRQ_SPY, col ="TRQ_SPY"),lwd=1) + ylim(0,4) +
ylab("Realized Quarticity") + xlab("Seconds") + 
scale_x_continuous(breaks = c(secs[1],secs[6:9]), labels = c("1", "60", "300", "900", "1800"))

p2 <- ggplot() + geom_smooth(aes(x=secs, y= RQ_TLT, col = "RQ_TLT"), span = 0.3,  se=FALSE) + 
geom_smooth(aes(secs, RQ_SPY, col ="RQ_SPY"), span = 0.3,  se=FALSE) + 
geom_smooth(aes(secs, TRQ_SPY, col ="TRQ_SPY"), span =0.3,  se=FALSE) + 
geom_smooth(aes(secs, TRQ_TLT, col ="TRQ_TLT"), span = 0.3,  se=FALSE) + 
ylim(0,3.45) + ylab("Realized Quarticity") + xlab("Seconds") + 
scale_x_continuous(breaks = c(secs[1],secs[6:9]), labels = c("1", "60", "300", "900", "1800")) + 
theme(legend.title = element_blank(), legend.position = c(0.8, 0.8), legend.background = element_rect(fill="lightblue", 
size=0.5, linetype="solid", colour = 'darkblue'), 
legend.text = element_text(colour="black", size=9, face="bold"), legend.text.align = 0) + 
scale_fill_manual(name ='', values = c("blue", "red")) 




#########################################################################################################
#
#
#					THE DRD-HAR MODEL WITHOUT USING STANDARDIZED RESIDUALS
#
#
#########################################################################################################


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
	mini <-  sum(mini)

	lOut <- list()
	lOut[["mini"]] <- mini 
	lOut[["minis"]] <- minis

	return(lOut)
}


ineqconstraint <- function(vPar, ...){

	dA <- vPar[1]  
	dB <- vPar[2]
	dC <- vPar[3]

	return(dA+dB+dC)

}

min.RSS <- function(vPar, data){

	rss <- minimizingfunc(data, vPar)$mini

	return(rss)

}


EstimatecorrHAR <- function(variances, correlation = NULL, proxy, trace=1, ineqfun = ineqconstraint, ineqLB = 0.00, ineqUB = 0.9999){
  #mEta is standardized residuals from univariate models!
  #variances should be produced by the univariate models and not come from realized measures!
  #correlation should be found via eg. five min samples but adhere to the same frequency as variances from HAR models. 

  dA = 0.01
  dB  = 0.001185
  dC = 0.2
  ## vector of starting parameters
  vPar = c(dA, dB, dC)

  correlation <- as.matrix(correlation)
  corrday <- correlation
  corrweek <- rowMeans(cbind(correlation, mlag(correlation,4,mean(correlation))))
  corrmonth <- rowMeans(cbind(correlation, mlag(correlation,21,mean(correlation)))) 

  data <- list(proxy[23:2516], corrday[22:2515], corrweek[22:2515], corrmonth[22:2515])
  optimizer = solnp(vPar, fun = min.RSS, data = data, 
                    ineqfun = ineqfun, #the inequality constraint
                    ineqLB  = ineqLB, ## the inequality lower bound
                    ineqUB = ineqUB, ## the inequality lower bound, i.e. 0.0 <= a + b + c < 0.9999
                    ## lower and upper bounds for all parameters
                    LB = c(0, 0, 0), UB = c(0.9999, 0.9999, 0.9999),
                    control = list(tol = 1e-22, outer.iter = 800, inner.iter = 1200, delta = 1e-8,
                    trace = trace))
                    

  params <- optimizer$pars

  min <- tail(optimizer$values, 1)

  hessian <- optimizer$hessian

  scores <- matrix(0L, nrow=(nrow(variances)), ncol = 3)

  step <- 1e-5 * vPar

  for(i in 1:length(step)){

	h <- step[i]
    delta <- rep(0, length(vPar))
    delta[i] <- h
																
	loglikeminus <- minimizingfunc(data, vPar-delta)$minis
	loglikeplus <- minimizingfunc(data, vPar+delta)$minis

	scores[,i] <- (loglikeplus - loglikeminus)/(2*h)

  }

  J <- (t(scores) %*% scores)/2516

  I <- optimizer$hessian/2516

  I <- solve(I)[-1 ,-1]

  vars <- (I * J * I)/2516
  
  rse <- sqrt(diag(vars))

  t.stat <- vPar/rse


  #calculating covariances: 
  hhatcorrHAR <- cbind(corrday[22:2515], corrweek[22:2515], corrmonth[22:2515]) %*% matrix(params)
  

  covs <- hhatcorrHAR * sqrt(variances[,1]) * sqrt(variances[,2])

  vSigma2 <- array(0L, dim = c(2,2, (nrow(variances))))

  for(i in 1:(nrow(variances))){

	vSigma2[,,i] <- matrix(c(variances[i,1], covs[i], covs[i], variances[i,2]), ncol=2, nrow=2) 

  }

	  #R-squared

	Rsquared <- 1-var(correlation[23:2516] - hhatcorrHAR)/var(correlation[23:2516])



  lOut <- list()

  lOut[["vPar"]] <- params
  lOut[["vSigma2"]] <- vSigma2
  lOut[["R2"]] <- Rsquared
  lOut[["estcor"]] <- hhatcorrHAR
  lOut[["MSE"]] <- min/2516 
  lOut[["rse"]] <- rse
  lOut[["hessian"]] <- hessian
  lOut[["Tstats"]] <- t.stat

  return(lOut)

 }



#######################################################################################################################
#
#
#										In-sample estimations for DRD-HAR model 
# 
#
#######################################################################################################################

ones <- rep(1, 2516-22)
ones <- as.matrix(ones)

fiveminvol_TLT <- matrix((calccov[[1]][[7]][1,1,])) * 10000  #sqrt
volday_TLT <- fiveminvol_TLT
volweek_TLT <- rowMeans(cbind(fiveminvol_TLT, mlag(fiveminvol_TLT,4,mean(fiveminvol_TLT))))
volmonth_TLT <- rowMeans(cbind(fiveminvol_TLT, mlag(fiveminvol_TLT,21,mean(fiveminvol_TLT))))

#RV
fiveminvol_SPY <- matrix((calccov[[1]][[7]][2,2,]*10000)) #sqrt
volday_SPY <- fiveminvol_SPY
volweek_SPY <- rowMeans(cbind(fiveminvol_SPY, mlag(fiveminvol_SPY,4,mean(fiveminvol_SPY)))) #sqrt
volmonth_SPY <- rowMeans(cbind(fiveminvol_SPY, mlag(fiveminvol_SPY,21,mean(fiveminvol_SPY)))) #sqrt


HAR_TLT <- lm(fiveminvol_TLT[23:2516] ~ volday_TLT[22:(2516-1)] + volweek_TLT[22:(2516-1)] + volmonth_TLT[22:(2516-1)])
hhat_HAR_TLT <- cbind(ones, volday_TLT[22:2515], volweek_TLT[22:2515], volmonth_TLT[22:2515])  %*% matrix(coef(HAR_TLT))

HAR_SPY <- lm(fiveminvol_SPY[23:2516] ~ volday_SPY[22:(2516-1)] + volweek_SPY[22:(2516-1)] + volmonth_SPY[22:(2516-1)])
hhat_HAR_SPY <- cbind(ones, volday_SPY[22:2515], volweek_SPY[22:2515], volmonth_SPY[22:2515])  %*% matrix(coef(HAR_SPY))

variances <- cbind(hhat_HAR_TLT , hhat_HAR_SPY)

RC5min <- calccov[[1]][[7]][2,1,]/sqrt(calccov[[1]][[7]][1,1,] * calccov[[1]][[7]][2,2,])


#i've used estimation scheme based on realized correlations and not standardized residuals, which is, the same. 
tt2 <- EstimatecorrHAR(variances, correlation = RC5min, proxy = RC5min, 1)


RC1min <- calccov[[1]][[6]][2,1,]/sqrt(calccov[[1]][[6]][1,1,] * calccov[[1]][[6]][2,2,])
RC15min <- calccov[[1]][[8]][2,1,]/sqrt(calccov[[1]][[8]][1,1,] * calccov[[1]][[8]][2,2,])
RC30min <- calccov[[1]][[9]][2,1,]/sqrt(calccov[[1]][[9]][1,1,] * calccov[[1]][[9]][2,2,])
RCdaily <- calccov[[1]][[10]][2,1,]/sqrt(calccov[[1]][[10]][1,1,] * calccov[[1]][[10]][2,2,])


fivemincorr <- list(NULL,NULL,NULL,NULL,NULL,as.matrix(RC1min), as.matrix(RC5min), as.matrix(RC15min), as.matrix(RC30min), as.matrix(RCdaily))

HARCor_freq <- list()

fiveminvol_TLT <- matrix((calccov[[1]][[7]][1,1,])) * 10000
fiveminvol_SPY <- matrix((calccov[[1]][[7]][2,2,]*10000)) #sqrt

for(i in 6:9){
	vol_TLT <- matrix((calccov[[1]][[i]][1,1,])) * 10000  #sqrt
	volday_TLT <- vol_TLT
	volweek_TLT <- rowMeans(cbind(vol_TLT, mlag(vol_TLT,4,mean(vol_TLT))))
	volmonth_TLT <- rowMeans(cbind(vol_TLT, mlag(vol_TLT,21,mean(vol_TLT))))

	#RV
	vol_SPY <- matrix((calccov[[1]][[i]][2,2,]*10000)) #sqrt
	volday_SPY <- vol_SPY
	volweek_SPY <- rowMeans(cbind(vol_SPY, mlag(vol_SPY,4,mean(vol_SPY)))) #sqrt
	volmonth_SPY <- rowMeans(cbind(vol_SPY, mlag(vol_SPY,21,mean(vol_SPY)))) #sqrt


	HAR_TLT <- lm(fiveminvol_TLT[23:2516] ~ volday_TLT[22:(2516-1)] + volweek_TLT[22:(2516-1)] + volmonth_TLT[22:(2516-1)])
	hhat_HAR_TLT <- cbind(ones, volday_TLT[22:2515], volweek_TLT[22:2515], volmonth_TLT[22:2515])  %*% matrix(coef(HAR_TLT))

	HAR_SPY <- lm(fiveminvol_SPY[23:2516] ~ volday_SPY[22:(2516-1)] + volweek_SPY[22:(2516-1)] + volmonth_SPY[22:(2516-1)])
	hhat_HAR_SPY <- cbind(ones, volday_SPY[22:2515], volweek_SPY[22:2515], volmonth_SPY[22:2515])  %*% matrix(coef(HAR_SPY))

	variances <- cbind(hhat_HAR_TLT, hhat_HAR_SPY)


	HARCor_freq[[i]] <- EstimatecorrHAR(variances, correlation = fivemincorr[[i]],proxy = RC5min, 1)

}

params <- matrix(unlist(lapply(HARCor_freq, function(x) x$vPar)), ncol = 3, byrow = T)

apply(params, MARGIN = c(2), FUN =function(x) quantile(x, 0.1))
apply(params, MARGIN = c(2), FUN =function(x) quantile(x, 0.9))





qlikes <- matrix(0L, nrow = 2494, ncol = 4)

qlikes5min <- numeric()

for(i in 1:2494){
					#vsigma2 starts from 22 and ends at 2515
		qlikes5min[i] <- QLIKE(HARCor_freq[[7]]$vSigma2[,,i], calccov[[1]][[7]][,,22+i]*10000, 2)
	}


for(j in 1:4){
	for(i in 1:2494){
					#vsigma2 starts from 22 and ends at 2515
		qlikes[i, j] <- QLIKE(HARCor_freq[[j+5]]$vSigma2[,,i], calccov[[1]][[7]][,,22+i]*10000, 2)
	}
}

colMeans(qlikes)
mean(colMeans(qlikes))


mean(c(HARCor_freq[[6]]$R2, HARCor_freq[[7]]$R2, HARCor_freq[[8]]$R2, HARCor_freq[[9]]$R2))
HARCor_freq[[6]]$R2; HARCor_freq[[7]]$R2; HARCor_freq[[8]]$R2; HARCor_freq[[9]]$R2


#######################################################################################################################
#
#
#							BOOTSTRAPPING DRD-HAR model (TIME-SERIES REDUCES TO RESIDUALS)
# 
#
#######################################################################################################################

#REMEMBER THAT THE INTERDEPENDENCE OF THE UNIVARIATE MODELS IS COMPLETELY CAPTURED IN THE STANDARDIZED RESIDUALS. 
#TAKING THE COVARIANCE OF THE STANDARDIZED RESIDUALS SHOULD GET YOU THE CORRELATION (DUE TO NORMALIZING WITH STD.DEV)

library(boot)

library(np)
#uses the row to construct the block length. Therefore you should transpose it. 
ll <- tsboot(dailyretotc[,1], mean ,R = 1000, l = 97, sim = "fixed", endcorr = TRUE)
indexation <- boot.array(ll)


HARCor_Boot <- matrix(0L, nrow = 1000, ncol = 3)

for(i in 1:1000){

	RC5min_boot <- RC5min[indexation[i, ]]

	HARCor_Boot[i, ] <- EstimatecorrHAR2(variances, correlation = RC5min_boot, proxy = RC5min_boot, 1)$vPar

}


standarderror_HARCor <- apply(HARCor_Boot, MARGIN = c(2), FUN = function(x) sd(x))


ttt <- RC5min[indexation[1, ]]



#######################################################################################################################
#
#
#					QLIKE COMPARISON USING 5-MIN FREQUENCY AND CONSTRUCTING IT USING BAR-CHART
# 
#
#######################################################################################################################

#UNIVARIATE INSANITY FILTER!
volatility.insanity.filter <- function(forecasts, lb, ub, adjujsted_forecast){

	ut1 <- forecasts
	items <- sum(forecasts < lb) + sum( forecasts > ub)

	ut1[forecasts < lb] = adjujsted_forecast
	ut1[forecasts > ub] = adjujsted_forecast

	lOut <- list()
	lOut[["vol"]] <- ut1
	lOut[["items"]] <- items

	return(lOut)

}




#Get correlations for all measures on 5 minute scheme. 
#Theres 8 measures in total, with two being rscov+ and rscor-, which will be discarded, so 6. 

#its not fiveminutes!
Correlations_5min_measures <- matrix(0L, ncol=6, nrow = 2516)

#list of kernels based on frequency. 5 minute is list object 7. 
kernel.param <- readRDS("bandwidthH.rds")

for(i in 1:2516){

	Correlations_5min_measures[i,1] <- realCov(mergedfrequencies[[6]][[i]]*100, T)[2,1]
	Correlations_5min_measures[i,2] <-  preavthrCOV(mergedfrequencies[[6]][[i]]*100, F)[2,1] / (sqrt(preavthrCOV(mergedfrequencies[[6]][[i]]*100, F)[1,1] * preavthrCOV(mergedfrequencies[[6]][[i]]*100, F)[2,2]))
	Correlations_5min_measures[i,3] <-  preavBPCOV(mergedfrequencies[[6]][[i]]*100,F,T,F,1)[2,1] #bpcov
	Correlations_5min_measures[i,5] <-  preavCov(mergedfrequencies[[6]][[i]]*100,T,T,T,1)[2,1] #mrc
	Correlations_5min_measures[i,6] <-  rKernelCov(list(mergedfrequencies[[6]][[i]][,1]* 100, mergedfrequencies[[6]][[i]][,2]* 100), 
		cor =TRUE, makeReturns = FALSE, kernel.type = "Parzen", kernel.param  = kernel.param[[6]][i])[2,1]


}

#really slow, thats why its in its own loop. 
for(i in 1:2516){
	Correlations_5min_measures[i,4] <-  preavBPCOV(mergedfrequencies[[1]][[i]]*100,T,T,T,1)[2,1] #pbpcov
	print(sprintf("%s", i))
}



#saveRDS(Correlations_5min_measures, "Correlations_for_barplot.rds")




#------------------------------------------------univariate models


#RUN ALL 
{
Correlations_barplot <- readRDS("Correlations_for_barplot.rds")


rq_TLT_1min <- numeric() 
rq_SPY_1min <- numeric() 

trq_TLT_1min <- numeric()
trq_SPY_1min <- numeric()

#for pbpcov 
trq_TLT_1sec <- numeric()
trq_SPY_1sec <- numeric()

for(i in 1:length(mergedfrequencies[[6]])){
				 #realized quarticity
	rq_TLT_1min[i] <- rQuar(mergedfrequencies[[6]][[i]][,1] * 100)
	rq_SPY_1min[i] <- rQuar(mergedfrequencies[[6]][[i]][,2] * 100)
				#realized tripower quarticity
	trq_TLT_1min[i] <- RTQ(mergedfrequencies[[6]][[i]][,1] * 100)
	trq_SPY_1min[i] <- RTQ(mergedfrequencies[[6]][[i]][,2] * 100)

				#for pbpcov
	trq_TLT_1sec[i] <- RTQ(mergedfrequencies[[1]][[i]][,1] * 100)
	trq_SPY_1sec[i] <- RTQ(mergedfrequencies[[1]][[i]][,2] * 100)

}


rq_TLT_1min <- matrix(rq_TLT_1min)
rq_SPY_1min <- matrix(rq_SPY_1min)


rqSPYweek <- rowMeans(cbind(rq_SPY_1min, mlag(rq_SPY_1min,4,mean(rq_SPY_1min))))
sqrtrq_SPYweek <- sqrt(rqSPYweek) - mean(sqrt(rqSPYweek))

rqSPYmonth <- rowMeans(cbind(rq_SPY_1min, mlag(rq_SPY_1min,21,mean(rq_SPY_1min))))
sqrtrq_SPYmonth <- sqrt(rqSPYmonth) - mean(sqrt(rqSPYmonth))


rqTLTweek <- rowMeans(cbind(rq_TLT_1min, mlag(rq_TLT_1min,4,mean(rq_TLT_1min))))
sqrtrq_TLTweek <- sqrt(rqTLTweek) - mean(sqrt(rqTLTweek))

rqTLTmonth <- rowMeans(cbind(rq_TLT_1min, mlag(rq_TLT_1min,21,mean(rq_TLT_1min))))
sqrtrq_TLTmonth <- sqrt(rqTLTmonth) - mean(sqrt(rqTLTmonth))

#------------------HAR preparation------------------------
#RV
oneminvol_TLT <- matrix((calccov[[1]][[6]][1,1,])) * 10000  #sqrt
volday_TLT <- oneminvol_TLT
volweek_TLT <- rowMeans(cbind(oneminvol_TLT, mlag(oneminvol_TLT,4,mean(oneminvol_TLT))))
volmonth_TLT <- rowMeans(cbind(oneminvol_TLT, mlag(oneminvol_TLT,21,mean(oneminvol_TLT))))

oneminvol_SPY <- matrix((calccov[[1]][[6]][2,2,])) * 10000  #sqrt
volday_SPY <- oneminvol_SPY
volweek_SPY <- rowMeans(cbind(oneminvol_SPY, mlag(oneminvol_SPY,4,mean(oneminvol_SPY))))
volmonth_SPY <- rowMeans(cbind(oneminvol_SPY, mlag(oneminvol_SPY,21,mean(oneminvol_SPY))))

#RSpos + RSneg
voldayposvar_TLT <- matrix((calccov[[2]][[6]][1,1,])) * 10000 #sqrt
voldaynegvar_TLT <- matrix((calccov[[3]][[6]][1,1,])) * 10000 #sqrt
voldayposvar_SPY <- matrix((calccov[[2]][[6]][2,2,])) * 10000 #sqrt
voldaynegvar_SPY <- matrix((calccov[[3]][[6]][2,2,])) * 10000 #sqrt


#MRC
oneminmrcvol_TLT <- matrix((calccov[[7]][[6]][1,1,])) * 10000  #sqrt
mrcvolday_TLT <- oneminmrcvol_TLT 
mrcvolweek_TLT <- rowMeans(cbind(oneminmrcvol_TLT, mlag(oneminmrcvol_TLT,4,mean(oneminmrcvol_TLT))))
mrcvolmonth_TLT <- rowMeans(cbind(oneminmrcvol_TLT, mlag(oneminmrcvol_TLT,21,mean(oneminmrcvol_TLT))))

oneminmrcvol_SPY <- matrix((calccov[[7]][[6]][2,2,])) * 10000  #sqrt
mrcvolday_SPY <- oneminmrcvol_SPY 
mrcvolweek_SPY <- rowMeans(cbind(oneminmrcvol_SPY, mlag(oneminmrcvol_SPY,4,mean(oneminmrcvol_SPY))))
mrcvolmonth_SPY <- rowMeans(cbind(oneminmrcvol_SPY, mlag(oneminmrcvol_SPY,21,mean(oneminmrcvol_SPY))))

#MRK
oneminmrkvol_TLT <- matrix((calccov[[8]][[6]][1,1,])) * 10000  #sqrt
mrkvolday_TLT <- oneminmrkvol_TLT 
mrkvolweek_TLT <- rowMeans(cbind(oneminmrkvol_TLT, mlag(oneminmrkvol_TLT,4,mean(oneminmrkvol_TLT))))
mrkvolmonth_TLT <- rowMeans(cbind(oneminmrkvol_TLT, mlag(oneminmrkvol_TLT,21,mean(oneminmrkvol_TLT))))

oneminmrkvol_SPY <- matrix((calccov[[8]][[6]][2,2,])) * 10000  #sqrt
mrkvolday_SPY <- oneminmrkvol_SPY 
mrkvolweek_SPY <- rowMeans(cbind(oneminmrkvol_SPY, mlag(oneminmrkvol_SPY,4,mean(oneminmrkvol_SPY))))
mrkvolmonth_SPY <- rowMeans(cbind(oneminmrkvol_SPY, mlag(oneminmrkvol_SPY,21,mean(oneminmrkvol_SPY))))


#BPV
oneminbpvol_TLT <- matrix((calccov[[5]][[6]][1,1,])*10000) #sqrt
bpvolday_TLT <- oneminbpvol_TLT
bpvolweek_TLT <- rowMeans(cbind(oneminbpvol_TLT, mlag(oneminbpvol_TLT,4,mean(oneminbpvol_TLT))))
bpvolmonth_TLT <- rowMeans(cbind(oneminbpvol_TLT, mlag(oneminbpvol_TLT,21,mean(oneminbpvol_TLT))))

oneminbpvol_SPY <- matrix((calccov[[5]][[6]][2,2,])*10000) #sqrt
bpvolday_SPY <- oneminbpvol_SPY
bpvolweek_SPY <- rowMeans(cbind(oneminbpvol_SPY, mlag(oneminbpvol_SPY,4,mean(oneminbpvol_SPY))))
bpvolmonth_SPY <- rowMeans(cbind(oneminbpvol_SPY, mlag(oneminbpvol_SPY,21,mean(oneminbpvol_SPY))))

#PBPCOV
oneminpbpvol_TLT <- matrix((calccov[[6]][[1]][1,1,])*10000) #sqrt
pbpvolday_TLT <- oneminpbpvol_TLT
pbpvolweek_TLT <- rowMeans(cbind(oneminpbpvol_TLT, mlag(oneminpbpvol_TLT,4,mean(oneminpbpvol_TLT))))
pbpvolmonth_TLT <- rowMeans(cbind(oneminpbpvol_TLT, mlag(oneminpbpvol_TLT,21,mean(oneminpbpvol_TLT))))

oneminpbpvol_SPY <- matrix((calccov[[6]][[1]][2,2,])*10000) #sqrt
pbpvolday_SPY <- oneminbpvol_SPY
pbpvolweek_SPY <- rowMeans(cbind(oneminpbpvol_SPY, mlag(oneminpbpvol_SPY,4,mean(oneminpbpvol_SPY))))
pbpvolmonth_SPY <- rowMeans(cbind(oneminpbpvol_SPY, mlag(oneminpbpvol_SPY,21,mean(oneminpbpvol_SPY))))


#TCOV
onemintvvol_TLT <- matrix((calccov[[4]][[6]][1,1,])*10000) #sqrt
tvvolday_TLT <- onemintvvol_TLT
tvvolweek_TLT <- rowMeans(cbind(onemintvvol_TLT, mlag(onemintvvol_TLT,4,mean(onemintvvol_TLT))))
tvvolmonth_TLT <- rowMeans(cbind(onemintvvol_TLT, mlag(onemintvvol_TLT,21,mean(onemintvvol_TLT))))

onemintvvol_SPY <- matrix((calccov[[4]][[6]][2,2,])*10000) #sqrt
tvvolday_SPY <- onemintvvol_SPY
tvvolweek_SPY <- rowMeans(cbind(onemintvvol_SPY, mlag(onemintvvol_SPY,4,mean(onemintvvol_SPY))))
tvvolmonth_SPY <- rowMeans(cbind(onemintvvol_SPY, mlag(onemintvvol_SPY,21,mean(onemintvvol_SPY))))


#RQ and TRQ
sqrtrq_TLT_1min <- sqrt(rq_TLT_1min) - mean(sqrt(rq_TLT_1min))
sqrttrq_TLT_1min <- sqrt(trq_TLT_1min) - mean(sqrt(trq_TLT_1min))
sqrttrq_TLT_1sec <- sqrt(trq_TLT_1sec) - mean(sqrt(trq_TLT_1sec))

sqrtrq_SPY_1min <- sqrt(rq_SPY_1min) - mean(sqrt(rq_SPY_1min))
sqrttrq_SPY_1min <- sqrt(trq_SPY_1min) - mean(sqrt(trq_SPY_1min))
sqrttrq_SPY_1sec <- sqrt(trq_SPY_1sec) - mean(sqrt(trq_SPY_1sec))

#jumpparams
jumpparam_TLT <- ifelse(oneminvol_TLT - oneminbpvol_TLT>0, oneminvol_TLT - oneminbpvol_TLT, 0)
jumpparam_SPY <- ifelse(oneminvol_SPY - oneminbpvol_SPY>0, oneminvol_SPY - oneminbpvol_SPY, 0)

ones <- matrix(rep(1, 2516-22))


#-------------- PROXY

RC5min <- calccov[[1]][[7]][2,1,]/sqrt(calccov[[1]][[7]][1,1,] * calccov[[1]][[7]][2,2,])

RV5min_TLT <-  matrix(calccov[[1]][[7]][1,1,] * 10000)
RV5min_SPY <-  matrix(calccov[[1]][[7]][2,2,] * 10000)

#--------------


#FOLLOW THIS SYSTEM:
#REMEMBTER TO REGRESS ON RV5MIN!


#EVERYTHING HAVE BEEN CHECKED TO SEE IF THEY NEEDED INSANITY FILTER. ONLY FEW PLACES NEEDED INSANITY FILTER!

#---------------------DRD-HAR
#DRD-HAR-RV
HAR_TLT_RV <- lm(RV5min_TLT[23:2516] ~ volday_TLT[22:(2516-1)] + volweek_TLT[22:(2516-1)] + volmonth_TLT[22:(2516-1)])
hhat_HAR_TLT_RV_1min <- cbind(ones, volday_TLT[22:2515], volweek_TLT[22:2515], volmonth_TLT[22:2515])  %*% matrix(coef(HAR_TLT_RV))

HAR_SPY_RV <- lm(RV5min_SPY[23:2516] ~ volday_SPY[22:(2516-1)] + volweek_SPY[22:(2516-1)] + volmonth_SPY[22:(2516-1)])
hhat_HAR_SPY_RV_1min <- cbind(ones, volday_SPY[22:2515], volweek_SPY[22:2515], volmonth_SPY[22:2515])  %*% matrix(coef(HAR_SPY_RV))

DRD_HAR_RV_1min <- EstimatecorrHAR(cbind(hhat_HAR_TLT_RV_1min, hhat_HAR_SPY_RV_1min), correlation = Correlations_barplot[,1], proxy = RC5min, 1)

#DRD-HAR-MRC
HAR_TLT_MRC <- lm(RV5min_TLT[23:2516] ~ mrcvolday_TLT[22:(2516-1)] + mrcvolweek_TLT[22:(2516-1)] + mrcvolmonth_TLT[22:(2516-1)])
hhat_HAR_TLT_MRC_1min <- cbind(ones, mrcvolday_TLT[22:2515], mrcvolweek_TLT[22:2515], mrcvolmonth_TLT[22:2515])  %*% matrix(coef(HAR_TLT_MRC))

HAR_SPY_MRC <- lm(RV5min_SPY[23:2516] ~ mrcvolday_SPY[22:(2516-1)] + mrcvolweek_SPY[22:(2516-1)] + mrcvolmonth_SPY[22:(2516-1)])
hhat_HAR_SPY_MRC_1min <- cbind(ones, mrcvolday_SPY[22:2515], mrcvolweek_SPY[22:2515], mrcvolmonth_SPY[22:2515])  %*% matrix(coef(HAR_SPY_MRC))

DRD_HAR_MRC_1min <- EstimatecorrHAR(cbind(hhat_HAR_TLT_MRC_1min, hhat_HAR_SPY_MRC_1min), correlation = Correlations_barplot[,5], proxy = RC5min, 1)

#DRD-HAR-MRK
HAR_TLT_mrk <- lm(RV5min_TLT[23:2516] ~ mrkvolday_TLT[22:(2516-1)] + mrkvolweek_TLT[22:(2516-1)] + mrkvolmonth_TLT[22:(2516-1)])
hhat_HAR_TLT_mrk_1min <- cbind(ones, mrkvolday_TLT[22:2515], mrkvolweek_TLT[22:2515], mrkvolmonth_TLT[22:2515])  %*% matrix(coef(HAR_TLT_mrk))

HAR_SPY_mrk <- lm(RV5min_SPY[23:2516] ~ mrkvolday_SPY[22:(2516-1)] + mrkvolweek_SPY[22:(2516-1)] + mrkvolmonth_SPY[22:(2516-1)])
hhat_HAR_SPY_mrk_1min <- cbind(ones, mrkvolday_SPY[22:2515], mrkvolweek_SPY[22:2515], mrkvolmonth_SPY[22:2515])  %*% matrix(coef(HAR_SPY_mrk))

DRD_HAR_mrk_1min <- EstimatecorrHAR(cbind(hhat_HAR_TLT_mrk_1min, hhat_HAR_SPY_mrk_1min), correlation = Correlations_barplot[,6], proxy = RC5min, 1)
#----------------------------
#---------------------DRD-HARQ

#DRD-HARQ-RV
HARQ_TLT_RV <- lm(RV5min_TLT[23:2516] ~ volday_TLT[22:(2516-1)] +  volweek_TLT[22:(2516-1)] + volmonth_TLT[22:(2516-1)] + I(volday_TLT[22:(2516-1)] * sqrtrq_TLT_1min[22:(2516-1)]))
hhat_HARQ_TLT_RV <- cbind(ones, volday_TLT[22:2515], volweek_TLT[22:2515], volmonth_TLT[22:2515],volday_TLT[22:(2516-1)] * sqrtrq_TLT_1min[22:(2516-1)])  %*% matrix(coef(HARQ_TLT_RV))

HARQ_SPY_RV <- lm(RV5min_SPY[23:2516] ~ volday_SPY[22:(2516-1)] +  volweek_SPY[22:(2516-1)] + volmonth_SPY[22:(2516-1)] + I(volday_SPY[22:(2516-1)] * sqrtrq_SPY_1min[22:(2516-1)]))
hhat_HARQ_SPY_RV <- cbind(ones, volday_SPY[22:2515], volweek_SPY[22:2515], volmonth_SPY[22:2515],volday_SPY[22:(2516-1)] * sqrtrq_SPY_1min[22:(2516-1)])  %*% matrix(coef(HARQ_SPY_RV))

#SAME CORRELATION AS HAR MODELS. ONLY THING THAT CHANGED IS THE UNIVARIATE VOLS. 
DRD_HARQ_RV_1min <- EstimatecorrHAR(cbind(hhat_HARQ_TLT_RV, hhat_HARQ_SPY_RV), correlation = Correlations_barplot[,1], proxy = RC5min, 1)

#DRD-HARQ-MRC
HARQ_TLT_MRC <- lm(RV5min_TLT[23:2516] ~ mrcvolday_TLT[22:(2516-1)] +  mrcvolweek_TLT[22:(2516-1)] + mrcvolmonth_TLT[22:(2516-1)] + I(mrcvolday_TLT[22:(2516-1)] * sqrtrq_TLT_1min[22:(2516-1)]))
hhat_HARQ_TLT_MRC <- cbind(ones, mrcvolday_TLT[22:2515], mrcvolweek_TLT[22:2515], mrcvolmonth_TLT[22:2515], mrcvolday_TLT[22:(2516-1)] * sqrtrq_TLT_1min[22:(2516-1)])  %*% matrix(coef(HARQ_TLT_MRC))

HARQ_SPY_MRC <- lm(RV5min_SPY[23:2516] ~ mrcvolday_SPY[22:(2516-1)] +  mrcvolweek_SPY[22:(2516-1)] + mrcvolmonth_SPY[22:(2516-1)] + I(mrcvolday_SPY[22:(2516-1)] * sqrtrq_SPY_1min[22:(2516-1)]))
hhat_HARQ_SPY_MRC <- cbind(ones, mrcvolday_SPY[22:2515], mrcvolweek_SPY[22:2515], mrcvolmonth_SPY[22:2515], mrcvolday_SPY[22:(2516-1)] * sqrtrq_SPY_1min[22:(2516-1)])  %*% matrix(coef(HARQ_SPY_MRC))

#SAME CORRELATION AS HAR MODELS. ONLY THING THAT CHANGED IS THE UNIVARIATE VOLS. 
DRD_HARQ_MRC_1min <- EstimatecorrHAR(cbind(hhat_HARQ_TLT_MRC, hhat_HARQ_SPY_MRC), correlation = Correlations_barplot[,5], proxy = RC5min, 1)

#DRD-HARQ-MRK
HARQ_TLT_mrk <- lm(RV5min_TLT[23:2516] ~ mrkvolday_TLT[22:(2516-1)] +  mrkvolweek_TLT[22:(2516-1)] + mrkvolmonth_TLT[22:(2516-1)] + I(mrkvolday_TLT[22:(2516-1)] * sqrtrq_TLT_1min[22:(2516-1)]))
hhat_HARQ_TLT_mrk <- cbind(ones, mrkvolday_TLT[22:2515], mrkvolweek_TLT[22:2515], mrkvolmonth_TLT[22:2515], mrkvolday_TLT[22:(2516-1)] * sqrtrq_TLT_1min[22:(2516-1)])  %*% matrix(coef(HARQ_TLT_mrk))

HARQ_SPY_mrk <- lm(RV5min_SPY[23:2516] ~ mrkvolday_SPY[22:(2516-1)] +  mrkvolweek_SPY[22:(2516-1)] + mrkvolmonth_SPY[22:(2516-1)] + I(mrkvolday_SPY[22:(2516-1)] * sqrtrq_SPY_1min[22:(2516-1)]))
hhat_HARQ_SPY_mrk <- cbind(ones, mrkvolday_SPY[22:2515], mrkvolweek_SPY[22:2515], mrkvolmonth_SPY[22:2515], mrkvolday_SPY[22:(2516-1)] * sqrtrq_SPY_1min[22:(2516-1)])  %*% matrix(coef(HARQ_SPY_mrk))

#SAME CORRELATION AS HAR MODELS. ONLY THING THAT CHANGED IS THE UNIVARIATE VOLS. 
DRD_HARQ_mrk_1min <- EstimatecorrHAR(cbind(hhat_HARQ_TLT_mrk, hhat_HARQ_SPY_mrk), correlation = Correlations_barplot[,6], proxy = RC5min, 1)
#--------------------------------
#--------------------------------------DRD-HARQ-F

#DRD-HARQ-F-RV
HARQF_TLT_RV <- lm(RV5min_TLT[23:2516] ~ volday_TLT[22:(2516-1)] +  volweek_TLT[22:(2516-1)] + 
volmonth_TLT[22:(2516-1)] + I(volday_TLT[22:(2516-1)] * sqrtrq_TLT_1min[22:(2516-1)]) + 
I(volweek_TLT[22:(2516-1)] * sqrtrq_TLTweek[22:(2516-1)]) + I(volmonth_TLT[22:(2516-1)] * sqrtrq_TLTmonth[22:(2516-1)]))

#does not need insanity filter!
hhat_HARQF_TLT_RV <- cbind(ones, volday_TLT[22:2515], volweek_TLT[22:2515], volmonth_TLT[22:2515],
volday_TLT[22:(2516-1)] * sqrtrq_TLT_1min[22:(2516-1)], volweek_TLT[22:(2516-1)] * sqrtrq_TLTweek[22:(2516-1)], volmonth_TLT[22:(2516-1)] * sqrtrq_TLTmonth[22:(2516-1)])  %*% matrix(coef(HARQF_TLT_RV))


HARQF_SPY_RV <- lm(RV5min_SPY[23:2516] ~ volday_SPY[22:(2516-1)] +  volweek_SPY[22:(2516-1)] + 
volmonth_SPY[22:(2516-1)] + I(volday_SPY[22:(2516-1)] * sqrtrq_SPY_1min[22:(2516-1)]) + 
I(volweek_SPY[22:(2516-1)] * sqrtrq_SPYweek[22:(2516-1)]) + I(volmonth_SPY[22:(2516-1)] * sqrtrq_SPYmonth[22:(2516-1)]))

hhat_HARQF_SPY_RV <- cbind(ones, volday_SPY[22:2515], volweek_SPY[22:2515], volmonth_SPY[22:2515],
volday_SPY[22:(2516-1)] * sqrtrq_SPY_1min[22:(2516-1)], volweek_SPY[22:(2516-1)] * sqrtrq_SPYweek[22:(2516-1)], volmonth_SPY[22:(2516-1)] * sqrtrq_SPYmonth[22:(2516-1)])  %*% matrix(coef(HARQF_SPY_RV))

#needs insanity filter!
hhat_HARQF_SPY_RV <- volatility.insanity.filter(hhat_HARQF_SPY_RV, min(RV5min_SPY), max(RV5min_SPY), mean(RV5min_SPY))

DRD_HARQF_RV_1min <- EstimatecorrHAR(cbind(hhat_HARQF_TLT_RV, hhat_HARQF_SPY_RV$vol), correlation = Correlations_barplot[,1], proxy = RC5min, 1)



#DRD-HARQ-F-MRC
HARQF_TLT_MRC <- lm(RV5min_TLT[23:2516] ~ mrcvolday_TLT[22:(2516-1)] +  mrcvolweek_TLT[22:(2516-1)] + 
mrcvolmonth_TLT[22:(2516-1)] + I(mrcvolday_TLT[22:(2516-1)] * sqrtrq_TLT_1min[22:(2516-1)]) + 
I(mrcvolweek_TLT[22:(2516-1)] * sqrtrq_TLTweek[22:(2516-1)]) + I(mrcvolmonth_TLT[22:(2516-1)] * sqrtrq_TLTmonth[22:(2516-1)]))

hhat_HARQF_TLT_MRC <- cbind(ones, mrcvolday_TLT[22:2515], mrcvolweek_TLT[22:2515], mrcvolmonth_TLT[22:2515],
mrcvolday_TLT[22:(2516-1)] * sqrtrq_TLT_1min[22:(2516-1)], mrcvolweek_TLT[22:(2516-1)] * sqrtrq_TLTweek[22:(2516-1)], mrcvolmonth_TLT[22:(2516-1)] * sqrtrq_TLTmonth[22:(2516-1)])  %*% matrix(coef(HARQF_TLT_MRC))


HARQF_SPY_MRC <- lm(RV5min_SPY[23:2516] ~ mrcvolday_SPY[22:(2516-1)] +  mrcvolweek_SPY[22:(2516-1)] + 
mrcvolmonth_SPY[22:(2516-1)] + I(mrcvolday_SPY[22:(2516-1)] * sqrtrq_SPY_1min[22:(2516-1)]) + 
I(mrcvolweek_SPY[22:(2516-1)] * sqrtrq_SPYweek[22:(2516-1)]) + I(mrcvolmonth_SPY[22:(2516-1)] * sqrtrq_SPYmonth[22:(2516-1)]))

hhat_HARQF_SPY_MRC <- cbind(ones, mrcvolday_SPY[22:2515], mrcvolweek_SPY[22:2515], mrcvolmonth_SPY[22:2515],
mrcvolday_SPY[22:(2516-1)] * sqrtrq_SPY_1min[22:(2516-1)], mrcvolweek_SPY[22:(2516-1)] * sqrtrq_SPYweek[22:(2516-1)], mrcvolmonth_SPY[22:(2516-1)] * sqrtrq_SPYmonth[22:(2516-1)])  %*% matrix(coef(HARQF_SPY_MRC))

hhat_HARQF_TLT_MRC <- volatility.insanity.filter(hhat_HARQF_TLT_MRC, min(RV5min_TLT), max(RV5min_TLT), mean(RV5min_TLT))

hhat_HARQF_SPY_MRC <- volatility.insanity.filter(hhat_HARQF_SPY_MRC, min(RV5min_SPY), max(RV5min_SPY), mean(RV5min_SPY))

DRD_HARQF_MRC_1min <- EstimatecorrHAR(cbind(hhat_HARQF_TLT_MRC$vol, hhat_HARQF_SPY_MRC$vol), correlation = Correlations_barplot[,5], proxy = RC5min, 1)


#DRD-HARQ-F-MRK

HARQF_TLT_mrk <- lm(RV5min_TLT[23:2516] ~ mrkvolday_TLT[22:(2516-1)] +  mrkvolweek_TLT[22:(2516-1)] + 
mrkvolmonth_TLT[22:(2516-1)] + I(mrkvolday_TLT[22:(2516-1)] * sqrtrq_TLT_1min[22:(2516-1)]) + 
I(mrkvolweek_TLT[22:(2516-1)] * sqrtrq_TLTweek[22:(2516-1)]) + I(mrkvolmonth_TLT[22:(2516-1)] * sqrtrq_TLTmonth[22:(2516-1)]))

hhat_HARQF_TLT_mrk <- cbind(ones, mrkvolday_TLT[22:2515], mrkvolweek_TLT[22:2515], mrkvolmonth_TLT[22:2515],
mrkvolday_TLT[22:(2516-1)] * sqrtrq_TLT_1min[22:(2516-1)], mrkvolweek_TLT[22:(2516-1)] * sqrtrq_TLTweek[22:(2516-1)], mrkvolmonth_TLT[22:(2516-1)] * sqrtrq_TLTmonth[22:(2516-1)])  %*% matrix(coef(HARQF_TLT_mrk))


HARQF_SPY_mrk <- lm(RV5min_SPY[23:2516] ~ mrkvolday_SPY[22:(2516-1)] +  mrkvolweek_SPY[22:(2516-1)] + 
mrkvolmonth_SPY[22:(2516-1)] + I(mrkvolday_SPY[22:(2516-1)] * sqrtrq_SPY_1min[22:(2516-1)]) + 
I(mrkvolweek_SPY[22:(2516-1)] * sqrtrq_SPYweek[22:(2516-1)]) + I(mrkvolmonth_SPY[22:(2516-1)] * sqrtrq_SPYmonth[22:(2516-1)]))

hhat_HARQF_SPY_mrk <- cbind(ones, mrkvolday_SPY[22:2515], mrkvolweek_SPY[22:2515], mrkvolmonth_SPY[22:2515],
mrkvolday_SPY[22:(2516-1)] * sqrtrq_SPY_1min[22:(2516-1)], mrkvolweek_SPY[22:(2516-1)] * sqrtrq_SPYweek[22:(2516-1)], mrkvolmonth_SPY[22:(2516-1)] * sqrtrq_SPYmonth[22:(2516-1)])  %*% matrix(coef(HARQF_SPY_mrk))


hhat_HARQF_SPY_mrk <- volatility.insanity.filter(hhat_HARQF_SPY_mrk, min(RV5min_SPY), max(RV5min_SPY), mean(RV5min_SPY))

DRD_HARQF_mrk_1min <- EstimatecorrHAR(cbind(hhat_HARQF_TLT_mrk, hhat_HARQF_SPY_mrk$vol), correlation = Correlations_barplot[,6], proxy = RC5min, 1)

#-------------------------------------------------------------
#---------------------------DRD-HARJ

#DRD-HARJ-RV
HARJ_TLT_RV <- lm(RV5min_TLT[23:2516] ~ volday_TLT[22:(2516-1)] + volweek_TLT[22:(2516-1)] + volmonth_TLT[22:(2516-1)] + jumpparam_TLT[22:(2516-1)])
hhat_HARJ_TLT_RV <- cbind(ones, volday_TLT[22:2515], volweek_TLT[22:2515], volmonth_TLT[22:2515],jumpparam_TLT[22:(2516-1)])  %*% matrix(coef(HARJ_TLT_RV))

HARJ_SPY_RV <- lm(RV5min_SPY[23:2516] ~ volday_SPY[22:(2516-1)] + volweek_SPY[22:(2516-1)] + volmonth_SPY[22:(2516-1)] + jumpparam_SPY[22:(2516-1)])
hhat_HARJ_SPY_RV <- cbind(ones, volday_SPY[22:2515], volweek_SPY[22:2515], volmonth_SPY[22:2515],jumpparam_SPY[22:(2516-1)])  %*% matrix(coef(HARJ_SPY_RV))

DRD_HARJ_RV_1min <- EstimatecorrHAR(cbind(hhat_HARJ_TLT_RV, hhat_HARJ_SPY_RV), correlation = Correlations_barplot[,1], proxy = RC5min, 1)


#DRD-HARJ-MRC
HARJ_TLT_MRC <- lm(RV5min_TLT[23:2516] ~ mrcvolday_TLT[22:(2516-1)] + mrcvolweek_TLT[22:(2516-1)] + mrcvolmonth_TLT[22:(2516-1)] + jumpparam_TLT[22:(2516-1)])
hhat_HARJ_TLT_MRC <- cbind(ones, mrcvolday_TLT[22:2515], mrcvolweek_TLT[22:2515], mrcvolmonth_TLT[22:2515],jumpparam_TLT[22:(2516-1)])  %*% matrix(coef(HARJ_TLT_MRC))

HARJ_SPY_MRC <- lm(RV5min_SPY[23:2516] ~ mrcvolday_SPY[22:(2516-1)] + mrcvolweek_SPY[22:(2516-1)] + mrcvolmonth_SPY[22:(2516-1)] + jumpparam_SPY[22:(2516-1)])
hhat_HARJ_SPY_MRC <- cbind(ones, mrcvolday_SPY[22:2515], mrcvolweek_SPY[22:2515], mrcvolmonth_SPY[22:2515],jumpparam_SPY[22:(2516-1)])  %*% matrix(coef(HARJ_SPY_MRC))

DRD_HARJ_MRC_1min <- EstimatecorrHAR(cbind(hhat_HARJ_TLT_MRC, hhat_HARJ_SPY_MRC), correlation = Correlations_barplot[,5], proxy = RC5min, 1)

#DRD-HARJ-MRK
HARJ_TLT_mrk <- lm(RV5min_TLT[23:2516] ~ mrkvolday_TLT[22:(2516-1)] + mrkvolweek_TLT[22:(2516-1)] + mrkvolmonth_TLT[22:(2516-1)] + jumpparam_TLT[22:(2516-1)])
hhat_HARJ_TLT_mrk <- cbind(ones, mrkvolday_TLT[22:2515], mrkvolweek_TLT[22:2515], mrkvolmonth_TLT[22:2515],jumpparam_TLT[22:(2516-1)])  %*% matrix(coef(HARJ_TLT_mrk))

HARJ_SPY_mrk <- lm(RV5min_SPY[23:2516] ~ mrkvolday_SPY[22:(2516-1)] + mrkvolweek_SPY[22:(2516-1)] + mrkvolmonth_SPY[22:(2516-1)] + jumpparam_SPY[22:(2516-1)])
hhat_HARJ_SPY_mrk <- cbind(ones, mrkvolday_SPY[22:2515], mrkvolweek_SPY[22:2515], mrkvolmonth_SPY[22:2515],jumpparam_SPY[22:(2516-1)])  %*% matrix(coef(HARJ_SPY_mrk))

DRD_HARJ_mrk_1min <- EstimatecorrHAR(cbind(hhat_HARJ_TLT_mrk, hhat_HARJ_SPY_mrk), correlation = Correlations_barplot[,6], proxy = RC5min, 1)

#-----------------------------------------------------------------------
#---------------------------DRD-HARQJ

#DRD-HARQJ-RV
HARQJ_TLT_RV <- lm(RV5min_TLT[23:2516] ~ volday_TLT[22:(2516-1)] +  volweek_TLT[22:(2516-1)] + volmonth_TLT[22:(2516-1)] + I(volday_TLT[22:(2516-1)] * sqrtrq_TLT_1min[22:(2516-1)]) + jumpparam_TLT[22:(2516-1)])
hhat_HARQJ_TLT_RV <- cbind(ones, volday_TLT[22:2515], volweek_TLT[22:2515], volmonth_TLT[22:2515],volday_TLT[22:(2516-1)] * sqrtrq_TLT_1min[22:(2516-1)], jumpparam_TLT[22:(2516-1)])  %*% matrix(coef(HARQJ_TLT_RV))

HARQJ_SPY_RV <- lm(RV5min_SPY[23:2516] ~ volday_SPY[22:(2516-1)] +  volweek_SPY[22:(2516-1)] + volmonth_SPY[22:(2516-1)] + I(volday_SPY[22:(2516-1)] * sqrtrq_SPY_1min[22:(2516-1)]) + jumpparam_SPY[22:(2516-1)])
hhat_HARQJ_SPY_RV <- cbind(ones, volday_SPY[22:2515], volweek_SPY[22:2515], volmonth_SPY[22:2515],volday_SPY[22:(2516-1)] * sqrtrq_SPY_1min[22:(2516-1)], jumpparam_SPY[22:(2516-1)])  %*% matrix(coef(HARQJ_SPY_RV))

DRD_HARQJ_RV_1min <- EstimatecorrHAR(cbind(hhat_HARQJ_TLT_RV, hhat_HARQJ_SPY_RV), correlation = Correlations_barplot[,1], proxy = RC5min, 1)

#DRD-HARQJ-MRC
HARQJ_TLT_MRC <- lm(RV5min_TLT[23:2516] ~ mrcvolday_TLT[22:(2516-1)] +  mrcvolweek_TLT[22:(2516-1)] + mrcvolmonth_TLT[22:(2516-1)] + I(mrcvolday_TLT[22:(2516-1)] * sqrtrq_TLT_1min[22:(2516-1)]) + jumpparam_TLT[22:(2516-1)])
hhat_HARQJ_TLT_MRC <- cbind(ones, mrcvolday_TLT[22:2515], mrcvolweek_TLT[22:2515], mrcvolmonth_TLT[22:2515],mrcvolday_TLT[22:(2516-1)] * sqrtrq_TLT_1min[22:(2516-1)], jumpparam_TLT[22:(2516-1)])  %*% matrix(coef(HARQJ_TLT_MRC))

HARQJ_SPY_MRC <- lm(RV5min_SPY[23:2516] ~ mrcvolday_SPY[22:(2516-1)] +  mrcvolweek_SPY[22:(2516-1)] + mrcvolmonth_SPY[22:(2516-1)] + I(mrcvolday_SPY[22:(2516-1)] * sqrtrq_SPY_1min[22:(2516-1)]) + jumpparam_SPY[22:(2516-1)])
hhat_HARQJ_SPY_MRC <- cbind(ones, mrcvolday_SPY[22:2515], mrcvolweek_SPY[22:2515], mrcvolmonth_SPY[22:2515],mrcvolday_SPY[22:(2516-1)] * sqrtrq_SPY_1min[22:(2516-1)], jumpparam_SPY[22:(2516-1)])  %*% matrix(coef(HARQJ_SPY_MRC))

DRD_HARQJ_MRC_1min <- EstimatecorrHAR(cbind(hhat_HARQJ_TLT_MRC, hhat_HARQJ_SPY_MRC), correlation = Correlations_barplot[,5], proxy = RC5min, 1)

#DRD-HARQJ-mrk
HARQJ_TLT_mrk <- lm(RV5min_TLT[23:2516] ~ mrkvolday_TLT[22:(2516-1)] +  mrkvolweek_TLT[22:(2516-1)] + mrkvolmonth_TLT[22:(2516-1)] + I(mrkvolday_TLT[22:(2516-1)] * sqrtrq_TLT_1min[22:(2516-1)]) + jumpparam_TLT[22:(2516-1)])
hhat_HARQJ_TLT_mrk <- cbind(ones, mrkvolday_TLT[22:2515], mrkvolweek_TLT[22:2515], mrkvolmonth_TLT[22:2515],mrkvolday_TLT[22:(2516-1)] * sqrtrq_TLT_1min[22:(2516-1)], jumpparam_TLT[22:(2516-1)])  %*% matrix(coef(HARQJ_TLT_mrk))

HARQJ_SPY_mrk <- lm(RV5min_SPY[23:2516] ~ mrkvolday_SPY[22:(2516-1)] +  mrkvolweek_SPY[22:(2516-1)] + mrkvolmonth_SPY[22:(2516-1)] + I(mrkvolday_SPY[22:(2516-1)] * sqrtrq_SPY_1min[22:(2516-1)]) + jumpparam_SPY[22:(2516-1)])
hhat_HARQJ_SPY_mrk <- cbind(ones, mrkvolday_SPY[22:2515], mrkvolweek_SPY[22:2515], mrkvolmonth_SPY[22:2515],mrkvolday_SPY[22:(2516-1)] * sqrtrq_SPY_1min[22:(2516-1)], jumpparam_SPY[22:(2516-1)])  %*% matrix(coef(HARQJ_SPY_mrk))

DRD_HARQJ_mrk_1min <- EstimatecorrHAR(cbind(hhat_HARQJ_TLT_mrk, hhat_HARQJ_SPY_mrk), correlation = Correlations_barplot[,6], proxy = RC5min, 1)


#-----------------------------------------------------------------
#--------------------------DRD-SHAR----

#DRD-SHAR

SHAR_TLT  <- lm(RV5min_TLT[23:2516] ~ voldayposvar_TLT[22:(2516-1)] + voldaynegvar_TLT[22:(2516-1)] + volweek_TLT[22:(2516-1)] + volmonth_TLT[22:(2516-1)])
hhat_SHAR_TLT <- cbind(ones, voldayposvar_TLT[22:(2516-1)], voldaynegvar_TLT[22:(2516-1)], volweek_TLT[22:2515], volmonth_TLT[22:2515])  %*% matrix(coef(SHAR_TLT))

SHAR_SPY  <- lm(RV5min_SPY[23:2516] ~ voldayposvar_SPY[22:(2516-1)] + voldaynegvar_SPY[22:(2516-1)] + volweek_SPY[22:(2516-1)] + volmonth_SPY[22:(2516-1)])
hhat_SHAR_SPY <- cbind(ones, voldayposvar_SPY[22:(2516-1)], voldaynegvar_SPY[22:(2516-1)], volweek_SPY[22:2515], volmonth_SPY[22:2515])  %*% matrix(coef(SHAR_SPY))

hhat_SHAR_SPY <- volatility.insanity.filter(hhat_SHAR_SPY, min(RV5min_SPY), max(RV5min_SPY), mean(RV5min_SPY))


DRD_SHAR_1min <- EstimatecorrHAR(cbind(hhat_SHAR_TLT, hhat_SHAR_SPY$vol), correlation = Correlations_barplot[,1], proxy = RC5min, 1)


#----------------------------------------------------------------------------------------
#---------------------------DRD-CHAR-----------

#DRD-CHAR-TV   (TCov)
CHAR_TLT_tv <- lm(RV5min_TLT[23:2516] ~ tvvolday_TLT[22:(2516-1)] + tvvolweek_TLT[22:(2516-1)] + tvvolmonth_TLT[22:(2516-1)])
hhat_CHAR_TLT_tv <- cbind(ones, tvvolday_TLT[22:2515], tvvolweek_TLT[22:2515], tvvolmonth_TLT[22:2515])  %*% matrix(coef(CHAR_TLT_tv))

CHAR_SPY_tv <- lm(RV5min_SPY[23:2516] ~ tvvolday_SPY[22:(2516-1)] + tvvolweek_SPY[22:(2516-1)] + tvvolmonth_SPY[22:(2516-1)])
hhat_CHAR_SPY_tv <- cbind(ones, tvvolday_SPY[22:2515], tvvolweek_SPY[22:2515], tvvolmonth_SPY[22:2515])  %*% matrix(coef(CHAR_SPY_tv))

DRD_CHAR_tv_1min <- EstimatecorrHAR(cbind(hhat_CHAR_TLT_tv, hhat_CHAR_SPY_tv), correlation = Correlations_barplot[,2], proxy = RC5min, 1)


#DRD-CHAR-BPV
CHAR_TLT_BPV <- lm(RV5min_TLT[23:2516] ~ bpvolday_TLT[22:(2516-1)] + bpvolweek_TLT[22:(2516-1)] + bpvolmonth_TLT[22:(2516-1)])
hhat_CHAR_TLT_BPV <- cbind(ones, bpvolday_TLT[22:2515], bpvolweek_TLT[22:2515], bpvolmonth_TLT[22:2515])  %*% matrix(coef(CHAR_TLT_BPV))

CHAR_SPY_BPV <- lm(RV5min_SPY[23:2516] ~ bpvolday_SPY[22:(2516-1)] + bpvolweek_SPY[22:(2516-1)] + bpvolmonth_SPY[22:(2516-1)])
hhat_CHAR_SPY_BPV <- cbind(ones, bpvolday_SPY[22:2515], bpvolweek_SPY[22:2515], bpvolmonth_SPY[22:2515])  %*% matrix(coef(CHAR_SPY_BPV))

DRD_CHAR_BPV_1min <- EstimatecorrHAR(cbind(hhat_CHAR_TLT_BPV, hhat_CHAR_SPY_BPV), correlation = Correlations_barplot[,3], proxy = RC5min, 1)


#DRD-CHAR-PBPV
CHAR_TLT_pbp <- lm(RV5min_TLT[23:2516] ~ pbpvolday_TLT[22:(2516-1)] + pbpvolweek_TLT[22:(2516-1)] + pbpvolmonth_TLT[22:(2516-1)])
hhat_CHAR_TLT_pbp <- cbind(ones, pbpvolday_TLT[22:2515], pbpvolweek_TLT[22:2515], pbpvolmonth_TLT[22:2515])  %*% matrix(coef(CHAR_TLT_pbp))

CHAR_SPY_pbp <- lm(RV5min_SPY[23:2516] ~ pbpvolday_SPY[22:(2516-1)] + pbpvolweek_SPY[22:(2516-1)] + pbpvolmonth_SPY[22:(2516-1)])
hhat_CHAR_SPY_pbp <- cbind(ones, pbpvolday_SPY[22:2515], pbpvolweek_SPY[22:2515], pbpvolmonth_SPY[22:2515])  %*% matrix(coef(CHAR_SPY_pbp))

DRD_CHAR_pbp_1sec <- EstimatecorrHAR(cbind(hhat_CHAR_TLT_pbp, hhat_CHAR_SPY_pbp), correlation = Correlations_barplot[,4], proxy = RC5min, 1)


#---------------------------------------------------------------
#--------------------------------DRD-CHARQ

#DRD-CHARQ-TV
CHARQ_TLT_tv <- lm(RV5min_TLT[23:2516] ~ tvvolday_TLT[22:(2516-1)] +  tvvolweek_TLT[22:(2516-1)] + tvvolmonth_TLT[22:(2516-1)] + I(tvvolday_TLT[22:(2516-1)] * sqrttrq_TLT_1min[22:(2516-1)]))
hhat_CHARQ_TLT_tv <- cbind(ones, tvvolday_TLT[22:2515], tvvolweek_TLT[22:2515], tvvolmonth_TLT[22:2515], tvvolday_TLT[22:(2516-1)] * sqrttrq_TLT_1min[22:(2516-1)])  %*% matrix(coef(CHARQ_TLT_tv))

CHARQ_SPY_tv <- lm(RV5min_SPY[23:2516] ~ tvvolday_SPY[22:(2516-1)] +  tvvolweek_SPY[22:(2516-1)] + tvvolmonth_SPY[22:(2516-1)] + I(tvvolday_SPY[22:(2516-1)] * sqrttrq_SPY_1min[22:(2516-1)]))
hhat_CHARQ_SPY_tv <- cbind(ones, tvvolday_SPY[22:2515], tvvolweek_SPY[22:2515], tvvolmonth_SPY[22:2515], tvvolday_SPY[22:(2516-1)] * sqrttrq_SPY_1min[22:(2516-1)])  %*% matrix(coef(CHARQ_SPY_tv))

DRD_CHARQ_tv_1min <- EstimatecorrHAR(cbind(hhat_CHARQ_TLT_tv, hhat_CHARQ_SPY_tv), correlation = Correlations_barplot[,2], proxy = RC5min, 1)


#DRD-CHARQ-BPV
CHARQ_TLT_BPV <- lm(RV5min_TLT[23:2516] ~ bpvolday_TLT[22:(2516-1)] +  bpvolweek_TLT[22:(2516-1)] + bpvolmonth_TLT[22:(2516-1)] + I(bpvolday_TLT[22:(2516-1)] * sqrttrq_TLT_1min[22:(2516-1)]))
hhat_CHARQ_TLT_BPV <- cbind(ones, bpvolday_TLT[22:2515], bpvolweek_TLT[22:2515], bpvolmonth_TLT[22:2515], bpvolday_TLT[22:(2516-1)] * sqrttrq_TLT_1min[22:(2516-1)])  %*% matrix(coef(CHARQ_TLT_BPV))

CHARQ_SPY_BPV <- lm(RV5min_SPY[23:2516] ~ bpvolday_SPY[22:(2516-1)] +  bpvolweek_SPY[22:(2516-1)] + bpvolmonth_SPY[22:(2516-1)] + I(bpvolday_SPY[22:(2516-1)] * sqrttrq_SPY_1min[22:(2516-1)]))
hhat_CHARQ_SPY_BPV <- cbind(ones, bpvolday_SPY[22:2515], bpvolweek_SPY[22:2515], bpvolmonth_SPY[22:2515], bpvolday_SPY[22:(2516-1)] * sqrttrq_SPY_1min[22:(2516-1)])  %*% matrix(coef(CHARQ_SPY_BPV))

DRD_CHARQ_BPV_1min <- EstimatecorrHAR(cbind(hhat_CHARQ_TLT_BPV, hhat_CHARQ_SPY_BPV), correlation = Correlations_barplot[,3], proxy = RC5min, 1)

#DRD-CHARQ-PBP
CHARQ_TLT_pbp <- lm(RV5min_TLT[23:2516] ~ pbpvolday_TLT[22:(2516-1)] +  pbpvolweek_TLT[22:(2516-1)] + pbpvolmonth_TLT[22:(2516-1)] + I(pbpvolday_TLT[22:(2516-1)] * sqrttrq_TLT_1sec[22:(2516-1)]))
hhat_CHARQ_TLT_pbp <- cbind(ones, pbpvolday_TLT[22:2515], pbpvolweek_TLT[22:2515], pbpvolmonth_TLT[22:2515], pbpvolday_TLT[22:(2516-1)] * sqrttrq_TLT_1sec[22:(2516-1)])  %*% matrix(coef(CHARQ_TLT_pbp))

CHARQ_SPY_pbp <- lm(RV5min_SPY[23:2516] ~ pbpvolday_SPY[22:(2516-1)] +  pbpvolweek_SPY[22:(2516-1)] + pbpvolmonth_SPY[22:(2516-1)] + I(pbpvolday_SPY[22:(2516-1)] * sqrttrq_SPY_1sec[22:(2516-1)]))
hhat_CHARQ_SPY_pbp <- cbind(ones, pbpvolday_SPY[22:2515], pbpvolweek_SPY[22:2515], pbpvolmonth_SPY[22:2515], pbpvolday_SPY[22:(2516-1)] * sqrttrq_SPY_1sec[22:(2516-1)])  %*% matrix(coef(CHARQ_SPY_pbp))

DRD_CHARQ_pbp_1sec <- EstimatecorrHAR(cbind(hhat_CHARQ_TLT_pbp, hhat_CHARQ_SPY_pbp), correlation = Correlations_barplot[,4], proxy = RC5min, 1)

} #END OF RUN ALL

#you can compute QLIKES IN THE END, WHEN YOU HAVE ALL THE ESTIMATES FOR DRD-TYPE MODELS. 

Qlikes_barplot <- matrix(0L, nrow = 2494, ncol = 22)

for(i in 1:2494){

	Qlikes_barplot[i, 1] <- QLIKE(DRD_HAR_RV_1min$vSigma2[,,i], calccov[[1]][[7]][,,22+i]*10000, 2)
	Qlikes_barplot[i, 2] <- QLIKE(DRD_HAR_MRC_1min$vSigma2[,,i], calccov[[1]][[7]][,,22+i]*10000, 2)
	Qlikes_barplot[i, 3] <- QLIKE(DRD_HAR_mrk_1min$vSigma2[,,i], calccov[[1]][[7]][,,22+i]*10000, 2)
	Qlikes_barplot[i, 4] <- QLIKE(DRD_HARQ_RV_1min$vSigma2[,,i], calccov[[1]][[7]][,,22+i]*10000, 2)
	Qlikes_barplot[i, 5] <- QLIKE(DRD_HARQ_MRC_1min$vSigma2[,,i], calccov[[1]][[7]][,,22+i]*10000, 2)
	Qlikes_barplot[i, 6] <- QLIKE(DRD_HARQ_mrk_1min$vSigma2[,,i], calccov[[1]][[7]][,,22+i]*10000, 2)
	Qlikes_barplot[i, 7] <- QLIKE(DRD_HARQF_RV_1min$vSigma2[,,i], calccov[[1]][[7]][,,22+i]*10000, 2)
	Qlikes_barplot[i, 8] <- QLIKE(DRD_HARQF_MRC_1min$vSigma2[,,i], calccov[[1]][[7]][,,22+i]*10000, 2)
	Qlikes_barplot[i, 9] <- QLIKE(DRD_HARQF_mrk_1min$vSigma2[,,i], calccov[[1]][[7]][,,22+i]*10000, 2)
	Qlikes_barplot[i, 10] <- QLIKE(DRD_HARJ_RV_1min$vSigma2[,,i], calccov[[1]][[7]][,,22+i]*10000, 2)
	Qlikes_barplot[i, 11] <- QLIKE(DRD_HARJ_MRC_1min$vSigma2[,,i], calccov[[1]][[7]][,,22+i]*10000, 2)
	Qlikes_barplot[i, 12] <- QLIKE(DRD_HARJ_mrk_1min$vSigma2[,,i], calccov[[1]][[7]][,,22+i]*10000, 2)
	Qlikes_barplot[i, 13] <- QLIKE(DRD_HARQJ_RV_1min$vSigma2[,,i], calccov[[1]][[7]][,,22+i]*10000, 2)
	Qlikes_barplot[i, 14] <- QLIKE(DRD_HARQJ_MRC_1min$vSigma2[,,i], calccov[[1]][[7]][,,22+i]*10000, 2)
	Qlikes_barplot[i, 15] <- QLIKE(DRD_HARQJ_mrk_1min$vSigma2[,,i], calccov[[1]][[7]][,,22+i]*10000, 2)
	Qlikes_barplot[i, 16] <- QLIKE(DRD_SHAR_1min$vSigma2[,,i], calccov[[1]][[7]][,,22+i]*10000, 2)
	Qlikes_barplot[i, 17] <- QLIKE(DRD_CHAR_tv_1min$vSigma2[,,i], calccov[[1]][[7]][,,22+i]*10000, 2)
	Qlikes_barplot[i, 18] <- QLIKE(DRD_CHAR_BPV_1min$vSigma2[,,i], calccov[[1]][[7]][,,22+i]*10000, 2)
	Qlikes_barplot[i, 19] <- QLIKE(DRD_CHAR_pbp_1sec$vSigma2[,,i], calccov[[1]][[7]][,,22+i]*10000, 2)
	Qlikes_barplot[i, 20] <- QLIKE(DRD_CHARQ_tv_1min$vSigma2[,,i], calccov[[1]][[7]][,,22+i]*10000, 2)
	Qlikes_barplot[i, 21] <- QLIKE(DRD_CHARQ_BPV_1min$vSigma2[,,i], calccov[[1]][[7]][,,22+i]*10000, 2)
	Qlikes_barplot[i, 22] <- QLIKE(DRD_CHARQ_pbp_1sec$vSigma2[,,i], calccov[[1]][[7]][,,22+i]*10000, 2)



}

colMeans(Qlikes_barplot, na.rm = T)

barplot_names <- c("DRD-HAR-RCov", "DRD-HAR-MRC", "DRD-HAR-MRK", "DRD-HARQ-RCov", "DRD-HARQ-MRC", "DRD-HARQ-MRK",
	 "DRD-HARQF-RCov", "DRD-HARQF-MRC", "DRD-HARQF-MRK",  "DRD-HARJ-RCov", "DRD-HARJ-MRC", "DRD-HARJ-MRK",
	"DRD-HARQJ-RCov", "DRD-HARQJ-MRC", "DRD-HARQJ-MRK", "DRD-SHAR", "DRD-CHAR-TCov", "DRD-CHAR-BPCov", 
	"DRD-CHAR-PBPCov", "DRD-CHARQ-TCov", "DRD-CHARQ-BPCov", "DRD-CHARQ-PBPCov")


#You have rescaled everything so the loss at max goes to 1.5. The correct values are in another table. 

testset <- data.frame(realized = barplot_names, measures = colMeans(Qlikes_barplot, na.rm = T))

library(reshape2)
testset <- melt(testset, id.vars = "realized")
#testset <- testset[c(1:5, 26:30, 51:55, 76:80, 101:105), ] 
#testset <-testset[c(1:3, 13:15, 25:27), ]



Position <- factor(testset$realized, level = barplot_names, labels = barplot_names)


p3 <- ggplot(data = testset, aes(x = Position, y=value)) + geom_bar( aes(fill = Position), stat = "identity", position = "dodge") +
ylab("QLIKE Loss") + theme(plot.title = element_text(hjust = 0.5, face = "bold")) + 
geom_text(aes(label=round(testset$value,3)), vjust=1.6, color="white", size=3) + scale_fill_discrete(name = "DRD-type models") + 
theme(legend.title = element_text(size = 14), legend.text = element_text(size = 12), 
	axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())  



#Combining all 3 plots:
library(gridExtra)
layout <- rbind(c(1,1), c(2,3))

p4 <- grid.arrange(p3, p2, p1, ncol=2, layout_matrix = layout)

ggsave(plot = p4, "barplot_signatureplot_autocorrelation.eps", device = "eps")


##############################################################################################################################
#
#
#													RiskMetrics on 5-minute frequency
#
#
##############################################################################################################################


ewma.filter.realized <- function(cov, half.life, correlation = F, lambda = NULL){

	iT <- length(cov[1,1,])

	covar <- array(0L, dim = c(2,2, iT))

	#unconditional covariance for first month.
	sampleR <- matrix(
	c(mean(cov[1,1,1:100], na.rm = T),
	  mean(cov[2,1,1:100], na.rm = T),
	  mean(cov[2,1,1:100], na.rm = T),
	  mean(cov[2,2,1:100], na.rm = T)),
	  ncol=2, nrow=2) 
	
	covar[,,1] <- sampleR


	#calculation for half-life: lmd^t = 0.5 (half-life def) <=> t = ln(0.5)/ln(lmd) <=> lmd = exp(ln(0.5)/t)
	if(is.null(lambda)){

		lambda <- exp(log(0.5)/half.life)

	}

	print(sprintf("Lambda: %s", lambda))

	if(is.null(half.life)){

		half.life <- log(0.5)/log(lambda)

	}
	
	print(sprintf("half life: %s", half.life))

	d <- 2  # asset dimension

	for(i in 2:iT){

	covar[,,i] <-  lambda * covar[,,i-1] + (1-lambda) * cov[,,i-1]

	} 

	if(correlation){

		corr <- array(0L, dim = c(2,2, iT))
		#corr[,,1] <- cor(x[1:100, ])

		for(i in 1:iT){

		d <- sqrt(diag(covar[,,i]))
		d <- diag(d, ncol(covar), ncol(covar))

		corr[,,i] <- inv(d) %*% covar[,,i] %*% inv(d) 
		}

		return(corr)
	}

return(covar)

}


test2 <- ewma.filter.realized(calccov[[1]][[7]], 12)

min.qlike.riskmetrics <- function(cov, half.life, proxy, correlation = F, lambda){

	filter <- ewma.filter.realized(cov, half.life,F, lambda)

	qlikes.riskmetrics <- numeric()
	for(i in 1:2515){

		qlikes.riskmetrics[i] <- QLIKE(filter[,,i], proxy[,,i+1], 2)

	}

 return(mean(qlikes.riskmetrics))

}



qlikes.riskmetrics <- numeric()
for(i in 1:10){
	qlikes.riskmetrics[i] <- min.qlike.riskmetrics(calccov[[1]][[i]], NULL, calccov[[1]][[7]], F, 0.94, T)
}



min.qlike.riskmetrics(calccov[[1]][[7]], NULL, calccov[[1]][[7]], correlation = F, 0.94)

##############################################################################################################################
#
#
#													DISCARDED STUFF
#
#
##############################################################################################################################

EstimatecorrHAR2 <- function(returns, variances, proxy, trace=1, ineqfun = ineqconstraint, ineqLB = 0.00, ineqUB = 0.9999){
  #mEta is standardized residuals from univariate models!
  #variances should be produced by the univariate models and not come from realized measures!
  #correlation should be found via eg. five min samples but adhere to the same frequency as variances from HAR models. 
  mEta <- numeric()
  for(i in 1:nrow(variances)){

  	mEta2 <- standardized.res.intraday(returns[[21+i]], variances[i, ])
  	mEta[i] <- realCov(mEta2, T)[2,1]

  }
  dA = 0.45 
  dB  = 0.09185
  dC = 0.26
  ## vector of starting parameters
  vPar = c(dA, dB, dC)


  mEta <- as.matrix(mEta)
  corrday <- mEta
  corrweek <- rowMeans(cbind(mEta, mlag(mEta,4,mean(mEta))))
  corrmonth <- rowMeans(cbind(mEta, mlag(mEta,21,mean(mEta)))) 

  data <- list(proxy[23:2516], corrday, corrweek, corrmonth)


  optimizer = solnp(vPar, fun = min.RSS, data = data, 
                    ineqfun = ineqfun, #the inequality constraint
                    ineqLB  = ineqLB, ## the inequality lower bound
                    ineqUB = ineqUB, ## the inequality lower bound, i.e. 0.0 <= a + b + c < 0.9999
                    ## lower and upper bounds for all parameters
                    LB = c(0, 0, 0), UB = c(0.9999, 0.9999, 0.9999),
                    control = list(tol = 1e-18, outer.iter = 800, inner.iter = 1200, delta = 1e-7,
                    trace = trace))
                    

  params <- optimizer$pars

  min <- tail(optimizer$values, 1)

  hessian <- optimizer$hessian

  scores <- matrix(0L, nrow=(nrow(variances)), ncol = 3)

  step <- 1e-5 * vPar

  for(i in 1:length(step)){

	h <- step[i]
    delta <- rep(0, length(vPar))
    delta[i] <- h
																
	loglikeminus <- minimizingfunc(mEta, vPar-delta)$minis
	loglikeplus <- minimizingfunc(mEta, vPar+delta)$minis

	scores[,i] <- (loglikeplus - loglikeminus)/(2*h)

  }

  J <- (t(scores) %*% scores)/length(mEta[[1]])

  I <- optimizer$hessian/length(mEta[[1]])

  I <- solve(I)[-1 ,-1]

  vars <- (I * J * I)/length(mEta[[1]])
  
  rse <- sqrt(diag(vars))

  t.stat <- vPar/rse


  #calculating covariances: 
  
  hhatcorrHAR <- cbind(corrday[22:2515], corrweek[22:2515], corrmonth[22:2515]) %*% matrix(params)

  covs <- hhatcorrHAR * sqrt(variances[,1]) * sqrt(variances[,2])

  vSigma2 <- array(0L, dim = c(2,2, (nrow(variances))))

  for(i in 1:(nrow(variances))){

    vSigma2[,,i] <- matrix(c(variances[i,1], covs[i], covs[i], variances[i,2]), ncol=2, nrow=2) 

  }

	  #R-squared

  Rsquared <- 1-var(proxy[23:2516] - hhatcorrHAR)/var(proxy[23:2516])
  



  lOut <- list()

  lOut[["vPar"]] <- params
  lOut[["vSigma2"]] <- vSigma2
  lOut[["R2"]] <- Rsquared
  lOut[["estcor"]] <- hhatcorrHAR
  lOut[["MSE"]] <- min/2516 
  lOut[["rse"]] <- rse
  lOut[["hessian"]] <- hessian
  lOut[["mEta"]] <- mEta
  lOut[["Tstats"]] <- t.stat

  return(lOut)

 }


HARQ_TLT_RV <- lm(RV5min_TLT[23:2516] ~ volday_TLT[22:(2516-1)] +  volweek_TLT[22:(2516-1)] + volmonth_TLT[22:(2516-1)] + I(volday_TLT[22:(2516-1)] * sqrtrq_TLT[22:(2516-1)]))
hhat_HARQ_TLT_RV <- cbind(ones, volday_TLT[22:2515], volweek_TLT[22:2515], volmonth_TLT[22:2515],volday_TLT[22:(2516-1)] * sqrtrq_TLT_1min[22:(2516-1)])  %*% matrix(coef(HARQ_TLT_RV))

HARQ_SPY_RV <- lm(RV5min_SPY[23:2516] ~ volday_SPY[22:(2516-1)] +  volweek_SPY[22:(2516-1)] + volmonth_SPY[22:(2516-1)] + I(volday_SPY[22:(2516-1)] * sqrtrq_SPY[22:(2516-1)]))
hhat_HARQ_SPY_RV <- cbind(ones, volday_SPY[22:2515], volweek_SPY[22:2515], volmonth_SPY[22:2515],volday_SPY[22:(2516-1)] * sqrtrq_SPY_1min[22:(2516-1)])  %*% matrix(coef(HARQ_SPY_RV))

DRD_HARQ_RV_1min <- EstimatecorrHAR2(mergedfrequencies[[6]], cbind(hhat_HARQ_TLT_RV, hhat_HARQ_SPY_RV), proxy = RC5min, 1)





#-------------------------------correlation HAR across minute frequencies -------------------------------

HARCor_freq <- list()

fiveminvol_TLT <- matrix((calccov[[1]][[7]][1,1,])) * 10000
fiveminvol_SPY <- matrix((calccov[[1]][[7]][2,2,]*10000)) #sqrt

for(i in 6:9){
	vol_TLT <- matrix((calccov[[1]][[i]][1,1,])) * 10000  #sqrt
	volday_TLT <- vol_TLT
	volweek_TLT <- rowMeans(cbind(vol_TLT, mlag(vol_TLT,4,mean(vol_TLT))))
	volmonth_TLT <- rowMeans(cbind(vol_TLT, mlag(vol_TLT,21,mean(vol_TLT))))

	#RV
	vol_SPY <- matrix((calccov[[1]][[i]][2,2,]*10000)) #sqrt
	volday_SPY <- vol_SPY
	volweek_SPY <- rowMeans(cbind(vol_SPY, mlag(vol_SPY,4,mean(vol_SPY)))) #sqrt
	volmonth_SPY <- rowMeans(cbind(vol_SPY, mlag(vol_SPY,21,mean(vol_SPY)))) #sqrt


	HAR_TLT <- lm(fiveminvol_TLT[23:2516] ~ volday_TLT[22:(2516-1)] + volweek_TLT[22:(2516-1)] + volmonth_TLT[22:(2516-1)])
	hhat_HAR_TLT <- cbind(ones, volday_TLT[22:2515], volweek_TLT[22:2515], volmonth_TLT[22:2515])  %*% matrix(coef(HAR_TLT))

	HAR_SPY <- lm(fiveminvol_SPY[23:2516] ~ volday_SPY[22:(2516-1)] + volweek_SPY[22:(2516-1)] + volmonth_SPY[22:(2516-1)])
	hhat_HAR_SPY <- cbind(ones, volday_SPY[22:2515], volweek_SPY[22:2515], volmonth_SPY[22:2515])  %*% matrix(coef(HAR_SPY))

	variances <- cbind(hhat_HAR_TLT, hhat_HAR_SPY)


	HARCor_freq[[i]] <- EstimatecorrHAR(dailyretotc[22:2515]*100, variances, fivemincorr_allfreq[[i]][2,1,], 1)

}

options(digits = 22)
vol_TLT <- matrix((calccov[[1]][[7]][1,1,])) * 10000  #sqrt
	volday_TLT <- vol_TLT
	volweek_TLT <- rowMeans(cbind(vol_TLT, mlag(vol_TLT,4,mean(vol_TLT))))
	volmonth_TLT <- rowMeans(cbind(vol_TLT, mlag(vol_TLT,21,mean(vol_TLT))))

	#RV
	vol_SPY <- matrix((calccov[[1]][[7]][2,2,]*10000)) #sqrt
	volday_SPY <- vol_SPY
	volweek_SPY <- rowMeans(cbind(vol_SPY, mlag(vol_SPY,4,mean(vol_SPY)))) #sqrt
	volmonth_SPY <- rowMeans(cbind(vol_SPY, mlag(vol_SPY,21,mean(vol_SPY)))) #sqrt


	HAR_TLT <- lm(fiveminvol_TLT[23:2516] ~ volday_TLT[22:(2516-1)] + volweek_TLT[22:(2516-1)] + volmonth_TLT[22:(2516-1)])
	hhat_HAR_TLT <- cbind(ones, volday_TLT[22:2515], volweek_TLT[22:2515], volmonth_TLT[22:2515])  %*% matrix(coef(HAR_TLT))

	#ttt2 <- hhat_HAR_TLT

	HAR_SPY <- lm(fiveminvol_SPY[23:2516] ~ volday_SPY[22:(2516-1)] + volweek_SPY[22:(2516-1)] + volmonth_SPY[22:(2516-1)])
	hhat_HAR_SPY <- cbind(ones, volday_SPY[22:2515], volweek_SPY[22:2515], volmonth_SPY[22:2515])  %*% matrix(coef(HAR_SPY))

	variances <- cbind(hhat_HAR_TLT*100, hhat_HAR_SPY*0.02)


	HARCor_freq <- EstimatecorrHAR(dailyretotc[22:2515]*100, sqrt(variances), fivemincorr_allfreq[[i]][2,1,], 1)


	ts.plot(sqrt(ttt2) - sqrt(hhat_HAR_TLT))



#-----------------test on how different estimates will be when removing first 5 min each day--------------

data_SPY <- readRDS("dataSPY.rds")
#adhoc deletion

for(i in 1:2516){

data_SPY[[i]] <- data_SPY[[i]][-c(1:280), ]

}


fiveminpricetest <- list()

for(i in 1:2516){

	fiveminpricetest[[i]] <- aggregatets(data_SPY[[i]], on="minutes", k=5)
}

logretSPY <- lapply(fiveminpricetest, function(x) (diff(log(x)))[-1] * 100)


rvarspy <- numeric()
for(i in 1:2516){
rvarspy[i] <- sqrt(rCov(logretSPY[[i]]))
}


vv <- (matrix(rvarspy))
vv2 <- rowMeans(cbind(vv, mlag(vv,4,mean(vv))))
vv3 <- rowMeans(cbind(vv, mlag(vv,21,mean(vv))))

hartest <- lm(rvarspy[23:2516] ~ vv[22:2515] + vv2[22:2515] +vv3[22:2515])
hhat_HAR_SPY <- cbind(ones, vv[22:2515], vv2[22:2515], vv3[22:2515])  %*% matrix(coef(hartest))
QLIKE_HAR_SPY <- colMeans(QLIKE(hhat_HAR_SPY, rvarspy[23:2516], 1))

#----------------Below is former validity check of own HAR models with package:-------------------------------------------

#HAREstimate is from HArmodel package.

library(HARModel)

HAREstimate(RM, BPV = NULL, RQ = NULL, periods = c(1,5,22),
periodsJ = NULL, periodsRQ = NULL, type = "HAR",
insanityFilter = TRUE, h = 1)


#TLT
HAR_5min_TLT <- HAREstimate(sqrt(calccov[[1]][[7]][1,1,])*100)
HARQ_5min_TLT <- HAREstimate(sqrt(calccov[[1]][[7]][1,1,])*100, type = "HARQ", RQ = rq_TLT, periodsRQ = 1)
HARQF_5min_TLT <- HAREstimate(sqrt(calccov[[1]][[7]][1,1,])*100, type = "HARQ", RQ = rq_TLT, periodsRQ = c(1, 5, 22))
HARJ_5min_TLT <-HAREstimate(sqrt(calccov[[1]][[7]][1,1,])*100, BPV = sqrt(calccov[[5]][[7]][1,1,]) * 100,
	type ="HARJ", periodsJ = 1)
HARQJ_5min_TLT <- HAREstimate(sqrt(calccov[[1]][[7]][1,1,])*100, BPV = sqrt(calccov[[5]][[7]][1,1,]) * 100,
	RQ = rq_TLT, type ="HARQ-J", periodsJ = 1, periodsRQ = 1)
CHAR_5min_TLT <- HAREstimate(sqrt(calccov[[1]][[7]][1,1,])*100, BPV = sqrt(calccov[[5]][[7]][1,1,]) * 100,
	type ="CHAR")
CHARQ_5min_TLT <- HAREstimate(RM = sqrt(calccov[[1]][[7]][1,1,])*100, BPV = sqrt(calccov[[5]][[7]][1,1,]) * 100,
	RQ = trq_TLT, type ="CHARQ", periodsRQ = 1)

#SHAR:
fiveminvol_TLT <- matrix(sqrt(calccov[[1]][[7]][1,1,])*100)
voldayposvar <- matrix(sqrt(calccov[[2]][[7]][1,1,])) * 100
voldaynegvar <- matrix(sqrt(calccov[[3]][[7]][1,1,])) * 100
volweek <- rowMeans(cbind(fiveminvol_TLT, mlag(fiveminvol_TLT,4,mean(fiveminvol_TLT))))
volmonth <- rowMeans(cbind(fiveminvol_TLT, mlag(fiveminvol_TLT,21,mean(fiveminvol_TLT))))

SHAR_5min_TLT <- lm(fiveminvol_TLT[23:2516] ~ voldayposvar[22:(2516-1)] + voldaynegvar[22:(2516-1)] + volweek[22:(2516-1)] + volmonth[22:(2516-1)])



#SPY
HAR_5min_SPY <- HAREstimate(sqrt(calccov[[1]][[7]][2,2,]*10000)
HARQ_5min_SPY <- HAREstimate(sqrt(calccov[[1]][[7]][2,2,]*10000), type = "HARQ", RQ = rq_SPY, periodsRQ = 1)
HARQF_5min_SPY <- HAREstimate(sqrt(calccov[[1]][[7]][2,2,])*100, type = "HARQ", RQ = rq_SPY, periodsRQ = c(1, 5, 22))
HARJ_5min_SPY <-HAREstimate(sqrt(calccov[[1]][[7]][2,2,])*100, BPV = sqrt(calccov[[5]][[7]][2,2,]) * 100,
	type ="HARJ", periodsJ = 1)
HARQJ_5min_SPY <- HAREstimate(sqrt(calccov[[1]][[7]][2,2,])*100, BPV = sqrt(calccov[[5]][[7]][2,2,]) * 100,
	RQ = rq_SPY, type ="HARQ-J", periodsJ = 1, periodsRQ = 1)
CHAR_5min_SPY <- HAREstimate(sqrt(calccov[[1]][[7]][2,2,])*100, BPV = sqrt(calccov[[5]][[7]][2,2,]) * 100,
	type ="CHAR")
CHARQ_5min_SPY <- HAREstimate(RM = sqrt(calccov[[1]][[7]][2,2,])*100, BPV = sqrt(calccov[[5]][[7]][2,2,]) * 100,
	RQ = trq_SPY, type ="CHARQ", periodsRQ = 1)
#SHAR:
fiveminvol_SPY <- matrix(sqrt(calccov[[1]][[7]][2,2,])*100)
voldayposvar_SPY <- matrix(sqrt(calccov[[2]][[7]][2,2,])) * 100
voldaynegvar_SPY <- matrix(sqrt(calccov[[3]][[7]][2,2,])) * 100
volweek <- rowMeans(cbind(fiveminvol_SPY, mlag(fiveminvol_SPY,4,mean(fiveminvol_SPY))))
volmonth <- rowMeans(cbind(fiveminvol_SPY, mlag(fiveminvol_SPY,21,mean(fiveminvol_SPY))))

SHAR_5min_SPY <- lm(fiveminvol_SPY[23:2516] ~ voldayposvar_SPY[22:(2516-1)] + voldaynegvar_SPY[22:(2516-1)] + volweek[22:(2516-1)] + volmonth[22:(2516-1)])





fiveminvol <- matrix(sqrt(calccov[[1]][[7]][2,2,])) * 100
volday <- fiveminvol
volweek <- rowMeans(cbind(fiveminvol, mlag(fiveminvol,4,mean(fiveminvol))))
volmonth <- rowMeans(cbind(fiveminvol, mlag(fiveminvol,21,mean(fiveminvol))))


#needs 23 data points to initialize. This is in agreement with HARestimate from HARmodel library when removing intercept part.
testerCHAR <- lm(sqrt(calccov[[1]][[7]][1,1,23:2516])*100 ~ 1+
	volday[22:(2516-1)] + volweek[22:(2516-1)] + volmonth[22:(2516-1)])

summary(testerCHAR)

sqrtrqtlt <- sqrt(rq_SPY) - mean(sqrt(rq_SPY))

testerHARQ <- lm(sqrt(calccov[[1]][[7]][2,2,23:2516])*100 ~ 1+
	volday[22:(2516-1)] + volweek[22:(2516-1)] + volmonth[22:(2516-1)] + I(volday[22:(2516-1)]*sqrtrqtlt[22:(2516-1)]))


