setwd('C:/Users/Emil/Desktop/Correlation check')

#TLT <- read.csv('daily_TLT.csv', header = T) 
#SPY <- read.csv('daily_SPY.csv', header = T) 
#TIP <- read.csv('daily_TIP.csv', header = T) 
#LQD2 <- read.csv('daily_LQD.csv', header = T) 


#Functions:

library(matlib)
realCov <- function(matrix, correlation = FALSE){
  # Realized Covariance estimator for 1 day. 
  #
  # Args: 
  #
  #   matrix: A matrix of returns. Nxk 
  #
  #	  N: number of intraday returns, k: number of asset dimension, T: number of days. 
  realcovariance <- matrix()
  	
  realcovariance <- t(matrix) %*% (matrix)

  if(correlation){

		d <- sqrt(diag(realcovariance))
		d <- diag(d, ncol(realcovariance), ncol(realcovariance))

		corr <- inv(d) %*% realcovariance %*% inv(d) 

		return(corr)
	}
  
 return( realcovariance ) 
}

#half-life: A 1 year estimation window with a 3 month half-life indicate that half.life = 1/4, maybe? 
ewma.filter <- function(x, half.life, correlation = F, lambda = NULL){

	iT <- nrow(x)

	covar <- array(0L, dim = c(ncol(x),ncol(x), iT))

	#unconditional covariance for first month.
	covar[,,1] <- cov(x[1:100, ])


	#calculation for half-life: lmd^t = 0.5 (half-life def) <=> t = ln(0.5)/ln(lmd) <=> lmd = exp(ln(0.5)/t)
	if(is.null(lambda)){

		lambda <- exp(log(0.5)/half.life)

	}

	print(sprintf("Lambda: %s", lambda))
	half.life <- log(0.5)/log(lambda)
	print(sprintf("half life: %s", half.life))

	for(i in 2:iT){


	covar[,,i] <-  lambda * covar[,,i-1] + (1-lambda) * t(x[i-1,])%*%x[i-1,]

	} 

	if(correlation){

		corr <- array(0L, dim = c(ncol(x),ncol(x), iT))
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


ewma.filter2006 <- function(data, correlation = F, tau1 = NULL, rho = NULL){

	#Reference values found in RiskMetrics 2006:
	#Rho and tau1 controls the half-life construction in the 1993 model. 
    #set Rho =1 and vary tau1 to vary the half-life. 
	if(is.null(tau1)){
    	tau1 <- 4
    }
    
    if(is.null(rho)){
    	rho <- sqrt(2)
	}

	tau0 <- 1560
    taumax <- 512
    if(is.null(tau1) && is.null(rho)){

    	kmax   <- round(log(taumax/tau1)/log(rho))

    }
    else{kmax <- 14 }
	#The one business day horizon at which
	#the data are available is denoted by dt.
	dt <- 1
	tau <- numeric(kmax)
	mu <- numeric(kmax)
	iT <- nrow(data)

	vol <- array(0L, dim = c(ncol(data), ncol(data), iT))
	volk <- list()
	weights <- numeric(kmax)


	for(k in 1:kmax){
		#now it is initialized using sample covariance
		vol[,,1] <- cov(data[1:100, ])
		tau[k] <- tau1 * rho^(k-1)  
		mu[k] <- exp(-dt/tau[k])
		weights[k] <- 1-log(tau[k]/tau0)
		#Calculate recursion of kth EWMA.

		#need to initialize the vol formula. He initializes it using a backcast??

		for(t in 2:iT){

			vol[,,t] <- mu[k] * vol[,,t-1] + (1 - mu[k]) *  t(data[t-1,])%*%data[t-1, ]

		}

		#saves kth recursion in list. 
		volk[[k]] <- vol 
		

	}
	#we pre-calculated the weights and made them sum to 1 by normalizing by the sum of the weights:
	weigths <- weights / sum(weights)
	
	#efficient vol is now a weighted sum of the k EWMA:

	effvol <- array(0L, dim=c(ncol(data), ncol(data), iT))

	for(k in 1:kmax){

		effvol <- effvol + weights[k] * volk[[k]]

	}

	if(correlation){

		corr <- array(0L, dim = c(ncol(data),ncol(data), iT))
		#corr[,,1] <- cor(x[1:100, ])

		for(i in 1:iT){

		d <- sqrt(diag(effvol[,,i]))
		d <- diag(d, ncol(effvol), ncol(effvol))

		corr[,,i] <- inv(d) %*% effvol[,,i] %*% inv(d) 
		}

		return(corr)
	}

return(effvol)

}


#data get and preparation: 
library(alphavantager)

av_api_key('0WXHEIY0O87LX4A3')

TLT <- as.data.frame(av_get(symbol = "TLT", av_fun = "TIME_SERIES_DAILY", outputsize = "full"))
rownames(TLT) <- TLT$timestamp
SPY <- as.data.frame(av_get(symbol = "SPY", av_fun = "TIME_SERIES_DAILY", outputsize = "full"))
rownames(SPY) <- SPY$timestamp
IAU <- as.data.frame(av_get(symbol = "IAU", av_fun = "TIME_SERIES_DAILY", outputsize = "full"))
rownames(IAU) <- IAU$timestamp
LQD <- as.data.frame(av_get(symbol = "LQD", av_fun = "TIME_SERIES_DAILY", outputsize = "full"))
rownames(LQD) <- LQD$timestamp
#Real estate
IYR <- as.data.frame(av_get(symbol = "IYR", av_fun = "TIME_SERIES_DAILY", outputsize = "full"))
rownames(LQD) <- IYR$timestamp


#returns

library(xts)
returns_TLT <- as.xts(diff(log(TLT[,4])), order.by = as.Date(TLT[,1], format='%d/%m/%Y')[-1])
returns_SPY <- as.xts(diff(log(SPY[,4])), order.by = as.Date(SPY[,1])[-1])
returns_IAU <- as.xts(diff(log(IAU[,4])), order.by = as.Date(IAU[,1])[-1])
returns_LQD <- as.xts(diff(log(LQD[,4])), order.by = as.Date(LQD[,1])[-1])
returns_IYR <- as.xts(diff(log(IYR[,4])), order.by = as.Date(IYR[,1])[-1])

#Sequencing them into proper time series. 
returns_TLT <- returns_TLT[seq(from= as.Date('2010-01-02'), to = as.Date('2020-01-02'), by=1), ]
returns_SPY <- returns_SPY[seq(from= as.Date('2010-01-02'), to = as.Date('2020-01-02'), by=1), ]
returns_LQD <- returns_LQD[seq(from= as.Date('2010-01-02'), to = as.Date('2020-01-02'), by=1), ]
returns_IAU <- returns_IAU[seq(from= as.Date('2010-01-02'), to = as.Date('2020-01-02'), by=1), ]
returns_IYR <- returns_IYR[seq(from= as.Date('2010-01-02'), to = as.Date('2020-01-02'), by=1), ]

returns_collected <- cbind(returns_TLT, returns_SPY, returns_IAU, returns_LQD, returns_IYR)


#Analysis
asset_cov <- array(0L, dim =  c(5,5,length(returns_collected[,1])))

for(i in 1:length(returns_collected[,1])){


	asset_cov[,,i] <- realCov(returns_collected[i,], F)

}


realcor_quarterly <- na.omit(rollapply(returns_collected, 60, function(x) cor(x), by.column = F, align = 'left'))

#they will always be different, since 2006 uses the average of 14 different ewmas. 
ewma <- ewma.filter2006(returns_collected, T, 120, 1) 
ewma1993 <- ewma.filter(returns_collected, 60, T) #60

#plot(ewma1993[2,1,], type ='l', col = 'red') + lines(ewma[2,1,])

arraycorrelation <- array(0L, dim = c(5,5,length(realcor_quarterly[,1])))

for(i in 1:length(realcor_quarterly[,1])){

	arraycorrelation[,,i] <- realcor_quarterly[i,]

}



colnames(arraycorrelation) <- c('TLT', 'SPY', 'IAU', 'LQD', 'IYR')
rownames(arraycorrelation) <- c('TLT', 'SPY', 'IAU', 'LQD', 'IYR')

#check. It works!
cor(returns_collected[1:60,])

arraycorrelation[,,1]


library(ggplot2)


SMA <- data.frame(arraycorrelation[2,1,], arraycorrelation[3,1,], arraycorrelation[4,1,], arraycorrelation[5,1,],
	arraycorrelation[3,2,], arraycorrelation[4,2,], arraycorrelation[5,2,], arraycorrelation[4,3,], arraycorrelation[5,3, ], arraycorrelation[5,4, ])

EWMA <- data.frame(ewma[2,1,], ewma[3,1,], ewma[4,1,], ewma[5,1,], ewma[3,2,], ewma[4,2,], ewma[5,2,], ewma[4,3,], ewma[5,3, ], ewma[5,4, ])

colnames(SMA) <- c('SPY/TLT', 'IAU/TLT', 'LQD/TLT', 'IAU/SPY', 'LQD/SPY', 'LQD/IAU')


p1 <- ggplot() + geom_line(aes(index(returns_collected[-c(1:59)]), SMA[,1], col='SMA: SPY/TLT')) + 
	ylab('Correlation') + xlab('Year') + geom_hline(aes(yintercept=0)) + geom_line(aes(index(returns_collected[-c(1:59)]), EWMA[-c(1:59),1], col='EWMA: SPY/TLT'))
p2 <- ggplot() + geom_line(aes(index(returns_collected[-c(1:59)]), SMA[,2], col='SMA: IAU/TLT')) + 
	ylab('Correlation') + xlab('Year') + geom_hline(aes(yintercept=0)) + geom_line(aes(index(returns_collected[-c(1:59)]), EWMA[-c(1:59),2], col='EWMA: IAU/TLT'))
p3 <- ggplot() + geom_line(aes(index(returns_collected[-c(1:59)]), SMA[,3], col='SMA: LQD/TLT')) + 
	ylab('Correlation') + xlab('Year') + geom_hline(aes(yintercept=0)) + geom_line(aes(index(returns_collected[-c(1:59)]), EWMA[-c(1:59),3], col='EWMA: LQD/TLT'))
p4 <- ggplot() + geom_line(aes(index(returns_collected[-c(1:59)]), SMA[,4], col='SMA: IYR/TLT')) + 
	ylab('Correlation') + xlab('Year') + geom_hline(aes(yintercept=0)) + geom_line(aes(index(returns_collected[-c(1:59)]), EWMA[-c(1:59),4], col='EWMA: IYR/TLT'))
p5 <- ggplot() + geom_line(aes(index(returns_collected[-c(1:59)]), SMA[,5], col='SMA: IAU/SPY')) + 
	ylab('Correlation') + xlab('Year') + geom_hline(aes(yintercept=0)) + geom_line(aes(index(returns_collected[-c(1:59)]), EWMA[-c(1:59),5], col='EWMA: IAU/SPY'))
p6 <- ggplot() + geom_line(aes(index(returns_collected[-c(1:59)]), SMA[,6], col='SMA: LQD/SPY')) + 
	ylab('Correlation') + xlab('Year') + geom_hline(aes(yintercept=0)) + geom_line(aes(index(returns_collected[-c(1:59)]), EWMA[-c(1:59),6], col='EWMA: LQD/SPY'))
p7 <- ggplot() + geom_line(aes(index(returns_collected[-c(1:59)]), SMA[,7], col='SMA: IYR/SPY')) + 
	ylab('Correlation') + xlab('Year') + geom_hline(aes(yintercept=0)) + geom_line(aes(index(returns_collected[-c(1:59)]), EWMA[-c(1:59),7], col='EWMA: IYR/SPY'))
p8 <- ggplot() + geom_line(aes(index(returns_collected[-c(1:59)]), SMA[,8], col='SMA: LQD/IAU')) + 
	ylab('Correlation') + xlab('Year') + geom_hline(aes(yintercept=0)) + geom_line(aes(index(returns_collected[-c(1:59)]), EWMA[-c(1:59),8], col='EWMA: LQD/IAU'))
p9 <- ggplot() + geom_line(aes(index(returns_collected[-c(1:59)]), SMA[,9], col='SMA: IYR/IAU')) + 
	ylab('Correlation') + xlab('Year') + geom_hline(aes(yintercept=0)) + geom_line(aes(index(returns_collected[-c(1:59)]), EWMA[-c(1:59),9], col='EWMA: IYR/IAU'))
p10 <- ggplot() + geom_line(aes(index(returns_collected[-c(1:59)]), SMA[,10], col='SMA: IYR/LQD')) + 
	ylab('Correlation') + xlab('Year') + geom_hline(aes(yintercept=0)) + geom_line(aes(index(returns_collected[-c(1:59)]), EWMA[-c(1:59),10], col='EWMA: IYR/LQD'))
library(gridExtra)

#grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10, ncol=1)

grid.arrange(p1,p2)

grid.arrange(p3,p4)

grid.arrange(p5,p6)

grid.arrange(p7,p8)

grid.arrange(p9,p10)

#storng correlation
grid.arrange(p3, p7)

#Recode the EWMA filter using RM-2005 with reference values. Then put it into the correlation graphs together with
#the simple moving average. Code IYR into the correlation parameter also. 

#Code the eigenvalue decomposition method for verifying the fact that your covariance matrices satisfies regularity conditions.
#Remember that covariance matrix needs to satisfy positive semi-definiteness.  
#for large covariance matrices, the matrix might become singular. If we have an N dim vector with M samples of the N dim vector, then
#when M<N the covariance matrix of this vector will become singular with M-N eigenvalues equal to zero. 



#------------------------------------------------------------------------------------------------------------
library(ggplot2)
library(RColorBrewer)

ewma.SPY <- ewma.filter2006(returns_SPY, F, 22) 
ewma.SPY <- xts(as.vector(ewma.SPY), order.by = as.Date(index(returns_SPY)))

retplot <- ggplot() +  geom_line(aes(index(returns_SPY),returns_SPY,col='Daily log-returns')) +
xlab('Year') + ylab('Returns') + theme(plot.title = element_text(hjust = 0.5)) +   
geom_line(aes(index(returns_SPY), ewma.SPY, col ="EWMA(22)"), size = 1) +
scale_y_continuous(sec.axis = sec_axis(~.*1, name = "Variance")) +
scale_color_manual(values=c("#56B4E9","darkslategrey")) + 
theme(legend.title=element_blank(), legend.justification=c(0,1), legend.position=c(0.7,0.97),
	legend.background = element_rect(fill="lightblue",
    size=0.5, linetype="solid", 
    colour ="darkblue"), 
    legend.text = element_text(colour="black", size=10, face="bold"))


retplot




#-------------ACF plot---------------------------------------
#logreturns:

#fiveminreturns:
mergedfrequencies <- readRDS("mergedfrequencies.rds")
fiveminreturns <- do.call.rbind(mergedfrequencies[[7]])[,2]



bacf <- acf(returns_SPY, plot = FALSE)
bacfdf <- with(bacf, data.frame(lag, acf))
ciline <- qnorm((1 - 0.95)/2)/sqrt(length(returns_SPY))

#squared logreturns:
bacf2 <- acf(returns_SPY^2, plot = FALSE)
bacfdf2 <- with(bacf2, data.frame(lag, acf))

#squared logreturns intraday

calccov <- readRDS("calculatedcovariances.rds")


#fiveminrealizedvariances:
fiveminrv <- calccov[[1]][[7]][2,2,]
bacf2 <- acf(fiveminrv, plot = FALSE)
bacfdf2 <- with(bacf2, data.frame(lag, acf))



data2 <- data.frame(bacfdf, bacfdf2)

q <- ggplot(data = data2, mapping = aes(x = lag, y = acf)) + geom_hline(aes(yintercept = 0)) + 
geom_hline(aes(yintercept = ciline),col = 'black', linetype = "dashed", size = 1) + 
geom_hline(aes(yintercept = -ciline), col = 'black', linetype = "dashed", size = 1) + 
geom_segment(aes(xend =lag, yend = 0, linetype = 'Log-returns'), col = 'red', size=1) + xlab('Lag') + ylab('ACF') + 
geom_segment(aes(xend=lag, yend=0, y = bacfdf2[,2], linetype = 'Squared log-returns'), col = 'red', size = 1, alpha = 0.6) +
theme(legend.title = element_blank(), legend.position = c(0.8, 0.8), legend.background = element_rect(fill="lightblue", 
size=0.5, linetype="solid", colour = 'darkblue'), 
legend.text = element_text(colour="black", size=10, face="bold")) + 
scale_fill_manual(name ='', values = c("blue", "red")) 


acf <- q + geom_point(aes(bacfdf[,1], bacfdf[,2]), col='black') + geom_point(aes(bacfdf2[,1], bacfdf2[,2]), col='red')
acf


library(gridExtra)
tt <- grid.arrange(retplot, acf, ncol=2)

ggsave(tt, file="introductiongraph.pdf", device = "pdf")


#-----------------------previous tick vs interpolation ----------------

tickmethod <- aggregatets(log(timeseriesSPY[[1]]), FUN = 'previoustick', on = 'seconds', k = 5)

interpolationmethod <- aggregatets(log(timeseriesSPY[[1]]), FUN = 'mean', on = 'seconds', k = 5)

ggplot() + geom_step(aes((1:50)*5, tickmethod[1:50], col = 'Tick method'), size = 1) + 
geom_line(aes((1:50)*5, interpolationmethod[1:50], col = 'Interpolation'), size = 1) + 
ylab('Log-price') + xlab('Time in seconds') +
theme(legend.title = element_blank(), legend.position = c(0.8, 0.8), legend.background = element_rect(fill="lightblue", 
size=0.5, linetype="solid", colour = 'darkblue'), 
legend.text = element_text(colour="black", size=10, face="bold"))


#Using the rckernel: Only works for intraday returns with a MxN list for n being asset dim and number of list elements
# and M intraday returns for list i.
#therefore it needs to get looped for all days in 10 years. 

data(lltc.xts)
data(sbux.xts)

#assets should be in a list per day. So this list should get looped for every day, thus creating 252 checks, in a year. 

check <- list(lltc.xts, sbux.xts)


#original parameters, if you don't specify yourself. 


#rKernelCov <- function(rData, cor = FALSE,  alignBy = "seconds", alignPeriod = 1,
                      # makeReturns = FALSE, kernelType = "rectangular", kernelParam = 1,
                       #kernelDOFadj = TRUE)

rcKernel <- rKernelCov(rdata = check, makeReturns = FALSE, kernel.type = "Parzen")




#----------------------------------------realized semicov graph----------------------------------------------

realsemicov <- function(matrix,type, correlation = FALSE){
  # Realized SemiCovariance estimator for 1 day. 
  #
  # Args: 
  #
  #   matrix: A matrix of returns. Nxk 
  #
  #	  N: number of intraday returns, k: number of asset dimension, T: number of days. 
  realcovariance <- matrix()
  
  	positive <- matrix * (matrix >= 0)
  	negative <- matrix * (matrix <= 0)

  if(type =='P'){
  	realcovariance <- t(positive) %*% (positive)
  }
  if(type == 'N'){
  	realcovariance <- t(negative) %*% (negative)
  }
  if(type == 'M'){
  	realcovariance <- t(positive) %*% negative + t(negative) %*% positive
  }


  if(correlation){

		d <- sqrt(diag(realcovariance))
		d <- diag(d, ncol(realcovariance), ncol(realcovariance))

		corr <- inv(d) %*% realcovariance %*% inv(d) 

		return(corr)
	}
  
 return( realcovariance ) 
}


realsemicov(returns_collected,'N',F)


 rollingP <- na.omit(rollapply(returns_collected[,1:2], 90, function(x) realsemicov(x,'P',F), by.column = F, align = 'left'))
 rollingN <- na.omit(rollapply(returns_collected[,1:2], 90, function(x) realsemicov(x,'N',F), by.column = F, align = 'left'))
 rollingM <- na.omit(rollapply(returns_collected[,1:2], 90, function(x) realsemicov(x,'M',F), by.column = F, align = 'left'))


arraysrollingP <- array(0L, dim = c(2,2,length(rollingP[,1])))
arraysrollingN <- array(0L, dim = c(2,2,length(rollingN[,1])))
arraysrollingM <- array(0L, dim = c(2,2,length(rollingM[,1])))



for(i in 1:length(rollingP[,1])){


	arraysrollingP[,,i] <- rollingP[i,]
	arraysrollingN[,,i] <- rollingN[i,]
	arraysrollingM[,,i] <- rollingM[i,]


}


#----trying to compute leverage effect by cor(v,X).

covariancenumb <- integer(length(returns_collected[,2]))

for(i in 1:length(returns_collected[,2])){

	covariancenumb[i] <- realCov(returns_collected[i,1:2],T)[1,2]

}


rollingleverage <- na.omit(rollapply(returns_collected[,1:2], 90, function(x) realCov(x,T), by.column = F, align = 'left'))


r1 <- ggplot() + geom_line(aes(index(returns_collected)[-c(1:89)], (arraysrollingP[1,2,]+arraysrollingN[1,2,]), col ='Concordant'), lwd=1) + 
		geom_line(aes(index(returns_collected)[-c(1:89)], arraysrollingM[1,2,], col ='Mixed'), lwd = 1) + 
			geom_line(aes(index(returns_collected)[-c(1:89)], rollingleverage[,2]/80, col="Correlation"), alpha = 0.5, lwd = 1)+
			scale_y_continuous(sec.axis = sec_axis(~.*80, name = "Correlation")) + ylab('Semicovariance') + xlab('Dates') +
				theme(legend.title=element_blank(), legend.justification=c(0,1), legend.position=c(0.7,0.47),
	legend.background = element_rect(fill="lightblue",
    size=0.5, linetype="solid", 
    colour ="darkblue"), 
    legend.text = element_text(colour="black", size=10, face="bold"))





