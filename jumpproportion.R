#jump proportion table. 



source("Previous functions/projectfunctions.R")
source("functions.R")

library(xts)
library(highfrequency)
library(matlib)

dataTLT <- readRDS("dataTLT.rds")
dataSPY <- readRDS("dataSPY.rds")



#you only need open-to-close

secs <-c(1,5,15,20,30,60,300,900,1800)


frequenciesTLT <- list()
frequenciesSPY <- list()

tempTLT <- list()
tempSPY <- list()



for(j in 1:length(secs)){
	for(i in 1:length(dataTLT)){

		tempTLT[[i]] <- aggregatets(dataTLT[[i]], on="seconds", k=secs[j])
		tempSPY[[i]] <- aggregatets(dataSPY[[i]], on="seconds", k=secs[j])
	}

	frequenciesTLT[[j]] <- tempTLT
	frequenciesSPY[[j]] <- tempSPY
	print(sprintf("%s",j))
}

#log-returns

for(j in 1:length(secs)){
	for(i in 1:length(dataTLT)){

		frequenciesTLT[[j]][[i]] <- diff(log(frequenciesTLT[[j]][[i]]))[-1]
		frequenciesSPY[[j]][[i]] <- diff(log(frequenciesSPY[[j]][[i]]))[-1]

	}
}


#-----------------preparing open-to-close returns. Everything is calculated separately. 

getDates <- unlist(lapply(dataTLT, function(x) as.character(index(x[1]))))

for(i in 1:length(getDates)){
	getDates[i] <- strsplit(getDates, " ")[[i]][1]
}

opentocloseTLT <- list()
opentocloseSPY <- list()

for(i in 1:length(dataTLT)){

	opentocloseTLT[[i]] <- dataTLT[[i]][c(1,length(dataTLT[[i]]))]
	opentocloseSPY[[i]] <- dataSPY[[i]][c(1,length(dataSPY[[i]]))]
}

for(i in 1:length(dataTLT)){

	opentocloseTLT[[i]] <- diff(log(opentocloseTLT[[i]]))[-1]
	opentocloseSPY[[i]] <- diff(log(opentocloseSPY[[i]]))[-1]

}

mergedopentoclose <- list()

for(i in 1:length(dataTLT)){

	mergedopentoclose[[i]] <- cbind(opentocloseTLT[[i]], opentocloseSPY[[i]])

}



#mergedopentoclose <- lapply(mergedopentoclose, function(x) x[c(1,nrow(x))])
mergedopentoclose <- lapply(mergedopentoclose, function(x) colSums(x, na.rm = T))

mergedopentoclose <- lapply(mergedopentoclose, function(x) matrix(x, nrow=1, ncol=2, byrow = T))


for(i in 1:length(dataTLT)){


	mergedopentoclose[[i]] <- xts(mergedopentoclose[[i]], order.by = as.Date(getDates[i]))

}


dailyacross <- array(unlist(mergedopentoclose), c(2,2,2516))


#constructing an array instead of list!
for(i in 1:2516){

	dailyacross[,,i] <- t(dailyacross[,,i])

}




#opentocloseTLTnonneg <- lapply(opentocloseTLT, function(x) x[c(which(x>0)[1], max(which(x>0)))])
#opentocloseSPYnonneg <- lapply(opentocloseSPY, function(x) x[c(which(x>0)[1], max(which(x>0)))])



#opentocloseTLT <- lapply(opentocloseTLT, function(x) x[c(1,length(x))])
#opentocloseSPY <- lapply(opentocloseSPY, function(x) x[c(1,length(x))])




#--------------------jump proportion calculation-------------------------------


#ad-hoc method to computing the daily jump proportion
JPdaily <- function(data){
	Rcov <- realCov(data) 
	TLTRV <- Rcov[1,1]
	SPYRV <- Rcov[2,2]

	BPcov <- preavBPCOV(data,F,F,F)
	TLTBV <- BPcov[1,1]
	SPYBV <- BPcov[2,2]

	JPTLT <- (TLTRV - TLTBV)/TLTRV
	JPSPY <- (SPYRV - SPYBV)/SPYRV

	lout <- list(JPTLT,TLTRV, TLTBV, JPSPY, SPYRV, SPYBV)

	names(lout) <- c("JPTLT", "TLTRV", "TLTBV", "JPSPY", "SPYRV", "SPYBV")

	return(lout)
}

#minute, no pre-averaging
JPminute <-function(data){

	RV <- realCov(data)
	BV <- preavBPCOV(data,F,F,F)

	JP <- (RV - BV)/RV

	lout <- list(RV, BV, JP)

	names(lout) <- c("RV", "BV", "JP")

	return(lout)
}

#pre-average version. 

JPsecond <-function(data, theta=1){

	RV <- preavCov(data,T,T,F,theta)
	BV <- preavBPCOV(data,T,F,T,theta)

	JP <- (RV - BV)/RV

	lout <- list(RV, BV, JP)

	names(lout) <- c("RV", "BV", "JP")

	return(lout)
}


#first analysis for open-to-close:

jpdaily <- list()

for(i in 1:length(dataTLT)){

	jpdaily[[i]] <- JPdaily(dailyacross[,,i])

}


#forcing jp interval between [0,1] due to extremely high measurement error
jptlt_daily <- sapply(jpdaily, function(x) x$JPTLT)
jptlt_daily[jptlt_daily < 0] <- 0 

rvtlt_daily <- sapply(jpdaily, function(x) x$TLTRV)
bvtlt_daily <- sapply(jpdaily, function(x) x$TLTBV)

sqrt(mean(bvtlt_daily)*252)*100
sqrt(mean(rvtlt_daily)*252)*100
mean(jptlt_daily)*100


 #(mean(rvtlt_daily+0.0000118) 

jpspy_daily <- sapply(jpdaily, function(x) x$JPSPY)
jpspy_daily[jpspy_daily < 0] <- 0
rvspy_daily <- sapply(jpdaily, function(x) x$SPYRV)
bvspy_daily <- sapply(jpdaily, function(x) x$SPYBV)

sqrt(mean(bvspy_daily)*252)*100
sqrt(mean(rvspy_daily)*252)*100
mean(jpspy_daily)*100

#mean(rvspy_daily+0.00002344)
#Do the analysis. First five needs jpsecond, last five needs jpminute. 

JPTLT <- list()
JPSPY <- list()

RVTLT <- list()
RVSPY <- list()

BVTLT <- list()
BVSPY <- list()

tempTLT <- list()
tempSPY <- list()

tempTLT2 <- list()
tempSPY2 <- list()

#Splitting up second and minute calculations due to second calc being incredibly slow. 

for(j in 1:4){
	for(i in 1:length(dataTLT)){
		tempTLT2[[i]] <- JPminute(frequenciesTLT[[j+5]][[i]])
		tempSPY2[[i]] <- JPminute(frequenciesSPY[[j+5]][[i]])
	}
	print(sprintf("%s",j))
	JPTLT[[j+5]] <- unlist(lapply(tempTLT2, function(x) x$JP))
	RVTLT[[j+5]] <- unlist(lapply(tempTLT2, function(x) x$RV))
	BVTLT[[j+5]] <- unlist(lapply(tempTLT2, function(x) x$BV))

	JPSPY[[j+5]] <- unlist(lapply(tempSPY2, function(x) x$JP))
	RVSPY[[j+5]] <- unlist(lapply(tempSPY2, function(x) x$RV))
	BVSPY[[j+5]] <- unlist(lapply(tempSPY2, function(x) x$BV))


}


for(j in 1:5)){
	for(i in 1:length(dataTLT)){
		tempTLT[[i]] <- JPsecond(frequenciesTLT[[j]][[i]])

		tempSPY[[i]] <- JPsecond(frequenciesSPY[[j]][[i]])

		print(sprintf("Second frequency of %s over day %s",secs[j],i))

	}
	
	JPTLT[[j]] <- unlist(lapply(tempTLT, function(x) x$JP))

	RVTLT[[j]] <- unlist(lapply(tempTLT, function(x) x$RV))

	BVTLT[[j]] <- unlist(lapply(tempTLT, function(x) x$BV))

	JPSPY[[j]] <- unlist(lapply(tempSPY, function(x) x$JP))

	RVSPY[[j]] <- unlist(lapply(tempSPY, function(x) x$RV))

	BVSPY[[j]] <- unlist(lapply(tempSPY, function(x) x$BV))

}


#saveRDS(list(JPTLT, RVTLT, BVTLT, JPSPY, RVSPY, BVSPY), file = "jumpproportiondata.rds") 



#Estimated jp can go outside the interval [0,1]. One solution would be to set all jp<0 equal to 0. 


for(i in 1:length(JPTLT)){
	JPTLT[[i]][JPTLT[[i]]<0] <- 0
	RVTLT[[i]][RVTLT[[i]]<0] <- 0
	BVTLT[[i]][BVTLT[[i]]<0] <- 0

	JPSPY[[i]][JPSPY[[i]]<0] <- 0
	RVSPY[[i]][RVSPY[[i]]<0] <- 0
	BVSPY[[i]][BVSPY[[i]]<0] <- 0
}



meanJPTLT <- unlist(lapply(JPTLT, function(x) mean(x, na.rm = T)))
meanRVTLT <- unlist(lapply(RVTLT, function(x) mean(x, na.rm = T)))
meanBVTLT <- unlist(lapply(BVTLT, function(x) mean(x, na.rm = T)))

meanJPSPY <- unlist(lapply(JPSPY, function(x) mean(x, na.rm = T)))
meanRVSPY <- unlist(lapply(RVSPY, function(x) mean(x, na.rm = T)))
meanBVSPY <- unlist(lapply(BVSPY, function(x) mean(x, na.rm = T)))


meanJPTLT*100 #Percentage
sqrt(252*meanRVTLT)*100 #Percentage annualized standard deviation.
sqrt(252*meanBVTLT)*100


meanJPSPY*100
sqrt(252*meanRVSPY)*100
sqrt(252*meanBVSPY)*100


#avg
colMeans(rbind(sqrt(252*meanRVTLT)*100,sqrt(252*meanRVSPY)*100))

colMeans(rbind(sqrt(252*meanBVTLT)*100,sqrt(252*meanBVSPY)*100))

colMeans(rbind(meanJPTLT*100,meanJPSPY*100))




#--------------------------------visualising Jump variation----------------
jp <- readRDS("jumpproportiondata.rds")


dailyacross <- array(unlist(mergedopentoclose), c(2,2,2516))

for(i in 1:2516){

	dailyacross[,,i] <- t(dailyacross[,,i])

}


Jpdaily <- list()

for(i in 1:length(dataTLT)){

	Jpdaily[[i]] <- JPdaily(dailyacross[,,i])

}


rvtlt_daily <- sapply(Jpdaily, function(x) x$TLTRV)
bvtlt_daily <- sapply(Jpdaily, function(x) x$TLTBV)

rvspy_daily <- sapply(Jpdaily, function(x) x$SPYRV)
bvspy_daily <- sapply(Jpdaily, function(x) x$SPYBV)



#5min

jvtlt_5min <- numeric()

jvspy_5min <- numeric()


jvtlt_daily <- numeric()

jvspy_daily <- numeric()


for(i in 1:length(jp[[2]][[7]])){

	jvtlt_5min[i] <- max(sqrt(jp[[2]][[7]][i]*252) - sqrt(jp[[3]][[7]][i]*252),0)*100
	jvspy_5min[i] <- max(sqrt(jp[[5]][[7]][i]*252) - sqrt(jp[[6]][[7]][i]*252),0)*100

	jvtlt_daily[i] <- max(sqrt(rvtlt_daily[i]*252) - sqrt(bvtlt_daily[i]*252),0)*100
	jvspy_daily[i] <- max(sqrt(rvspy_daily[i]*252) - sqrt(bvspy_daily[i]*252),0)*100
}


library(ggplot2)
library(gridExtra)
library(ggpubr)



p1 <- ggplot() + geom_line(aes(as.Date(getDates), jvtlt_5min), col ="cyan4") + xlab("") + ylab("JV (%)") + 
ggtitle("TLT") + theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
p2 <- ggplot() + geom_line(aes(as.Date(getDates), jvtlt_daily),col ="slateblue4")+ xlab("") + ylab("JV (%)") +
ggtitle("TLT") + theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
p3 <- ggplot() + geom_line(aes(as.Date(getDates), jvspy_5min),col ="cyan4") + xlab("") + ylab("JV (%)") +
ggtitle("SPY") + theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
p4 <- ggplot() + geom_line(aes(as.Date(getDates), jvspy_daily),col ="slateblue4")+ xlab("") + ylab("JV (%)") +
ggtitle("SPY") + theme(plot.title = element_text(hjust = 0.5), legend.position = "none")



ggarrange(p1,p2,p3,p4,ncol=2, nrow=2)

#ggsave("JVdiff.eps", device = "eps")





#theta plot (done on in 2019):

#THIS IS NOT USEFUL. DISCARDED.

theta <- seq(0.1,2,0.1)

JPTLT_plot <- list()
JPSPY_plot <- list()

tempTLT <- matrix(0L, ncol = length(days2019), nrow=length(theta))
tempSPY <- matrix(0L, ncol = length(days2019), nrow=length(theta))

days2019 <- seq(2265, 2516,1)

for(j in 1:(length(frequenciesTLT)/2)){ #all second frequencies
	for(i in 1:length(theta)){ #amount of thetas to test
		for(k in 1:length(days2019)){ #first 10 days
			tempTLT[i,k] <- c(JPsecond(frequenciesTLT[[j]][[days2019[k]]],theta[i])$JP)
			tempSPY[i,k] <- c(JPsecond(frequenciesSPY[[j]][[days2019[k]]],theta[i])$JP)

			print(sprintf("Theta: %s, Frequency: %s, Day: %s",i,j,k))
		}
	}
	JPTLT_plot[[j]] <- tempTLT
	JPSPY_plot[[j]] <- tempSPY
}

saveRDS(list(JPTLT_plot, JPSPY_plot), file = "jumppropforplot.rds")

library(ggplot2)
library(zoo)

ggplot() + geom_line(aes(theta, rowMeans(JPSPY_plot[[3]])))

tester <- JPSPY_plot[[1]]

tester[tester<0] <- NaN

ts.plot(rowMeans(tester, na.rm = T))

tt <- rowMeans(rollapply(tester, 10, function(x) rowMeans(x), by.column = F, align = 'center'))

lines(theta, rowMeans(JPTLT_plot[[6]]))


ts.plot(tt)