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


opentocloseTLT <- list()
opentocloseSPY <- list()

for(i in 1:length(dataTLT)){

	opentocloseTLT[[i]] <- diff(log(dataTLT[[i]]))[-1]
	opentocloseSPY[[i]] <- diff(log(dataSPY[[i]]))[-1]
}


#find first non-negative return and last non-negative return. 
#Reason: Each day we have two returns, one of them will be zero making BV zero and we will get negative jump proportion. 

opentocloseTLTnonneg <- lapply(opentocloseTLT, function(x) x[c(which(x>0)[1], max(which(x>0)))])
opentocloseSPYnonneg <- lapply(opentocloseSPY, function(x) x[c(which(x>0)[1], max(which(x>0)))])



opentocloseTLT <- lapply(opentocloseTLT, function(x) x[c(1,length(x))])
opentocloseSPY <- lapply(opentocloseSPY, function(x) x[c(1,length(x))])


frequenciesSPY[[10]] <- opentocloseSPY
frequenciesTLT[[10]] <- opentocloseTLT


#--------------------jump proportion calculation-------------------------------

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

JPsecond <-function(data){

	RV <- preavCov(data,T,T,F,1)
	BV <- preavBPCOV(data,T,F,T,1)

	JP <- (RV - BV)/RV

	lout <- list(RV, BV, JP)

	names(lout) <- c("RV", "BV", "JP")

	return(lout)
}


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

for(j in 1:(length(frequenciesTLT)/2)){
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


for(j in 1:(length(frequenciesTLT)/2)){
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


saveRDS(list(JPTLT, RVTLT, BVTLT, JPSPY, RVSPY, BVSPY), file = "jumpproportiondata.rds") 



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
