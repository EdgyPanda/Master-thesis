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







tester <- list() #opentoclose
tester2 <- list() #opentoclose first non negative return. 
for(i in 1:length(frequenciesTLT[[10]])){

	tester[[i]] <- JPminute(frequenciesTLT[[10]][[i]])
	tester2[[i]] <- JPminute(opentocloseTLTnonneg[[i]])
}

#In essence, just use normal open-to-close log-returns. 




test <- lapply(tester2, function(x) x$JP)
test <- unlist(test)

JPminute(opentocloseTLTnonneg[[8]])