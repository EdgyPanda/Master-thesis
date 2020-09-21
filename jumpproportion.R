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


#opentocloseTLTnonneg <- lapply(opentocloseTLT, function(x) x[c(which(x>0)[1], max(which(x>0)))])
#opentocloseSPYnonneg <- lapply(opentocloseSPY, function(x) x[c(which(x>0)[1], max(which(x>0)))])



#opentocloseTLT <- lapply(opentocloseTLT, function(x) x[c(1,length(x))])
#opentocloseSPY <- lapply(opentocloseSPY, function(x) x[c(1,length(x))])


frequenciesSPY[[10]] <- opentocloseSPY
frequenciesTLT[[10]] <- opentocloseTLT


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

	return(c(JPTLT, JPSPY))
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