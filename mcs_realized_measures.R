#MCS REALIZED MEASURES & summary statistics for both assets.  


source("Previous functions/projectfunctions.R")
source("functions.R")

library(xts)
library(highfrequency)
library(matlib)
library(MCS)

dataTLT <- readRDS("dataTLT.rds")
dataSPY <- readRDS("dataSPY.rds")







bandwidthH <- function(list, sparsedata){

	n <- as.vector(sapply(list,length))

	#n <- as.vector(sapply(list, length))

	c <- ((12)^2/0.269)^(1/5)

	#w <- -1/(n-1) * t(data[1:(n-1)]) %*% data[2:n]
	w <- numeric()
	IV <- numeric()
	for(i in 1:length(list)){

		w[i] <- 1/(2*n[i]) * realCov(list[[i]])
		IV[i] <- realCov(sparsedata[[i]])

	}

	e <- w/(IV/n)

	#e <- noisetosignal(list)

	H <- c * e^(4/5) * n^(3/5)

	return(H)

}





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

mergedopentoclose <- list()

for(i in 1:length(dataTLT)){

	mergedopentoclose[[i]] <- na.omit(cbind(opentocloseTLT[[i]], opentocloseSPY[[i]]))

}

opentocloseTLT <- lapply(opentocloseTLT, function(x) x[c(1,length(x))])
opentocloseSPY <- lapply(opentocloseSPY, function(x) x[c(1,length(x))])
mergedopentoclose <- lapply(mergedopentoclose, function(x) x[c(1,nrow(x))])


#merging list elements. 

mergedfrequencies <- list()
tempmerged <- list()

for(j in 1:length(secs)){
	for(i in 1:length(dataTLT)){

		tempmerged[[i]] <- na.omit(cbind(frequenciesTLT[[j]][[i]], frequenciesSPY[[j]][[i]]))

	}
	mergedfrequencies[[j]] <- tempmerged
}



#--------------------------------------------------SUMMARY STATISTICS------------------------------------------

#We do it on preferred sampling schemes. 


RCov_5min <- lapply(mergedfrequencies[[7]], function(x) realCov(x))
RCov_daily <- lapply(mergedopentoclose, function(x) realCov(x))

MRC_30sec <- lapply(mergedfrequencies[[4]], function(x) preavCov(x, T, T, F, 1))
MRK_30sec