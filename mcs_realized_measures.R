#MCS REALIZED MEASURES & summary statistics for both assets.  


source("Previous functions/projectfunctions.R")
source("functions.R")

library(xts)
library(highfrequency)
library(matlib)
library(MCS)

dataTLT <- readRDS("dataTLT.rds")
dataSPY <- readRDS("dataSPY.rds")






#need sparse sampled data for bandwidth selection. 

getDates <- unlist(lapply(dataTLT, function(x) as.character(index(x[1]))))

for(i in 1:length(getDates)){
	getDates[i] <- strsplit(getDates, " ")[[i]][1]
}


sparseTLT20min <- list()
sparseSPY20min <- list()
opentocloseTLT <- list()
opentocloseSPY <- list()

for(i in 1:length(dataTLT)){

	sparseTLT20min[[i]] <- aggregatets(dataTLT[[i]], on = "minutes", k = 20)
	sparseSPY20min[[i]] <- aggregatets(dataSPY[[i]], on = "minutes", k = 20)

}

for(i in 1:length(dataTLT)){

	sparseTLT20min[[i]] <- diff(log(sparseTLT20min[[i]]))[-1]
	sparseSPY20min[[i]] <- diff(log(sparseSPY20min[[i]]))[-1]

}

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

mergedopentoclose <- lapply(mergedopentoclose, function(x) colSums(x, na.rm = T))

mergedopentoclose <- lapply(mergedopentoclose, function(x) matrix(x, nrow=1, ncol=2, byrow = T))


for(i in 1:length(dataTLT)){


	mergedopentoclose[[i]] <- xts(mergedopentoclose[[i]], order.by = as.Date(getDates[i]))

}



#log-returns and not percentage log-returns.
mergedfrequencies <- readRDS("mergedfrequencies.rds")

mergedfrequencies[[10]] <- mergedopentoclose

#----------------------Finding optimal bandwidth for all frequencies: 

#This is loaded into the bandwidthH.rds for better access. Below takes 15 mins to run. 



#H <- list()

#frequenciesTLT <- lapply(mergedfrequencies, function(x) sapply(x, function(z) z[,1]))

#frequenciesSPY <- lapply(mergedfrequencies, function(x) sapply(x, function(z) z[,2]))



#for(i in 1:length(mergedfrequencies)){

#	temp <- cbind(bandwidthH(frequenciesTLT[[i]],sparseTLT20min), 
#		bandwidthH(frequenciesSPY[[i]],sparseSPY20min))
#
#	H[[i]] <- rowMeans(temp)
#	print(sprintf("%s", i))
#}

#saveRDS(H,"bandwidthH.rds")

H <- readRDS("bandwidthH.rds")


# -------------------------------------------Calculating realized measures  across frequencies -------------------
#
#
#


Rcov_frequencies <- list()

tempRcov <- array(0L, dim = c(2,2,length(dataTLT)))

Rcovpos_frequencies <- list()

tempRcovpos <- array(0L, dim = c(2,2,length(dataTLT)))

Rcovneg_frequencies <- list()

tempRcovneg <- array(0L, dim = c(2,2,length(dataTLT)))

Tcov_frequencies <- list()

tempTcov <- array(0L, dim = c(2,2,length(dataTLT)))

for(i in 1:length(mergedfrequencies)){
	for(j in 1:length(dataTLT)){

		tempRcov[,,j] <- realCov(mergedfrequencies[[i]][[j]])
		tempRcovpos[,,j] <- realsemicov(mergedfrequencies[[i]][[j]], "P")
		tempRcovneg[,,j] <- realsemicov(mergedfrequencies[[i]][[j]], "N")
		tempTcov[,,j] <- preavthrCOV(mergedfrequencies[[i]][[j]])

		print(sprintf("frequency: %s, day: %s", i,j))
	}
	Rcov_frequencies[[i]] <- tempRcov	
	Rcovpos_frequencies[[i]] <- tempRcovpos
	Rcovneg_frequencies[[i]] <- tempRcovneg
	Tcov_frequencies[[i]] <- tempTcov
}


BPcov_frequencies <- list()

tempBPcov <- array(0L, dim = c(2,2,length(dataTLT)))

MRC_frequencies <- list()

tempMRC <- array(0L, dim = c(2,2,length(dataTLT)))

PBPcov_frequencies <- list()

tempPBPcov <- array(0L, dim = c(2,2,length(dataTLT)))

MRK_frequencies <- list()

tempMRK <- array(0L, dim = c(2,2,length(dataTLT)))


for(i in 1:(length(mergedfrequencies)-1)){
	for(j in 1:length(dataTLT)){

		tempBPcov[,,j] <- preavBPCOV(mergedfrequencies[[i]][[j]],F,F,F)
		tempPBPcov[,,j] <- preavBPCOV(mergedfrequencies[[i]][[j]],T,F,T,1)
		tempMRC[,,j] <- preavCov(mergedfrequencies[[i]][[j]], T, T, F, 1)
		tempMRK[,,j] <- rKernelCov(list(mergedfrequencies[[i]][[j]][,1], mergedfrequencies[[i]][[j]][,2]), 
		makeReturns = FALSE, kernel.type = "Parzen", kernel.param  = H[[i]][j])
		print(sprintf("frequency: %s, day: %s", i,j))

	}

	BPcov_frequencies[[i]] <- tempBPcov
	PBPcov_frequencies[[i]] <- tempPBPcov
	MRC_frequencies[[i]] <- tempMRC
	MRK_frequencies[[i]] <- tempMRK
}


#estimators that doesn't work on daily data: BPCov (sampling across days works), PBPCov, MRC.

saveRDS(list(Rcov_frequencies, Rcovpos_frequencies, Rcovneg_frequencies, Tcov_frequencies, BPcov_frequencies, 
	PBPcov_frequencies, MRC_frequencies, MRK_frequencies), file = "calculatedcovariances.rds")