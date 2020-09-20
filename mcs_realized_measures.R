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


sparseTLT20min <- list()
sparseSPY20min <- list()


for(i in 1:length(dataTLT)){

	sparseTLT20min[[i]] <- aggregatets(dataTLT[[i]], on = "minutes", k = 20)
	sparseSPY20min[[i]] <- aggregatets(dataSPY[[i]], on = "minutes", k = 20)

}

for(i in 1:length(dataTLT)){

	sparseTLT20min[[i]] <- diff(log(sparseTLT20min[[i]]))[-1]
	sparseSPY20min[[i]] <- diff(log(sparseSPY20min[[i]]))[-1]

}
#--------------------------------------------------

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

#finding optimal bandwidth for MRK:

H <- cbind(bandwidthH(frequenciesTLT[[3]],sparseTLT20min), bandwidthH(frequenciesSPY[[3]],sparseSPY20min))

H <- rowMeans(H)

getDates <- unlist(lapply(dataTLT, function(x) as.character(index(x[1]))))

for(i in 1:length(getDates)){
	getDates[i] <- strsplit(getDates, " ")[[i]][1]
}

library(ggplot2)

ggplot() + geom_line(aes(as.Date(getDates), H, group = 1)) + geom_hline(yintercept = mean(H), col="steelblue", lwd=1)

#We do it on preferred sampling schemes. 


RCov_5min <- lapply(mergedfrequencies[[7]], function(x) realCov(x))
RCov_daily <- lapply(mergedopentoclose, function(x) realCov(x))

MRC_30sec <- lapply(mergedfrequencies[[4]], function(x) preavCov(x, T, T, F, 1))
MRK_30sec <- 