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

	opentocloseTLT[[i]] <- diff(log(opentocloseTLT[[i]]))[-1] * 100
	opentocloseSPY[[i]] <- diff(log(opentocloseSPY[[i]]))[-1] * 100

}


mergedopentoclose <- list()

for(i in 1:length(dataTLT)){

	mergedopentoclose[[i]] <- cbind(opentocloseTLT[[i]], opentocloseSPY[[i]])

}


#log-returns and not percentage log-returns.
mergedfrequencies <- readRDS("mergedfrequencies.rds")



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

