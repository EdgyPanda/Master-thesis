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


mergedfrequencies <- readRDS("mergedfrequencies.rds")

#Finding optimal bandwidth for all frequencies: 

