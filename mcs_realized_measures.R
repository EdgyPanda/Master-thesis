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


library(ggplot2)

ggplot() + geom_line(aes(as.Date(getDates), H, group = 1)) + geom_hline(yintercept = mean(H), col="steelblue", lwd=1)

#We do it on preferred sampling schemes. 


tt <- lapply(mergedfrequencies[[7]], function(x) x*100)

RCov_5min <- lapply(mergedfrequencies[[7]], function(x) realCov(x))
BPCov_5min <- lapply(mergedfrequencies[[7]], function(x) preavBPCOV(x, F, F, F))
TCov_5min <- lapply(mergedfrequencies[[7]], function(x) preavthrCOV(x,F))
RSCovpos_5min <- lapply(mergedfrequencies[[7]], function(x) realsemicov(x,"P"))
RSCovneg_5min <- lapply(mergedfrequencies[[7]], function(x) realsemicov(x,"N"))
RCov_daily <- lapply(mergedopentoclose, function(x) realCov(x))

MRC_30sec <- lapply(mergedfrequencies[[5]], function(x) preavCov(x, T, T, F, 1))
PBPCov_30sec <- lapply(mergedfrequencies[[7]], function(x) preavBPCOV(x, T, F, T, 1))

MRK_15sec <- list()

for(i in 1:length(dataTLT)){

	MRK_15sec[[i]] <-  rKernelCov(list(mergedfrequencies[[3]][[i]][,1], mergedfrequencies[[3]][[i]][,2]), makeReturns = FALSE, kernel.type = "Parzen", kernel.param  = H[i])


}


RCov_5min <- matrix(unlist(lapply(RCov_5min, function(x) cbind(x[1,1], x[2,2], x[2,1]/(sqrt(x[2,2])*sqrt(x[1,1]))))), 
	nrow = 2516, ncol=3, byrow=T)

colnames(RCov_5min) <- c("TLT", "SPY", "Correlation")

BPCov_5min <- matrix(unlist(lapply(BPCov_5min, function(x) cbind(x[1,1], x[2,2], x[2,1]/(sqrt(x[2,2])*sqrt(x[1,1]))))), 
	nrow = 2516, ncol=3, byrow=T)

colnames(BPCov_5min) <- c("TLT", "SPY", "Correlation")

TCov_5min <- matrix(unlist(lapply(TCov_5min, function(x) cbind(x[1,1], x[2,2], x[2,1]/(sqrt(x[2,2])*sqrt(x[1,1]))))), 
	nrow = 2516, ncol=3, byrow=T)

colnames(TCov_5min) <- c("TLT", "SPY", "Correlation")

RSCovpos_5min <- matrix(unlist(lapply(RSCovpos_5min, function(x) cbind(x[1,1], x[2,2], x[2,1]/(sqrt(x[2,2])*sqrt(x[1,1]))))), 
	nrow = 2516, ncol=3, byrow=T)

colnames(RSCovpos_5min) <- c("TLT", "SPY", "Correlation")

RSCovneg_5min <- matrix(unlist(lapply(RSCovneg_5min, function(x) cbind(x[1,1], x[2,2], x[2,1]/(sqrt(x[2,2])*sqrt(x[1,1]))))), 
	nrow = 2516, ncol=3, byrow=T)

colnames(RSCovneg_5min) <- c("TLT", "SPY", "Correlation")

RCov_daily <- matrix(unlist(lapply(RCov_daily, function(x) cbind(x[1,1], x[2,2], x[2,1]/(sqrt(x[2,2])*sqrt(x[1,1]))))), 
	nrow = 2516, ncol=3, byrow=T)

colnames(RCov_daily) <- c("TLT", "SPY", "Correlation")

MRC_30sec <- matrix(unlist(lapply(MRC_30sec, function(x) cbind(x[1,1], x[2,2], x[2,1]/(sqrt(x[2,2])*sqrt(x[1,1]))))), 
	nrow = 2516, ncol=3, byrow=T)

colnames(MRC_30sec) <- c("TLT", "SPY", "Correlation")

PBPCov_30sec <- matrix(unlist(lapply(PBPCov_30sec, function(x) cbind(x[1,1], x[2,2], x[2,1]/(sqrt(x[2,2])*sqrt(x[1,1]))))), 
	nrow = 2516, ncol=3, byrow=T)

colnames(PBPCov_30sec) <- c("TLT", "SPY", "Correlation")

MRK_15sec <- matrix(unlist(lapply(MRK_15sec, function(x) cbind(x[1,1], x[2,2], x[2,1]/(sqrt(x[2,2])*sqrt(x[1,1]))))), 
	nrow = 2516, ncol=3, byrow=T)

colnames(MRK_15sec) <- c("TLT", "SPY", "Correlation")


daily_SPY <- read.csv("daily_SPY.csv", header = T)
daily_SPY <- diff(log(daily_SPY$close))[-1]

daily_TLT <- read.csv("daily_TLT.csv", header = T)
daily_TLT <- diff(log(daily_TLT$close))[-1]

merged <- cbind(daily_TLT, daily_SPY)

#HOW SHOULD YOU TRANSFORM IT?!?
RCov_5min[,1:2] <- RCov_5min[,1:2]*(24/6.5) * 100

BPCov_5min[,1:2] <- 252*BPCov_5min[,1:2]*100*(24/6.5)

TCov_5min[,1:2] <- 252*TCov_5min[,1:2]*100*(24/6.5)

RSCovpos_5min[,1:2] <- 252*RSCovpos_5min[,1:2]*100*(24/6.5)

RSCovneg_5min[,1:2] <- 252*RSCovneg_5min[,1:2]*100*(24/6.5)

RCov_daily[,1:2] <- 252*RCov_daily[,1:2]*(24/6.5)

MRC_30sec[,1:2] <- 252*MRC_30sec[,1:2]*100*(24/6.5)

PBPCov_30sec[,1:2] <- 252*PBPCov_30sec[,1:2]*100*(24/6.5)

MRK_15sec[,1:2] <- 252*MRK_15sec[,1:2]*100*(24/6.5)

library(PerformanceAnalytics)

#using performanceAnalytics package you can easily get your summary statistics. 

#TLT
table.Stats(RCov_5min[,1]*1e5)
table.Autocorrelation(RCov_5min[,1])

table.Stats(BPCov_5min[,1])
table.Autocorrelation(BPCov_5min[,1])

table.Stats(TCov_5min[,1])
table.Autocorrelation(TCov_5min[,1])

table.Stats(RSCovpos_5min[,1])
table.Autocorrelation(RSCovpos_5min[,1])

table.Stats(RSCovneg_5min[,1])
table.Autocorrelation(RSCovneg_5min[,1])

table.Stats(RCov_daily[,1])
table.Autocorrelation(RCov_daily[,1])

table.Stats(MRC_30sec[,1])
table.Autocorrelation(MRC_30sec[,1])

table.Stats(PBPCov_30sec[,1])
table.Autocorrelation(PBPCov_30sec[,1])

table.Stats(MRK_15sec[,1])
table.Autocorrelation(MRK_15sec[,1])

#SPY
table.Stats(RCov_5min[,2])
table.Autocorrelation(RCov_5min[,2])

table.Stats(BPCov_5min[,2])
table.Autocorrelation(BPCov_5min[,2])

table.Stats(TCov_5min[,2])
table.Autocorrelation(TCov_5min[,2])

table.Stats(RSCovpos_5min[,2])
table.Autocorrelation(RSCovpos_5min[,2])

table.Stats(RSCovneg_5min[,2])
table.Autocorrelation(RSCovneg_5min[,2])

table.Stats(RCov_daily[,2])
table.Autocorrelation(RCov_daily[,2])

table.Stats(MRC_30sec[,2])
table.Autocorrelation(MRC_30sec[,2])

table.Stats(PBPCov_30sec[,2])
table.Autocorrelation(PBPCov_30sec[,2])

table.Stats(MRK_15sec[,2])
table.Autocorrelation(MRK_15sec[,2])

#Correlation
table.Stats(RCov_5min[,3])
table.Autocorrelation(RCov_5min[,3])

table.Stats(BPCov_5min[,3])
table.Autocorrelation(BPCov_5min[,3])

table.Stats(TCov_5min[,3])
table.Autocorrelation(TCov_5min[,3])

table.Stats(RSCovpos_5min[,3])
table.Autocorrelation(RSCovpos_5min[,3])

table.Stats(RSCovneg_5min[,3])
table.Autocorrelation(RSCovneg_5min[,3])

table.Stats(RCov_daily[,3])
table.Autocorrelation(RCov_daily[,3])

table.Stats(MRC_30sec[,3])
table.Autocorrelation(MRC_30sec[,3])

table.Stats(PBPCov_30sec[,3])
table.Autocorrelation(PBPCov_30sec[,3])

table.Stats(MRK_15sec[,3])
table.Autocorrelation(MRK_15sec[,3])

#introductory graph

meanpriceTLT <- unlist(lapply(dataTLT, function(x) mean(x)))

p1 <- ggplot() + geom_line(aes(as.Date(getDates), RCov_5min[,1]/(24/6.5)), col = "firebrick") + ylab("Ann. variance (%)") + 
xlab("Dates") + ggtitle("TLT") + theme(plot.title = element_text(hjust = 0.5)) 


p2 <- ggplot() + geom_line(aes(as.Date(getDates), RCov_5min[,2]/(24/6.5)), col = "deepskyblue1") + ylab("Ann. variance (%)") + 
xlab("Dates") + ggtitle("SPY") + theme(plot.title = element_text(hjust = 0.5))


p3 <- ggplot() + geom_line(aes(as.Date(getDates), RCov_5min[,3]), col = "maroon") + ylab("Correlation") + 
xlab("Dates") + theme(plot.title = element_text(hjust = 0.5)) + 
geom_hline(yintercept = 0, linetype = "dashed")


library(gridExtra)

pf <- grid.arrange(p1,p2,p3, ncol=1)


ggsave(pf, file="varianceassets.eps", device="eps")