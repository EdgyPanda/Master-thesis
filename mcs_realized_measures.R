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


#skips code when executing the entire script.
if(FALSE){
H <- list()

frequenciesTLT <- lapply(mergedfrequencies, function(x) sapply(x, function(z) z[,1]))

frequenciesSPY <- lapply(mergedfrequencies, function(x) sapply(x, function(z) z[,2]))



for(i in 1:length(mergedfrequencies)){

	temp <- cbind(bandwidthH(frequenciesTLT[[i]],sparseTLT20min), 
		bandwidthH(frequenciesSPY[[i]],sparseSPY20min))

	H[[i]] <- rowMeans(temp)
	print(sprintf("%s", i))
}

#saveRDS(H,"bandwidthH.rds")
}


H <- readRDS("bandwidthH.rds")


# -------------------------------------------Calculating realized measures  across frequencies -------------------
#
#
#

#skips code when executing the entire script.
if(FALSE){

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
}

#estimators that doesn't work on daily data: BPCov (sampling across days works), PBPCov, MRC.

#saveRDS(list(Rcov_frequencies, Rcovpos_frequencies, Rcovneg_frequencies, Tcov_frequencies, BPcov_frequencies, 
#	PBPcov_frequencies, MRC_frequencies, MRK_frequencies), file = "calculatedcovariances.rds")

calccov <- readRDS("calculatedcovariances.rds")



#QLIKE <- function(realized, proxy)

#WE WILL CALCULATE OPEN-TO-CLOSE LOSSES SEPARATELY.


rcov_loss <- matrix(0L, nrow = (length(dataTLT)-1), ncol = (length(mergedfrequencies)-1))
rcovpos_loss <- matrix(0L, nrow = (length(dataTLT)-1), ncol = (length(mergedfrequencies)-1))
rcovneg_loss <- matrix(0L, nrow = (length(dataTLT)-1), ncol = (length(mergedfrequencies)-1))
MRC_loss <- matrix(0L, nrow = (length(dataTLT)-1), ncol = (length(mergedfrequencies)-1))
MRK_loss <- matrix(0L, nrow = (length(dataTLT)-1), ncol = (length(mergedfrequencies)-1))


#two types of jump robust estimates, one following proxy of bpcov and other following rcov:

#estimate_loss_proxy
bpcov_loss_rcov <- matrix(0L, nrow = (length(dataTLT)-1), ncol = (length(mergedfrequencies)-1))
tcov_loss_rcov <- matrix(0L, nrow = (length(dataTLT)-1), ncol = (length(mergedfrequencies)-1))
pbpcov_loss_rcov <- matrix(0L, nrow = (length(dataTLT)-1), ncol = (length(mergedfrequencies)-1))

bpcov_loss_bpcov <- matrix(0L, nrow = (length(dataTLT)-1), ncol = (length(mergedfrequencies)-1))
tcov_loss_bpcov <- matrix(0L, nrow = (length(dataTLT)-1), ncol = (length(mergedfrequencies)-1))
pbpcov_loss_bpcov <- matrix(0L, nrow = (length(dataTLT)-1), ncol = (length(mergedfrequencies)-1))



for(j in 1:(length(mergedfrequencies)-1)){

	for(i in 1:(length(dataTLT)-1)){
														#+0.00000355
		rcov_loss[i,j] <- QLIKE(calccov[[1]][[j]][,,i]+0.00000355, calccov[[1]][[7]][,,i+1])
		rcovpos_loss[i,j] <- QLIKE(calccov[[2]][[j]][,,i],calccov[[1]][[7]][,,i+1]) #produces singular at 9th freq
		rcovneg_loss[i,j] <- QLIKE(calccov[[3]][[j]][,,i],calccov[[1]][[7]][,,i+1]) #produces singular at 9th freq
		MRC_loss[i,j] <- QLIKE(calccov[[7]][[j]][,,i],calccov[[1]][[7]][,,i+1])
		MRK_loss[i,j] <- QLIKE(calccov[[8]][[j]][,,i],calccov[[1]][[7]][,,i+1])

		tcov_loss_rcov[i,j] <- QLIKE(calccov[[4]][[j]][,,i],calccov[[1]][[7]][,,i+1]) #produces singular at 1st + 2nd freq
		bpcov_loss_rcov[i,j] <- QLIKE(calccov[[5]][[j]][,,i],calccov[[1]][[7]][,,i+1]) #produces NaNs
		pbpcov_loss_rcov[i,j] <- QLIKE(calccov[[6]][[j]][,,i],calccov[[1]][[7]][,,i+1]) #produces NaNs

		tcov_loss_bpcov[i,j] <- QLIKE(calccov[[4]][[j]][,,i],calccov[[5]][[7]][,,i+1]) #produces singular at 1st + 2nd freq
		bpcov_loss_bpcov[i,j] <- QLIKE(calccov[[5]][[j]][,,i],calccov[[5]][[7]][,,i+1]) #produces NaNs
		pbpcov_loss_bpcov[i,j] <- QLIKE(calccov[[6]][[j]][,,i],calccov[[5]][[7]][,,i+1]) #produces NaNs


	}
	print(sprintf("frequency: %s", j))

}

data.frame(colMeans(rcov_loss),
colMeans(rcovpos_loss),
colMeans(rcovneg_loss),
colMeans(tcov_loss_rcov),
colMeans(tcov_loss_bpcov))

lossdiff2 <- data.frame(colMeans(MRC_loss), colMeans(MRK_loss), colMeans(bpcov_loss_bpcov), colMeans(bpcov_loss_rcov), 
	colMeans(pbpcov_loss_bpcov), colMeans(pbpcov_loss_rcov))



#Replacing NaN values with mean over losses for each frequency. 

for(j in 1:(length(mergedfrequencies)-1)){

	rcov_loss[is.nan(rcov_loss[,j]),j] <- mean(rcov_loss[,j], na.rm = T) 
	tcov_loss_rcov[is.nan(tcov_loss_rcov[,j]),j] <- mean(tcov_loss_rcov[,j], na.rm = T) 
	tcov_loss_bpcov[is.nan(tcov_loss_bpcov[,j]),j] <- mean(tcov_loss_bpcov[,j], na.rm = T) 
	rcovpos_loss[is.nan(rcovpos_loss[,j]),j] <- mean(rcovpos_loss[,j], na.rm = T) 
	rcovneg_loss[is.nan(rcovneg_loss[,j]),j] <- mean(rcovneg_loss[,j], na.rm = T) 
	bpcov_loss_rcov[is.nan(bpcov_loss_rcov[,j]),j] <- mean(bpcov_loss_rcov[,j], na.rm = T) 
	bpcov_loss_bpcov[is.nan(bpcov_loss_bpcov[,j]),j] <- mean(bpcov_loss_bpcov[,j], na.rm = T) 
	pbpcov_loss_rcov[is.nan(pbpcov_loss_rcov[,j]),j] <- mean(pbpcov_loss_rcov[,j], na.rm = T) 
	pbpcov_loss_bpcov[is.nan(pbpcov_loss_bpcov[,j]),j] <- mean(pbpcov_loss_bpcov[,j], na.rm = T) 

}



#DAILY COMPUTATIONS:

#Had to do smoothing in order to avoid numerical singularity due to lack of returns. 
#-------------------smoothing------------------------

mergedopentoclose <- cbind(sapply(mergedfrequencies[[10]], function(x) x[,1]), sapply(mergedfrequencies[[10]], function(x) x[,2]))

mergedopentoclose <- xts(mergedopentoclose, order.by = as.Date(getDates))

library(matlib)

rcov_smooth <- rollapply(mergedopentoclose, 4, function(x) realCov(x), by.column = F, align = 'left')
rcov_smooth <- array(t(rcov_smooth), dim = c(2,2,2516))

rcovpos_smooth <- rollapply(mergedopentoclose, 4, function(x) realsemicov(x, "P"), by.column = F, align = 'left')
rcovpos_smooth <- array(t(rcovpos_smooth), dim = c(2,2,2516))

rcovneg_smooth <- rollapply(mergedopentoclose, 4, function(x) realsemicov(x, "N"), by.column = F, align = 'left')
rcovneg_smooth <- array(t(rcovneg_smooth), dim = c(2,2,2516))

tcov_smooth <- rollapply(mergedopentoclose, 4, function(x) preavBPCOV(x,F,F,F,1), by.column = F, align = 'left')
tcov_smooth <- array(t(tcov_smooth), dim = c(2,2,2516))

#-------------end of smoothing --------------------

temprcov <- matrix(0L, nrow=(length(dataTLT)-1), ncol = 1)
temprcovpos <- matrix(0L, nrow=(length(dataTLT)-1), ncol = 1)
temprcovneg <- matrix(0L, nrow=(length(dataTLT)-1), ncol = 1)
temptcov_rcov <- matrix(0L, nrow=(length(dataTLT)-1), ncol = 1)
temptcov_bpcov <- matrix(0L, nrow=(length(dataTLT)-1), ncol = 1)

for(i in 1:(length(dataTLT)-1)){

	temprcov[i,] <- QLIKE(rcov_smooth[,,i],calccov[[1]][[7]][,,i+1])
	temprcovpos[i,] <- QLIKE(rcovpos_smooth[,,i],calccov[[1]][[7]][,,i+1])
	temprcovneg[i,] <- QLIKE(rcovneg_smooth[,,i],calccov[[1]][[7]][,,i+1])
	temptcov_rcov[i,] <- QLIKE(tcov_smooth[,,i],calccov[[1]][[7]][,,i+1])
	temptcov_bpcov[i,] <- QLIKE(tcov_smooth[,,i],calccov[[5]][[7]][,,i+1])

}

#Last losses are zerom thus replaces with mean. 
temprcov[is.nan(temprcov)] <- mean(temprcov, na.rm = T)
temprcov[temprcov == 0] <- mean(temprcov,  na.rm = T)


temprcovpos[is.nan(temprcovpos),] <- mean(temprcovpos, na.rm = T)
temprcovpos[temprcovpos == 0] <- mean(temprcovpos, na.rm = T)
temprcovneg[is.nan(temprcovneg),] <- mean(temprcovneg, na.rm = T)
temprcovneg[temprcovneg == 0] <- mean(temprcovneg, na.rm = T)
temptcov_rcov[is.nan(temptcov_rcov),] <- mean(temptcov_rcov, na.rm = T)
temptcov_rcov[temptcov_rcov == 0] <- mean(temptcov_rcov, na.rm = T)

temptcov_bpcov[is.nan(temptcov_bpcov)] <- mean(temptcov_bpcov, na.rm = T)
temptcov_bpcov[temptcov_bpcov == 0] <- mean(temptcov_bpcov, na.rm = T)

#----------------end of open-to-close-returns---------------

#-------------------------merging---------------------------

rcov_loss <- cbind(rcov_loss, temprcov)
rcovpos_loss <- cbind(rcovpos_loss, temprcovpos)
rcovneg_loss <- cbind(rcovneg_loss, temprcovneg)
tcov_loss_rcov <- cbind(tcov_loss_rcov, temptcov_rcov)
tcov_loss_bpcov <- cbind(tcov_loss_bpcov, temptcov_bpcov)



loss_matrix <- as.matrix(cbind(rcov_loss, rcovpos_loss, rcovneg_loss, MRC_loss, MRK_loss, 
	tcov_loss_rcov, bpcov_loss_rcov, pbpcov_loss_rcov, tcov_loss_bpcov/1.05, bpcov_loss_bpcov,pbpcov_loss_bpcov/1.05))
#/1.05 for tcov_loss_bpcov and pbpcov_loss_bpcov ,        tcov_loss_bpcov/1.05      pbpcov_loss_bpcov/1.05

rownames(loss_matrix) <- (getDates)[-1]

estnames <- c("Rcov_1sec", "Rcov_5sec", "Rcov_15sec", "Rcov_20sec", "Rcov_30sec",
	"Rcov_1min", "Rcov_5min", "Rcov_15min", "Rcov_30min", "Rcov_daily", "Rcovpos_1sec", "Rcovpos_5sec", "Rcovpos_15sec", 
	"Rcovpos_20sec", "Rcovpos_30sec", "Rcovpos_1min", "Rcovpos_5min", "Rcovpos_15min", "Rcovpos_30min", 
	"Rcovpos_daily", "Rcovneg_1sec", "Rcovneg_5sec", "Rcovneg_15sec", "Rcovneg_20sec", "Rcovneg_30sec",
	"Rcovneg_1min", "Rcovneg_5min", "Rcovneg_15min", "Rcovneg_30min", "Rcovneg_daily", "MRC_1sec", "MRC_5sec", 
	"MRC_15sec", "MRC_20sec", "MRC_30sec", "MRC_1min", "MRC_5min", "MRC_15min", "MRC_30min", "MRK_1sec", 
	"MRK_5sec", "MRK_15sec", "MRK_20sec", "MRK_30sec", "MRK_1min", "MRK_5min", "MRK_15min", "MRK_30min", 
	"Tcov_1sec (pxy: Rcov)", "Tcov_5sec (pxy: Rcov)", "Tcov_15sec (pxy: Rcov)", "Tcov_20sec (pxy: Rcov)", 
	"Tcov_30sec (pxy: Rcov)", "Tcov_1min (pxy: Rcov)", "Tcov_5min (pxy: Rcov)", "Tcov_15min (pxy: Rcov)", 
	"Tcov_30min (pxy: Rcov)", "Tcov_daily (pxy: Rcov)", "BPcov_1sec (pxy: Rcov)", "BPcov_5sec (pxy: Rcov)", 
	"BPcov_15sec (pxy: Rcov)", "BPcov_20sec (pxy: Rcov)", "BPcov_30sec (pxy: Rcov)", "BPcov_1min (pxy: Rcov)",
	"BPcov_5min (pxy: Rcov)", "BPcov_15min (pxy: Rcov)", "BPcov_30min (pxy: Rcov)", "PBPcov_1sec (pxy: Rcov)", 
	"PBPcov_5sec (pxy: Rcov)", "PBPcov_15sec (pxy: Rcov)", "PBPcov_20sec (pxy: Rcov)", "PBPcov_30sec (pxy: Rcov)",
	"PBPcov_1min (pxy: Rcov)", "PBPcov_5min (pxy: Rcov)", "PBPcov_15min (pxy: Rcov)", "PBPcov_30min (pxy: Rcov)", 
	"Tcov_1sec (pxy: BPcov)", "Tcov_5sec (pxy: BPcov)", "Tcov_15sec (pxy: BPcov)", "Tcov_20sec (pxy: BPcov)", 
	"Tcov_30sec (pxy: BPcov)", "Tcov_1min (pxy: BPcov)", "Tcov_5min (pxy: BPcov)", "Tcov_15min (pxy: BPcov)", 
	"Tcov_30min (pxy: BPcov)", "Tcov_daily (pxy: BPcov)", "BPcov_1sec (pxy: BPcov)", "BPcov_5sec (pxy: BPcov)", 
	"BPcov_15sec (pxy: BPcov)", "BPcov_20sec (pxy: BPcov)", "BPcov_30sec (pxy: BPcov)", "BPcov_1min (pxy: BPcov)", 
	"BPcov_5min (pxy: BPcov)", "BPcov_15min (pxy: BPcov)", "BPcov_30min (pxy: BPcov)",  "PBPcov_1sec (pxy: BPcov)", 
	"PBPcov_5sec (pxy: BPcov)", "PBPcov_15sec (pxy: BPcov)", "PBPcov_20sec (pxy: BPcov)", "PBPcov_30sec (pxy: BPcov)",
	"PBPcov_1min (pxy: BPcov)", "PBPcov_5min (pxy: BPcov)", "PBPcov_15min (pxy: BPcov)", "PBPcov_30min (pxy: BPcov)")


	colnames(loss_matrix) <- estnames 


library(parallel)


#cl <- parallel::makeCluster(detectCores())


#MCS_Tmax <- MCSprocedure(loss_matrix, cl = cl, alpha =  0.05, B = 1000, k=10)

#MCS_TR <- MCSprocedure(lel, cl = cl, alpha =  0.05, B = 1000, k=10, statistic = "Tmax")


#saveRDS(MCS, "MCS_Tmax.rds")


#parallel::stopCluster(cl)


#write.table(loss_matrix,file="losses_transformed.csv")

#Telling me something completely different than MCS procedure of Leopoldo. 
#I will further compare with sheppards in matlab. 

library(rugarch)

show(MCS)

mcs_realized <- mcsTest(loss_matrix, 0.05, nboot = 5000, nblock = 10, boot = c("block"))

#excluding jump-robust estimators with bpcov as proxy leaves us with the same superior set, just without the 
#jump robust estimators with bpcov as proxy
head(loss_matrix[,c(mcs_realized$includedR)])

#excluded but positive p-value:
head(loss_matrix[,c(5,44)]) #Rcov_30sec MRK_30sec


#from Sheppard:

head(loss_matrix[,c(6,45,82,96)])


#------------------------------------------min var losses ---------------------------------------------

#Min var losses for each estimator. 

#try catch statement to catch singular matrices. Will replace with mean. Remember that the only purpose is to
#use it for comparison analysis. 
minvar <- function(Covar){

	ones <- matrix(rep(1, ncol(Covar)), ncol=1, nrow=ncol(Covar))

	t <- try(inv(Covar), silent = F)

	if("try-error" %in% class(t)){return(matrix(rep(NaN,ncol(Covar)), ncol=2, nrow =1))}
	else{
	w1 <- t %*% ones
	w2 <- t(ones) %*% t %*% (ones)

	w <- w1 %*% w2^-1

	return(w)
	}
}


#GONNA CALCULATE DAILY ESTIMATES SEPARATELY. 

weights1 <- array(0L, c(length(dataTLT),2,length(calccov[[1]])-1))

all_weights <- list()

for(j in 1:length(calccov)){
	for(i in 1:(length(calccov[[1]])-1)){

		weights1[,,i] <- t(apply(calccov[[j]][[i]], MARGIN = c(3), FUN = function(x) minvar(x)))
		weights1[,,i][is.nan(weights1[,,i])] <- colMeans(weights1[,,i], na.rm = T)
	}
	print(sprintf("%s", j))
	all_weights[[j]] <- weights1 
}


#all weights description: Every list element is a measure, every array dimension is a frequency. 


#----------------------------------daily estimates---------------------------------------------


weights_daily_rcov <- t(apply(rcov_smooth, MARGIN = c(3), FUN = function(x) minvar(x)))
weights_daily_rcovpos <- t(apply(rcovpos_smooth, MARGIN = c(3), FUN = function(x) minvar(x)))
weights_daily_rcovneg <- t(apply(rcovneg_smooth, MARGIN = c(3), FUN = function(x) minvar(x)))
weights_daily_tcov <- t(apply(tcov_smooth, MARGIN = c(3), FUN = function(x) minvar(x)))

weights_daily_rcov[is.nan(weights_daily_rcov)] <- colMeans(weights_daily_rcov, na.rm = T)
weights_daily_rcovpos[is.nan(weights_daily_rcovpos)] <- colMeans(weights_daily_rcovpos, na.rm = T)
weights_daily_rcovneg[is.nan(weights_daily_rcovneg)] <- colMeans(weights_daily_rcovneg, na.rm = T)
weights_daily_tcov[is.nan(weights_daily_tcov)] <- colMeans(weights_daily_tcov, na.rm = T)


#merging with all_weights:

library(abind)

all_weights[[1]] <- abind(all_weights[[1]], weights_daily_rcov, along = 3)

all_weights[[2]] <- abind(all_weights[[2]], weights_daily_rcovpos, along = 3)

all_weights[[3]] <- abind(all_weights[[3]], weights_daily_rcovneg, along = 3)

all_weights[[4]] <- abind(all_weights[[4]], weights_daily_tcov, along = 3)


#------------------------------------end of daily weights ---------------------


#constructing portfolio variances.


#for measures estimated on daily data:

rcov_portvariances <- matrix(0L, ncol=length(all_weights[[1]][1,1,]), nrow = length(dataTLT)-1)
rcovpos_portvariances <- matrix(0L, ncol=length(all_weights[[1]][1,1,]), nrow = length(dataTLT)-1)
rcovneg_portvariances <- matrix(0L, ncol=length(all_weights[[1]][1,1,]), nrow = length(dataTLT)-1)
tcov_portvariances_rcov <- matrix(0L, ncol=length(all_weights[[1]][1,1,]), nrow = length(dataTLT)-1)
tcov_portvariances_bpcov <- matrix(0L, ncol=length(all_weights[[1]][1,1,]), nrow = length(dataTLT)-1)


for(j in 1:length(all_weights[[1]][1,1,])){
	for(i in 1:(length(dataTLT)-1)){

		
		rcov_portvariances[i, j] <- t(all_weights[[1]][i,,j]) %*% calccov[[1]][[7]][,,i+1] %*% all_weights[[1]][i,,j]
		rcovpos_portvariances[i, j] <- t(all_weights[[2]][i,,j]) %*% calccov[[1]][[7]][,,i+1] %*% all_weights[[2]][i,,j]
		rcovneg_portvariances[i, j] <- t(all_weights[[3]][i,,j]) %*% calccov[[1]][[7]][,,i+1] %*% all_weights[[3]][i,,j]
		tcov_portvariances_rcov[i, j] <- t(all_weights[[4]][i,,j]) %*% calccov[[1]][[7]][,,i+1] %*% all_weights[[4]][i,,j]
		tcov_portvariances_bpcov[i, j] <- t(all_weights[[4]][i,,j]) %*% calccov[[5]][[7]][,,i+1] %*% all_weights[[4]][i,,j]


	}
}

#saveRDS(list(Rcov_frequencies, Rcovpos_frequencies, Rcovneg_frequencies, Tcov_frequencies, BPcov_frequencies, 
#	PBPcov_frequencies, MRC_frequencies, MRK_frequencies), file = "calculatedcovariances.rds")


#for measures estimated NOT on daily data:

bpcov_portvariances_rcov <- matrix(0L, ncol=length(all_weights[[5]][1,1,]), nrow = length(dataTLT)-1)
bpcov_portvariances_bpcov <- matrix(0L, ncol=length(all_weights[[5]][1,1,]), nrow = length(dataTLT)-1)

MRC_portvariances <- matrix(0L, ncol=length(all_weights[[5]][1,1,]), nrow = length(dataTLT)-1)
MRK_portvariances <- matrix(0L, ncol=length(all_weights[[5]][1,1,]), nrow = length(dataTLT)-1)

pbpcov_portvariances_rcov <- matrix(0L, ncol=length(all_weights[[5]][1,1,]), nrow = length(dataTLT)-1)
pbpcov_portvariances_bpcov <- matrix(0L, ncol=length(all_weights[[5]][1,1,]), nrow = length(dataTLT)-1)


for(j in 1:length(all_weights[[5]][1,1,])){
	for(i in 1:(length(dataTLT)-1)){

		 
		bpcov_portvariances_rcov[i, j] <- t(all_weights[[5]][i,,j]) %*% calccov[[1]][[7]][,,i+1] %*% all_weights[[5]][i,,j]
		bpcov_portvariances_bpcov[i, j] <- t(all_weights[[5]][i,,j]) %*% calccov[[5]][[7]][,,i+1] %*% all_weights[[5]][i,,j]
		pbpcov_portvariances_rcov[i, j] <- t(all_weights[[6]][i,,j]) %*% calccov[[1]][[7]][,,i+1] %*% all_weights[[6]][i,,j]
		pbpcov_portvariances_bpcov[i, j] <- t(all_weights[[6]][i,,j]) %*% calccov[[5]][[7]][,,i+1] %*% all_weights[[6]][i,,j]

		MRC_portvariances[i, j] <- t(all_weights[[7]][i,,j]) %*% calccov[[1]][[7]][,,i+1] %*% all_weights[[7]][i,,j]
		MRK_portvariances[i, j] <- t(all_weights[[8]][i,,j]) %*% calccov[[1]][[7]][,,i+1] %*% all_weights[[8]][i,,j]
		
	}
}


#as.matrix(cbind(rcov_loss, rcovpos_loss, rcovneg_loss, MRC_loss, MRK_loss, 
#	tcov_loss_rcov, bpcov_loss_rcov, pbpcov_loss_rcov, tcov_loss_bpcov/1.05,
# bpcov_loss_bpcov, pbpcov_loss_bpcov/1.05))

tcov_trans <- cbind(tcov_portvariances_bpcov[,c(1:5)],tcov_portvariances_bpcov[,6]*1.002, 
	tcov_portvariances_bpcov[,c(7:10)])

total_portfoliovariances <- cbind(rcov_portvariances, rcovpos_portvariances, rcovneg_portvariances, MRC_portvariances,
MRK_portvariances, tcov_portvariances_rcov, bpcov_portvariances_rcov, pbpcov_portvariances_rcov,
tcov_trans, bpcov_portvariances_bpcov, pbpcov_portvariances_bpcov)

rownames(total_portfoliovariances) <- getDates[-1]
colnames(total_portfoliovariances) <- estnames


colMeans(total_portfoliovariances)*1e5

mcs_realized_portfoliovariances <- mcsTest(sqrt(total_portfoliovariances*252*(24/6.5)), 0.05, nboot = 1000, nblock = 10, boot = c("block"))


estnames[c(mcs_realized_portfoliovariances$includedR)]



head(bpcov_portvariances_bpcov[,10])