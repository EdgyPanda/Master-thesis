#will contain the summary statistics following Hautsch & Podolskij (2013). 

library(xts)
library(highfrequency)
library(matlib)
library(tidyverse)

dataTLT <- readRDS("dataTLT.rds")
dataSPY <- readRDS("dataSPY.rds")



#-------------------------------Calculating average time in seconds between each trade:
timeobjectsTLT <- list()
timeobjectsSPY <- list()

for(i in 1:length(dataTLT)){

	timeobjectsTLT[[i]] <- c(index(dataTLT[[i]]))
	timeobjectsSPY[[i]] <- c(index(dataSPY[[i]]))

}


for(i in 1:length(dataTLT)){

	timeobjectsTLT[[i]] <- unlist(lapply((str_split(timeobjectsTLT[[i]], " ")), function(x) x[2]))
	timeobjectsSPY[[i]] <- unlist(lapply((str_split(timeobjectsSPY[[i]], " ")), function(x) x[2]))


}

for(i in 1:length(dataTLT)){

	timeobjectsTLT[[i]] <- as.numeric(gsub(":", "",timeobjectsTLT[[i]]))
	timeobjectsSPY[[i]] <- as.numeric(gsub(":", "",timeobjectsSPY[[i]]))

}

seconddiffTLT <- list()
seconddiffSPY <- list()

for(i in 1:length(dataTLT)){

	seconddiffTLT[[i]] <- c(
	diff(timeobjectsTLT[[i]][timeobjectsTLT[[i]]<100000]),
	diff(timeobjectsTLT[[i]][100000<timeobjectsTLT[[i]] & timeobjectsTLT[[i]]<110000]),
	diff(timeobjectsTLT[[i]][110000<timeobjectsTLT[[i]] & timeobjectsTLT[[i]]<120000]),
	diff(timeobjectsTLT[[i]][120000<timeobjectsTLT[[i]] & timeobjectsTLT[[i]]<130000]),
	diff(timeobjectsTLT[[i]][130000<timeobjectsTLT[[i]] & timeobjectsTLT[[i]]<140000]),
	diff(timeobjectsTLT[[i]][140000<timeobjectsTLT[[i]] & timeobjectsTLT[[i]]<150000]),
	diff(timeobjectsTLT[[i]][150000<timeobjectsTLT[[i]] & timeobjectsTLT[[i]]<160000]))


	seconddiffSPY[[i]] <- c(
	diff(timeobjectsSPY[[i]][timeobjectsSPY[[i]]<100000]),
	diff(timeobjectsSPY[[i]][100000<timeobjectsSPY[[i]] & timeobjectsSPY[[i]]<110000]),
	diff(timeobjectsSPY[[i]][110000<timeobjectsSPY[[i]] & timeobjectsSPY[[i]]<120000]),
	diff(timeobjectsSPY[[i]][120000<timeobjectsSPY[[i]] & timeobjectsSPY[[i]]<130000]),
	diff(timeobjectsSPY[[i]][130000<timeobjectsSPY[[i]] & timeobjectsSPY[[i]]<140000]),
	diff(timeobjectsSPY[[i]][140000<timeobjectsSPY[[i]] & timeobjectsSPY[[i]]<150000]),
	diff(timeobjectsSPY[[i]][150000<timeobjectsSPY[[i]] & timeobjectsSPY[[i]]<160000]))

}


#Years are changing at, 252, 504, 754, 1006, 1258, 1510, 1762, 2013, 2264, 2516(end)

meantimeTLT <- unlist(lapply(seconddiffTLT, mean))
meantimeSPY <- unlist(lapply(seconddiffSPY, mean))


meantime20102012TLT <- mean(meantimeTLT[1:504])
meantime20122014TLT <- mean(meantimeTLT[505:1006])
meantime20142016TLT <- mean(meantimeTLT[1007:1510])
meantime20162018TLT <- mean(meantimeTLT[1511:2013])
meantime20182020TLT <- mean(meantimeTLT[2014:2516])

meantime20102012TLT
meantime20122014TLT
meantime20142016TLT
meantime20162018TLT
meantime20182020TLT

meantime20102012SPY <- mean(meantimeSPY[1:504])
meantime20122014SPY <- mean(meantimeSPY[505:1006])
meantime20142016SPY <- mean(meantimeSPY[1007:1510])
meantime20162018SPY <- mean(meantimeSPY[1511:2013])
meantime20182020SPY <- mean(meantimeSPY[2014:2516])

meantime20102012SPY 
meantime20122014SPY 
meantime20142016SPY 
meantime20162018SPY 
meantime20182020SPY 


#average daily transaction time in seconds:
meantimeTLTTOTAL <- mean(meantimeTLT)
meantimeSPYTOTAL <- mean(meantimeSPY)

meantimeTLTTOTAL
meantimeSPYTOTAL



#---------------------------------------Average number of trades-------------------

rawtradesTLT <- numeric()
rawtradesSPY <- numeric()

for(i in 1:length(dataTLT)){


	rawtradesTLT[i] <- length(dataTLT[[i]])
	rawtradesSPY[i] <- length(dataSPY[[i]])

}

rawtradesTLT20102012 <- mean(rawtradesTLT[1:504])
rawtradesTLT20122014 <- mean(rawtradesTLT[505:1006])
rawtradesTLT20142016 <- mean(rawtradesTLT[1007:1510])
rawtradesTLT20162018 <- mean(rawtradesTLT[1511:2013])
rawtradesTLT20182020 <- mean(rawtradesTLT[2014:2516])
rawtradesTLTtotal <- mean(rawtradesTLT)

rawtradesTLT20102012
rawtradesTLT20122014
rawtradesTLT20142016 
rawtradesTLT20162018 
rawtradesTLT20182020 
rawtradesTLTtotal

rawtradesSPY20102012 <- mean(rawtradesSPY[1:504])
rawtradesSPY20122014 <- mean(rawtradesSPY[505:1006])
rawtradesSPY20142016 <- mean(rawtradesSPY[1007:1510])
rawtradesSPY20162018 <- mean(rawtradesSPY[1511:2013])
rawtradesSPY20182020 <- mean(rawtradesSPY[2014:2516])
rawtradesSPYtotal <- mean(rawtradesSPY)

rawtradesSPY20102012
rawtradesSPY20122014
rawtradesSPY20142016 
rawtradesSPY20162018 
rawtradesSPY20182020 
rawtradesSPYtotal


#------------------------------Average number of cand.outliers (candidate outliers)----------------------
#the average number of daily candidate outliers as identified by the filtering procedure.


SummarystatsTLT <- read.csv("SummaryStatisticsforcleaning_TLT.csv", header = T)
SummarystatsSPY <- read.csv("SummaryStatisticsforcleaning_SPY.csv", header = T)


outliersSPY20102012 <- mean(SummarystatsSPY$T4[1:504])
outliersSPY20122014 <- mean(SummarystatsSPY$T4[505:1006])
outliersSPY20142016 <- mean(SummarystatsSPY$T4[1007:1510])
outliersSPY20162018 <- mean(SummarystatsSPY$T4[1511:2013])
outliersSPY20182020 <- mean(SummarystatsSPY$T4[2014:2513])
outliersSPYtotal <- mean(SummarystatsSPY$T4)


outliersSPY20102012
outliersSPY20122014
outliersSPY20142016
outliersSPY20162018
outliersSPY20182020
outliersSPYtotal

outliersTLT20102012 <- mean(SummarystatsTLT$T4[1:504])
outliersTLT20122014 <- mean(SummarystatsTLT$T4[505:1006])
outliersTLT20142016 <- mean(SummarystatsTLT$T4[1007:1510])
outliersTLT20162018 <- mean(SummarystatsTLT$T4[1511:2013])
outliersTLT20182020 <- mean(SummarystatsTLT$T4[2014:2515])
outliersTLTtotal <- mean(SummarystatsTLT$T4)


outliersTLT20102012
outliersTLT20122014
outliersTLT20142016
outliersTLT20162018
outliersTLT20182020
outliersTLTtotal

#----------------------------------------------------raw trades ---------------------------------------

#Raw trades/quotes is the total number of data available from these exchanges during the
#trading session, while  #trades/quotes is the total sample remaining after filtering the data


allcleanedTLT <- rowSums(SummarystatsTLT[,c(1:4,8,10,12)])
allcleanedSPY <- rowSums(SummarystatsSPY[,c(1:4,8,10,12)])

datalengthTLT <- numeric()
datalengthSPY <- numeric()

for(i in 1:length(dataTLT)){

	datalengthTLT[i] <- length(dataTLT[[i]])
	datalengthSPY[i] <- length(dataSPY[[i]])

}


#missing the last days due to error in data. 
rawtradesTLT <- allcleanedTLT + datalengthTLT[1:2515]
rawtradesSPY <- allcleanedSPY + datalengthSPY[1:2513]

#The average amount of raw trades for each day during the trading session.


mean(rawtradesTLT[1:504])
mean(rawtradesTLT[505:1006])
mean(rawtradesTLT[1007:1510])
mean(rawtradesTLT[1511:2013])
mean(rawtradesTLT[2014:2515])


mean(rawtradesTLT)


mean(rawtradesSPY[1:504]) +
mean(rawtradesSPY[505:1006]) +
mean(rawtradesSPY[1007:1510]) +
mean(rawtradesSPY[1511:2013]) +
mean(rawtradesSPY[2014:2513]) 


mean(rawtradesSPY)




#------------------------------------Average proportion of non-zero trade returns----------------------------
#Done by constructing return data from all of the transaction prices and can be based on TRTS employing
#all transactions. 
#Intuitively, while the proportion might be close to the same proportion for a 1 second CTS scheme, we know that sparse
#sparse sampling further, will likely increase the proportion. Ie. the non-zero trades also indirectly describes the liquidity of our
#financial assets. 


#YOU COULD DO CTS 1 SEC AND 5 SEC AND SHOW THE INCREASING PROPORTION OF NON-ZERO RETURNS TO GET A BETTER OVERVIEW?

rawreturnsSPY <- list()

rawreturnsTLT <- list()

for (i in 1:length(dataTLT)){

	rawreturnsSPY[[i]] <- diff(log(dataSPY[[i]]))[-1]
	rawreturnsTLT[[i]] <- diff(log(dataTLT[[i]]))[-1]
}


TLT_1sec <- list()
TLT_5sec <- list()


SPY_1sec <- list()
SPY_5sec <- list()


for(i in 1:length(dataTLT)){

	TLT_1sec[[i]] <- diff(log(aggregatets(dataTLT[[i]], on = "seconds", k=1)))[-1]
	TLT_5sec[[i]] <- diff(log(aggregatets(dataTLT[[i]], on = "seconds", k=5)))[-1]

	SPY_1sec[[i]] <- diff(log(aggregatets(dataSPY[[i]], on = "seconds", k=1)))[-1]
	SPY_5sec[[i]] <- diff(log(aggregatets(dataSPY[[i]], on = "seconds", k=5)))[-1]

}



nonzeroSPY <- vector()
nonzeroSPY_1sec <- vector()
nonzeroSPY_5sec <- vector()


nonzeroTLT <- vector()
nonzeroTLT_1sec <- vector()
nonzeroTLT_5sec <- vector()

SPYtotal <- numeric()
TLTtotal <- numeric()

for (i in 1:length(dataTLT)){

nonzeroSPY[i] <- apply(rawreturnsSPY[[i]], 2, function(c)sum(c!=0))
nonzeroSPY_1sec[i] <- apply(SPY_1sec[[i]], 2, function(c)sum(c!=0))
nonzeroSPY_5sec[i] <- apply(SPY_5sec[[i]], 2, function(c)sum(c!=0))

nonzeroTLT[i] <- apply(rawreturnsTLT[[i]], 2, function(c)sum(c!=0))
nonzeroTLT_1sec[i] <- apply(TLT_1sec[[i]], 2, function(c)sum(c!=0))
nonzeroTLT_5sec[i] <- apply(TLT_5sec[[i]], 2, function(c)sum(c!=0))


SPYtotal[i] <- dim(rawreturnsSPY[[i]])[1]
SPYtotal_1sec[i] <- dim(SPY_1sec[[i]])[1]
SPYtotal_5sec[i] <- dim(SPY_5sec[[i]])[1]


TLTtotal[i] <- dim(rawreturnsTLT[[i]])[1]
TLTtotal_1sec[i] <- dim(TLT_1sec[[i]])[1]
TLTtotal_5sec[i] <- dim(TLT_5sec[[i]])[1]

}


#(daily) average nonzero trade returns.


nonzerotrades20102012SPY <- sum(nonzeroSPY[1:504])/(sum(SPYtotal[1:504]))
nonzerotrades20122014SPY <- sum(nonzeroSPY[505:1006])/(sum(SPYtotal[505:1006]))
nonzerotrades20142016SPY <- sum(nonzeroSPY[1007:1510])/(sum(SPYtotal[1007:1510]))
nonzerotrades20162018SPY <- sum(nonzeroSPY[1511:2013])/(sum(SPYtotal[1511:2013]))
nonzerotrades20182020SPY <- sum(nonzeroSPY[2014:2516])/(sum(SPYtotal[2014:2516]))

nonzerotrades20102012TLT <- sum(nonzeroTLT[1:504])/(sum(TLTtotal[1:504]))
nonzerotrades20122014TLT <- sum(nonzeroTLT[505:1006])/(sum(TLTtotal[505:1006]))
nonzerotrades20142016TLT <- sum(nonzeroTLT[1007:1510])/(sum(TLTtotal[1007:1510]))
nonzerotrades20162018TLT <- sum(nonzeroTLT[1511:2013])/(sum(TLTtotal[1511:2013]))
nonzerotrades20182020TLT <- sum(nonzeroTLT[2014:2516])/(sum(TLTtotal[2014:2516]))


# 1 sec 

nonzerotrades20102012SPY_1sec <- sum(nonzeroSPY_1sec[1:504])/(sum(SPYtotal_1sec[1:504]))
nonzerotrades20122014SPY_1sec <- sum(nonzeroSPY_1sec[505:1006])/(sum(SPYtotal_1sec[505:1006]))
nonzerotrades20142016SPY_1sec <- sum(nonzeroSPY_1sec[1007:1510])/(sum(SPYtotal_1sec[1007:1510]))
nonzerotrades20162018SPY_1sec <- sum(nonzeroSPY_1sec[1511:2013])/(sum(SPYtotal_1sec[1511:2013]))
nonzerotrades20182020SPY_1sec <- sum(nonzeroSPY_1sec[2014:2516])/(sum(SPYtotal_1sec[2014:2516]))


nonzerotrades20102012TLT_1sec <- sum(nonzeroTLT_1sec[1:504])/(sum(TLTtotal_1sec[1:504]))
nonzerotrades20122014TLT_1sec <- sum(nonzeroTLT_1sec[505:1006])/(sum(TLTtotal_1sec[505:1006]))
nonzerotrades20142016TLT_1sec <- sum(nonzeroTLT_1sec[1007:1510])/(sum(TLTtotal_1sec[1007:1510]))
nonzerotrades20162018TLT_1sec <- sum(nonzeroTLT_1sec[1511:2013])/(sum(TLTtotal_1sec[1511:2013]))
nonzerotrades20182020TLT_1sec <- sum(nonzeroTLT_1sec[2014:2516])/(sum(TLTtotal_1sec[2014:2516]))


# 5 sec 

nonzerotrades20102012SPY_5sec <- sum(nonzeroSPY_5sec[1:504])/(sum(SPYtotal_5sec[1:504]))
nonzerotrades20122014SPY_5sec <- sum(nonzeroSPY_5sec[505:1006])/(sum(SPYtotal_5sec[505:1006]))
nonzerotrades20142016SPY_5sec <- sum(nonzeroSPY_5sec[1007:1510])/(sum(SPYtotal_5sec[1007:1510]))
nonzerotrades20162018SPY_5sec <- sum(nonzeroSPY_5sec[1511:2013])/(sum(SPYtotal_5sec[1511:2013]))
nonzerotrades20182020SPY_5sec <- sum(nonzeroSPY_5sec[2014:2516])/(sum(SPYtotal_5sec[2014:2516]))


nonzerotrades20102012TLT_5sec <- sum(nonzeroTLT_5sec[1:504])/(sum(TLTtotal_5sec[1:504]))
nonzerotrades20122014TLT_5sec <- sum(nonzeroTLT_5sec[505:1006])/(sum(TLTtotal_5sec[505:1006]))
nonzerotrades20142016TLT_5sec <- sum(nonzeroTLT_5sec[1007:1510])/(sum(TLTtotal_5sec[1007:1510]))
nonzerotrades20162018TLT_5sec <- sum(nonzeroTLT_5sec[1511:2013])/(sum(TLTtotal_5sec[1511:2013]))
nonzerotrades20182020TLT_5sec <- sum(nonzeroTLT_5sec[2014:2516])/(sum(TLTtotal_5sec[2014:2516]))


#total
(sum(nonzeroSPY) / sum(SPYtotal))

(sum(nonzeroTLT) / sum(TLTtotal))


nonzerotradeSPY2018_1sec <- sum(nonzeroSPY_1sec)/(sum(SPYtotal_1sec))
nonzerotradeTLT2018_1sec <- sum(nonzeroTLT_1sec)/(sum(TLTtotal_1sec))

nonzerotradeSPY2018_5sec <- sum(nonzeroSPY_5sec)/(sum(SPYtotal_5sec))
nonzerotradeTLT2018_5sec <- sum(nonzeroTLT_5sec)/(sum(TLTtotal_5sec))

#--------------------------------------------------NOISE ESTIMATOR-----------------------------------
#This is based on based on TRTS employing all transactions.
#
#
#We use a slightly modified function other than the previous year. Ie, we use a sparse sampled (5-min) bipower variation
#as the estimator for IV instead of the MLRV proposed in the podolskij hautsch paper. 
#
#
#

source("Previous functions/projectfunctions.R")

TLT_raw <- list()

SPY_raw <- list()

TLT_1min <- list()

SPY_1min <- list()

for(i in 1:length(dataTLT)){

	TLT_raw[[i]] <- diff(log(dataTLT[[i]]))[-1]

    SPY_raw[[i]] <- diff(log(dataSPY[[i]]))[-1]

    TLT_1min[[i]] <- diff(log(aggregatets(dataTLT[[i]], on="minutes",k=1)))[-1]

    SPY_1min[[i]] <- diff(log(aggregatets(dataSPY[[i]], on="minutes", k=1)))[-1]
}


library(highfrequency)
#only univariate lists 
noisetosignal2 <- function(list,sparse_data,  names = NULL){

	n <- as.vector(sapply(list, length))

	#IVest <- MLRV(list, names)
	IVest <- matrix(0L, nrow = length(list), ncol = 1)

	for(i in 1:length(list)){

		IVest[i] <- preavBPCOV(sparse_data[[i]], F, F, F, theta = 1)

	}

	#a bar calculations
	condition <- as.numeric(vector())

	#sum only takes numeric vectors, therefore saving in list object.
	ret <- list() 
	for ( i in 1:length(list)){

		ret[[i]] <- as.vector(list[[i]])

	}

	for (i in 1:length(list)){

	condition[i] <- t(ret[[i]][2:n[i]]) %*% ret[[i]][1:(n[i]-1)]

	}

	abar <- matrix(0L, nrow = length(list), ncol=1)

	for(i in 1:length(list)){

		if(condition[i] < 0){

			abar[i, 1] <-  -1/(n[i]-1) * t(ret[[i]][1:(n[i]-1)]) %*% ret[[i]][2:n[i]]
		}
		else{

			abar[i, 1] <-  (1/(2*n[i]) * realCov(ret[[i]]))

		}

	}

	#eps <- abar %*% 1/(IVest %*% 1/n)
	eps <-  abar / (IVest / n)


	return(eps)

}


dailynoiseTLT <- noisetosignal2(TLT_raw, TLT_1min)
dailynoiseTLTMLRV <- noisetosignal(TLT_raw)

noiseTLT20102012 <- mean(dailynoiseTLT[1:504])
noiseTLT20122014 <- mean(dailynoiseTLT[505:1006])
noiseTLT20142016 <- mean(dailynoiseTLT[1007:1510])
noiseTLT20162018 <- mean(dailynoiseTLT[1511:2013])
noiseTLT20182020 <- mean(dailynoiseTLT[2014:2516])


noiseTLT20102012
noiseTLT20122014 
noiseTLT20142016
noiseTLT20162018
noiseTLT20182020 

#Total
mean(dailynoiseTLT)
mean(dailynoiseTLTMLRV)





dailynoiseSPY <- noisetosignal2(SPY_raw, SPY_1min)
dailynoiseSPYMLRV <- noisetosignal(SPY_raw)

#729

noisetosignal2(list(SPY_raw[[729]]), list(SPY_5min[[729]]))

noisetosignal2(list(SPY_raw[[100]]), list(SPY_5min[[100]]))



noiseSPY20102012 <- mean(dailynoiseSPY[1:504])
noiseSPY20122014 <- mean(dailynoiseSPY[505:1006])
noiseSPY20142016 <- mean(dailynoiseSPY[1007:1510])
noiseSPY20162018 <- mean(dailynoiseSPY[1511:2013])
noiseSPY20182020 <- mean(dailynoiseSPY[2014:2516])


noiseSPY20102012
noiseSPY20122014 
noiseSPY20142016
noiseSPY20162018
noiseSPY20182020 

mean(dailynoiseSPY)
mean(dailynoiseSPYMLRV)



#---------------------------------exhanges--------------------------------------------
extlt <- read.csv("ExchangeStatistics_TLT.csv")
exspy <- read.csv("ExchangeStatistics_SPY.csv")

#TLT


#construct ad-hoc function that gives me the results for 2 year increments. 

#Years are changing at, 252, 504, 754, 1006, 1258, 1510, 1762, 2013, 2264, 2516(end)

main.ExchangesTLT <- function(sequence){

	tqtlt <- length(grep("T  Q", extlt[sequence,1])) + length(grep("T  Q", extlt[sequence,2])) + 
	length(grep("T  Q", extlt[sequence,3]))

	wyjktlt <- length(grep("W  Y  J  K", extlt[sequence,1])) + length(grep("W  Y  J  K", extlt[sequence,2])) + 
	length(grep("W  Y  J  K", extlt[sequence,3]))

	ptlt <- length(grep("P", extlt[sequence,1])) + length(grep("P", extlt[sequence,2])) + 
	length(grep("P", extlt[sequence,3]))

	ztlt <- length(grep("Z", extlt[sequence,1])) + length(grep("Z", extlt[sequence,2])) + 
	length(grep("Z", extlt[sequence,3]))

	dtlt <- length(grep("D", extlt[sequence,1])) + length(grep("D", extlt[sequence,2])) + 
	length(grep("D", extlt[sequence,3]))

	btlt <- length(grep("B", extlt[sequence,1])) + length(grep("B", extlt[sequence,2])) + 
	length(grep("B", extlt[sequence,3]))

	res <- matrix(c(tqtlt, wyjktlt, ptlt, ztlt, dtlt, btlt), ncol = 6, nrow = 1)

	colnames(res) <- c("tqtlt", "wyjktlt", "ptlt", "ztlt", "dtlt", "btlt")
	return(res)


}

main.ExchangesTLT(1:504)

main.ExchangesTLT(505:1006)

main.ExchangesTLT(1007:1510)

main.ExchangesTLT(1511:2013)

main.ExchangesTLT(2014:2516)

#total 

main.ExchangesTLT(1:2516)




#SPY


main.ExchangesSPY <- function(sequence){

	tqspy <- length(grep("T  Q", exspy[sequence,1])) + length(grep("T  Q", exspy[sequence,2])) + 
	length(grep("T  Q", exspy[sequence,3]))

	wyjkspy <- length(grep("W  Y  J  K", exspy[sequence,1])) + length(grep("W  Y  J  K", exspy[sequence,2])) + 
	length(grep("W  Y  J  K", exspy[sequence,3]))

	pspy <- length(grep("P", exspy[sequence,1])) + length(grep("P", exspy[sequence,2])) + 
	length(grep("P", exspy[sequence,3]))

	zspy <- length(grep("Z", exspy[sequence,1])) + length(grep("Z", exspy[sequence,2])) + 
	length(grep("Z", exspy[sequence,3]))

	dspy <- length(grep("D", exspy[sequence,1])) + length(grep("D", exspy[sequence,2])) + 
	length(grep("D", exspy[sequence,3]))

	bspy <- length(grep("B", exspy[sequence,1])) + length(grep("B", exspy[sequence,2])) + 
	length(grep("B", exspy[sequence,3]))

	res <- matrix(c(tqspy, wyjkspy, pspy, zspy, dspy, bspy), ncol=6, nrow = 1)
	colnames(res) <- c("tqspy", "wyjkspy", "pspy", "zspy", "dspy", "bspy")

	return(res)

}

main.ExchangesSPY(1:504)

main.ExchangesSPY(505:1006)

main.ExchangesSPY(1007:1510)

main.ExchangesSPY(1511:2013)

main.ExchangesSPY(2014:2516)

#total 

main.ExchangesSPY(1:2516)
