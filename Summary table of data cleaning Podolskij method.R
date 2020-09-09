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
#the average number of candidate outliers  as identified by the filtering procedure.


SummarystatsTLT <- read.csv("SummaryStatisticsforcleaning_TLT.csv", header = T)
SummarystatsSPY <- read.csv("SummaryStatisticsforcleaning_SPY.csv", header = T)


outliersSPY20102012 <- sum(SummarystatsSPY$T4[1:504])
outliersSPY20122014 <- sum(SummarystatsSPY$T4[505:1006])
outliersSPY20142016 <- sum(SummarystatsSPY$T4[1007:1510])
outliersSPY20162018 <- sum(SummarystatsSPY$T4[1511:2013])
outliersSPY20182020 <- sum(SummarystatsSPY$T4[2014:2513])
outliersSPYtotal <- sum(SummarystatsSPY$T4)


outliersSPY20102012
outliersSPY20122014
outliersSPY20142016
outliersSPY20162018
outliersSPY20182020
outliersSPYtotal

outliersTLT20102012 <- sum(SummarystatsTLT$T4[1:504])
outliersTLT20122014 <- sum(SummarystatsTLT$T4[505:1006])
outliersTLT20142016 <- sum(SummarystatsTLT$T4[1007:1510])
outliersTLT20162018 <- sum(SummarystatsTLT$T4[1511:2013])
outliersTLT20182020 <- sum(SummarystatsTLT$T4[2014:2515])
outliersTLTtotal <- sum(SummarystatsTLT$T4)


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
mean(rawtradesTLT)
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



nonzeroSPY <- vector()

nonzeroTLT <- vector()


SPYtotal <- numeric()
TLTtotal <- numeric()

for (i in 1:length(dataTLT)){

nonzeroSPY[i] <- apply(rawreturnsSPY[[i]], 2, function(c)sum(c!=0))

nonzeroTLT[i] <- apply(rawreturnsTLT[[i]], 2, function(c)sum(c!=0))

SPYtotal[i] <- dim(rawreturnsSPY[[i]])[1]

TLTtotal[i] <- dim(rawreturnsTLT[[i]])[1]

}


#(daily) average nonzero trade returns.


nonzerotrades20102012SPY <- sum(nonzeroSPY[1:504])/(sum(SPYtotal[1:504]))
nonzerotrades20122014SPY <- sum(nonzeroSPY[505:1006])/(sum(SPYtotal[505:1006]))
nonzerotrades20142016SPY <- sum(nonzeroSPY[1007:1510])/(sum(SPYtotal[1007:1510]))
nonzerotrades20162018SPY <- sum(nonzeroSPY[1511:2013])/(sum(SPYtotal[1511:2013]))
nonzerotrades20182020SPY <- sum(nonzeroSPY[2014:2516])/(sum(SPYtotal[2014:2516]))

nonzerotrades20102012TLT <- sum(nonzeroSPY[1:504])/(sum(SPYtotal[1:504]))
nonzerotrades20122014TLT <- sum(nonzeroSPY[505:1006])/(sum(SPYtotal[505:1006]))
nonzerotrades20142016TLT <- sum(nonzeroSPY[1007:1510])/(sum(SPYtotal[1007:1510]))
nonzerotrades20162018TLT <- sum(nonzeroSPY[1511:2013])/(sum(SPYtotal[1511:2013]))
nonzerotrades20182020TLT <- sum(nonzeroSPY[2014:2516])/(sum(SPYtotal[2014:2516]))

#As a test, here's 2018: (CORRECT)

nonzerotradeSPY2018 <- sum(nonzeroSPY[2014:2264])/(sum(SPYtotal[2014:2264]))
nonzerotradeTLT2018 <- sum(nonzeroTLT[2014:2263])/(sum(TLTtotal[2014:2263]))



#total
(sum(nonzeroSPY) / sum(SPYtotal))


(sum(nonzeroTLT) / sum(TLTtotal))


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


dailynoiseTLT <- noisetosignal2(TLT_raw, TLT_5min)
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

tqtlt <- length(grep("T  Q", extlt[,1])) + length(grep("T  Q", extlt[,2])) + length(grep("T  Q", extlt[,3]))

wyjktlt <- length(grep("W  Y  J  K", extlt[,1])) + length(grep("W  Y  J  K", extlt[,2])) + length(grep("W  Y  J  K", extlt[,3]))

ptlt <- length(grep("P", extlt[,1])) + length(grep("P", extlt[,2])) + length(grep("P", extlt[,3]))

ztlt <- length(grep("Z", extlt[,1])) + length(grep("Z", extlt[,2])) + length(grep("Z", extlt[,3]))

dtlt <- length(grep("D", extlt[,1])) + length(grep("D", extlt[,2])) + length(grep("D", extlt[,3]))

btlt <- length(grep("B", extlt[,1])) + length(grep("B", extlt[,2])) + length(grep("B", extlt[,3]))


c(tqtlt, wyjktlt, ptlt, ztlt, dtlt, btlt)

#SPY

tqspy <- length(grep("T  Q", exspy[,1])) + length(grep("T  Q", exspy[,2])) + length(grep("T  Q", exspy[,3]))

wyjkspy <- length(grep("W  Y  J  K", exspy[,1])) + length(grep("W  Y  J  K", exspy[,2])) + length(grep("W  Y  J  K", exspy[,3]))

pspy <- length(grep("P", exspy[,1])) + length(grep("P", exspy[,2])) + length(grep("P", exspy[,3]))

zspy <- length(grep("Z", exspy[,1])) + length(grep("Z", exspy[,2])) + length(grep("Z", exspy[,3]))

dspy <- length(grep("D", exspy[,1])) + length(grep("D", exspy[,2])) + length(grep("D", exspy[,3]))

bspy <- length(grep("B", exspy[,1])) + length(grep("B", exspy[,2])) + length(grep("B", exspy[,3]))

c(tqspy, wyjkspy, pspy, zspy, dspy, bspy)