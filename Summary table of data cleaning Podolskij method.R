#will contain the summary statistics following Hautsch & Podolskij (2013). 

library(xts)
library(highfrequency)
library(matlib)
library(tidyverse)

dataTLT <- readRDS("dataTLT.rds")
dataSPY <- readRDS("dataSPY.rds")



#Calculating average time in seconds between each trade:
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
