#Automated datacleaning using constructed functions from last years project. 
#
#
# Save cleaned data in .rds format. 
#
#
#
library(dataCleanFunc)
library(highfrequency)
setwd("C:/Users/emil7/Dropbox/Uni - Mathematics and Economics/Speciale - Master Thesis/Data cleaned/TLT/Merged files")

#TLT until 20120515 (due to corrupt file)

data1 <- read.csv("MasterData - until 20120515 (not included).csv", header = T)

cleanData1 <- nanRemover(data1)

dates1 <- getDates(data1)

timestamps1 <- timestamp(cleanData1, dates1)

timeseriesTLT1 <- timeseriesList(data1, cleanData1, timestamps1)

#free up some RAM.
rm(data1, dates1, timestamps1)

#TLT 20120515 - 20150910 (due to corrupt file) - have to limit data2 manually to remove NA columns produced by error.

data2 <- read.csv("MasterData - from 20120515 till 20150910 (not included).csv", header = T)[,1:1670]

cleanData2 <- nanRemover(data2)

dates2 <- getDates(data2)

timestamps2 <- timestamp(cleanData2, dates2)

timeseriesTLT2 <- timeseriesList(data2, cleanData2, timestamps2)

rm(data2, dates2, timestamps2)

#TLT 20150910 - 20191231

data3 <- read.csv("MasterData.csv", header = T)

cleanData3 <- nanRemover(data3)

dates3 <- getDates(data3)

timestamps3 <- timestamp(cleanData3, dates3)

timeseriesTLT3 <- timeseriesList(data3, cleanData3, timestamps3)

rm(data3, dates3, timestamps3)



#--------------collecting all of the list into one--------------------------------


t1 <- c(timeseriesTLT1, timeseriesTLT2)


dataTLT <- c(t1, timeseriesTLT3)

saveRDS(dataTLT, "dataTLT.rds")







#-------------------------------SPY-----------------------------------------------

setwd("C:/Users/emil7/Dropbox/Uni - Mathematics and Economics/Speciale - Master Thesis/Data cleaned/SPY/Merged files")

#SPY until 20150512 (due to corrupt file)

data1 <- read.csv("MasterData until 20150512 (not included).csv", header = T)

cleanData1 <- nanRemover(data1)

dates1 <- getDates(data1)

timestamps1 <- timestamp(cleanData1, dates1)

timeseriesSPY1 <- timeseriesList(data1, cleanData1, timestamps1)

#free up some RAM.
rm(data1, dates1, timestamps1)



data2 <- read.csv("MasterData from 20150512 till 20161130 (not included).csv", header = T)

cleanData2 <- nanRemover(data2)

dates2 <- getDates(data2)

timestamps2 <- timestamp(cleanData2, dates2)

timeseriesSPY2 <- timeseriesList(data2, cleanData2, timestamps2)

rm(data2, dates2, timestamps2)




data3 <- read.csv("MasterData from 20161130 until 20190517 (not included).csv", header = T)

cleanData3 <- nanRemover(data3)

dates3 <- getDates(data3)

timestamps3 <- timestamp(cleanData3, dates3)

timeseriesSPY3 <- timeseriesList(data3, cleanData3, timestamps3)

rm(data3, dates3, timestamps3)



data4 <- read.csv("MasterData.csv", header = T)

cleanData4 <- nanRemover(data4)

dates4 <- getDates(data4)

timestamps4 <- timestamp(cleanData4, dates4)

timeseriesSPY4 <- timeseriesList(data4, cleanData4, timestamps4)

rm(data4, dates4, timestamps4)



S1 <- c(timeseriesSPY1, timeseriesSPY2, timeseriesSPY3, timeseriesSPY4)


saveRDS(S1, "dataSPY.rds")