#--------------------------pre-averaging stability theta---------------------------------
#
#
#
#
#
source("Previous functions/projectfunctions.R")


library(xts)
library(highfrequency)
library(matlib)

dataTLT <- readRDS("dataTLT.rds")
dataSPY <- readRDS("dataSPY.rds")

theta <- seq(0.1,2,0.1)

TLT_1sec <- list()
TLT_5sec <- list()

SPY_1sec <- list()
SPY_5sec <- list()

for(i in 1:length(dataTLT)){

	TLT_1sec[[i]] <- diff(log(aggregatets(dataTLT[[i]], on="seconds",k=1)))[-1]
	TLT_5sec[[i]] <- diff(log(aggregatets(dataTLT[[i]], on="seconds",k=5)))[-1]

	SPY_1sec[[i]] <- diff(log(aggregatets(dataSPY[[i]], on="seconds", k=1)))[-1]
    SPY_5sec[[i]] <- diff(log(aggregatets(dataSPY[[i]], on="seconds", k=5)))[-1]
}

Merged_1sec <- list()
Merged_5sec <- list()
for (i  in 1:length(dataTLT)){

  Merged_1sec[[i]] <- na.omit(cbind(TLT_1sec[[i]], SPY_1sec[[i]]))
  Merged_5sec[[i]] <- na.omit(cbind(TLT_5sec[[i]], SPY_5sec[[i]]))
  print(sprintf("%s out of %s",i,length(dataTLT)))
}


#free up ram:
rm(dataTLT, dataSPY)

#Testing for all days requires too much time, therefore testing for 10 different days, for each year.

MRC_1sec_TLT <- matrix(0L, ncol = 100, nrow = length(theta))
MRC_5sec_TLT <- matrix(0L, ncol = 100, nrow = length(theta))
BPMRC_1sec_TLT <- matrix(0L, ncol = 100, nrow = length(theta))
BPMRC_5sec_TLT <- matrix(0L, ncol = 100, nrow = length(theta))

MRC_1sec_SPY <- matrix(0L, ncol = 100, nrow = length(theta))
MRC_5sec_SPY <- matrix(0L, ncol = 100, nrow = length(theta))
BPMRC_1sec_SPY <- matrix(0L, ncol = 100, nrow = length(theta))
BPMRC_5sec_SPY <- matrix(0L, ncol = 100, nrow = length(theta))

MRC_1sec_COR <- matrix(0L, ncol = 100, nrow = length(theta))
MRC_5sec_COR <- matrix(0L, ncol = 100, nrow = length(theta))
BPMRC_1sec_COR <- matrix(0L, ncol = 100, nrow = length(theta))
BPMRC_5sec_COR <- matrix(0L, ncol = 100, nrow = length(theta))


#Constructing a sequence of 10 day intervals in each year (each year varies with trading days). 
set.seed(1)
iT <- c(sample(1:250, 1), sample(251:500,1), sample(505:754,1), sample(755:1006,1),
	sample(1007:1258,1), sample(1259:1510,1), sample(1511:1762,1), sample(1763:2013,1), 
	sample(2014:2264,1), sample(2265:2516,1))

iT <- c(seq(iT[1],iT[1]+9),seq(iT[2],iT[2]+9),seq(iT[3],iT[3]+9),seq(iT[4],iT[4]+9),
	seq(iT[5],iT[5]+9),seq(iT[6],iT[6]+9),seq(iT[7],iT[7]+9),seq(iT[8],iT[8]+9),
	seq(iT[9],iT[9]+9),seq(iT[10],iT[10]+9))


#We are testing for stability in theta. Therefore days are independent and does not matter. 

#trying to optimize a bit:

temp1 <- array(0L, dim =c(2,2,length(theta)))
temp2 <- array(0L, dim =c(2,2,length(theta)))
temp3 <- array(0L, dim =c(2,2,length(theta)))
temp4 <- array(0L, dim =c(2,2,length(theta)))
MRC_1sec <- list() 
MRC_5sec <- list()
BPMRC_1sec <- list()
BPMRC_5sec <- list()

library(tictoc)
#Do not run, takes approx 7 hours and 30 mins. 
tic()
for(i in 1:100){
	for(j in 1:length(theta)){
		temp1[,,j] <- preavCov(Merged_1sec[[iT[i]]], T, T, F, theta = theta[j])
		temp2[,,j] <- preavCov(Merged_5sec[[iT[i]]], T, T, F, theta = theta[j])
	    temp3[,,j] <- preavBPCOV(Merged_1sec[[iT[i]]], TRUE, FALSE, TRUE, theta = theta[j])
	    temp4[,,j] <- preavBPCOV(Merged_5sec[[iT[i]]], TRUE, FALSE, TRUE, theta = theta[j])

		print(sprintf("Theta iteration %s, for the current day %s and %s", j,iT[i], i))

	}
	MRC_1sec[[i]] <- temp1
	MRC_5sec[[i]] <- temp2
	BPMRC_1sec[[i]] <- temp3
	BPMRC_5sec[[i]] <- temp4
}
toc()


saveRDS(MRC_1sec, "MRC_1sec.rds")
saveRDS(MRC_5sec, "MRC_5sec.rds")
saveRDS(BPMRC_1sec, "BPMRC_1sec.rds")
saveRDS(BPMRC_5sec, "BPMRC_5sec.rds")


MRC_1sec_TLT_IV <- matrix(0L, ncol = 100, nrow = length(theta))
MRC_5sec_TLT_IV  <- matrix(0L, ncol = 100, nrow = length(theta))
BPMRC_1sec_TLT_IV  <- matrix(0L, ncol = 100, nrow = length(theta))
BPMRC_5sec_TLT_IV  <- matrix(0L, ncol = 100, nrow = length(theta))


MRC_1sec_SPY_IV <- matrix(0L, ncol = 100, nrow = length(theta))
MRC_5sec_SPY_IV  <- matrix(0L, ncol = 100, nrow = length(theta))
BPMRC_1sec_SPY_IV  <- matrix(0L, ncol = 100, nrow = length(theta))
BPMRC_5sec_SPY_IV  <- matrix(0L, ncol = 100, nrow = length(theta))

MRC_1sec_COR <- matrix(0L, ncol = 100, nrow = length(theta))
MRC_5sec_COR  <- matrix(0L, ncol = 100, nrow = length(theta))
BPMRC_1sec_COR <- matrix(0L, ncol = 100, nrow = length(theta))
BPMRC_5sec_COR  <- matrix(0L, ncol = 100, nrow = length(theta))


for(i in 1:100){
	for(j in 1:length(theta)){

		#TLT
		MRC_1sec_TLT_IV[j,i] <- MRC_1sec[[i]][1,1,j]
		MRC_5sec_TLT_IV[j,i] <- MRC_5sec[[i]][1,1,j]
		BPMRC_1sec_TLT_IV[j,i] <- BPMRC_1sec[[i]][1,1,j]
		BPMRC_5sec_TLT_IV[j,i] <- BPMRC_5sec[[i]][1,1,j]

		#SPY
		MRC_1sec_SPY_IV[j,i] <- MRC_1sec[[i]][2,2,j]
		MRC_5sec_SPY_IV[j,i] <- MRC_5sec[[i]][2,2,j]
		BPMRC_1sec_SPY_IV[j,i] <- BPMRC_1sec[[i]][2,2,j]
		BPMRC_5sec_SPY_IV[j,i] <- BPMRC_5sec[[i]][2,2,j]

		#COV
		MRC_1sec_COR[j,i] <- MRC_1sec[[i]][2,1,j]
		MRC_5sec_COR[j,i] <- MRC_5sec[[i]][2,1,j]
		BPMRC_1sec_COR[j,i] <- BPMRC_1sec[[i]][2,1,j]
		BPMRC_5sec_COR[j,i] <- BPMRC_5sec[[i]][2,1,j]

	}

}

#Get correlation
for(i in 1:100){
	for(j in 1:length(theta)){
		MRC_1sec_COR[j,i] <- MRC_1sec_COR[j,i]/(sqrt(MRC_1sec_TLT_IV[j,i])*sqrt(MRC_1sec_SPY_IV[j,i]))
		MRC_5sec_COR[j,i] <- MRC_5sec_COR[j,i]/(sqrt(MRC_5sec_TLT_IV[j,i])*sqrt(MRC_5sec_SPY_IV[j,i]))
		BPMRC_1sec_COR[j,i] <- BPMRC_1sec_COR[j,i]/(sqrt(BPMRC_1sec_TLT_IV[j,i])*sqrt(BPMRC_1sec_SPY_IV[j,i]))
		BPMRC_5sec_COR[j,i] <- BPMRC_5sec_COR[j,i]/(sqrt(BPMRC_5sec_TLT_IV[j,i])*sqrt(BPMRC_5sec_SPY_IV[j,i]))
	}
}

library(stats)
library(matlib)


MRC_1sec_TLT_mean <- rowMeans(MRC_1sec_TLT_IV)*1e5
MRC_5sec_TLT_mean <- rowMeans(MRC_5sec_TLT_IV)*1e5
BPMRC_1sec_TLT_mean <- rowMeans(BPMRC_1sec_TLT_IV)*1e5
BPMRC_5sec_TLT_mean <- rowMeans(BPMRC_5sec_TLT_IV)*1e5

MRC_1sec_SPY_mean <- rowMeans(MRC_1sec_SPY_IV)*1e5
MRC_5sec_SPY_mean <- rowMeans(MRC_5sec_SPY_IV)*1e5
BPMRC_1sec_SPY_mean <- rowMeans(BPMRC_1sec_SPY_IV)*1e5
BPMRC_5sec_SPY_mean <- rowMeans(BPMRC_5sec_SPY_IV)*1e5

MRC_1sec_COR_mean <- rowMeans(MRC_1sec_COR)
MRC_5sec_COR_mean <- rowMeans(MRC_5sec_COR)
BPMRC_1sec_COR_mean <- rowMeans(BPMRC_1sec_COR)
BPMRC_5sec_COR_mean <- rowMeans(BPMRC_5sec_COR)


library(ggplot2)
library(latex2exp)
library(ggthemes)
library(RColorBrewer)


#HEX CODE FOR CORRELATIONS: #fb9a99, #e31a1c. 

#note to self, scale colour manual should be on the last graph together with the theme. 

#for TLT
p1 <- ggplot() + 
geom_line( aes(x = theta, y= c(MRC_1sec_TLT_mean), color="MRC_1sec_TLT"), lwd = 1.2) +
geom_line( aes(x = theta, y= c(MRC_5sec_TLT_mean), color="MRC_5sec_TLT"), lwd = 1.2) + 
geom_line( aes(x = theta, y= c(BPMRC_1sec_TLT_mean), color="BPMRC_1sec_TLT"), lwd = 1.2) + 
geom_line( aes(x = theta, y= c(BPMRC_5sec_TLT_mean), color="BPMRC_5sec_TLT"), lwd = 1.2) +
theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"),
	legend.title = element_blank(),  axis.title=element_text(size=12)) +
labs(title = "TLT",
       x = "Theta",
       y = "QV estimates",
       colour = "Estimators") + ylim(0,5) 


#for SPY
p2 <- ggplot() + 
geom_line( aes(x = theta, y= c(MRC_1sec_SPY_mean), color="MRC_1sec_SPY"), lwd = 1.2) +
geom_line( aes(x = theta, y= c(MRC_5sec_SPY_mean), color="MRC_5sec_SPY"), lwd = 1.2) + 
geom_line( aes(x = theta, y= c(BPMRC_1sec_SPY_mean), color="BPMRC_1sec_SPY"), lwd = 1.2) + 
geom_line( aes(x = theta, y= c(BPMRC_5sec_SPY_mean), color="BPMRC_5sec_SPY"), lwd = 1.2) +
theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"),
	legend.title = element_blank(),  axis.title=element_text(size=12)) +
labs(title = "SPY",
       x = "Theta",
       y = "QV estimates",
       colour = "Estimators") + ylim(0,5) 

#for Correlation
p3 <- ggplot() + 
geom_line( aes(x = theta, y= c(MRC_1sec_COR_mean), color="MRC_1sec"), lwd = 1.2) +
geom_line( aes(x = theta, y= c(MRC_5sec_COR_mean), color='MRC_5sec'), lwd = 1.2) + 
geom_line( aes(x = theta, y= c(BPMRC_1sec_COR_mean), color='BPMRC_1sec'), lwd = 1.2) + 
geom_line( aes(x = theta, y= c(BPMRC_5sec_COR_mean), color='BPMRC_5sec'), lwd = 1.2) +
theme_grey() +
theme(legend.position = c(0.70, 0.23), legend.background = element_rect(fill="lightblue", size=0.5,
 linetype="solid"), 
	plot.title = element_text(hjust = 0.5, face = "bold"),  axis.title=element_text(size=12)) +
labs(title = "TLT & SPY",
       x = "Theta",
       y = "Correlation estimates",
       colour = "Estimators") + ylim(-0.8, 0)

library(gridExtra)

 p4 <- grid.arrange(p1,p2,p3, ncol= 3)


ggsave(file="thetaplot.eps", plot = p4)

