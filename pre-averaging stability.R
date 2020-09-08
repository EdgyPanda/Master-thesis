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




saveRDS(MRC_1sec, "MRC_1sec_100days2015.rds")
saveRDS(MRC_5sec, "MRC_5sec_100days2015.rds")
saveRDS(BPMRC_1sec, "BPMRC_1sec_100days2015.rds")
saveRDS(BPMRC_5sec, "BPMRC_5sec_100days2015.rds")
#saveRDS(MRC_1sec, "MRC_1sec.rds")
#saveRDS(MRC_5sec, "MRC_5sec.rds")
#saveRDS(BPMRC_1sec, "BPMRC_1sec.rds")
#saveRDS(BPMRC_5sec, "BPMRC_5sec.rds")


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
geom_line( aes(x = theta, y= c(MRC_1sec_TLT_mean), color="#a6cee3"), lwd = 1.2) +
geom_line( aes(x = theta, y= c(MRC_5sec_TLT_mean), color="#1f78b4"), lwd = 1.2) + 
geom_line( aes(x = theta, y= c(BPMRC_1sec_TLT_mean), color="#b2df8a"), lwd = 1.2) + 
geom_line( aes(x = theta, y= c(BPMRC_5sec_TLT_mean), color="#33a02c"), lwd = 1.2) +
theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"),
	legend.title = element_blank(),  axis.title=element_text(size=12)) +
labs(title = "TLT",
       x = "Theta",
       y = "IV estimates",
       colour = "Estimators") + ylim(0,5) +
scale_color_manual(values = c('#a6cee3' = '#a6cee3', '#1f78b4' = '#1f78b4', '#b2df8a'='#b2df8a', '#33a02c'='#33a02c'),
	labels = unname(TeX(c("$MRC_{t}^{1 sec}",  "$MRC_{t}^{5 sec}", "PBPCov_{t}^{1 sec}", "PBPCov_{t}^{5 sec}"))))


#for SPY
p2 <- ggplot() + 
geom_line( aes(x = theta, y= c(MRC_1sec_SPY_mean), color="#a6cee3"), lwd = 1.2) +
geom_line( aes(x = theta, y= c(MRC_5sec_SPY_mean), color="#1f78b4"), lwd = 1.2) + 
geom_line( aes(x = theta, y= c(BPMRC_1sec_SPY_mean), color="#b2df8a"), lwd = 1.2) + 
geom_line( aes(x = theta, y= c(BPMRC_5sec_SPY_mean), color="#33a02c"), lwd = 1.2) +
theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"),
	legend.title = element_blank(),  axis.title=element_text(size=12)) +
labs(title = "SPY",
       x = "Theta",
       y = "IV estimates",
       colour = "Estimators") +ylim(0,8) +
scale_color_manual(values = c('#a6cee3' = '#a6cee3', '#1f78b4' = '#1f78b4', '#b2df8a'='#b2df8a', '#33a02c'='#33a02c'),
	labels = unname(TeX(c("$MRC_{t}^{1 sec}",  "$MRC_{t}^{5 sec}", "PBPCov_{t}^{1 sec}", "PBPCov_{t}^{5 sec}"))))

#for Correlation
p3 <- ggplot() + 
geom_line( aes(x = theta, y= c(MRC_1sec_COR_mean), color="#a6cee3"), lwd = 1.2) +
geom_line( aes(x = theta, y= c(MRC_5sec_COR_mean), color='#1f78b4'), lwd = 1.2) + 
geom_line( aes(x = theta, y= c(BPMRC_1sec_COR_mean), color='#b2df8a'), lwd = 1.2) + 
geom_line( aes(x = theta, y= c(BPMRC_5sec_COR_mean), color='#33a02c'), lwd = 1.2) +
scale_color_manual(values = c('#a6cee3' = '#a6cee3', '#1f78b4' = '#1f78b4', '#b2df8a'='#b2df8a', '#33a02c'='#33a02c'),
	labels = unname(TeX(c("$MRC_{t}^{1 sec}", "$MRC_{t}^{5 sec}", "PBPCov_{t}^{1 sec}", "PBPCov_{t}^{5 sec}")))) +
theme_grey() +
theme(legend.position = c(0.70, 0.23), legend.background = element_rect(fill="lightblue", size=0.5,
 linetype="solid"), 
	plot.title = element_text(hjust = 0.5, face = "bold"),  axis.title=element_text(size=12)) +
labs(title = "TLT & SPY",
       x = "Theta",
       y = "Correlation estimates",
       colour = "Estimators") + ylim(-0.8, 0)

library(gridExtra)

grid.arrange(p1,p2,p3, ncol= 3)






#BELOW CONTAINS CODE FROM LAST YEAR

####################################The Choice of theta######################################

#cts5secmerged contains one NA value. 


theta <- seq(0.1, 2, 0.04)

estimates1sec <- matrix(0L, ncol=50, nrow=length(theta))

estimates1secBP <- matrix(0L, ncol=50, nrow=length(theta))

estimates5sec <- matrix(0L, ncol=50, nrow=length(theta))

estimates5secBP <- matrix(0L, ncol=50, nrow=length(theta))

estimatesTLT <- matrix(0L, ncol = 50, nrow = length(theta))

estimatesBPTLT <- matrix(0L, ncol=50, nrow=length(theta))

estimatesSPY <- matrix(0L, ncol = 50, nrow = length(theta))

estimatesBPSPY <- matrix(0L, ncol=50, nrow=length(theta))

fivesecestimateTLT <- matrix(0L, ncol=50, nrow=length(theta))
fivesecestimateBPTLT <- matrix(0L, ncol=50, nrow=length(theta))

fivesecestimateSPY <- matrix(0L, ncol=50, nrow=length(theta))
fivesecestimateBPSPY <- matrix(0L, ncol=50, nrow=length(theta))

estimate_rt_TLTsigma <- matrix(0L, ncol=50, nrow=length(theta))
estimate_rt_SPYsigma <- matrix(0L, ncol=50, nrow=length(theta))
estimate_rt_rho <- matrix(0L, ncol=50, nrow=length(theta))


estimate_rt_BPTLTsigma <- matrix(0L, ncol=50, nrow=length(theta))
estimate_rt_BPSPYsigma <- matrix(0L, ncol=50, nrow=length(theta))
estimate_rt_BPrho <- matrix(0L, ncol=50, nrow=length(theta))

#50 days will take approx 13 hours for correlation calculations. 



for (i in 1:50){

	for (j in 1:length(theta)){

		#TLT
		estimatesTLT[j,i] <- preavCov(cts1secTLT[[i]], FALSE, FALSE, FALSE, theta = theta[j])
		estimatesBPTLT[j,i] <- preavBPCOV(cts1secTLT[[i]], TRUE, FALSE, TRUE, theta = theta[j])

		fivesecestimateTLT[j,i] <- preavCov(cts5secTLT[[i]], TRUE, TRUE, FALSE, theta = theta[j])
		fivesecestimateBPTLT[j,i] <- preavBPCOV(cts5secTLT[[i]], TRUE, FALSE, TRUE, theta = theta[j])

		#SPY
		estimatesSPY[j,i] <- preavCov(cts1secSPY[[i]], FALSE, FALSE, FALSE, theta = theta[j])
		estimatesBPSPY[j,i] <- preavBPCOV(cts1secSPY[[i]], TRUE, FALSE, TRUE, theta = theta[j])

		fivesecestimateSPY[j,i] <- preavCov(cts5secSPY[[i]], TRUE, TRUE, FALSE, theta = theta[j])
		fivesecestimateBPSPY[j,i] <- preavBPCOV(cts5secSPY[[i]], TRUE, FALSE, TRUE, theta = theta[j])

		estimates1sec[j,i] <- preavCov(cts1secmerged[[i]], TRUE, TRUE, TRUE,  theta[j])[2,1]

		estimates1secBP[j,i] <- preavBPCOV(cts1secmerged[[i]], TRUE, TRUE, TRUE, theta[j])[2,1]

		estimates5sec[j,i] <- preavCov(cts5secmerged[[i]], TRUE, TRUE, TRUE,  theta[j])[2,1]

		estimates5secBP[j,i] <- preavBPCOV(cts5secmerged[[i]], TRUE, TRUE, TRUE, theta[j])[2,1]

		estimate_rt_TLTsigma[j,i] <- preavCov(rt_both[[i]], TRUE, TRUE, FALSE, theta[j])[1,1]
		estimate_rt_SPYsigma[j,i] <- preavCov(rt_both[[i]], TRUE, TRUE, FALSE, theta[j])[2,2]
		estimate_rt_rho[j,i] <- preavCov(rt_both[[i]], TRUE, TRUE, TRUE, theta[j])[2,1]

		estimate_rt_BPTLTsigma[j,i] <- preavBPCOV(rt_both[[i]], TRUE, FALSE, TRUE, theta[j])[1,1]
		estimate_rt_BPSPYsigma[j,i] <- preavBPCOV(rt_both[[i]], TRUE, FALSE, TRUE, theta[j])[2,2]
		estimate_rt_BPrho[j,i] <- preavBPCOV(rt_both[[i]], TRUE, TRUE, TRUE, theta[j])[2,1]



		print(sprintf("Theta iteration %s, for the current day %s of %s", j,i, length(cts1secTLT)))
		
	}
}


meanestimate_rt_TLTsigma <- rowMeans(estimate_rt_TLTsigma)
meanestimate_rt_SPYsigma <- rowMeans(estimate_rt_SPYsigma)
meanestimate_rt_rho <- rowMeans(estimate_rt_rho)

#BPCov
meanestimate_rt_BPTLTsigma <- rowMeans(estimate_rt_BPTLTsigma)
meanestimate_rt_BPSPYsigma <- rowMeans(estimate_rt_BPSPYsigma)
meanestimate_rt_BPrho <- rowMeans(estimate_rt_BPrho)

#correlation 
estimate1secmean <- rowMeans(estimates1sec)
estimate1secBPmean <- rowMeans(estimates1secBP)

estimates5secmean <- rowMeans(estimates5sec)
estimates5secBPmean <- rowMeans(estimates5secBP)

plot(estimate1secmean)


#sqrt() everything to get the volatility...

#5 sec department:
fivesecmeanTLT <- rowMeans(fivesecestimateTLT)
fivesecmeanBPTLT <- rowMeans(fivesecestimateBPTLT)

fivesecmeanSPY <- rowMeans(fivesecestimateSPY)
fivesecmeanBPSPY <- rowMeans(fivesecestimateBPSPY)



#1 sec department:
estmeanTLT <- rowMeans(estimatesTLT)

estmeanBPTLT <- rowMeans(estimatesBPTLT)

estmeanSPY <- rowMeans(estimatesSPY)

estmeanBPSPY <- rowMeans(estimatesBPSPY)

plot(estmeanBPTLT)

library(ggplot2)
library(ggthemes)
library(RColorBrewer)

nb.cols <- 18
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

sumdataTLT <- data.frame(cbind(theta, MRC_1sec_TLT_mean, BPMRC_1sec_TLT_mean, MRC_5sec_TLT_mean,
 BPMRC_5sec_TLT_mean)) #1e5*meanestimate_rt_TLTsigma, 1e5*meanestimate_rt_BPTLTsigma))
colnames(sumdataTLT) <- c("Theta", "RCOVpa_1sec_TLT", "BCOVpa_1sec_TLT", "RCOVpa_5sec_TLT", "BCOVpa_5sec_TLT", "RCOV_rt_TLT", "BCOV_rt_TLT")


sumdataSPY <- data.frame(cbind(theta, MRC_1sec_SPY_mean, BPMRC_1sec_SPY_mean , MRC_5sec_SPY_mean, BPMRC_5sec_SPY_mean)) # 1e5*meanestimate_rt_SPYsigma, 1e5*meanestimate_rt_BPSPYsigma))
colnames(sumdataSPY) <- c("Theta", "RVpa_1sec_SPY", "BVpa_1sec_SPY", "RVpa_5sec_SPY", "BVpa_5sec_SPY", "RCOV_rt_SPY", "BCOV_rt_SPY")

sumdatarho <- data.frame(cbind(theta, estimate1secmean, estimate1secBPmean, estimates5secmean, estimates5secBPmean,meanestimate_rt_rho , meanestimate_rt_BPrho))
colnames(sumdatarho) <- c("Theta", "RCOVpa_1sec", "BCOVpa_1sec", "RCOVpa_5sec", "BCOVpa_5sec", "RCOV_rt", "BCOV_rt")


p1 <- ggplot() + 
geom_line(data = sumdataTLT, aes(x = sumdataTLT[,1], y=c(sumdataTLT[,2]), color='RCovpa_1sec_TLT'), lwd = 1.2) + 
geom_line(data = sumdataTLT, aes(x = sumdataTLT[,1], y=c(sumdataTLT[,3]), color='BPCovpa_1sec_TLT'), lwd = 1.2) + 
geom_line(data = sumdataTLT, aes(x = sumdataTLT[,1], y=c(sumdataTLT[,4]), color='RCovpa_5sec_TLT'), lwd = 1.2) + 
geom_line(data = sumdataTLT, aes(x = sumdataTLT[,1], y=c(sumdataTLT[,5]), color='BCovpa_5sec_TLT'), lwd = 1.2) +
scale_color_brewer(palette = "Set2") + 
theme( plot.title = element_text(hjust = 0.5, face = "bold"),legend.title = element_blank(),  axis.title=element_text(size=12)) + ylim(0, 4) +
labs(title = "TLT",
       x = "Theta",
       y = "IV estimates",
       colour = "Estimators") 

p2 <- ggplot() + 
geom_line(data = sumdataSPY, aes(x = sumdataSPY[,1], y=sumdataSPY[,2], color='RCovpa_1sec_SPY'), lwd = 1.2) + 
geom_line(data = sumdataSPY, aes(x = sumdataSPY[,1], y=sumdataSPY[,3], color='BPCovpa_1sec_SPY'), lwd = 1.2) + 
geom_line(data = sumdataSPY, aes(x = sumdataSPY[,1], y=sumdataSPY[,4], color='RCovpa_5sec_SPY'), lwd = 1.2) + 
geom_line(data = sumdataSPY, aes(x = sumdataSPY[,1], y=sumdataSPY[,5], color='BPCovpa_5sec_SPY'), lwd = 1.2) +
theme(legend.title = element_blank(),  plot.title = element_text(hjust = 0.5, face = "bold"),
  axis.title=element_text(size=12)) + ylim(0, 13) + 
scale_color_brewer(palette = "Set2") + 
labs(title = "SPY",
       x = "Theta",
       y = "IV estimates",
       colour = "Estimators") 



p3 <- ggplot(data = sumdatarho) + 
geom_line(aes(x = sumdatarho[,1], y=sumdatarho[,2], color = 'RCovpa_1sec'), lwd = 1.2) + 
geom_line(aes(x = sumdatarho[,1], y=sumdatarho[,3], color='BPCovpa_1sec'), lwd = 1.2) + 
geom_line(aes(x = sumdatarho[,1], y=sumdatarho[,4], color='RCovpa_5sec'), lwd = 1.2) + 
geom_line(aes(x = sumdatarho[,1], y=sumdatarho[,5], color='BPCovpa_5sec'), lwd = 1.2) +
scale_color_brewer(palette = "Set2") + 
theme_grey() +
theme(legend.position = c(0.70, 0.23), legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid"), plot.title = element_text(hjust = 0.5, face = "bold"),  axis.title=element_text(size=12)) +
ylim(-0.8, 0) + 
labs(title = "TLT & SPY",
       x = "Theta",
       y = "IV estimates",
       colour = "Estimators") 



library(gridExtra)
grid.arrange(p1, p2, p3, nrow = 1)