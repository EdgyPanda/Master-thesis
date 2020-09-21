#semicovariance plot. 


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



sparseTLT5min <- list()
sparseSPY5min <- list()


for(i in 1:length(dataTLT)){

	sparseTLT5min[[i]] <- aggregatets(dataTLT[[i]], on = "minutes", k = 5)
	sparseSPY5min[[i]] <- aggregatets(dataSPY[[i]], on = "minutes", k = 5)

}

for(i in 1:length(dataTLT)){

	sparseTLT5min[[i]] <- diff(log(sparseTLT5min[[i]]))[-1]
	sparseSPY5min[[i]] <- diff(log(sparseSPY5min[[i]]))[-1]

}


merged5min <- list()

for(i in 1:length(dataTLT)){

	merged5min[[i]] <- na.omit(cbind(sparseTLT5min[[i]], sparseSPY5min[[i]]))

}


library(ggplot2)

RSCovpos_5min <- lapply(merged5min, function(x) realsemicov(x,"P"))
RSCovneg_5min <- lapply(merged5min, function(x) realsemicov(x,"N"))
RSCovmixed_5min <- lapply(merged5min, function(x) realsemicov(x,"M"))

RSCovpos_5min <- matrix(unlist(lapply(RSCovpos_5min, function(x) cbind(x[1,1], x[2,2], x[2,1]))), 
	nrow = 2516, ncol=3, byrow=T)

colnames(RSCovpos_5min) <- c("TLT", "SPY", "Correlation")

RSCovneg_5min <- matrix(unlist(lapply(RSCovneg_5min, function(x) cbind(x[1,1], x[2,2], x[2,1]))), 
	nrow = 2516, ncol=3, byrow=T)

colnames(RSCovneg_5min) <- c("TLT", "SPY", "Correlation")

RSCovmixed_5min <- matrix(unlist(lapply(RSCovmixed_5min, function(x) cbind(x[1,1], x[2,2], x[2,1]))), 
	nrow = 2516, ncol=3, byrow=T)

colnames(RSCovmixed_5min) <- c("TLT", "SPY", "Covariance")

p1 <- ggplot() + geom_line(aes(as.Date(getDates), RSCovmixed_5min[,3], col = "Mixed")) + 
geom_line(aes(as.Date(getDates), RSCovpos_5min[,3]+RSCovneg_5min[,3], col = "Concordant")) + 
ylab("Semicovariances") + xlab("Dates") +
theme_grey() + ylim(-5e-04,9.53611e-05) + 
theme( legend.position = c(0.80, 0.23), legend.background = element_rect(fill="lightblue", size=0.8,
 linetype="solid"), 
	plot.title = element_text(hjust = 0.5, face = "bold"),  axis.title=element_text(size=12)) 

ggsave(p1, file = "Semicovariances.eps", device = "eps")


max(abs(RSCovmixed_5min[,3])) > max(RSCovpos_5min[,3]+RSCovneg_5min[,3])