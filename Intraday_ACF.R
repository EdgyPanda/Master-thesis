#--------------------------------------------------Intra day ACF using old 2018 cleaned SPY returns from ATP topic course:


dataTLT <- readRDS("dataTLT.rds")
dataSPY <- readRDS("dataSPY.rds")
library(highfrequency)




#Constructing a list for all second frequencies, which contains list-elements of all days of each asset for specified
#second frequency. Then I'll average out acf for all days and average out across frequencies. 

secs <- seq(5,60,5)

acrossfreqspy <- list()
acrossfreqtlt <- list()

aggregatespy <- list()
aggregatetlt <- list()

for(j in 1:length(secs)){
	for(i in 1:length(dataSPY)){

		aggregatespy[[i]] <- aggregatets(dataSPY[[i]], on="seconds",k=secs[j])
		aggregatetlt[[i]] <- aggregatets(dataTLT[[i]], on="seconds",k=secs[j])

	}
	acrossfreqspy[[j]] <- aggregatespy
	acrossfreqtlt[[j]] <- aggregatetlt
}



for(j in 1:length(secs)){
	for(i in 1:length(dataSPY)){

		acrossfreqspy[[j]][[i]] <- diff(log(acrossfreqspy[[j]][[i]]))[-1]
		acrossfreqtlt[[j]][[i]] <- diff(log(acrossfreqtlt[[j]][[i]]))[-1]

	}
	print(sprintf("%s", j))
}



#SPY GETS NEGATIVE 1ST AUTOCORRELATION AT K=15 SECONDS (-2.62872e-03), K=20 seconds (-0.005257), K=25 seconds (-2.2531e-03)
# K = 45 seconds (-6.175382e-03), k = 55 SECONDS (SMALLEST: -7.004245e-03)


library(ggplot2)
library(forecast)



acfnumbersspy <- array(0L, dim = c(31,2516,length(secs)))
acfnumberstlt <- array(0L, dim = c(31,2516,length(secs)))

for(j in 1:length(secs)){
	for(i in 1:length(aggregatespy)){

	acfnumbersspy[, i, j] <- acf(acrossfreqspy[[j]][[i]], lag.max  = 30, plot =  F)$acf
	acfnumberstlt[, i, j] <- acf(acrossfreqtlt[[j]][[i]], lag.max  = 30, plot =  F)$acf


	}
}

#HERE GIVES YOU THE NEGATIVE AUTOCORRELATION FOR THE SECOND LAG. 
#RowMeans for arrays takes mean of all rows in each matrix for all j. So, it already averages across days and frequencies!
avgacfspy <- rowMeans(acfnumbersspy)
avgacftlt <- rowMeans(acfnumberstlt)


data2 <- data.frame(avgacfspy, avgacftlt)[-1,]
lag <- 1:30
ciline <- qnorm((1 - 0.95)/2)/sqrt(length(dataSPY))


q <- ggplot(data = data2, mapping = aes(x = lag, y = data2[,1])) + geom_hline(aes(yintercept = 0)) + 
geom_hline(aes(yintercept = ciline),col = 'black', linetype = "dashed", size = 1) + 
geom_hline(aes(yintercept = -ciline), col = 'black', linetype = "dashed", size = 1) + 
geom_segment(aes(xend =lag, yend = 0, linetype = 'Log-returns SPY'), col = 'red', size=1) + xlab('Lag') + ylab('ACF') + 
geom_segment(aes(xend=lag, yend=0, y = data2[,2], linetype = 'Log-returns TLT'), col = 'red', size = 1, alpha = 0.6) +
theme(legend.title = element_blank(), legend.position = c(0.8, 0.8), legend.background = element_rect(fill="lightblue", 
size=0.5, linetype="solid", colour = 'darkblue'), 
legend.text = element_text(colour="black", size=9, face="bold"), legend.text.align = 1) + 
scale_fill_manual(name ='', values = c("blue", "red")) 


acf <- q + geom_point(aes(lag, data2[,1]), col='black') + geom_point(aes(lag, data2[,2]), col='steelblue')
acf

grid.arrange(retplot, acf, ncol=2)