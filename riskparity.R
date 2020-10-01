#risk parity portfolio preliminaries 


source("functions.R")
library(riskParityPortfolio)
library(matlib)
library(xts)
library(highfrequency)
library(matlib)
library(MCS)

dtb3 <- read.csv("DTB3.csv", header =  F)

dtb3[,2] <- as.numeric(levels(dtb3[,2]))[dtb3[,2]]

dtb3 <- dtb3[!is.na(dtb3[,2]), ]

ts.plot(dtb3[,2])

logriskfreerate <- log(1 + dtb3[,2]/100) #daily



riskparity_2dim <- function(matrix, risktarget, rt = F){
  
  bonds <- matrix[1,1]
  stocks <- matrix[2,2]
  
  w_1 <- sqrt(bonds)^-1 / (sqrt(stocks)^-1 + sqrt(bonds)^-1)
  w_2 <- sqrt(stocks)^-1 / (sqrt(stocks)^-1 + sqrt(bonds)^-1)

  w <- matrix(c(w_1, w_2), ncol=1, nrow=2) #dimnames = list(c(), c("Bond", "Stock"))
  
  portrisk <- as.numeric(sqrt(t(w) %*% (matrix) %*% w)) 
  
  riskcont <-  (w * (matrix %*% w)) / portrisk

  relativeriskcont <- (w * (matrix %*% w)) / portrisk^2

  if(rt){

  	alpha <- risktarget / portrisk

  	w_new <- w %*% alpha

  	w_new <- matrix(c(w_new[1], w_new[2]), ncol=1, nrow=2)

  	portrisk <-  as.numeric(sqrt(t(w_new) %*% matrix %*% w_new))  #gives marginal risk for each asset. 

  	riskcont <- (w_new * (matrix %*% w_new)) / portrisk

  	w_riskfree <- uniroot(function(x) colSums(w_new)+x-1, interval = c(-100,100))$root

  	w_new <- matrix(c(w_new, w_riskfree), ncol=1, nrow=3)
  	rownames(w_new) <- c("TLT", "SPY", "riskfree")

  	lout <- list(w_new, portrisk, riskcont)
  	
  	names(lout) <- c("w", "portrisk", "riskcont")

 	return(lout)
  }

  lout <- list(w, portrisk, riskcont, relativeriskcont)
  names(lout) <- c("w", "portrisk", "riskcont", "relativeriskcont")

  return(lout)
  
}



calccov <- readRDS("calculatedcovariances.rds")

mergedfrequencies <- readRDS("mergedfrequencies.rds")


#-------------------------------portfolio volatility's sensitivity to correlation--------------------------

#getting average vol for tlt and spy.

meanvolTLT <- mean(sqrt(calccov[[1]][[7]][1,1,]*252*(24/6.5)))

meanvolSPY <-  mean(sqrt(calccov[[1]][[7]][2,2,]*252*(24/6.5)))

correlations <- seq(-0.9,0.9,0.01)

sensitivity <- numeric()

sensitivity2 <- numeric()

for(i in 1:length(correlations)){

	newcovariance <- matrix(c(meanvolTLT^2, meanvolTLT*meanvolSPY*correlations[i],
		meanvolTLT*meanvolSPY*correlations[i], meanvolSPY^2), ncol = 2, nrow = 2)

	w2 <- riskparity_2dim(newcovariance,0,F)$w

	portrisk2 <- riskparity_2dim(newcovariance,0,F)$portrisk

	sensitivity2[i] <-  (w2[1]*w2[2]*meanvolTLT*meanvolSPY)  /  (portrisk2)


}


library(ggplot2)


#Shows the volatility's sensitivity to correlation for the unlevered risk-parity portfolio. In essence, 
#upscaling and downscaling the portfolios using a leverage parameter only shifts the graph.

ggplot() + geom_line(aes(correlations, sensitivity2*100), col = "red", lwd = 1) + 
  scale_x_continuous(breaks = round(seq(-0.9,0.9, by = 0.1),1)) + ylab("portfolio volatility (%)")

#-------------------------------weight distribution dependent on excess returns---------------------------------------

 #done for returns until you have control over a risk-free asset.

 #fix bond return at 3%

source("Previous functions/projectfunctions.R")

stockvol <- seq(0,0.15,0.001)*1e-5
bondvol <- rep(0.03, length(stockvol))*1e-5

covs <- array(0L, c(2,2,length(stockvol)))

for(i in 1:length(stockvol)){

	covs[,,i] <- matrix(c(bondvol[i]^2, 0, 0, stockvol[i]^2), ncol=2, nrow=2)


}

weightsfordistribution <- matrix(0L, ncol=2, nrow=length(stockvol))

for(i in 1:length(stockvol)){

	weightsfordistribution[i, ] <- riskparity_2dim(covs[,,i])$w


}

library(PerformanceAnalytics)

rownames(weightsfordistribution) <- stockvol * 1e5



ggplot() + geom_line(aes(stockvol*1e5, weightsfordistribution[,2]))


#------------------------------------------trying effcient frontier and "risk-parity line"--------------------------

#starts at minvar portfolio and then goes to 100% stocks. 
#
#
# Be aware that TLT did a better job than SPY therefore you can only construct this where you go
# 100% into TLT instead of SPY. 


#data get and preparation: 
library(alphavantager)
library(xts)

av_api_key('0WXHEIY0O87LX4A3')

TLT <- as.data.frame(av_get(symbol = "TLT", av_fun = "TIME_SERIES_DAILY", outputsize = "full"))
rownames(TLT) <- TLT$timestamp
SPY <- as.data.frame(av_get(symbol = "SPY", av_fun = "TIME_SERIES_DAILY", outputsize = "full"))
rownames(SPY) <- SPY$timestamp

returns_TLT <- as.xts(diff(log(TLT[,4])), order.by = as.Date(TLT[,1], format='%d/%m/%Y')[-1])
returns_SPY <- as.xts(diff(log(SPY[,4])), order.by = as.Date(SPY[,1])[-1])

returns_TLT <-  returns_TLT[seq(from= as.Date('2010-01-02'), to = as.Date('2019-12-31'), by=1), ]
returns_SPY <- returns_SPY[seq(from= as.Date('2010-01-02'), to = as.Date('2019-12-31'), by=1), ]

merged_ret <- cbind(returns_TLT, returns_SPY)


ggplot() + geom_line(aes(index(returns_TLT), 1+cumsum(returns_TLT), col="TLT")) + 
geom_line(aes(index(returns_TLT), 1+cumsum(returns_SPY), col="SPY")) 

covlol <- realCov(merged_ret)


colnames(covlol) <- c("TLT", "SPY")

minvarweights <- minvar(covlol)

minvarret <- merged_ret %*% minvarweights

#SPY
w2 <- seq(minvarweights[2],1,0.001)
#SPY
w1 <- 1-w2

w_synthetic <- matrix(cbind(w1,w2), ncol = 2, nrow = length(w1))

portfolios <- matrix(0L, ncol = length(w1), nrow = length(daily_logret[,2]))

for(i in 1:length(w1)){

	portfolios[,i] <- (merged_ret) %*% w_synthetic[i, ]

}


expectedreturns <- colMeans(portfolios)*252*100

portdev <- apply(portfolios, MARGIN = c(2), FUN = function(x) sd(x))* sqrt(252)*100

#expectedreturns <- sort(expectedreturns[1:58], decreasing = T)
#portdev <- sort(portdev[1:58], decreasing = T)




ggplot() + geom_line(aes(portdev, expectedreturns))  +
geom_point(aes(sd(minvarret)*sqrt(252)*100,colMeans(minvarret)*252*100)) +
geom_point(aes(portdev[36],expectedreturns[36])) +
geom_text(aes(sd(minvarret)*sqrt(252)*100,colMeans(minvarret)*252*100, 
	label="Minimum variance portfolio"),hjust=-.05, vjust=0) + 
geom_text(aes(portdev[36],expectedreturns[36], 
	label="60/40 bond/equity"),hjust=-0.05, vjust=0)


#60/40 check

w1t <- 0.6
w2t <- 0.4

six4ret <- daily_logret %*% matrix(c(w1t, w2t), ncol=1, nrow=2) 

sd(six4ret)*sqrt(252)*100

mean(six4ret)*252*100





test <- seq(0.1,1,0.01)

portrettest <- matrix(0L, ncol = length(test), nrow = length(daily_logret[,2]))

for(i in 1:length(test)){

	portrettest[,i] <- daily_logret[,2] * test[i]

}

 devtest <- apply(portrettest, MARGIN = c(2), FUN = function(x) sd(x)) * sqrt(252)*100

meantest <- colMeans(portrettest) * 252 * 100


ts.plot(1+cumsum(daily_logret[,1]), col = "red") + lines(1+cumsum(daily_logret[,2]))

















#0.5 risk free asset? 

#assume 1% each year. Then with a notional on 100, we have scale it out everyday:

dailygains = (0.01*100)/252


priceriskfree <- numeric(2515)

priceriskfree[1] <- 100

#this is not completely correct
for(i in 2:2516){

	priceriskfree[i] <- priceriskfree[i-1] + (0.01*priceriskfree[i-1])/252


}


riskfreereturns <- diff(log(priceriskfree))[-1]


returns2 <-  cbind(lel[-c(1,2),], riskfreereturns)


portret <- t(t(riskparity_2dim(newcovariance, 0.1,T)$w) %*% t(returns2))

portret2 <- t(t(riskparity_2dim(newcovariance, 0.1,T)$w[1:2]) %*% t(lel))

ts.plot(1+cumsum(portret)) + lines(1+cumsum(portret2), col="red")


lel <- t(sapply(mergedfrequencies[[10]], function(x) cbind(x[,1], x[,2])))