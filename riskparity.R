#risk parity portfolio preliminaries 


source("functions.R")
source("Previous functions/projectfunctions.R")
library(riskParityPortfolio)
library(matlib)
library(xts)
library(highfrequency)
library(matlib)
library(MCS)



riskparity_2dim <- function(matrix, risktarget, rt = F){
  #
  #
  #
  #NOTE TO YOURSELF: Palomar uses non-sqrt portrisk, which gives 
  #reasonable risk for the unlevered risk-parity portfolio. Using sqrt
  #
  #REMEMBER TO USE COV(),SD()  etc. FOR CLOSE-TO-CLOSE RETURNS, SINCE REALCOV OVERESTIMATES
  #THE COVARIANCE FOR CLOSE-TO-CLOSE RETURNS

  bonds <- matrix[1,1]
  stocks <- matrix[2,2]
  
  w_1 <- sqrt(bonds)^-1 / (sqrt(stocks)^-1 + sqrt(bonds)^-1)
  w_2 <- sqrt(stocks)^-1 / (sqrt(stocks)^-1 + sqrt(bonds)^-1)

  w <- matrix(c(w_1, w_2), ncol=1, nrow=2) #dimnames = list(c(), c("Bond", "Stock"))
  
  #Palomar uses portfolio variance as portfolio risk, thus no sqrt. 
  portrisk <- as.numeric(sqrt(t(w) %*% (matrix) %*% w))
  
  riskcont <-  w * (matrix %*% w)/portrisk

  relativeriskcont <- (w * (matrix %*% w)) / portrisk^2

  if(rt){

  	alpha <- risktarget / portrisk

  	w_new <- w %*% alpha

  	w_new <- matrix(c(w_new[1], w_new[2]), ncol=1, nrow=2)
  				#here is sqrt, while palomar uses variance, no sqrt.
  	portrisk <-  as.numeric(sqrt((t(w_new) %*% matrix %*% w_new)))  #gives marginal risk for each asset. 

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

riskparity_2dim(cov(merged_ret))

riskparity_2dim(calccov[[5]][[6]][,,1])$portrisk*sqrt(252)*100

calccov <- readRDS("calculatedcovariances.rds")

mergedfrequencies <- readRDS("mergedfrequencies.rds")


#-------------------------------portfolio volatility's sensitivity to correlation--------------------------

#getting average vol for tlt and spy.

meanvolTLT <- mean(sqrt(calccov[[1]][[7]][1,1,]*252))

meanvolSPY <-  mean(sqrt(calccov[[1]][[7]][2,2,]*252))

correlations <- seq(-0.9,0.9,0.01)

sensitivity <- numeric()

sensitivity2 <- numeric()

for(i in 1:length(correlations)){

	newcovariance <- matrix(c(meanvolTLT^2, meanvolTLT*meanvolSPY*correlations[i],
		meanvolTLT*meanvolSPY*correlations[i], meanvolSPY^2), ncol = 2, nrow = 2)

	w2 <- riskparity_2dim(newcovariance,0,F)$w

	portrisk2 <- (riskparity_2dim(newcovariance,0,F)$portrisk)

	sensitivity2[i] <-  (w2[1]*w2[2]*meanvolTLT*meanvolSPY)  /  (portrisk2)


}


library(ggplot2)


#Shows the volatility's sensitivity to correlation for the unlevered risk-parity portfolio. In essence, 
#upscaling and downscaling the portfolios using a leverage parameter only shifts the graph.

ggplot() + geom_line(aes(correlations, sensitivity2), col = "red", lwd = 1) + 
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

	weightsfordistribution[i, ] <- riskparity_2dim(covs[,,i])$w[1:2]


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

#
#
#
#USING CLOSE-TO-CLOSE RETURNS --> COV(), SD() AND NOT INTRADAY MEASURES (since they dont have proper scaling). 
library(alphavantager)
library(xts)

library(ggplot2)
source("functions.R")
source("APIKEY.R")

av_api_key(apikey)

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

covlol <- cov(merged_ret)

colnames(covlol) <- c("TLT", "SPY")

minvarweights <- minvar(covlol)

minvarret <- merged_ret %*% minvarweights

#SPY
w2 <- seq(minvarweights[2],1,0.001)
#SPY
w1 <- 1-w2

w_synthetic <- matrix(cbind(w1,w2), ncol = 2, nrow = length(w1))

portfolios <- matrix(0L, ncol = length(w1), nrow = length(merged_ret[,2]))

for(i in 1:length(w1)){

	portfolios[,i] <- (merged_ret) %*% w_synthetic[i, ]

}

#ADDED 2 AS 2% FOR THE RISKFREE ASSET. IT IS ALL THE WAY DOWN TO THE GGPLOT GRAPH 

expectedreturns <- colMeans(portfolios)*252*100  + 2

portdev <- apply(portfolios, MARGIN = c(2), FUN = function(x) sd(x))* sqrt(252)*100

#expectedreturns <- sort(expectedreturns[1:58], decreasing = T)
#portdev <- sort(portdev[1:58], decreasing = T)

#unlevered risk-parity:



rpunlevered <- riskparity_2dim(covlol)$w

retrpunlevered <- merged_ret %*% rpunlevered 

meanrpunlevered <- mean(retrpunlevered) * 252 *100  + 2

sdrpunlevered <- sd(retrpunlevered) * sqrt(252) * 100

#constructing risk-parity leverage line. 
alpha <- seq(0.95,1.7,0.01)

leverageline <- rpunlevered %*% alpha

leveragelineret <- merged_ret %*% leverageline 

leveragelinemeans <- colMeans(leveragelineret) * 252 * 100 +2
leveragelinestds <- apply(leveragelineret, MARGIN = c(2), FUN = function(x) sd(x))* sqrt(252)*100

rplevered <- riskparity_2dim(covlol, 0.0846, T)$w[1:2]

retrplevered <- merged_ret %*% rplevered 

meanrplevered <- mean(retrplevered) * 252 *100 

sdrplevered <- sd(retrplevered) * sqrt(252) * 100

#Finding the levered risk parity portfolio with same standard deviation as 80/20 portfolio.
root <- uniroot(function(x)  sd(merged_ret %*% (rpunlevered %*% x))*sqrt(252)*100 - portdev[366], interval = c(0,100))$root

rplevered <- rpunlevered %*%  root

retrplevered <- merged_ret %*% rplevered 

meanrplevered <- mean(retrplevered) * 252 *100 + 2

sdrplevered <- sd(retrplevered) * sqrt(252) * 100


ggplot() + geom_line(aes(portdev, expectedreturns, col = "Efficient frontier"), lwd=1)  +
geom_line(aes(leveragelinestds, leveragelinemeans, col ="Leverage line"), lwd=1) +
geom_point(aes(sd(minvarret)*sqrt(252)*100,colMeans(minvarret)*252*100+2)) +
geom_point(aes(portdev[166],expectedreturns[166])) +
geom_point(aes(portdev[366],expectedreturns[366])) +
geom_point(aes(sdrpunlevered,meanrpunlevered)) + 
geom_text(aes(sd(minvarret)*sqrt(252)*100,colMeans(minvarret)*252*100+2, 
	label="Minimum variance portfolio"),hjust=-.05, vjust=0) + 
geom_text(aes(portdev[166],expectedreturns[166], 
	label="60/40 equity/bond"),hjust=-0.05, vjust=0.5) + 
geom_text(aes(portdev[366],expectedreturns[366], 
	label="80/20 equity/bond"),hjust=-0.05, vjust=0.5) + 
geom_text(aes(sdrpunlevered,meanrpunlevered, 
	label="Risk-parity unlevered"),hjust=-0.09, vjust=-0.5) +
geom_point(aes(sdrplevered,meanrplevered)) + 
geom_text(aes(sdrplevered,meanrplevered, 
	label="Risk-parity levered"),hjust=-0.05, vjust=0.5) + ylab("Annualized expected returns (%)")+
xlab("Annualized risk (%)") + 
theme(legend.title = NULL,legend.position = c(0.70, 0.23), legend.background = element_rect(fill="lightblue", size=0.5,
 linetype="solid"), 
	plot.title = element_text(hjust = 0.5, face = "bold"),  axis.title=element_text(size=12))




#------------------------------------CALCULATING VOL-SCALED PORTFOLIOS-------------------------------------

TLT <- as.data.frame(av_get(symbol = "TLT", av_fun = "TIME_SERIES_DAILY", outputsize = "full"))
rownames(TLT) <- TLT$timestamp
SPY <- as.data.frame(av_get(symbol = "SPY", av_fun = "TIME_SERIES_DAILY", outputsize = "full"))
rownames(SPY) <- SPY$timestamp

returns_TLT <- as.xts(diff(log(TLT[,4])), order.by = as.Date(TLT[,1], format='%d/%m/%Y')[-1])
returns_SPY <- as.xts(diff(log(SPY[,4])), order.by = as.Date(SPY[,1])[-1])


returns_TLT <-  returns_TLT[seq(from= as.Date('2010-01-04'), to = as.Date('2019-12-31'), by=1), ]
returns_SPY <- returns_SPY[seq(from= as.Date('2010-01-04'), to = as.Date('2019-12-31'), by=1), ]



merged_ret <- cbind(returns_TLT, returns_SPY)


#KEY NOTE: IT IS IMPORTANT THAT BOTH ROLLING FORECASTS START WITH SAME FIXED WINDOW TO ACHIEVE CONSISTENCY.
#OTHERWISE ROLLING VARS WILL BE UNDER/OVER ESTIMATED. 

varsmerged <- na.omit(rollapply(merged_ret, 60, function(x) realCov(x), by.column = T, align = 'right', by = 1)) 

stdTLT <- sqrt(varsmerged[,1]) #* sqrt(252)
stdSPY <- sqrt(varsmerged[,2]) #* sqrt(252)
#target is in terms of var. That implies that you need to think: For what variance do I get 10% vol. 
#Ie. sqrt(0.01)=0.1.NOPE

#decimalnumbers

#finds when first rolling sd.dev is calculated
start <- length(returns_TLT) - length(varsmerged[,1])

target <- 0.10

scaledTLT <- (target/stdTLT)  *  returns_TLT[-c(1:start), ] #* 252

scaledSPY <- (target/stdSPY) * returns_SPY[-c(1:start), ] #* 252

riskfreealloc <- (target/stdTLT)  + (target/stdSPY)


#60/40 portfolio:



portret6040 <- 1*scaledSPY + 0*scaledTLT


rollingdevportret6040 <- sqrt(na.omit(rollapply(portret6040, 60, function(x) realCov(x), by.column = F, align = 'right')))

ggplot() + geom_line(aes(1:length(rollingdevportret6040), rollingdevportret6040)) + geom_hline(yintercept = target)




#intraday tryout:

calccov <- readRDS("calculatedcovariances.rds")

mergedfrequencies <- readRDS("mergedfrequencies.rds")

total30minreturns <- do.call.rbind(mergedfrequencies[[9]])

varsmerged2 <- na.omit(rollapply(total30minreturns, 12*60, function(x) (realCov(x)), by.column = T, 
	align = 'right'))


thirtymintlt <- sqrt(varsmerged2[,1]) 
thirtyminspy <- sqrt(varsmerged2[,2]) 

start <- length(total30minreturns[,1]) - length(varsmerged2[,1])

scaled30mintlt <- (target/thirtymintlt)  *  total30minreturns[-c(1:start), 1] 
scaled30minspy <- (target/thirtyminspy)  *  total30minreturns[-c(1:start), 2] 

thirtyminport <- 1*scaled30minspy + scaled30mintlt * 0

rollingdevportret6040thirtymin <- (na.omit(rollapply(thirtyminport, 12*60, function(x) sqrt(realCov(x)), by.column = F, 
	align = 'right', by = 12)))

ggplot() + geom_line(aes(1:length(rollingdevportret6040thirtymin), rollingdevportret6040thirtymin, col="30 min")) + 
geom_hline(yintercept = target) + 
geom_line(aes(1:length(rollingdevportret6040), rollingdevportret6040, col="Daily"))



#5min


total5minreturns <- do.call.rbind(mergedfrequencies[[7]])

varsmerged5min <- na.omit(rollapply(total5minreturns, 78*60, function(x) (realCov(x)), by.column = T, 
	align = 'right'))


fivemintlt <- sqrt(varsmerged5min[,1]) 
fiveminspy <- sqrt(varsmerged5min[,2]) 

start <- length(total5minreturns[,1]) - length(varsmerged5min[,1])

scaled5mintlt <- (target/fivemintlt)  *  total5minreturns[-c(1:start), 1] 
scaled5minspy <- (target/fiveminspy)  *  total5minreturns[-c(1:start), 2] 

fiveminport <- 1*scaled5minspy + scaled5mintlt * 0

rollingdevportret6040fivemin <- (na.omit(rollapply(fiveminport, 78*60, function(x) sqrt(realCov(x)), by.column = F, 
	align = 'right', by = 78)))



mergedfrequencies[[7]]


#------------------------------------------TRYING WITH SIM DATA---------------------------------

#do ewma to figure it out. 

# do variances instead. 

brownian <- function(J, N, sigma, alpha, beta, gamma,rho ,time){ 


simulations <- brownian(2,23400*2,10,0, 0, 1,0, 1)
logretsims <- diff(log(simulations[,3]))[-1]
logretsims2 <- diff(log(simulations[,4]))[-1]
logretsims <- cbind(logretsims, logretsims2)


logretsims <- xts(logretsims, order.by = as.Date(1:length(logretsims[,1])))

sqrt(realCov(logretsims)) 

sdlogret3 <- sqrt(na.omit(rollapply(logretsims[,1], 23400*2, function(x) (realCov(x)), by.column = F, align = 'left'))) #* sqrt(252)
sdlogret4 <- (na.omit(rollapply(logretsims[,2], 23400*2, function(x) sqrt(realCov(x)), by.column = F, align = 'left'))) #* sqrt(252)

#sdlogret3 <- ewma.filter2006((logretsims),F,T)

#sdlogret32 <- apply(sdlogret3, MARGIN = c(3), FUN = function(x) x[1,1])

#sdlogret3 <- sdlogret3 * sqrt(252)
#sdlogret4 <- sdlogret4 * sqrt(252)


scaledret3 <- logretsims[,1] * (0.1/sdlogret3) 
sclaedret4 <- logretsims[-c(1:46799),2] * (10/sdlogret4) 


scaledport <- scaledret3 * 1 

sdport <-  ewma.filter2006((scaledport),F,F)  #sqrt(na.omit(rollapply(scaledport, 38400, function(x) (realCov(x)), by.column = F, align = 'left')))




ts.plot(sdport)