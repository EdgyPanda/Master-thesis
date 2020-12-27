#risk parity portfolio preliminaries 


source("functions.R")
source("Previous functions/projectfunctions.R")
library(riskParityPortfolio)
library(matlib)
library(xts)
library(highfrequency)
library(matlib)
library(MCS)
library(PerformanceAnalytics)


#####################################################################################################
#
#
#                                  GETTING CLOSE-TO-CLOSE RETURNS
#
#
#####################################################################################################
library(alphavantager)
library(ggplot2)
#source("functions.R")
source("APIKEY.R")

av_api_key(apikey)

TLT <- as.data.frame(av_get(symbol = "TLT", av_fun = "TIME_SERIES_DAILY", outputsize = "full"))
rownames(TLT) <- TLT$timestamp
SPY <- as.data.frame(av_get(symbol = "SPY", av_fun = "TIME_SERIES_DAILY", outputsize = "full"))
rownames(SPY) <- SPY$timestamp

returns_TLT <- as.xts(diff(log(TLT[,4])), order.by = as.Date(TLT[,1], format='%d/%m/%Y')[-1])
returns_SPY <- as.xts(diff(log(SPY[,4])), order.by = as.Date(SPY[,1])[-1])

#####################################################################################################
#
#
#                     Calculating the risk-free rate from a 3-month T-bill:
#
#
#####################################################################################################


#Calculating the risk-free rate from a 3-month T-bill:
#3month T-bill. Has no coupons 
dtb3 <- read.csv("DTB3.csv", header =  T)

dtb3[,2] <- as.numeric(levels(dtb3[,2]))[dtb3[,2]]

dtb3[,2] <- na.approx(dtb3[,2])       #dtb3[!is.na(dtb3[,2]), ]



#this is how you should do it:
logriskfreerate <- log(1 + dtb3[,2]/(100*365)) 

logriskfreerate <- xts(logriskfreerate, order.by = as.Date(dtb3[,1]))


test <- (1+dtb3[,2]/100)^(1/(365))-1

ggplot() + geom_line(aes(1:length(test), test, col = "test")) + geom_line(aes(1:length(test), logriskfreerate, col="RF")) 

ts.plot(logriskfreerate*100, col = "red")


intersectMulti <- function(x=list()){
 for(i in 2:length(x)){
    if(i==2) foo <- x[[i-1]]
    foo <- intersect(foo,x[[i]]) #find intersection between ith and previous
 }
 return(x[[1]][match(foo, x[[1]])]) #get original to retain format
}

indexes <- intersectMulti(list(index(logriskfreerate), index(returns_TLT)))

logriskfreerate <- logriskfreerate[indexes]


#These are now excess returns.
returns_TLT <-  returns_TLT[seq(from= as.Date('2010-01-02'), to = as.Date('2019-12-31'), by=1), ] - logriskfreerate
returns_SPY <- returns_SPY[seq(from= as.Date('2010-01-02'), to = as.Date('2019-12-31'), by=1), ] - logriskfreerate

merged_ret <- cbind(returns_TLT, returns_SPY)



#####################################################################################################
#
#
#                                     PORTFOLIO ANALYSIS:
#
#
#####################################################################################################
library(CVXR)

portolioMaxSharpeRatio <- function(mu, Sigma) {
  w_ <- Variable(nrow(Sigma))
  prob <- Problem(Minimize(quad_form(w_, Sigma)),
                  constraints = list(w_ >= 0, t(mu) %*% w_ == 1))
  result <- solve(prob)
  return(as.vector(result$getValue(w_)/sum(result$getValue(w_))))
}




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

calccov <- readRDS("calculatedcovariances.rds")

mergedfrequencies <- readRDS("mergedfrequencies.rds")

riskparity_2dim(cov(merged_ret))

riskparity_2dim(calccov[[5]][[6]][,,1])$portrisk*sqrt(252)*100



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

p1 <- ggplot() + geom_line(aes(correlations, sensitivity2*100), col = "red", lwd = 1) + 
  scale_x_continuous(breaks = round(seq(-0.9,0.9, by = 0.1),1)) + ylab("portfolio volatility (%)") + xlab("Correlation")

#-------------------------------weight distribution dependent on excess returns---------------------------------------

 #done for returns until you have control over a risk-free asset.


#source("Previous functions/projectfunctions.R")

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
#
#
#
#
###################################################


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


expectedreturns <- colMeans(portfolios)*252*100  

portdev <- apply(portfolios, MARGIN = c(2), FUN = function(x) sd(x))* sqrt(252)*100

#expectedreturns <- sort(expectedreturns[1:58], decreasing = T)
#portdev <- sort(portdev[1:58], decreasing = T)

#unlevered risk-parity:



rpunlevered <- riskparity_2dim(covlol)$w

retrpunlevered <- merged_ret %*% rpunlevered 

meanrpunlevered <- mean(retrpunlevered) * 252 *100  

sdrpunlevered <- sd(retrpunlevered) * sqrt(252) * 100

#constructing risk-parity leverage line. For alpha = 0, then it will obviously be in origo. 
#it is not the capital market line. 
#0.95
alpha <- seq(0.95,1.7,0.01)

leverageline <- rpunlevered %*% alpha

riskfreeassetcont <- numeric()
for(i in 1:length(alpha)){

	riskfreeassetcont[i] <- uniroot(function(x) rowSums(t(leverageline))[i] + x - 1, interval = c(-100,100))$root

}

leverageline <- cbind(t(leverageline), riskfreeassetcont)
leverageline <- t(leverageline)


leveragelineret <- cbind(merged_ret, logriskfreerate) %*% leverageline 

leveragelinemeans <- colMeans(leveragelineret) * 252 * 100 
leveragelinestds <- apply(leveragelineret, MARGIN = c(2), FUN = function(x) sd(x))* sqrt(252)*100

rplevered <- riskparity_2dim(covlol, 0.0846, T)$w[1:2]

retrplevered <- merged_ret %*% rplevered 

meanrplevered <- mean(retrplevered) * 252 *100 

sdrplevered <- sd(retrplevered) * sqrt(252) * 100

#Finding the levered risk parity portfolio with same standard deviation as 80/20 portfolio.
root <- uniroot(function(x)  sd(merged_ret %*% (rpunlevered %*% x))*sqrt(252)*100 - portdev[366], interval = c(0,100))$root

rplevered <- rpunlevered %*%  root

retrplevered <- cbind(merged_ret,logriskfreerate) %*% t(cbind(t(rplevered), -0.4647607))

meanrplevered <- mean(retrplevered) * 252 *100

sdrplevered <- sd(retrplevered) * sqrt(252) * 100

#CML:

tangent <- portolioMaxSharpeRatio(colMeans(merged_ret*100)*252,covlol*10000)

tangentret <- merged_ret %*% tangent

tangentexpectedret <- mean(tangentret) * 100 * 252

tangentdeviation <- sd(tangentret) * 100 * sqrt(252)


sharpetangent <- (tangentexpectedret - mean(logriskfreerate) * 100 * 252)/tangentdeviation

#leverageline computes the weights of the risky-assets and riskfree. 

library(PerformanceAnalytics)

CML <-  mean(logriskfreerate) * 100 * 252 + sharpetangent * leveragelinestds


p2 <- ggplot() + geom_line(aes(portdev, expectedreturns, col = "Efficient frontier"), lwd=1)  +
geom_line(aes(leveragelinestds, leveragelinemeans, col ="Leverage line"), lwd=1) +
geom_point(aes(sd(minvarret)*sqrt(252)*100,colMeans(minvarret)*252*100)) +
geom_line(aes(leveragelinestds, CML, col ="CML"), lwd=1) +
geom_point(aes(portdev[166],expectedreturns[166])) +
geom_point(aes(portdev[366],expectedreturns[366])) +
geom_point(aes(tangentdeviation,tangentexpectedret)) +
geom_point(aes(sdrpunlevered,meanrpunlevered)) + 
geom_text(aes(sd(minvarret)*sqrt(252)*100,colMeans(minvarret)*252*100, 
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


library(gridExtra)

p3 <- grid.arrange(p1, p2, ncol=2)

ggsave(p3, file="portfoliosensandefficientfront.eps", device = "eps")


#------------------------------------CALCULATING VOL-SCALED PORTFOLIOS-------------------------------------

TLT <- as.data.frame(av_get(symbol = "TLT", av_fun = "TIME_SERIES_DAILY", outputsize = "full"))
rownames(TLT) <- TLT$timestamp
SPY <- as.data.frame(av_get(symbol = "SPY", av_fun = "TIME_SERIES_DAILY", outputsize = "full"))
rownames(SPY) <- SPY$timestamp

returns_TLT <- as.xts(diff(log(TLT[,4])), order.by = as.Date(TLT[,1], format='%d/%m/%Y')[-1])
returns_SPY <- as.xts(diff(log(SPY[,4])), order.by = as.Date(SPY[,1])[-1])


returns_TLT <-  returns_TLT[seq(from= as.Date('2010-01-04'), to = as.Date('2019-12-31'), by=1), ] - logriskfreerate
returns_SPY <- returns_SPY[seq(from= as.Date('2010-01-04'), to = as.Date('2019-12-31'), by=1), ] - logriskfreerate



merged_ret <- cbind(returns_TLT, returns_SPY)


#KEY NOTE: IT IS IMPORTANT THAT BOTH ROLLING FORECASTS START WITH SAME FIXED WINDOW TO ACHIEVE CONSISTENCY.
#OTHERWISE ROLLING VARS WILL BE UNDER/OVER ESTIMATED. 

varsmerged <- ewma.filter(merged_ret, 30, F, T)


stdTLT <- sqrt(varsmerged[1,1,]) #* sqrt(252)
stdSPY <- sqrt(varsmerged[2,2,]) #* sqrt(252)
#target is in terms of var. That implies that you need to think: For what variance do I get 10% vol. 
#Ie. sqrt(0.01)=0.1.NOPE

#decimalnumbers

#finds when first rolling sd.dev is calculated
start <- length(returns_TLT) - length(varsmerged[1,1,])

target <- 0.1

scaledTLT <- (target/stdTLT)  *  returns_TLT #* 252

scaledSPY <- (target/stdSPY) * returns_SPY #* 252

riskfreealloc <- (target/stdTLT)  + (target/stdSPY)


#60/40 portfolio:

portret6040 <- 0.6*scaledSPY + 0.4*scaledTLT


rollingdevportret6040 <- ewma.filter(portret6040, 30, F, T)
rollingdevportret6040 <- xts(sqrt(rollingdevportret6040[1,1,]), order.by = index(returns_TLT))


ggplot() + geom_line(aes(1:length(rollingdevportret6040), rollingdevportret6040)) + geom_hline(yintercept = target)


#leverage returns it to the risk-target


levparam <- 0.1 / rollingdevportret6040

portret6040lev <- levparam * 0.6 * scaledSPY + levparam * 0.4 * scaledTLT


rollingdevportret6040levered <- ewma.filter(portret6040lev, 30, F, T)

rollingdevportret6040levered <- xts(sqrt(rollingdevportret6040levered[1,1,]), order.by = index(returns_TLT))

ggplot() + geom_line(aes(index(rollingdevportret6040levered), rollingdevportret6040levered)) + geom_hline(yintercept = target)


#Scaled vs unscaled returns across assets: NOPE

unscaledacross <- rowMeans(merged_ret)

scaledacross <- rowMeans(cbind(scaledTLT, scaledSPY))

#cumulative frequency. 
ggplot() + geom_line(aes(index(merged_ret), 1+cumsum(unscaledacross), col = "Unscaled")) + 
geom_line(aes(index(merged_ret), 1+cumsum(scaledacross)*(1/10), col = "Scaled")) 


#1 year rolling returns: NOPE 

unscaledrolling <- na.omit(rollapply(unscaledacross, 252, function(x) mean(x), by.column = F, 
	align = 'right'))
scaledrolling <- na.omit(rollapply(scaledacross, 252, function(x) mean(x), by.column = F, 
	align = 'right'))

ggplot() + geom_line(aes(index(merged_ret)[252:2516], unscaledrolling, col = "Unscaled")) + 
geom_line(aes(index(merged_ret)[252:2516], scaledrolling, col = "Scaled"), alpha = 0.4) 



#leverage graph where you see for each weight, how much it undertargets the original risk target

weightSPY <- seq(0,1,0.01)
weightTLT <- sort(weightSPY, T)

portfolioreturnsfrontier <- matrix(0L, ncol = length(weightTLT), nrow= 2516)

for(i in 1:length(weightSPY)){

	portfolioreturnsfrontier[, i] <- weightSPY[i]*scaledSPY + weightTLT[i]*scaledTLT

}

rollingstd <- matrix(0L, ncol = length(weightTLT), nrow= 2516)

for(i in 1:length(weightSPY)){

	rollingstd[,i] <- ewma.filter(xts(portfolioreturnsfrontier[,i], order.by = as.Date(1:2516)), 30, F, T)

}

rollingstd2 <- apply(rollingstd, MARGIN = c(2), FUN = function(x) sqrt(x))


#average absolute deviation away from risk-target:

avgabsdev <- apply(rollingstd2, MARGIN = c(2), FUN = function(x) (0.1-mean(x)))


#leverage parameter:

leverageparams <- apply(rollingstd2, MARGIN = c(2), FUN = function(x) mean(0.1/x))


p4 <- ggplot() + geom_line(aes(weightTLT, avgabsdev, col = "Avg. dev from risk-target"), lwd = 1)  +
geom_line(aes(weightTLT, (leverageparams-1)/40, col = "Avg. leverage"), lwd = 1)  +
scale_y_continuous(sec.axis = sec_axis(~.*40+1, name = "Leverage")) + 
theme(legend.justification=c(0,1), legend.position=c(0.67,0.97),
	legend.background = element_rect(fill="lightblue",
    size=0.5, linetype="solid", 
    colour ="darkblue"), 
    legend.text = element_text(colour="black", size=8, face="bold")) + xlab("Weight") + ylab("Deviation")

ggsave(p4, file="weightandlev.eps", device = "eps")

