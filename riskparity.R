#risk parity portfolio preliminaries 


source("functions.R")
library(riskParityPortfolio)
library(matlib)
library(xts)
library(highfrequency)
library(matlib)
library(MCS)



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

ggplot() + geom_line(aes(correlations, sensitivity2), col = "red", lwd = 1) + 
  scale_x_continuous(breaks = round(seq(-0.9,0.9, by = 0.1),1)) + ylab("portfolio volatility (%)")


(riskparity_2dim(newcovariance, 0.1,T)$w) * (newcovariance %*% riskparity_2dim(newcovariance, 0.1,T)$w)/0.1

riskparity_2dim(newcovariance, 0.1,F)

riskparity_2dim(newcovariance, 0.1,T)



