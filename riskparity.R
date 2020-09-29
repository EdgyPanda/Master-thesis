#risk parity portfolio preliminaries 


source("functions.R")
library(riskParityPortfolio)
library(matlib)

riskparity_2dim <- function(matrix, risktarget, rt = F){
  
  bonds <- matrix[1,1]
  stocks <- matrix[2,2]
  
  w_1 <- sqrt(bonds)^-1 / (sqrt(stocks)^-1 + sqrt(bonds)^-1)
  w_2 <- sqrt(stocks)^-1 / (sqrt(stocks)^-1 + sqrt(bonds)^-1)

  w <- matrix(c(w_1, w_2), ncol=1, nrow=2) #dimnames = list(c(), c("Bond", "Stock"))
  
  portrisk <- sqrt(t(w) %*% (matrix) %*% w) 
  
  riskcont <- c(sqrt(bonds)*w_1, sqrt(stocks)*w_2)

  if(rt){

  	alpha <- risktarget / portrisk

  	w_new <- w %*% alpha

  	portrisk <-  w_new %*% portrisk  #gives marginal risk for each asset. 

  	lout <- list(w_new, portrisk, riskcont)
  	
  	names(lout) <- c("w", "portrisk", "riskcont")

 	return(lout)
  }

  lout <- list(w, portrisk, riskcont)
  names(lout) <- c("w", "portrisk", "riskcont")

  return(lout)
  
}



calccov <- readRDS("calculatedcovariances.rds")

weightsrp <- matrix(0L, nrow = 100, ncol = 2)

weightspalomar <-  matrix(0L, nrow = 100, ncol = 2)


for(i in 1:100){

weightsrp[i,] <- riskparity_2dim(calccov[[1]][[7]][,,i], 0.1, T)$portrisk

weightspalomar[i, ] <- riskParityPortfolio(calccov[[1]][[1]][,,i])$w

}

riskparity_2dim(calccov[[1]][[1]][,,1])$riskcont



lel <- matrix(c(diag(calccov[[1]][[1]][,,1])[1],0,0,diag(calccov[[1]][[1]][,,1])[2]), ncol = 2, nrow = 2)

minvar(lel)

