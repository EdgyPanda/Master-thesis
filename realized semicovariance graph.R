#-------------------------------------------------realized semicovariance graph ----------------------------
bivarcorr <- function(rho){

1/(2*pi) * (rho * acos(-rho) + sqrt(1-rho^2))

}


-2*bivarcorr(--0.75)
	bivarcorr(-0.75)


rhotest <- seq(-1,1,0.01)

bivarcorrdifferential <- function(x){

	1/(2*pi) * (acos(x) - (2*x)/(sqrt(1-x^2)) )


}

c1 <- ggplot() + geom_line(aes(rhotest, bivarcorrdifferential(rhotest)), col = "steelblue", lwd=1) + 
	geom_hline(yintercept = 0) + ylab('Differential') + xlab('Rho') + 
		scale_x_continuous(breaks = seq(-1,1,0.25))

grid.arrange(r1,c1, ncol = 2)

library(ggplot2)
ggplot() + geom_line(aes(rhotest, bivarcorr(rhotest))) + geom_line(aes(rhotest, -2*bivarcorr(-rhotest)))




