install.packages("edfun")
library(edfun)
library(MASS)

phi <- function(x){pnorm(x)}
phi_inv <- function(x){qnorm(x)}

cdfTran <- function(x){
	dist <- edfun(x, dfun=NULL)
	cdfTran <- dist$pfun(x)
	return(cdfTran)}
	
cdfInv <- function(x,y){
	dist <- edfun(x, dfun=NULL)
	cdfInv <- dist$qfun(y)
	return(cdfInv)}	


Ecdf<-cdfTran(E); hist(Ecdf)
Enorm <- phi_inv(Ecdf); hist(Enorm)
Enorm <- Enorm[Enorm != Inf]

Ccdf<-cdfTran(C1); hist(Ccdf)
Cnorm <- phi_inv(Ccdf); hist(Cnorm)
Cnorm <- Cnorm[Cnorm != Inf]

rhovec <- NA
for (i in 1:20){
	Cnorm <- sample(Cnorm, length(Cnorm))[1:length(Enorm)] #randomly cut Cnorm to be length of Enorm
	rhovec[i]<-cov(Enorm, Cnorm)
}
rho <- mean(rhovec)

Sigma <- matrix(c(1, rho, rho, 1), nrow=2)
ec <- mvrnorm(length(Enorm), mu=c(0,0), Sigma=Sigma)
cov(ec)
phie <- phi(ec[,1]); hist(phie)
phic <- phi(ec[,2]); hist(phic)

Epartial<-cdfInv(E, phie); hist(Epartial) #symmetric 
Cpartial<-cdfInv(C1, phic); hist(Cpartial)

plot(Epartial, Cpartial)
cov(Epartial, Cpartial)
plot(E,C1)
cov(E,C1)















