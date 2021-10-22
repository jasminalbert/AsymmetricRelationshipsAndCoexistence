install.packages("edfun") #package for empirical pdfs, cdfs...
library(edfun)
library(MASS)

#phi as in theory of main text
phi <- function(x){pnorm(x)} 	 #normal cdf
phi_inv <- function(x){qnorm(x)} #normal inverse cdf 

cdfTran <- function(x){ 		#transforms var by its cdf 
	dist <- edfun(x, dfun=NULL)	#makes uniform
	cdfTran <- dist$pfun(x)
	return(cdfTran)}
	
cdfInv <- function(x,y){		#tranform var by inverse cdf of 
	dist <- edfun(x, dfun=NULL)	#desired distribution
	cdfInv <- dist$qfun(y)
	return(cdfInv)}	


Ecdf<-cdfTran(E); hist(Ecdf)		#transform into uniform
Enorm <- phi_inv(Ecdf); hist(Enorm)	#make normal
Enorm <- Enorm[Enorm != Inf]		#take out infinity

Ccdf<-cdfTran(C1); hist(Ccdf)		#do again for C1
Cnorm <- phi_inv(Ccdf); hist(Cnorm)
Cnorm <- Cnorm[Cnorm != Inf]

#not same length so do multiple random cuts and caculate cov and get mean
rhovec <- NA
for (i in 1:20){
	Cnorm <- sample(Cnorm, length(Cnorm))[1:length(Enorm)] #randomly cut Cnorm to be length of Enorm
	rhovec[i]<-cov(Enorm, Cnorm)
}
rho <- mean(rhovec)

Sigma <- matrix(c(1, rho, rho, 1), nrow=2)	#rho = rho umlaut	
#make bivariate normal with rho umlaut
ec <- mvrnorm(length(Enorm), mu=c(0,0), Sigma=Sigma) 
cov(ec)
phie <- phi(ec[,1]); hist(phie)		#transform into uniform
phic <- phi(ec[,2]); hist(phic)

#transform into original E and C distribution
Epartial<-cdfInv(E, phie); hist(Epartial)  	
Cpartial<-cdfInv(C1, phic); hist(Cpartial)	#symmetric?

plot(Epartial, Cpartial)
cov(Epartial, Cpartial)
plot(E,C1)
cov(E,C1)















