source("align_ranks.R")
source("ForcedChemo_Chesson-C.R")
require(copula)
require(mvtnorm)

normcor <- function(X,Y){
	rank_X <- pobs(X)
	rank_Y <- pobs(Y)
	rho <- cor(qnorm(rank_X), qnorm(rank_Y))
	return(rho)
}

rho <- c(normcor(E, C1), normcor(E, C2))

dat1<- cbind(sort(E), sort(C1))
dat2<- cbind(sort(E), sort(C2))
Sigma <- list(	matrix(c(1, rho[1], rho[1], 1), nrow=2),
				matrix(c(1, rho[2], rho[2], 1), nrow=2))

reps <- 500
res1 <- array(NA, dim=c(dim(dat1),reps))
res2 <- res1
normcorSims1 <- NA
normcorSims2 <- NA
r1psharpsims <- NA
r2psharpsims <- NA


for (i in 1:reps){
	
	sims1<-rmvnorm(n=length(E), sigma=Sigma[[1]])
	res1[,,i]<-alignranks(dat1,sims1)
	normcorSims1[i] <- normcor(res1[,1,i],res1[,2,i])
	
	r1psharpsims[i] <- mean(r1C(res1[,1,i],res1[,2,i],parms))
	

	sims2<-rmvnorm(n=length(E), sigma=Sigma[[2]])
	res2[,,i]<-alignranks(dat2,sims2)
	normcorSims2[i] <- normcor(res2[,1,i],res2[,2,i])
	
	r2psharpsims[i] <- mean(r2C(res2[,1,i],res2[,2,i],parms))

	if (i%%50==0){print(i)}
}
#check
hist(normcorSims1, main='histogram of normcor1')
abline(v=rho[1], col='red') #good
hist(normcorSims2, main='histogram of normcor2')
abline(v=rho[2], col='red')


hist(r1psharpsims)
hist(r2psharpsims)

r1psharp <- median(r1psharpsims)
r2psharp <- median(r2psharpsims)

r1bar - r1sharp # = epsilon(EC)
epsilon1.ATA <- r1bar - r1psharp # = epsilon[EC]
epsilon1.cor <- r1psharp - r1sharp # = epsilon[E||C]
(r1bar - r1psharp) + (r1psharp - r1sharp) # = epsilon(EC)

r2bar - r2sharp # = epsilon(EC)
epsilon2.ATA <- r2bar - r2psharp # = epsilon[EC]
epsilon2.cor <- r2psharp - r2sharp # = epsilon[E||C]
(r2bar - r2psharp) + (r2psharp - r2sharp) # = epsilon(EC)


(r1bar - r1sharp) - (r2bar - r2sharp) #storage effect, Delta(EC)
epsilon1.ATA - epsilon2.ATA #Delta[EC], contribution of ATA
epsilon1.cor - epsilon2.cor #Delta[E||C], correlation per se
(epsilon1.ATA - epsilon2.ATA) + (epsilon1.cor - epsilon2.cor)








