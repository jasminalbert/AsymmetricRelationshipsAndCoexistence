
# E(t) -pobs-> r(E(t)) -> qnorm(r(E(t))) --> rho_umlaut

# rho_umlaut -rmvnorm (mvtnorm)-> e(t) -> pnorm(e) -aligned ranks-> 

# -> symmetric E and C
# compute correlations of E and C partials acorss differrent simulations, look at distributions across simulations
# compare to histograms of pereason and spearman correlation of original E and C
# write function called normcor (normalized correlation)
# should be same up to sampling variation
# normcor(x,y)
# return (cor(qnorm(pobs(x), qnorm(pobs(y)))))
# take norm cor of E and C and norm cor of E partial and C partial
# if norm cor is in the distribution proceede to look at peearson and spearman 
# use partial sharps to calculate partial sharp will get a distribution of IGR 
# fix notation first!
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

norm <- normcor(E, C1)
pear <- cor(E, C1, method='pearson')
kend <- cor(E, C1, method='kendall')
spear <- cor(E, C1, method='spearman')
rho <- c(norm, pear, kend, spear)

dat<- cbind(sort(E), sort(C1))
Sigma <- matrix(c(1, rho[1], rho[1], 1), nrow=2)

reps <- 500
res <- array(NA, dim=c(dim(dat),reps))
normcorSims <- NA
pearsonSims <- NA
spearmanSims <- NA

for (i in 1:reps){
	
	sims<-rmvnorm(n=length(E), sigma=Sigma)
	res[,,i]<-alignranks(dat,sims)
	normcorSims[i] <- normcor(res[,1,i],res[,2,i])
	pearsonSims[i] <- cor(res[,,i])[1,2]
	spearmanSims[i] <- cor(res[,,i], method = 'spearman')[1,2]

	if (i%%50==0){print(i)}
}
#check
hist(normcorSims, main='histogram of normcor')
abline(v=rho[1], col='red') #good

hist(spearmanSims, xlim=c(-1,-0.5),main='histogram of Scor')
abline(v=rho[4], col='red') #bad

hist(pearsonSims, xlim=c(-1,-0.5), main='histogram of Pcor')
abline(v=rho[2], col='red') #bad

r1psharpsims <- NA

for (i in 1:reps){
	r1psharpsims[i] <- mean(r1C(res[,1,i],res[,2,i],parms))
	#r2psharp[i] <- mean(r2C())
}

hist(r1psharpsims)

r1psharp <- median(r1psharpsims)













