###################################################
# FUNCTIONS REQUIRED FOR MAKING PARTIAL SHARPS VARS
# to parse pure correlation from tail asymmetries
###################################################
require(copula); require(mvtnorm)

### alignranks: 
## rearranges data so ranks match with sims
# part of various copula surrogate algorithms
# for each column j of simulation sims[,], the
# elements of dat[,j] are permuted so the ranks  
# of the result are aligned with those of sims[,j]
## ARGS
# dat		T by d matrix, each column sorted!!
# sims 		T by d matrix
## OUTPUT	T by d matrix
## note: algorithm assumes column sims[,j] has no ties
# for any j

alignranks <- function(dat, sims){
	simsrk <- apply(FUN=rank, MARGIN=2, X=sims)
	res <- array(NA, dim(sims))
	
	for (counter in 1:dim(dat)[2]){
		res[,counter] <- dat[simsrk[,counter], counter]
	}
	return(res)
}

### normcor:
## normalizes two variables then computes correlation
## ARGS
# X			vector length T
# Y 		vector length T
## OUTPUT	single value

normcor <- function(X, Y){
	rank_X <- pobs(X)
	rank_Y <- pobs(Y)
	rho <- cor(qnorm(rank_X), qnorm(rank_Y))
	return(rho)
}

### makePsharp:
## runs normcor and alignranks with bivariate standard 
## normal sims with Sigma matrix where rho = normcor(X, Y)
## and var(X) = var(Y) = 1
# generates rep number of sims to get rep=k number 
# of res[,,k] with columns res[,j,] as the partial
# sharped variables 
## ARGS
# X			vector length T
# Y 		vector length T
# rep		single value - number of repeats
## OUTPUT	T x 2 x rep array

makePsharp <- function(X, Y, reps){
	
	rho <- normcor(X, Y)
	
	X <- sort(X); Y <- sort(Y)
	
	Sigma <- matrix(c(1, rho, rho, 1), nrow=2)
	
	dat <- cbind(X, Y)
	
	res <- array(NA, dim=c(dim(dat), reps))
		
	for (i in 1:reps){
		sims <- rmvnorm(n=length(X), sigma=Sigma)
		res[,,i] <- alignranks(dat, sims)
	}
	return(res)
}






