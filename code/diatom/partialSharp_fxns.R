#functions required for making partial sharp distributions
#to parse pure correlation (correlation per se) from tail asymmetries
  #alignranks, normcor, makePsharp

##libraries used (invoked with ::): copule, stats, mvtnorm
  
### alignranks ###
## rearranges data so ranks match with sims
# part of various copula surrogate algorithms
# for each column j of simulation sims[,], the
# elements of dat[,j] are permuted so the ranks  
# of the result are aligned with those of sims[,j]
#ARGS:
  #dat		T by d matrix, each column sorted!!
  #sims 		T by d matrix
#OUT:
  #T by d matrix
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

### normcor ###
# normalizes two variables then computes correlation
#ARGS:
  #X			vector length T
  #Y 		vector length T
#OUT:
  #single value
normcor <- function(X, Y){
	rank_X <- copula::pobs(X)
	rank_Y <- copula::pobs(Y)
	rho <- stats::cor(stats::qnorm(rank_X), stats::qnorm(rank_Y))
	return(rho)
}

### makePsharp ###
# runs normcor and alignranks with bivariate standard 
# normal sims with Sigma matrix where rho = normcor(X, Y)
# and var(X) = var(Y) = 1
# generates rep number of sims to get rep=k number 
# of res[,,k] with columns res[,j,] as the partial
# sharped variables 
#ARGS:
  #X			vector length T
  #Y 		vector length T
  #rep		single value - number of repeats
#OUT:
  #T x 2 x rep array
makePsharp <- function(X, Y, reps){
	
	rho <- normcor(X, Y)
	
	X <- sort(X); Y <- sort(Y)
	
	Sigma <- matrix(c(1, rho, rho, 1), nrow=2)
	
	dat <- cbind(X, Y)
	
	res <- array(NA, dim=c(dim(dat), reps))
		
	for (i in 1:reps){
		sims <- mvtnorm::rmvnorm(n=length(X), sigma=Sigma)
		res[,,i] <- alignranks(dat, sims)
	}
	return(res)
}






