#function: co.periods
#computes lengths of periods of noticeable coexistence
#ARGS:
# pop - T x 6 df from popsim
# dom - dominance threshold, 0<dom<1
# N - total spaces
#OUTPUT: a list of three, co.periods for each noise regime

co.periods <- function(pop, dom, N){
	
	N1 <- pop[,1]
		
	coexist <- N1 < dom*N & N1 > (1-dom)*N
		
	RLE <- rle(coexist)
	
	return(RLE)
}




#test 
#source("./new_functions/popsim.R")

#testnoise <- matrix(rnorm(1200), ncol=12)

#pop <- popsim(exp(testnoise), 50, 25, 0.5)

#dom <- 0.90; N <- 50

#co.periods(pop, dom, N)

