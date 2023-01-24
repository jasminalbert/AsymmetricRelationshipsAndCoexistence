
numRes_loc <- "../results_numeric/"
betanoise_loc <- paste0(numRes_loc, "betanoise.RDS")
B <-readRDS(betanoise_loc)

### popsim ###
#ARGS:
#bnoise	2 column matrix of species noise
#N      total abundance (N1 + N2), numeric
#N1     initial abundance of species i, numeric
#up_i	upper bound of species i noise
#up_j 	upper bound of species j noise
#delta  death rate (0-1), numeric
#plot   T/F
#len	length x-axis of plot
#OUT:
#data.frame nrow M with columns N1 and N2 that are the abundance of sp1 and sp2 and other metrics over time
  
popsim <- function(N,N1,up_i,up_j,delta,bnoise,plot=TRUE, len=100) {
  
  #convert upper bound to sigma and mu
  sigma_i <- up_i 
  sigma_j <- up_j
  mu_i <- up_i/2 
  mu_j <- up_j/2
  
  #length
  M <- nrow(bnoise)
  B_i <- sigma_i * bnoise[,1]
  B_j <- sigma_j * bnoise[,2] 

    
  #species 2 initial abundance
  N2 <- N-N1
  tot_new_juvs_vec <- {}
    
  #simulate competitive lottery model to get abundances at each time step
  for (t in 1:M){
    #total new juveniles is sum of fecundity times abundance for each species
    tot_new_juvs <- (B_i[t]*N1[t]) + (B_j[t]*N2[t])
    tot_new_juvs_vec[t] <- tot_new_juvs
      
    #lottery model equation
    N1[t+1] <- (1-delta)*N1[t] + delta*N*((B_i[t]*N1[t])/tot_new_juvs)	
    N2[t+1] <- (1-delta)*N2[t] + delta*N*((B_j[t]*N2[t])/tot_new_juvs)
  }
  #
  if (plot==TRUE){
    plot(N1[1:len], type='l', ylim=c(0,N), xaxt='n', ylab='', yaxt='n',xlab='')
  }
  N2<-N2[-1000001]
  N1<-N1[-1000001]
  #lines(N2[1:len], col="red")
  res<- data.frame(N1=N1, N2=N2,B_i,B_i_juvs=B_i*N1,B_j, B_j_juvs=B_j*N2,
                   tot_new_juvs=tot_new_juvs_vec, 
                   i_new_adults=delta*N*((B_i*N1)/tot_new_juvs_vec), 
                    j_new_adults=delta*N*((B_j*N2)/tot_new_juvs_vec))
  res$i_wins <- res$i_new_adults>res$j_new_adults
  return(res)
}

#simulate pops with beta fecunsities and only difference being ATAs

### location to save results ###
popsim_loc <- paste0(numRes_loc, "betapopsim.RData")

N <- 50
N1 <- 40
up_i <- 1 #upper bound of i
up_j <- 1.05 
d <- 0.7

popL <- popsim(N,N1,up_i,up_j,d,B$l,len=500)
popS <- popsim(N,N1,up_i,up_j,d,B$s,len=500)
popR <- popsim(N,N1,up_i,up_j,d,B$r,len=500)
popsim <- list(l=popL, s=popS, r=popR)
save(popsim, file=popsim_loc)



