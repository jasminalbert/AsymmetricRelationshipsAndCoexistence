
numRes_loc <- "../results_numeric/"
betanoise_loc <- paste0(numRes_loc, "betanoise.RDS")
B <-readRDS(betanoise_loc)

b2p <- function(lb,ub){
  mu <- (lb+ub)/2
  sigma <- ub-lb
  return(c(mu=mu, sigma=sigma))
}

#load litle b's
# B_i = sigma_i * b_i + mu_i

### popsim ###
#ARGS:
#b1     species 1 environmental noise, vector length M
#b2     species 2 environmental noise, vector length M
#N      total abundance (N1 + N2), numeric
#N1     initial abundance of species 1, numeric
#delta  death rate (0-1), numeric
#M      length of noise, numeric
#OUT:
#data.frame 2xM+1 with columns N1 and N2 that are the abundance of sp1 and sp2 over time
  
popsim <- function(N,N1,lb_i,lb_j,up_i,up_j,delta,bnoise,plot=TRUE, len=100) {
  
  paramsi <- b2p(lb_i, up_i)
  paramsj <- b2p(lb_j, up_j)
  
  sigma_i <- paramsi["sigma"]; sigma_j <- paramsj["sigma"]
  mu_i <- paramsi["mu"]; mu_j <- paramsj["mu"]

  
  if (mu_i < (sigma_i/2)){
    stop("error: mu_i is not >= sigma_i/2")
  }
  if (mu_j < (sigma_j/2)){
    stop("error: mu_j is not >= sigma_j/2")
  }
  
  
  #length
  M <- nrow(bnoise)
  b_i <- bnoise[,"i"]
  b_j <- bnoise[,"j"]
  B_i <- sigma_i * b_i + mu_i
  B_j <- sigma_j * b_j + mu_j
  
    
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


popL <- popsim(50,1,0,0,1,1.05,0.4,b$l, len=500)
popR <- popsim(50,1,0,0,1,1.05,0.4,b$r, len=500)
popS <- popsim(50,1,0,0,1,1.05,0.4,b$um, len=500)

#head(popL)
#head(round(popR,3),100)


