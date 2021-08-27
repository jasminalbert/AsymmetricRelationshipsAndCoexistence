#pop sim function
popsim <- function(b1,b2,N,N1,delta,M){
  B1 <- exp(b1)
  B2 <- exp(b2)
  N2 <- N-N1
  
  for (t in 1:M){
    tot_new_juvs <- (B1[t]*N1[t]) + (B2[t]*N2[t])
    
    N1[t+1] <- (1-delta)*N1[t] + delta*N*((B1[t]*N1[t])/tot_new_juvs)	
    
    N2[t+1] <- (1-delta)*N2[t] + delta*N*((B2[t]*N2[t])/tot_new_juvs)
  }
  return(data.frame(N1=N1, N2=N2))
}


