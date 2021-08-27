library(MASS)
# make into function

# generate noise

makenoise <- function(M){
  
  #sigma and mu to make standard normal bivariate noise
  Sigma_norm <- matrix(c(1,0,0,1),nrow=2)
  mu_norm <- c(0,0)
  
  #standard normal
  b <- mvrnorm(M, mu_norm, Sigma_norm) #points (a1,a2)
  a1 <- b[,1]
  a2 <- b[,2]
  
  #A <- sample(ind, M*0.5) #50% chance
  A <- a1>0
  
  #left tail
  b[A,] <- rep(-abs(a1[A]))
  b[!A,] <- abs(b[!A,])
  b_l <- b
  #plot(b_l[,1],b_l[,2])
  
  #right tail
  b[A,] <- -b[A,]
  b[!A,] <- -b[!A,]
  b_r <- b
  #plot(b_r[,1],b_r[,2])
  
  #symmetric 
  Sigma_sym <- cor(b_l) 
  b_s <- mvrnorm(M, mu_norm, Sigma_sym)
  #plot(b_s[,1],b_s[,2])
  
  return(list(l=b_l,r=b_r,s=b_s))
}

#makenoise(100)





