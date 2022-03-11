library(MASS)
#makenoise function
#follows protocol described in SupMat section S3:Noise
#makes list of length 3 containing the three different types of bivariate noise: 
#left tail ATA, right tail ATA, symmetric
#ARGS:
  #M    numeric, length of noise
#OUT:
  #List length 3, each element is a matrix with dimensions 2 x M

makenoise <- function(M){
  
  #Sigma matrix and mu to make standard normal bivariate noise
  Sigma_norm <- matrix(c(1,0,0,1),nrow=2)
  mu_norm <- c(0,0)
  
  #bivariate standard normal
  b <- mvrnorm(M, mu_norm, Sigma_norm) #points 
  a1 <- b[,1]
  a2 <- b[,2]
  
  #half of points will be more than zero, 50%
  A <- a1>0 
  
  #left tail
  b[A,] <- (-abs(a1[A])) #50% b_tilde[,1] = b_tilde[,2] = -|a_1|
  b[!A,] <- abs(b[!A,]) #50% b_tilde[,1] = |a_1|, b_tilde[,2] = |a_2|
  b_l <- b
  
  #right tail - get by flipping left tail noise 
  b[A,] <- -b[A,] 
  b[!A,] <- -b[!A,]
  b_r <- b
  
  #symmetric 
  #with similar correlation as asymmetric noises
  Sigma_sym <- cor(b_l) 
  b_s <- mvrnorm(M, mu_norm, Sigma_sym)
  
  return(list(l=b_l,r=b_r,s=b_s))
}

#makenoise(100)





