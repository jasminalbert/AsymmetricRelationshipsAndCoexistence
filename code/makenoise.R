#This script makes sets of noise to be used in the lottery model analysis and saves their length, M, and correlation, rho
  #b: list of 4 different types of bivariate noise: left tail ATA, righ tail ATA, bivariate normal with cov (1,0,0,1) and symmetric
  #u: univariate standard normal 

##libraries used (invoked with ::): MASS, stats

### location to save results ###
numeric_results_loc <- "../results_numeric"

if(dir.exists(numeric_results_loc)==FALSE){
  dir.create(numeric_results_loc)
}
noise_loc <- paste0(numeric_results_loc, "/noise.RData")
M_loc <- paste0(numeric_results_loc, "/M.RDS")
rho_loc <- paste0(numeric_results_loc, "/rho.RDS")

### makenoise function ###
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
  b <- MASS::mvrnorm(M, mu_norm, Sigma_norm) #points
  colnames(b) <- c("i","j")
  b_ti <- b
  
  #half of points will be more than zero, 50%
  a1 <- b[,1]; a2 <- b[,2]
  A <- a1>0 
  
  #left tail
  b[A,] <- (-abs(a1[A])) #50% b[,1] = b[,2] = -|a_1|
  b[!A,] <- abs(b[!A,]) #50% b[,1] = |a_1|, b[,2] = |a_2|
  b_l <- b
  
  #right tail - get by flipping left tail noise 
  b[A,] <- -b[A,] 
  b[!A,] <- -b[!A,]
  b_r <- b
  
  #symmetric 
  #with similar correlation as asymmetric noises
  Sigma_sym <- stats::cor(b_l) 
  b_um <- MASS::mvrnorm(M, mu_norm, Sigma_sym)
  colnames(b_um) <- colnames(b)
  
  return(list(l=b_l,r=b_r,s=b_um,til=b_ti))
}

### make noise ###
set.seed(1360)

#noise length
M <- 1000000
#bivariate noise list
b <- makenoise(M)
#univariate standard normal
u <- stats::rnorm(M)

#get correlation
rho <- stats::cor(b_tilde$l)[1,2]
cat("M=",M, " rho=", rho)

### save results ###
save(b,u, file = noise_loc)
saveRDS(M, file=M_loc)
saveRDS(rho, file=rho_loc)

