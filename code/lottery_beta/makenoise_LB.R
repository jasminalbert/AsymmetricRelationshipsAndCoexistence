#This script makes sets of noise to be used in the lottery model analysis and saves their length, M, and correlation, rho
  #b_tilde: list of 3 different types of bivariate noise: left tail ATA, righ tail ATA, and symmetric
  #u_tilde: univariate standard normal 

##libraries used (invoked with ::): MASS, stats

### location to save results ###
numeric_results_loc <- "../results_numeric"

if(dir.exists(numeric_results_loc)==FALSE){
  dir.create(numeric_results_loc)
}
noise_loc <- paste(numeric_results_loc, "/noise.RData", sep = "")
M_loc <- paste(numeric_results_loc, "/M.RDS", sep = "")
rho_loc <- paste(numeric_results_loc, "/rho.RDS", sep = "")

### makenoise function ###
#follows protocol described in SupMat section 
#makes list of length 3 containing the three different types of bivariate noise: 
#left tail ATA, right tail ATA, symmetric, standard
#v_l, v_r, v_um (umlaut), v_ti (tilde)
#ARGS:
  #M    numeric, length of noise
#OUT:
  #List length 4, each element is a matrix with dimensions 2 x M
makenoise <- function(M){
  
  #Sigma matrix and mu to make standard normal bivariate noise
  Sigma_norm <- matrix(c(1,0,0,1),nrow=2)
  mu_norm <- c(0,0)
  
  #bivariate standard normal
  v <- MASS::mvrnorm(M, mu_norm, Sigma_norm) #points 
  colnames(v) <- c("i","j")
  v_ti <- v
  
  #half of points will be more than zero, 50%
  a1 <- v[,1]; a2 <- v[,2]
  
  A <- a1>0 
  
  #left tail
  v[A,] <- (-abs(a1[A])) #50% v[,1] = v[,2] = -|a_1|
  v[!A,] <- abs(v[!A,]) #50% v[,1] = |a_1|, v[,2] = |a_2|
  v_l <- v
  
  #right tail - get by flipping left tail noise 
  v[A,] <- -v[A,] 
  v[!A,] <- -v[!A,]
  v_r <- v
  
  #symmetric 
  #with similar correlation as asymmetric noises
  Sigma_sym <- stats::cor(v_l) 
  v_um <- MASS::mvrnorm(M, mu_norm, Sigma_sym)
  colnames(v_um) <- colnames(v)
  
  return(list(l=v_l,r=v_r,til=v_ti,um=v_um))
}

### make noise ###
set.seed(1360)

#noise length
M <- 1000000
#bivariate noise list
v <- makenoise(M)
#univariate standard normal
#u_tilde <- stats::rnorm(M) #still need??

#get correlation
rho <- stats::cor(v$l)[1,2]
cat("M=",M, " rho=", rho)

### save results ###
save(v, file = noise_loc)
saveRDS(M, file=M_loc)
saveRDS(rho, file=rho_loc)

