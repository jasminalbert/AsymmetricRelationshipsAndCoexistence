# transform noise by sigma and mu to make all the noise needed to compute DeltaI and related quantities

#b_tilde - 2 sets bivariate noise 
#u_tilde - univariate noise
#rho - common correlation of bivariate noise sets
#sigma - common sd of bivariate noise sets
#mudif - mu1-mu2

transform <- function(b_tilde, u_tilde, rho, sigma, mudif, b_s=FALSE){
  b_l1 <- sigma*b_tilde$l[,1] 
  b_l2 <- sigma*b_tilde$l[,2] - mudif
  
  b_r1 <- sigma*b_tilde$r[,1] 
  b_r2 <- sigma*b_tilde$r[,2] - mudif
  
  u_s <- sqrt(2*sigma^2 - 2*rho*sigma^2)*u_tilde + mudif 
  u_1 <- sqrt(2*sigma^2)*u_tilde + mudif
  u_2 <- sqrt(2*sigma^2)*u_tilde  
  
  if (b_s==FALSE){
    return(list(b_l1=b_l1, b_l2=b_l2, b_r1=b_r1, b_r2=b_r2, 
                u_s=u_s, u_1=u_1, u_2=u_2))
  } else{
    b_s1 <- sigma*b_tilde$s[,1] 
    b_s2 <- sigma*b_tilde$s[,2] - mudif
    
    return(list(b_l1=b_l1, b_l2=b_l2, b_r1=b_r1, b_r2=b_r2, 
                b_s1=b_s1, b_s2=b_s2, u_s=u_s, u_1=u_1, u_2=u_2))
  }
  
  
}


