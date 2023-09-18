#This script contains the function:
#	decompose
	#	for coexistence decomposition of lottery model

##libraries used (invoked with ::): stats

### decompose ###
# compute estimates for epsilons and Deltas of coexistence decomposition
# using two alternatives for qij
# see sup mat section 6: Efficient Computation
#ARGS:
  #mudif    mean difference; mu1-mu2
  #sigma    standard deviation of noise
  #delta    death rate
  #b_tilde  list length three of three types of bivariate noise with standard normal marginals
  #u        standard normal univariate noise
#OUT:
  #data.frame (7x8) containing all epsilon_i's epsilon_j's and 
  #deltas computed using the qij alternatives
decompose <- function(mudif,sigma,delta,b_tilde,u) {
  #notation follows from paper/sup mat
  
  #define different noise 
  bi_til<-b_tilde$l[,1] #left ATA for species i
  bj_til<-b_tilde$l[,2] #left ATA for species j
  bi_um <-b_tilde$s[,1] #symmetric for species i
  bj_um <-b_tilde$s[,2] #symmetric for species j
  
  #length
  M <- nrow(b_tilde$s) 
  
  #defining epsilons for rare species i (second list in S6)
  #1. epsilon^0_i (eq 39)
  e_0i <- log(1-delta+delta*exp(mudif))
  
  #2. epsilon^E_i bar (eq 42)
  #estimate
  e_Ei_hat <- mean(log(1-delta+delta*exp(sigma*u + mudif - (sigma^2)/2))) - e_0i
  #standard error
  e_Ei_se <- stats::sd(log(1-delta+delta*exp(sigma*u + mudif - (sigma^2)/2)))/sqrt(M)
  
  #3. epsilon^C_i bar (eq 45)
  #estimate
  e_Ci_hat <- mean(log(1-delta+delta*exp(-sigma*u + mudif + (sigma^2)/2))) - e_0i
  #standard error
  e_Ci_se <- stats::sd(log(1-delta+delta*exp(-sigma*u + mudif + (sigma^2)/2)))/sqrt(M)
  
  #4. epsilon^(E#C) bar (eq 51)
  #estimate
  e_ECsharpi_hat <- mean(log(1-delta+delta*exp(sigma*sqrt(2)*u + mudif))) - e_Ei_hat - e_Ci_hat - e_0i
  #standard error
  e_ECsharpi_se <- stats::sd( log(1-delta+delta*exp(sigma*sqrt(2)*u + mudif)) 
                       - log(1-delta+delta*exp(sigma*u + mudif - (sigma^2)/2)) 
                       - log(1-delta+delta*exp(-sigma*u + mudif + (sigma^2)/2)) )/sqrt(M)  
  
  #5. epsilon^[E||C] bar (eq 55)
  #estimate
  e_ECpipi_hat <- mean(log(1-delta+delta*exp(sigma*(bi_um - bj_um)+mudif)))-mean(log(1-delta+delta*exp(sigma*sqrt(2)*u+mudif)))
  #standard error
  s1 <- stats::sd(log(1-delta+delta*exp(sigma*(bi_um - bj_um)+mudif)))/sqrt(M)
  s2 <- stats::sd(log(1-delta+delta*exp(sigma*sqrt(2)*u+mudif)))/sqrt(M)
  e_ECpipi_se <- sqrt(s1^2 + s2^2)
  
  #5. epsilon^[EC] bar (eq 54)
  #estimate
  e_ECi_hat <- mean(log(1-delta+delta*exp(sigma*(bi_til-bj_til)+mudif))) - mean(log(1-delta+delta*exp(sigma*(bi_um-bj_um)+mudif)))   
  #standard error
  s3 <- stats::sd(log(1-delta+delta*exp(sigma*(bi_til-bj_til)+mudif)))/sqrt(M) #need sigma here?
  e_ECi_se <-sqrt(s3^2 + s1^2)
  
  
  #7. r_i\i bar (eq 36)
  #estimate 
  r_i_hat <- mean(log(1-delta+delta*exp(sigma*(bi_til-bj_til)+mudif)))
  #standard error 
  r_i_se <- s3
  
  #defining epsilions for common species j (third list in S6) 
  #1. epsilon^0_j (eq 67)
  e_0j <- 0
  
  #2. epsilon^E_j bar (eq 72)
  #estimate
  e_Ej_hat <- mean(log(1-delta+delta*exp(sigma*u - (sigma^2)/2))) - e_0j
  #standard error
  e_Ej_se <- stats::sd(log(1-delta+delta*exp(sigma*u - (sigma^2)/2)))/sqrt(M)
  
  #3. epsilon^C_j bar (eq 77)
  #estimate
  e_Cj_hat <- mean(log(1-delta+delta*exp(-sigma*u + (sigma^2)/2))) - e_0j
  #standard error
  e_Cj_se <- stats::sd(log(1-delta+delta*exp(-sigma*u + (sigma^2)/2)))/sqrt(M)
  
  #4. epsilon^(E#C)j bar (eq 85)
  #estimate
  e_ECsharpj_hat <- mean(log(1-delta+delta*exp(sigma*sqrt(2)*u))) - e_Ej_hat - e_Cj_hat - e_0j
  #standard error
  e_ECsharpj_se <- stats::sd( log(1-delta+delta*exp(sigma*sqrt(2)*u)) 
                       - log(1-delta+delta*exp(u - (sigma^2)/2)) #need sigma here?
                       - log(1-delta+delta*exp(-sigma*u + (sigma^2)/2)) )/sqrt(M)  
  
  #5. epsilon^[E||C]j bar (eq 93)
  #estimate
  e_ECpipj_hat <- mean(-log(1-delta+delta*exp(sigma*sqrt(2)*u)))
  #standard error
  e_ECpipj_se <- stats::sd(-log(1-delta+delta*exp(sigma*sqrt(2)*u)))/sqrt(M)
  
  #6. epsilon^[EC] bar (eq 91)
  e_ECj <- 0
  
  #7. r_j\i bar 
  r_j <- 0
  
  #Deltas (fourth list in S6)
  #q=1
  #Delta estimates:
  D0 <- e_0i                                  #Deltai_0
  DE <- e_Ei_hat - e_Ej_hat                   #Deltai_E
  DC <- e_Ci_hat - e_Cj_hat                   #Deltai_C
  DECsharp <- e_ECsharpi_hat - e_ECsharpj_hat #Deltai_(E#C)
  DEC <- e_ECi_hat - e_ECj                    #Deltai_[EC]
  DECpip <- e_ECpipi_hat - e_ECpipj_hat       #Deltai_[E||C]
  Dr <- r_i_hat 
  DrwoATA <- r_i_hat - DECpip
  
  #standard error of Delta estimates:
  DE_se <- stats::sd( log(1-delta+delta*exp(sigma*u + mudif - (sigma^2)/2)) - log(1-delta+delta*exp(sigma*u - (sigma^2)/2)) )/sqrt(M)
  DC_se <- stats::sd( log(1-delta+delta*exp(-sigma*u + mudif + (sigma^2)/2)) - log(1-delta+delta*exp(-sigma*u + (sigma^2)/2)) )/sqrt(M)
  DECsharp_se <- stats::sd( log(1-delta+delta*exp(sigma*sqrt(2)*u + mudif)) - log(1-delta+delta*exp(sigma*u + mudif - (sigma^2)/2))
                     - log(1-delta+delta*exp(-sigma*u + mudif + (sigma^2)/2)) 
                     - log(1-delta+delta*exp(sigma*sqrt(2)*u)) + log(1-delta+delta*exp(sigma*u - (sigma^2)/2))
                     + log(1-delta+delta*exp(-sigma*u + (sigma^2)/2)) )/sqrt(M)
  s1 <- stats::sd(log(1-delta+delta*exp(sigma*(bi_um - bj_um)+mudif)))/sqrt(M) 
  s2 <- stats::sd(log(1-delta+delta*exp(sigma*sqrt(2)*u+mudif)) - log(1-delta+delta*exp(sigma*sqrt(2)*u)))/sqrt(M) 
  DEC_se <- sqrt(s1^2 + s2^2)
  DECpip_se <- e_ECpipi_se 
  Dr_se <- r_i_se
  
  #q=...
  qij <- 1/((1-delta)*exp(-(mudif))+delta)
  
  #Delta estimates:
  DEq <- e_Ei_hat - qij*e_Ej_hat
  DCq <- e_Ci_hat - qij*e_Cj_hat
  DECsharpq <- e_ECsharpi_hat - qij*e_ECsharpj_hat
  DECq <- e_ECi_hat - qij*e_ECj
  DECpipq <- e_ECpipj_hat-qij*e_ECpipj_hat
  
  #standard error of Delta estimates:
  DEq_se <- stats::sd( log(1-delta+delta*exp(sigma*u + mudif - (sigma^2)/2)) - qij*log(1-delta+delta*exp(sigma*u - (sigma^2)/2)) )/sqrt(M)
  DCq_se <- stats::sd( log(1-delta+delta*exp(-sigma*u + mudif + (sigma^2)/2)) - qij*log(1-delta+delta*exp(-sigma*u + (sigma^2)/2)) )/sqrt(M)
  DECsharpq_se <- stats::sd( log(1-delta+delta*exp(sigma*sqrt(2)*u + mudif)) - log(1-delta+delta*exp(sigma*u + mudif - (sigma^2)/2))
                      - log(1-delta+delta*exp(-sigma*u + mudif + (sigma^2)/2)) 
                      - qij*log(1-delta+delta*exp(sigma*sqrt(2)*u)) + qij*log(1-delta+delta*exp(sigma*u - (sigma^2)/2))
                      + qij*log(1-delta+delta*exp(-sigma*u + (sigma^2)/2)) )/sqrt(M)
  s1 <- stats::sd(log(1-delta+delta*exp(sigma*(bi_um - bj_um)+mudif)))/sqrt(M) 
  s2 <- stats::sd(log(1-delta+delta*exp(sigma*sqrt(2)*u+mudif)) - qij*log(1-delta+delta*exp(sigma*sqrt(2)*u)))/sqrt(M) 
  DECq_se <- sqrt(s1^2 + s2^2)
  
  #return:
  #8 columns: ei, ei_se, ej, ej_se, D, D_se, Dq, Dq_se 
    # epsilon estimate and standard error for each species and 
    # Delta estimate and standard error for each qij alternative
  #7 rows: 0, E, C, (E#C), [EC], [E||C], r, rwoATA 
    # coexistence mechanims and gwr and gwr without ATA contribution 
  ei <- c(e_0i, e_Ei_hat, e_Ci_hat, e_ECsharpi_hat, e_ECi_hat, e_ECpipi_hat, r_i_hat, NA)
  ei_se <- c(0, e_Ei_se, e_Ci_se, e_ECsharpi_se, e_ECi_se, e_ECpipi_se, r_i_se, NA)
  ej <- c(e_0j, e_Ej_hat, e_Cj_hat, e_ECsharpj_hat, e_ECj, e_ECpipj_hat, r_j, NA)
  ej_se <- c(0, e_Ej_se, e_Cj_se, e_ECsharpj_se, 0, e_ECpipj_se, 0, NA)
  D <- c(D0, DE, DC, DECsharp, DEC, DECpip, Dr, DrwoATA)
  D_se <- c(0, DE_se, DC_se, DECsharp_se, DEC_se, DECpip_se, Dr_se, NA)
  Dq <- c(D0, DEq, DCq, DECsharpq, DECq, DECpipq, Dr, DrwoATA)
  Dq_se <- c(0, DEq_se, DCq_se, DECsharpq_se, DECq_se, DECpip_se, Dr_se, NA)
  
  res <- data.frame(ei=ei, ei_se=ei_se, ej=ej, ej_se=ej_se, D=D, D_se=D_se, Dq=Dq, Dq_se=Dq_se)
  row.names(res) <- c("0","E","C","(E#C)","[EC]","[E||C]","r", "rwoATA")
  
  return(res)
}

#dec <- decompose(-0.1,1.6,0.5,b_tilde,u_tilde)
#round(dec,3)




