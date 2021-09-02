decompose <- function(mudif,sigma,delta,b_tilde,u) {
  
  bi_til<-b_tilde$l[,1] #left ATA_i
  bj_til<-b_tilde$l[,2] #left ATA_j
  bi_um <-b_tilde$s[,1] #symmetric_i
  bj_um <-b_tilde$s[,2] #symmetric_j
  M <- nrow(b_tilde$s) 
  
  #invader i
  #epsilon^0_i
  e_0i <- log(1-delta+delta*exp(mudif))
  
  #epsilon^E_i bar (42)
  #hat
  e_Ei_hat <- mean(log(1-delta+delta*exp(sigma*u + mudif - (sigma^2)/2))) - e_0i
  #SE hat
  e_Ei_se <- sd(log(1-delta+delta*exp(sigma*u + mudif - (sigma^2)/2)))/sqrt(M)
  
  
  #epsilon^C_i bar (45)
  #hat
  e_Ci_hat <- mean(log(1-delta+delta*exp(-sigma*u + mudif + (sigma^2)/2))) - e_0i
  #SE hat
  e_Ci_se <- sd(log(1-delta+delta*exp(-sigma*u + mudif + (sigma^2)/2)))/sqrt(M)
  
  #epsilon^(E#C) bar (51)
  #hat
  e_ECsharpi_hat <- mean(log(1-delta+delta*exp(sigma*sqrt(2)*u + mudif))) - e_Ei_hat - e_Ci_hat - e_0i
  #SE hat
  e_ECsharpi_se <- sd( log(1-delta+delta*exp(sigma*sqrt(2)*u + mudif)) 
                       - log(1-delta+delta*exp(sigma*u + mudif - (sigma^2)/2)) 
                       - log(1-delta+delta*exp(-sigma*u + mudif + (sigma^2)/2)) )/sqrt(M)  
  
  #epsilon^[EC] bar (55)
  #hat
  e_ECi_hat <- mean(log(1-delta+delta*exp(sigma*(bi_um - bj_um)+mudif)))-mean(log(1-delta+delta*exp(sigma*sqrt(2)*u+mudif)))
  #SE hat
  s1 <- sd(log(1-delta+delta*exp(sigma*(bi_um - bj_um)+mudif)))/sqrt(M)
  s2 <- sd(log(1-delta+delta*exp(sigma*sqrt(2)*u+mudif)))/sqrt(M)
  e_ECi_se <- sqrt(s1^2 + s2^2)
  
  #epsilon^[E||C] bar (54)
  #hat
  e_ECpipi_hat <- mean(log(1-delta+delta*exp(sigma*(bi_til-bj_til)+mudif))) - mean(log(1-delta+delta*exp(sigma*(bi_um-bj_um)+mudif)))   
  #SE hat
  s3 <- sd(log(1-delta+delta*exp((bi_til-bj_til)+mudif)))/sqrt(M) #need sigma here?
  e_ECpipi_se <-sqrt(s3^2 + s1^2)
  
  #r_i\i bar (36)
  #hat
  r_i_hat <- mean(log(1-delta+delta*exp(sigma*(bi_til-bj_til)+mudif)))
  #SE hat
  r_i_se <- s3
  
  #resident j 
  #epsilon^0_j
  e_0j <- 0
  
  #epsilon^E_j bar (72)
  #hat
  e_Ej_hat <- mean(log(1-delta+delta*exp(sigma*u - (sigma^2)/2))) - e_0j
  #SE hat
  e_Ej_se <- sd(log(1-delta+delta*exp(sigma*u - (sigma^2)/2)))/sqrt(M)
  
  
  #epsilon^C_j bar (77)
  #hat
  e_Cj_hat <- mean(log(1-delta+delta*exp(-sigma*u + (sigma^2)/2))) - e_0j
  #SE hat
  e_Cj_se <- sd(log(1-delta+delta*exp(-sigma*u + (sigma^2)/2)))/sqrt(M)
  
  #epsilon^(E#C)j bar (85)
  #hat
  e_ECsharpj_hat <- mean(log(1-delta+delta*exp(sigma*sqrt(2)*u))) - e_Ej_hat - e_Cj_hat - e_0j
  #SE hat
  e_ECsharpj_se <- sd( log(1-delta+delta*exp(sigma*sqrt(2)*u)) 
                       - log(1-delta+delta*exp(u - (sigma^2)/2)) #need sigma here?
                       - log(1-delta+delta*exp(-sigma*u + (sigma^2)/2)) )/sqrt(M)  
  
  #epsilon^[EC]j bar (93)
  #hat
  e_ECj_hat <- mean(-log(1-delta+delta*exp(sigma*sqrt(2)*u)))
  #SE hat
  e_ECj_se <- sd(-log(1-delta+delta*exp(sigma*sqrt(2)*u)))
  
  #epsilon^[E||C] bar (91)
  e_ECpipj <- 0
  
  #r_j\i bar 36
  r_j <- 0
  
  #q=1
  #Delta hats:
  D0 <- e_0i  #Deltai_0
  DE <- e_Ei_hat - e_Ej_hat #Deltai_E
  DC <- e_Ci_hat - e_Cj_hat #Deltai_C
  DECsharp <- e_ECsharpi_hat - e_ECsharpj_hat #Deltai_(E#C)
  DEC <- e_ECi_hat - e_ECj_hat #Deltai_[EC]
  DECpip <- e_ECpipi_hat  #Deltai_[E||C]
  Dr <- r_i_hat 
  DrwoATA <- r_i_hat - DECpip
  
  #SE(Delta hats):
  DE_se <- sd( log(1-delta+delta*exp(sigma*u + mudif - (sigma^2)/2)) - log(1-delta+delta*exp(sigma*u - (sigma^2)/2)) )/sqrt(M)
  DC_se <- sd( log(1-delta+delta*exp(-sigma*u + mudif + (sigma^2)/2)) - log(1-delta+delta*exp(-sigma*u + (sigma^2)/2)) )/sqrt(M)
  DECsharp_se <- sd( log(1-delta+delta*exp(sigma*sqrt(2)*u + mudif)) - log(1-delta+delta*exp(sigma*u + mudif - (sigma^2)/2))
                     - log(1-delta+delta*exp(-sigma*u + mudif + (sigma^2)/2)) 
                     - log(1-delta+delta*exp(sigma*sqrt(2)*u)) + log(1-delta+delta*exp(sigma*u - (sigma^2)/2))
                     + log(1-delta+delta*exp(-sigma*u + (sigma^2)/2)) )/sqrt(M)
  s1 <- sd(log(1-delta+delta*exp(sigma*(bi_um - bj_um)+mudif)))/sqrt(M) #does this need a second term w q12?
  s2 <- sd(log(1-delta+delta*exp(sigma*sqrt(2)*u+mudif)) - log(1-delta+delta*exp(sigma*sqrt(2)*u)))/sqrt(M) 
  DEC_se <- sqrt(s1^2 + s2^2)
  DECpip_se <- e_ECpipi_se 
  Dr_se <- r_i_se
  
  
  #q=...
  qij <- 1/((1-delta)*exp(-(mudif))+delta)
  
  #Delta hats:
  DEq <- e_Ei_hat - qij*e_Ej_hat
  DCq <- e_Ci_hat - qij*e_Cj_hat
  DECsharpq <- e_ECsharpi_hat - qij*e_ECsharpj_hat
  DECq <- e_ECi_hat - qij*e_ECj_hat
  
  #SE(Delta hats):
  DEq_se <- sd( log(1-delta+delta*exp(sigma*u + mudif - (sigma^2)/2)) - qij*log(1-delta+delta*exp(sigma*u - (sigma^2)/2)) )/sqrt(M)
  DCq_se <- sd( log(1-delta+delta*exp(-sigma*u + mudif + (sigma^2)/2)) - qij*log(1-delta+delta*exp(-sigma*u + (sigma^2)/2)) )/sqrt(M)
  DECsharpq_se <- sd( log(1-delta+delta*exp(sigma*sqrt(2)*u + mudif)) - log(1-delta+delta*exp(sigma*u + mudif - (sigma^2)/2))
                      - log(1-delta+delta*exp(-sigma*u + mudif + (sigma^2)/2)) 
                      - qij*log(1-delta+delta*exp(sigma*sqrt(2)*u)) + qij*log(1-delta+delta*exp(sigma*u - (sigma^2)/2))
                      + qij*log(1-delta+delta*exp(-sigma*u + (sigma^2)/2)) )/sqrt(M)
  s1 <- sd(log(1-delta+delta*exp(sigma*(bi_um - bj_um)+mudif)))/sqrt(M) #does this need a second term w q12?
  s2 <- sd(log(1-delta+delta*exp(sigma*sqrt(2)*u+mudif)) - qij*log(1-delta+delta*exp(sigma*sqrt(2)*u)))/sqrt(M) 
  DECq_se <- sqrt(s1^2 + s2^2)
  
  #return 8 columns
  #8 columns: ei, ei_se, ej, ej_se, D, D_se, Dq, Dq_se 
  # rows: 0, E, C, (E#C), [EC], [E||C], r
  ei <- c(e_0i, e_Ei_hat, e_Ci_hat, e_ECsharpi_hat, e_ECi_hat, e_ECpipi_hat, r_i_hat, NA)
  ei_se <- c(0, e_Ei_se, e_Ci_se, e_ECsharpi_se, e_ECi_se, e_ECpipi_se, r_i_se, NA)
  ej <- c(e_0j, e_Ej_hat, e_Cj_hat, e_ECsharpj_hat, e_ECj_hat, e_ECpipj, r_j, NA)
  ej_se <- c(0, e_Ej_se, e_Cj_se, e_ECsharpj_se, e_ECj_se, 0, 0, NA)
  D <- c(D0, DE, DC, DECsharp, DEC, DECpip, Dr, DrwoATA)
  D_se <- c(0, DE_se, DC_se, DECsharp_se, DEC_se, DECpip_se, Dr_se, NA)
  Dq <- c(D0, DEq, DCq, DECsharpq, DECq, DECpip, Dr, DrwoATA)
  Dq_se <- c(0, DEq_se, DCq_se, DECsharpq_se, DECq_se, DECpip_se, Dr_se, NA)
  
  res <- data.frame(ei=ei, ei_se=ei_se, ej=ej, ej_se=ej_se, D=D, D_se=D_se, Dq=Dq, Dq_se=Dq_se)
  row.names(res) <- c("0","E","C","(E#C)","[EC]","[E||C]","r", "rwoATA")
  
  return(res)
}

#dec <- decompose(-0.1,1.6,0.5,b_tilde,u_tilde)
#round(dec,3)




