decompose <- function(mudif,sigma,delta,b_tilde,u) {
  
  bi_til<-b_tilde$l[,1]
  bj_til<-b_tilde$l[,2]
  bi_um <-b_tilde$s[,1]
  bj_um <-b_tilde$s[,2]
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


#decomposition function #need to check
decomp <- function(sigma, mu, delta, M){
  
  #exponentiate noise to get B
  load("./noise_etc.RData")
  rv <- transform(b_tilde, u_tilde, rho, sigma, mu, b_s=TRUE) #makes random vars defined by parameters
  Bl <- cbind(exp(rv$b_l1), exp(rv$b_l2))
  Br <- cbind(exp(rv$b_r1), exp(rv$b_r2))
  Bs <- cbind(exp(rv$b_s1), exp(rv$b_s2))
  
  #make B sharp
  bsharp_tilde <- makenoise(M)
  rho <- cor(bsharp_tilde$l)[1,2]
  rvsharp <- transform(bsharp_tilde, u_tilde, rho, sigma, mu, b_s=TRUE)
  Bl2sharp <- exp(rvsharp$b_l2)
  Br2sharp <- exp(rvsharp$b_r2)
  Bs2sharp <- exp(rvsharp$b_s2)
  
  #epsilon_0
  null_l <- log(1-delta+delta*(mean(Bl[,1])/mean(Bl[,2])))
  null_r <- log(1-delta+delta*(mean(Bl[,1])/mean(Bl[,2])))
  null_s <- log(1-delta+delta*(mean(Bl[,1])/mean(Bl[,2])))
  
  #epsilon_C
  epC1l <- mean(log(1-delta+delta*(mean(Bl[,1])/Bl[,2])))- null_l
  epC1r <- mean(log(1-delta+delta*(mean(Br[,1])/Br[,2])))- null_r
  epC1s <- mean(log(1-delta+delta*(mean(Bs[,1])/Bs[,2])))- null_s
  
  epC2l <- mean(log(1-delta+delta*(mean(Bl[,2])/Bl[,2])))
  epC2r <- mean(log(1-delta+delta*(mean(Br[,2])/Br[,2])))
  epC2s <- mean(log(1-delta+delta*(mean(Bs[,2])/Bs[,2])))
  
  #epsilon_E
  epE1l <- mean(log(1-delta+delta*(Bl[,1]/mean(Bl[,2]))))- null_l
  epE1r <- mean(log(1-delta+delta*(Br[,1]/mean(Br[,2]))))- null_r
  epE1s <- mean(log(1-delta+delta*(Bs[,1]/mean(Bs[,2]))))- null_s
  
  epE2l <- mean(log(1-delta+delta*(Bl[,2]/mean(Bl[,2]))))
  epE2r <- mean(log(1-delta+delta*(Br[,2]/mean(Br[,2]))))
  epE2s <- mean(log(1-delta+delta*(Bs[,2]/mean(Bs[,2]))))
  
  #epsilon_EC
  epEC1l <- mean(log(1-delta+delta*(Bl[,1]/Bl[,2])))- (null_l+ epC1l+ epE1l)
  epEC1r <- mean(log(1-delta+delta*(Br[,1]/Br[,2])))- (null_r+ epC1r+ epE1r)
  epEC1s <- mean(log(1-delta+delta*(Bs[,1]/Bs[,2])))- (null_s+ epC1s+ epE1s)
  
  epEC2l <- mean(log(1-delta+delta*(Bl[,2]/Bl[,2])))- (epC2l+ epE2l)
  epEC2r <- mean(log(1-delta+delta*(Br[,2]/Br[,2])))- (epC2r+ epE2r)
  epEC2s <- mean(log(1-delta+delta*(Bs[,2]/Bs[,2])))- (epC2s+ epE2s)
  
  #epsilon_(E#C)
  epECsharp1l <- mean(log(1-delta+delta*(Bl[,1]/Bl2sharp)))- (null_l+ epC1l+ epE1l)
  epECsharp1r <- mean(log(1-delta+delta*(Br[,1]/Br2sharp)))- (null_r+ epC1r+ epE1r)
  epECsharp1s <- mean(log(1-delta+delta*(Bs[,1]/Bs2sharp)))- (null_s+ epC1s+ epE1s)
  
  epECsharp2l <- mean(log(1-delta+delta*(Bl[,2]/Bl2sharp)))- (epC2l+ epE2l)
  epECsharp2r <- mean(log(1-delta+delta*(Br[,2]/Br2sharp)))- (epC2r+ epE2r)
  epECsharp2s <- mean(log(1-delta+delta*(Bs[,2]/Bs2sharp)))- (epC2s+ epE2s)
  
  #epsilon_(EC)
  epEC_1l <- epEC1l - epECsharp1l
  epEC_1r <- epEC1r - epECsharp1r
  epEC_1s <- epEC1s - epECsharp1s
  
  epEC_2l <- epEC2l - epECsharp2l
  epEC_2r <- epEC2r - epECsharp2r
  epEC_2s <- epEC2s - epECsharp2s
  
  #Delta_0
  Dnull_l <- null_l
  Dnull_r <- null_r
  Dnull_s <- null_s
  
  #Delta_E
  DE_l <- epE1l-epE2l
  DE_r <- epE1r-epE2r
  DE_s <- epE1s-epE2s
  
  #Delta_C
  DC_l <- epC1l-epC2l
  DC_r <- epC1r-epC2r
  DC_s <- epC1s-epC2s
  
  #Delta_(EC)
  DEC_l <- epEC_1l-epEC_2l
  DEC_r <- epEC_1r-epEC_2r
  DEC_s <- epEC_1s-epEC_2s
  
  #Delta_(E#C)
  DECsharpl <- epECsharp1l-epECsharp2l
  DECsharpr <- epECsharp1r-epECsharp2r
  DECsharpS <- epECsharp1s-epECsharp2s
  
  #format
  Dnull <- c(Dnull_l, Dnull_r, Dnull_s)
  
  DE <- c(DE_l, DE_r, DE_s)
  
  DC <- c(DC_l, DC_r, DC_s)
  
  DEC_ <- c(DEC_l, DEC_r, DEC_s)
  
  DECsharp <- c(DECsharpl, DECsharpr, DECsharpS)
  
  res <- data.frame(Dnull=Dnull, DE=DE, DC=DC, DEC_=DEC_, DECsharp=DECsharp,
                    row.names=c("left","right","sym"))
  return(res)
  
}

#decomp(1.6, c(0.3,0.4), 0.4, M)


