#numeric_results_loc <- "../results_numeric"
#betanoise_loc <- paste0(numeric_results_loc, "/betanoise.RDS")
#B<-readRDS(betanoise_loc)

decomp <- function(etai, etaj, delta, Blist, dir="LEFT", qij=FALSE){
	
	r <- function(E,C){log(1-delta+E/C)}
	rbar <- function(E,C){mean(r(E,C))}
	terms <- c("0","E", "C", "EsC", "EpsC", "ATA", "r")
	
	#storage for epsilons and standard error
	ei <- matrix(nrow=length(terms), ncol=1, dimnames=list(terms,"epsilon_i"))
	sei <- matrix(ei, dimnames=list(terms,"stanErr(e_i)"))
	ej <- matrix(ei, dimnames=list(terms,"epsilon_j"))
	sej <- matrix(sei, dimnames=list(terms,"stanErr(e_j)"))
	seD <- matrix(sei, dimnames=list(terms,"stanErr(D_i)"))
	
	#B variables
	if (dir=="RIGHT") {
		Bi <- Blist$r[,1]; Bj <- Blist$r[,2]
	} else{
		Bi <- Blist$l[,1]; Bj <- Blist$l[,2]
	}
	Bis <- Blist$til[,1]; Bjs <- Blist$til[,2]
	Bips <- Blist$s[,1]; Bjps <- Blist$s[,2]

	M <- length(Bi)
	
	## epsilon_i's (invader) ##
	ei["0",] <- r(etai,etaj/delta)
	sei["0",] <- NA
	
	ei["E",] <- rbar(etai*Bis, etaj/(delta*2))-ei["0",]
	sei["E",] <- sd(r(etai*Bis, etaj/(delta*2)))/sqrt(M)
	
	ei["C",] <- rbar(etai/2, (etaj*Bjs)/delta)-ei["0",]
	sei["C",] <- sd(r(etai/2, (etaj*Bjs)/delta))/sqrt(M)
	
	ei["EsC",] <- rbar(etai*Bis, (etaj/delta)*Bjs) - ei["0",] - ei["C",] - ei["E",]
	sei["EsC",] <- sd(r(etai*Bis, (etaj/delta)*Bjs)-r(etai/2,(etaj/delta)*Bjs) - r(etai*Bis,etaj/(2*delta)) )/sqrt(M)
	
	ei["EpsC",] <- rbar(etai*Bips, (etaj/delta)*Bjps) - rbar(etai*Bis, (etaj/delta)*Bjs)
	s1 <- sd(r(etai*Bips, (etaj/delta)*Bjps))/sqrt(M)
	s2 <- sd(r(etai*Bis, (etaj/delta)*Bjs))/sqrt(M)
	sei["EpsC",] <- sqrt(s1^2 + s2^2)

	ei["ATA",] <- rbar(etai*Bi, etaj*Bj/delta) - rbar(etai*Bips, etaj*Bjps/delta)
	s3 <- sd(r(etai*Bi, etaj*Bj/delta))/sqrt(M)
	sei["ATA",] <- sqrt(s1^2 + s3^2)

	ei["r",] <- rbar(etai*Bi, etaj*Bj/delta)
	sei["r",] <- s3

	## epsilon_j's (resident) ##
	ej["0",] <- 0; sej["0",] <- NA

	ej["E",] <- rbar(etaj*Bj, etaj/(2*delta))
	sej["E",] <- sd(r(etaj*Bj, etaj/(2*delta)))/sqrt(M)

	ej["C",] <- rbar(etaj/2, etaj*Bjs/delta)
	sej["C",] <- sd(r(etaj/2, etaj*Bjs/delta))/sqrt(M)

	ej["EsC",] <- rbar(etaj*Bj, etaj*Bjs/delta) - ej["C",]-ej["E",]
	sej["EsC",] <- sd(r(etaj*Bj, etaj*Bjs/delta)- r(etaj*Bj, etaj/(2*delta)) - r(etaj/2, etaj*Bjs/delta))/sqrt(M)

	ej["EpsC",] <- -rbar(etaj*Bj, etaj*Bjs/delta)
	sej["EpsC",] <- sd(r(etaj*Bj, etaj*Bjs/delta))/sqrt(M)

	ej["ATA",] <- 0
	sej["ATA",] <- NA

	ej["r",] <- 0
	sej["r",] <- NA

	## Delta_i's (invader epsilons - resident epsilons) ##
	
	if (qij==TRUE){
		qij <- (etai/etaj) / ((1-delta) + (etai/etaj)*delta)
	} else if (qij==FALSE) {
		qij <- 1
	} else {stop("qij arguement must be logical")}
		
	Di <- matrix(ei-qij*ej, dimnames=list(terms,"Delta_i"))
	
	seD["0",] <- NA
	
	s1 <- sd(r(etai*Bis, etaj/(2*delta)))/sqrt(M)
	s2 <- qij*sd(r(etaj*Bj, etaj/(2*delta)))/sqrt(M)
	seD["E",] <- sqrt(s1^2+s2^2)

	seD["C",] <- sd(r(etai/2, etaj*Bjs/delta) - qij*r(etaj/2, etaj*Bjs/delta))/sqrt(M)

	seD["EsC",] <- sd(r(etai*Bis, etaj*Bjs/delta) - r(etai*Bis, etaj/(2*delta)) - r(etai/2, etaj*Bjs/delta) - qij*r(etaj*Bj, etaj*Bjs/delta) + qij*r(etaj*Bj, etaj/(2*delta)) + qij*r(etaj/2, etaj*Bjs/delta))/sqrt(M)

	seD["EpsC",] <- sd(r(etai*Bips, etaj*Bjps/delta)-r(etai*Bis, etaj*Bjs/delta)+ qij*r(etaj*Bj, etaj*Bjs/delta))/sqrt(M)

	seD["ATA",] <- sei["ATA",]
	seD["r",] <- sei["r",]

	return(data.frame(cbind(ei,ej,Di,sei,sej,seD)))
}

b2p <- function(lb,ub){
	mu <- (lb+ub)/2
	sigma <- ub-lb
	return(c(mu=mu, sigma=sigma))
}

#load litle b's
# B_i = sigma_i * b_i + mu_i

decomposeB <- function(lb_i,lb_j,ub_i,ub_j, delta,blist,dir="LEFT",plot=FALSE) {
	
	paramsi <- b2p(lb_i, ub_i)
	paramsj <- b2p(lb_j, ub_j)
	
	sigma_i <- paramsi["sigma"]; sigma_j <- paramsj["sigma"]
	mu_i <- paramsi["mu"]; mu_j <- paramsj["mu"]
	
	  #convert upper bound to sigma and mu
  #sigma_i <- up_i 
  #sigma_j <- up_j
  #mu_i <- 
  #mu_j <- 
	
  #notation follows from paper/sup mat
  
  if (mu_i < (sigma_i/2)){
  	stop("error: mu_i is not >= sigma_i/2")
  }
    if (mu_j < (sigma_j/2)){
  	stop("error: mu_j is not >= sigma_j/2")
  }
  
  b <- data.frame(blist$l) #left ATA biv dis
  b_til <- data.frame(blist$til) #standard normal biv
  b_um <- data.frame(blist$um) #normal biv w same rho as ATA noise
  if (dir=="RIGHT"){
  	b <- data.frame(blist$r) #right ATA biv dis
  }
  
  #length
  M <- nrow(b$i) 
  
  #defining epsilons for rare species i (second list in S6)
  #1. epsilon^0_i (eq 39)
  e0i <- log(1-delta+delta*(mu_i/mu_j))
  
  #2. epsilon^E_i bar (eq 42)
  #estimate
  eEi_hat <- mean(log(1-delta+delta*((sigma_i*b$i + mu_i)/mu_j))) - e0i
  #standard error
  #eEi_se <- sd()/sqrt(M)
  
  #3. epsilon^C_i bar (eq 45)
  #estimate
  eCi_hat <- mean(log(1-delta+delta*(mu_i/(sigma_j*b$j + mu_j)))) - e0i
  #standard error
  #eCi_se <- sd()/sqrt(M)
  
  #4. epsilon^(E#C) bar (eq 51)
  #estimate
  eECshi_hat <- mean(log(1-delta+delta* ((sigma_i*b_til$i+mu_i)/(sigma_j*b_til$j+mu_j))))- eEi_hat-eCi_hat - e0i
  #standard error
  #eECshi_se <- sd()/sqrt(M)  
  
  #5. epsilon^[EC] bar (eq 55)
  #estimate
  eECi_hat <- mean(log(1-delta+delta* ((sigma_i*b$i+mu_i)/(sigma_j*b$j+mu_j)))) - mean(log(1-delta+delta* ((sigma_i*b_um$i+mu_i)/(sigma_j*b_um$j+mu_j)))) 
  #standard error
  #s1 <- sd()/sqrt(M)
  #s2 <- sd()/sqrt(M)
  #eECi_se <- sqrt(s1^2 + s2^2)
  
  #6. epsilon^[E||C] bar (eq 54)
  #estimate
  eECpsi_hat <- mean(log(1-delta+delta* ((sigma_i*b_um$i+mu_i)/(sigma_j*b_um$j+mu_j)))) - mean(log(1-delta+delta* ((sigma_i*b_til$i+mu_i)/(sigma_j*b_til$j+mu_j))))  
  #standard error
  #s3 <- sd()/sqrt(M) 
  #eECpipi_se <-sqrt(s3^2 + s1^2)
  
  #7. r_i\i bar (eq 36)
  #estimate 
  rii_hat <- mean(log(1-delta+delta*((sigma_i*b$i+mu_i)/(sigma_j*b$j+mu_j))))
  #standard error 
  #rii_se
  
  #defining epsilions for common species j (third list in S6) 
  #1. epsilon^0_j (eq 67)
  e0j <- 0
  
  #2. epsilon^E_j bar (eq 72)
  #estimate
  eEj_hat <- mean(log(1-delta+delta*((sigma_j*b$j + mu_j)/mu_j))) - e0j
  #standard error
  #eEj_se <- sd()/sqrt(M)
  
  #3. epsilon^C_j bar (eq 77)
  #estimate
  eCj_hat <- mean(log(1-delta+delta*(mu_j/(sigma_j*b$j + mu_j)))) - e0j
  #standard error
  #eCj_se <- sd()/sqrt(M)
  
  #4. epsilon^(E#C)j bar (eq 85)
  #estimate
  eECshj_hat <- mean(log(1-delta+delta* ((sigma_j*b_til$j+mu_j)/(sigma_j*b_til$i+mu_j)))) - eEj_hat - eCj_hat - e0j
  #standard error
  #eECsharpj_se <- sd( )/sqrt(M)  
  
  #5. epsilon^[EC]j bar (eq 93)
  #estimate
  eECj_hat <- 0
  #standard error
  #eECj_se <- sd()
  
  #6. epsilon^[E||C] bar (eq 91)
  eECpsj_hat <- -(mean(log(1-delta+delta* ((sigma_j*b_til$j+mu_j)/(sigma_j*b_til$i+mu_j)))))
  
  #7. r_j\i bar 
  r_ji <- 0
  
  #Deltas (fourth list in S6)
  #q=1
  #Delta estimates:
  D0 <- e0i                                  #Deltai_0
  DE <- eEi_hat - eEj_hat                   #Deltai_E
  DC <- eCi_hat - eCj_hat                   #Deltai_C
  DECsh <- eECshi_hat - eECshj_hat 			#Deltai_(E#C)
  DEC <- eECi_hat - eECj_hat                #Deltai_[EC]
  DECps <- eECpsi_hat - eECpsj_hat           #Deltai_[E||C]
  Dr <- rii_hat 
  DrwoATA <- rii_hat - DECps
  
  #standard error of Delta estimates:
  #DE_se <- sd()/sqrt(M)
  #DC_se <- sd( )/sqrt(M)
  #DECsharp_se <- sd()/sqrt(M)
  #s1 <- sd()/sqrt(M) 
  #s2 <- sd()/sqrt(M) 
  #DEC_se <- sqrt(s1^2 + s2^2)
  #DECpip_se <- e_ECpipi_se 
  #Dr_se <- r_i_se
  
  #q=...
  qij <- mu_i/(mu_j*(1-delta)+ mu_i*delta)
  
  #Delta estimates:
  DEq <- eEi_hat - qij*eEj_hat
  DCq <- eCi_hat - qij*eCj_hat
  DECshq <- eECshi_hat - qij*eECshj_hat
  DECq <- eECi_hat - qij*eECj_hat
  DECpsq <- eECpsi_hat - qij*eECpsj_hat 
  
  #standard error of Delta estimates:
  #DEq_se <- sd()/sqrt(M)
  #DCq_se <- sd()/sqrt(M)
  #DECsharpq_se <- sd()/sqrt(M) 
  #DECq_se <- sqrt(s1^2 + s2^2)
  
  #return:
  #8 columns: ei, ei_se, ej, ej_se, D, D_se, Dq, Dq_se 
    # epsilon estimate and standard error for each species and 
    # Delta estimate and standard error for each qij alternative
  #7 rows: 0, E, C, (E#C), [EC], [E||C], r, rwoATA 
    # coexistence mechanims and gwr and gwr without ATA contribution 
  ei <- c(e0i, eEi_hat, eCi_hat, eECshi_hat, eECi_hat, eECpsi_hat, rii_hat, NA)
  #ei_se <- c(0, eEi_se, eCi_se, eECshi_se, eECi_se, eECpsi_se, ri_se, NA)
  ej <- c(e0j, eEj_hat, eCj_hat, eECshj_hat, eECj_hat, eECpsj_hat, r_ji, NA)
  #ej_se <- c(0, eEj_se, eCj_se, eECshj_se, eECj_se, 0, 0, NA)
  D <- c(D0, DE, DC, DECsh, DEC, DECps, Dr, DrwoATA)
  #D_se <- c(0, DE_se, DC_se, DECsh_se, DEC_se, DECps_se, Dr_se, NA)
  Dq <- c(D0, DEq, DCq, DECshq, DECq, DECpsq, Dr, DrwoATA)
  #Dq_se <- c(0, DEq_se, DCq_se, DECshq_se, DECq_se, DECps_se, Dr_se, NA)
  
  #res <- data.frame(ei=ei, ei_se=ei_se, ej=ej, ej_se=ej_se, D=D, D_se=D_se, Dq=Dq, Dq_se=Dq_se)
  res <- data.frame(ei=ei, ej=ej, D=D, Dq=Dq)
  row.names(res) <- c("0","E","C","(E#C)","[EC]","[E||C]","r", "rwoATA")
    if (plot==TRUE){
  		par(mfrow=c(2,1))
  		hist(b$i*sigma_i+mu_i, main='', xlab='',sub=paste("mu_i=",mu_i,",mu_j=",mu_j,",sigma_i=",sigma_i,",sigma_j=",sigma_j), col="black", xlim=c(0,ifelse(ub_i>ub_j, ub_i, ub_j)))
  		hist(b$j*sigma_j+mu_j, add=TRUE)
  		barplot(res$D[1:7],names.arg=rownames(res[1:7,]), main=paste("delta=",delta))
  }

  return(res)
}




