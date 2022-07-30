#transforrming normal noise to beta
F <- function(noise,...){
	qbeta(pnorm(noise),...)-0.5
}

obj <- F(b[[1]], shape1=0.5, shape2=0.5)
hist(obj)

#making noise to be used
# A
## 1) (v_i_r, v_j_r), (v_i_l, v_j_l): 
	# left and right TA bivariate noise with standard normal marginals 
## 2) (v_i_tilde, v_j_tilde): 
	# bivaraite normal mean (0,0) and cov (1,0,0,1)
## 3) (v_i_umlaut, v_j_umlaut):
	# bivariate normal mean (0,0) and cov (1,rho,rho,1) 
		# where rho is cov between v_j and v_i
# B
## make beta b_l, b_r, b_umlaut, b_tilde 
	# by taking the F of the respective (v_i, v_j)'s 

source("./lottery_beta/makenoise_LB.R")
noise_loc <- paste(numeric_results_loc, "/noise.RData", sep = "")
load(noise_loc) #object is called v
hist(v$l[,1]); mean(v$l[,1]); sd(v$l[,1])
hist(v$l[,2]); mean(v$l[,2]); sd(v$l[,2])
hist(v$r[,1]); mean(v$r[,1]); sd(v$r[,1])
hist(v$r[,2]); mean(v$r[,2]); sd(v$r[,2])
hist(v$til[,1]); mean(v$til[,1]); sd(v$til[,1])
hist(v$til[,2]); mean(v$til[,2]); sd(v$til[,2])
hist(v$um[,1]); mean(v$um[,1]); sd(v$um[,1])
hist(v$um[,2]); mean(v$um[,2]); sd(v$um[,2])
plot(v$l[1:10000,1], v$l[1:10000,2]); cov(v$l)
plot(v$r[1:10000,1], v$r[1:10000,2]); cov(v$r)
plot(v$til[1:10000,1], v$til[1:10000,2]); cov(v$til)
plot(v$um[1:10000,1], v$um[1:10000,2]); cov(v$um)
b <- lapply(v, FUN=F, shape1=0.5, shape2=0.5)
#check
hist(b$l[,1]); mean(b$l[,1]); sd(b$l[,1])
hist(b$l[,2]); mean(b$l[,2]); sd(b$l[,2])
hist(b$r[,1]); mean(b$r[,1]); sd(b$r[,1])
hist(b$r[,2]); mean(b$r[,2]); sd(b$r[,2])
hist(b$til[,1]); mean(b$til[,1]); sd(b$til[,1])
hist(b$til[,2]); mean(b$til[,2]); sd(b$til[,2])
hist(b$um[,1]); mean(b$um[,1]); sd(b$um[,1])
hist(b$um[,2]); mean(b$um[,2]); sd(b$um[,2])
plot(b$l[1:10000,1], b$l[1:10000,2]); cov(b$l)
plot(b$r[1:10000,1], b$r[1:10000,2]); cov(b$r)
plot(b$til[1:10000,1], b$til[1:10000,2]); cov(b$til)
plot(b$um[1:10000,1], b$um[1:10000,2]); cov(b$um)
#good

#now can make Bi's and Bj's
#note: should only consider sigma and mu such that
	# mu = sigma/2; so that B is never negative 





lottery <- function(b, b_sh, mudif, sigma, delta, method){
	#b <- lapply(b, function(X) {cbind(sigma*X[,1],sigma*X[,2]+mudif)})
	#b_sh <- lapply(b_sh, function(X) {cbind(sigma*X[,1],sigma*X[,2]+mudif)})
	
	if (method=="beta"){
		b_uni <- lapply(b, FUN=pnorm)
		b_uni_sh <- lapply(b_sh, FUN=pnorm)
		B <- lapply(b_uni, FUN=qbeta, shape1=0.5, shape2=0.5)
		B_sh <- lapply(b_uni_sh, FUN=qbeta, shape1=0.5, shape2=0.5)
		
		B <- lapply(B, function(X) {cbind(sigma*X[,1],sigma*X[,2]+mudif)})
		B_sh <- lapply(B_sh, function(X) {cbind(sigma*X[,1],sigma*X[,2]+mudif)})
	}
	
	if (method=="exp"){
		b <- lapply(b, function(X) {cbind(sigma*X[,1],sigma*X[,2]+mudif)})
		b_sh <- lapply(b_sh, function(X) {cbind(sigma*X[,1],sigma*X[,2]+mudif)})	
		B <- lapply(b, FUN=exp)
		B_sh <- lapply(b_sh, FUN=exp)
	}
	
	B1 <- B$s[,1]; B2 <- B$s[,2]
	B1_sh <- B_sh$s[,1]; B2_sh <- B_sh$s[,2]
	
	B1a <- B$l[,1]; B2a <- B$l[,2]
	B1a_sh <- B_sh$l[,1]; B2a_sh <- B_sh$l[,2]	
	
	C1 <- C2 <- B2/delta
	C1a <- C2a <- B2a/delta

	r1 = log(1-delta + B1/C1);  
	r2 = log(1-delta + B2/C2); 
	r1a = log(1-delta + B1a/C1a);  
	r2a = log(1-delta + B2a/C2a);

	rsh1 = log(1-delta + B1_sh/C1)
	rsh2 = log(1-delta + B2_sh/C2) 
	rsh1a = log(1-delta + B1a_sh/C1a)
	rsh2a = log(1-delta + B2a_sh/C2a) 


	rbar1 = mean(r1); rshbar1 = mean(rsh1) 
	rbar2 = mean(r2); rshbar2 = mean(rsh2) 

	rbar1a = mean(r1a); rshbar1ata = mean(rsh1a) 
	rbar2a = mean(r2a); rshbar2ata = mean(rsh2a) 

	q12 = 1; 

	DeltaI = rbar1 - rshbar1 + q12*rshbar2;  

	DeltaIata = rbar1ata - rshbar1ata + q12*rshbar2ata;  

	ata <- DeltaI - DeltaIata
	
	return(c("ata"=ata, "rbar1"=rbar1, "rbar1ata"=rbar1ata))
}

b <- makenoise(1000000)[c(1,3)]
b_sh <- makenoise(1000000)[c(1,3)]

lottery(b, b_sh, mudif=0.2, sigma=2, delta=0.5, method="exp")
lottery(b, b_sh, mudif=0.2, delta=0.6, sigma=0.05, method="beta")


#mean(bnorm$l[,1]);mean(bnorm$s[,1]);mean(bnorm$l[,2]);mean(bnorm$s[,2])
b_exp<-

#hist(b_unis$l[,2])
#plot(b_unis$l[1:1000,1], b_unis$l[1:1000,2])
#plot(b_unis$s[1:1000,1], b_unis$s[1:1000,2])
#hist(b_unis$s[,1])
#hist(b_unis$s[,2])



b_betas <- lapply(b_unis, FUN=qbeta, shape1=0.5, shape2=0.5)
#hist(b_betas$l[,2])
#plot(b_betas$l[1:1000,1], b_betas$l[1:1000,2])
#plot(b_betas$s[1:1000,1], b_betas$s[1:1000,2])
#hist(b_betas$s[,1])
#hist(b_betas$s[,2])
b_betash <- lapply(b_unissh, FUN=qbeta, shape1=0.5, shape2=0.5)

B1 <- b_betas$s[,1]; B2 <- b_betas$s[,2]
B1ata <- b_betas$l[,1]; B2ata <- b_betas$l[,2]

B1sh <- b_betash$s[,1]; B2sh <- b_betash$s[,2]
B1atash <- b_betash$l[,1]; B2atash <- b_betash$l[,2]

B1norm <- bnorm$s[,1]; B2norm <- bnorm$s[,2]
B1normata <- bnorm$l[,1]; B2normata <- bnorm$l[,2]

B1shnorm <- bnormsh$s[,1]; B2shnorm <- bnormsh$s[,2]
B1atashnorm <- bnormsh$l[,1]; B2atashnorm <- bnormsh$l[,2]

## Step 3a: simulate to generate C1(t), C2(t), and r1(t).
## We don't actually need r2(t) because r2bar=0 when species 2 is resident. 
## 
## In the lottery model the resident species is always
## at population=N (total number of sites). The competition  
## parameters C1 and C2 both equal the ratio between the total number
## of competing larvae, and the number of open sites; so when species
## 2 is resident C1 = C2 = (B2*N)/(delta*N)=B2/delta. 
delta <- 0.4
C1 <- C2 <- B2/delta
C1ata <- C2ata <- B2ata/delta

C1norm <- C2norm <- B2norm/delta
C1atanorm <- C2atanorm <- B2normata/delta

r1 = log(1-delta + B1/C1);  
r2 = log(1-delta + B2/C2); 
r1ata = log(1-delta + B1ata/C1ata);  
r2ata = log(1-delta + B2ata/C2ata);

r1 = log(1-delta + B1norm/C1norm);  
r2 = log(1-delta + B2norm/C2norm); 
r1ata = log(1-delta + B1ata/C1ata);  
r2ata = log(1-delta + B2ata/C2ata);

## Step 3b: use C1(t) and C2(t) with the E-sharps to 
## calculate r1.sharp and r2.sharp 
rsh1 = log(1-delta + B1sh/C1)
rsh2 = log(1-delta + B2sh/C2) 
rsh1ata = log(1-delta + B1atash/C1ata)
rsh2ata = log(1-delta + B2atash/C2ata) 

## Step 4: compute the average growth rates
rbar1 = mean(r1); rshbar1 = mean(rsh1) 
rbar2 = mean(r2); rshbar2 = mean(rsh2) 

rbar1ata = mean(r1ata); rshbar1ata = mean(rsh1ata) 
rbar2ata = mean(r2); rshbar2ata = mean(rsh2ata) 

## Step 5: compute the scaling factors. For this model we know them. 
q12 = 1; 

## Step 5: calculate storage effect for species 1, Delta.Ib1
DeltaI = rbar1 - rshbar1 + q12*rshbar2;  

DeltaIata = rbar1ata - rshbar1ata + q12*rshbar2ata;  

DeltaI - DeltaIata

# Output 
cat("rbar.1 = ",round(rbar.1,digits=3),"\n"); 
cat("rsharp.1 = ",round(rsharp.1,digits=3),"\n"); 
cat("rsharp.2 = ",round(rsharp.2,digits=3),"\n"); 
cat("Delta.Ib1 = ", round(Delta.Ib1, digits=3),"\n");  
