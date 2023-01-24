#This script makes transforms b made in makenoise.R to beta distributed fecundities, B

### location of results ###
numeric_results_loc <- "../results_numeric"

noise_loc <- paste0(numeric_results_loc, "/noise.RData")
load(noise_loc) #list of noise, b, containing:

	## 1) (b_i_r, b_j_r), (b_i_l, b_j_l): 
		# left and right TA bivariate noise with standard normal 		marginals 
	## 2) (b_i_tilde, b_j_tilde): 
		# bivaraite normal mean (0,0) and cov (1,0,0,1)
	## 3) (b_i_s, b_j_s): (umlaut)
		# bivariate normal mean (0,0) and cov (1,rho,rho,1) 
			# where rho is cov between b_j_l and b_i_l

#now make B
## make beta b_l, b_r, b_s, b_tilde 
	# by taking the F of the respective (b_i, b_j)'s

#transforming normal noise to beta
F <- function(noise,...){
	qbeta(pnorm(noise),...)
}

#looking at b
hist(b$l[,1]); mean(b$l[,1]); sd(b$l[,1])
hist(b$l[,2]); mean(b$l[,2]); sd(b$l[,2])
hist(b$r[,1]); mean(b$r[,1]); sd(b$r[,1])
hist(b$r[,2]); mean(b$r[,2]); sd(b$r[,2])
hist(b$til[,1]); mean(b$til[,1]); sd(b$til[,1])
hist(b$til[,2]); mean(b$til[,2]); sd(b$til[,2])
hist(b$s[,1]); mean(b$s[,1]); sd(b$s[,1])
hist(b$s[,2]); mean(b$s[,2]); sd(b$s[,2])
plot(b$l[1:10000,1], b$l[1:10000,2]); cov(b$l)
plot(b$r[1:10000,1], b$r[1:10000,2]); cov(b$r)
plot(b$til[1:10000,1], b$til[1:10000,2]); cov(b$til)
plot(b$s[1:10000,1], b$s[1:10000,2]); cov(b$s)

#transform
betaB <- lapply(b, FUN=F, shape1=0.5, shape2=0.5)
#check
hist(betaB$l[,1]); mean(betaB$l[,1]); sd(betaB$l[,1])
hist(betaB$l[,2]); mean(betaB$l[,2]); sd(betaB$l[,2])
hist(betaB$r[,1]); mean(betaB$r[,1]); sd(betaB$r[,1])
hist(betaB$r[,2]); mean(betaB$r[,2]); sd(betaB$r[,2])
hist(betaB$til[,1]); mean(betaB$til[,1]); sd(betaB$til[,1])
hist(betaB$til[,2]); mean(betaB$til[,2]); sd(betaB$til[,2])
hist(betaB$s[,1]); mean(betaB$s[,1]); sd(betaB$s[,1])
hist(betaB$s[,2]); mean(betaB$s[,2]); sd(betaB$s[,2])
plot(betaB$l[1:10000,1], betaB$l[1:10000,2]); cov(betaB$l)
plot(betaB$r[1:10000,1], betaB$r[1:10000,2]); cov(betaB$r)
plot(betaB$til[1:10000,1], betaB$til[1:10000,2]); cov(betaB$til)
plot(betaB$s[1:10000,1], betaB$s[1:10000,2]); cov(betaB$s)
#good; save
betanoise_loc <- paste0(numeric_results_loc, "/betanoise.RDS")
saveRDS(betaB, file=betanoise_loc)

#now can make Bi's and Bj's






