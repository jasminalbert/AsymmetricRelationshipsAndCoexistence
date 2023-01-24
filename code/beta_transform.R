#transforming normal noise to beta
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

#source("./lottery_beta/makenoise_LB.R")
noise_loc <- paste(numeric_results_loc, "/noiseB.RData", sep = "")
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
betanoise_loc <- paste0(numeric_results_loc, "/betanoise.RDS")
saveRDS(b, file=betanoise_loc)

#now can make Bi's and Bj's
#note: should only consider sigma and mu such that
	# mu = sigma/2; so that B is never negative 





