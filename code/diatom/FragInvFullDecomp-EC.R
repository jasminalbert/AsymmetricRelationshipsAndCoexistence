source("partialSharp_fxns.R")
source("DataAndVKQfuns.R")
require(deSolve)


######################################
# Model 
######################################
forceChemo <- function(t,y,parms) {
	temp <- parms["Tbar"] + parms["a"]*sin(2*pi*t/parms["P"]); 
	V1 <- V1quad(temp); V2<- V2modfun(temp); K1<- K1flatfun(temp); K2 <-  K2flatfun(temp) 
	Q1 <- Q1fun(temp); Q2 <-  Q2fun(temp);
	R <- y[1]; x1<- y[2]; x2<- y[3]; 
	up1 <- V1*R/(K1 + R); 
	up2 <- V2*R/(K2 + R)  
	dR <- parms["D"]*(parms["S"]-R) - Q1*x1*up1 - Q2*x2*up2; 
 	D <- parms["D"]; names(D) <- NULL;
	dx1 <- x1*(up1-D); 
	dx2 <- x2*(up2-D); 
	names(V1) <- NULL; names(up1) <- NULL; names(up2) <- NULL; names(R) <- NULL; 
	return( list(dx = c(dR,dx1,dx2), temp = temp) ) 
}

rEC <-  function(E,C,parms) {E/C - parms["D"]} 


q12 <- 1
reps <- 500
# parameters #
parms <- c(Tbar=18, a=6, P=60, D=0.09, S=35)

########################################### 
# simulate with 1 as invader, 2 as resident 
###########################################
y0 <- c(R=.1,x1=0,x2=20)
times <- seq(0,3600,by=.1)
out1i <- ode(y0,times,func=forceChemo,parms=parms) 
e <- times > max(times-1200)
out1i <- out1i[e,] 
minus1 <- data.frame(out1i) 
names(minus1) <- c("time","R","x1","x2","temp") 
    
#### repeat with no fluctuations     
parms2 <- c(Tbar=18, a=0, P=60, D=0.09, S=35)    
out1i <- ode(y0,times,func=forceChemo,parms=parms2) 
e <- times > max(times-1200)
out1i <- out1i[e,]
minus1noVar <- data.frame(out1i)
names(minus1noVar) <- c("time","R","x1","x2","temp")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
# Calculations for 1 invading 2 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###############################
# E and C vars for 1 invading 2 
###############################
E1 <- V1quad(minus1$temp)
E10 <- mean(E1)

E1noVar <- V1quad(minus1noVar$temp)
E1star <- mean(E1noVar)

E2 <- V2modfun(minus1$temp)
E20 <- mean(E2)

E2noVar <- V2modfun(minus1noVar$temp)
E2star <- mean(E2noVar)

C1 <- (K1flatfun(minus1$temp)+minus1$R)/minus1$R
C10 <- mean(C1)

C1noVar <- (K1flatfun(minus1noVar$temp) + minus1noVar$R)/minus1noVar$R
C1star <- mean(C1noVar)

C2 <- (K2flatfun(minus1$temp)+minus1$R)/minus1$R
C20 <- mean(C2)

C2noVar <- (K2flatfun(minus1noVar$temp) + minus1noVar$R)/minus1noVar$R
C2star <- mean(C2noVar)

### "exact sharping" - all possible combinations 
EC1 <- expand.grid(E1, C1)
E1sharp <- EC1[,1]; C1sharp <- EC1[,2]
rm(EC1)

EC2 <- expand.grid(E2, C2)
E2sharp <- EC2[,1]; C2sharp <- EC2[,2]
rm(EC2)

### "partial sharping" - 
EC1psharp <- makePsharp(E1, C1, reps)
E1psharp <- EC1psharp[,1,]; C1psharp <- EC1psharp[,2,]

EC2psharp <- makePsharp(E2, C2, reps)
E2psharp <- EC2psharp[,1,]; C2psharp <- EC2psharp[,2,]

######################
# r's for 1 invading 2 
######################

# no variance baseline, r1(E*, C*)
r1star <- rEC(E1star, C1star, parms2) 
r2star <- rEC(E2star, C2star, parms2)

# mean (E,C) in presence of varying temperature, r(E0, C0)
r10 <- rEC(E10, C10, parms)
r20 <- rEC(E20, C20, parms)

# rbar(E#, C#)
r1sharp <- mean(rEC(E1sharp, C1sharp, parms))
r2sharp <- mean(rEC(E2sharp, C2sharp, parms))

# rbar(E||, C||)
r1psharpsims <- apply(EC1psharp, MARGIN=3, function(EC){ mean(rEC(EC[,1], EC[,2], parms)) } );
r1psharp <- median(r1psharpsims)

r2psharpsims <- apply(EC2psharp, MARGIN=3, function(EC){ mean(rEC(EC[,1], EC[,2], parms)) } );
r2psharp <- median(r2psharpsims)

# rbar(E, C)
r1bar <- mean(rEC(E1, C1, parms))
r2bar <- mean(rEC(E2, C2, parms))

###########################
# epsilons for 1 invading 2 
###########################

# epsilon{prime}
eps1 <- r10 - r1star
eps2 <- r20 - r2star

# epsilonE
eps1E <- mean(rEC(E1, C10, parms)) - r10
eps2E <- mean(rEC(E2, C20, parms)) - r20

# epsilonC
eps1C <- mean(rEC(E10, C1, parms)) - r10
eps2C <- mean(rEC(E20, C2, parms)) - r20

# epsilonEC
eps1EC <- r1bar - eps1E - eps1C - r10
eps2EC <- r2bar - eps2E - eps2C - r20

# epsilon(E#C)
eps1EsharpC <- r1sharp - eps1E - eps1C - r10
eps2EsharpC <- r2sharp - eps2E - eps2C - r20

# epsilon(EC)
eps1ECpar <- r1bar - r1sharp
eps2ECpar <- r2bar - r2sharp

# epsilon[EC]
eps1ECbrk <- r1bar - r1psharp
eps2ECbrk <- r2bar - r2psharp

# epsilon[E||C]
eps1EpsharpC <- r1psharp - r1sharp
eps2EpsharpC <- r2psharp - r2sharp

##########################
# Deltas for 1 invading 2
##########################

Delta1star <- r1star - q12*r2star
Delta10 <- eps1 - q12*eps2
Delta1E <- eps1E - q12*eps2E
Delta1C <- eps1C - q12*eps2C
Delta1EC <- eps1EC - q12*eps2EC

#Delta1EC = 
Delta1EsharpC <- eps1EsharpC - q12*eps2EsharpC
Delta1ECpar <- eps1ECpar - q12*eps2ECpar

#Delta1ECpar = |/
Delta1ECbrk <- eps1ECbrk - q12*eps2ECbrk
Delta1EpsharpC <- eps1EpsharpC - q12*eps2EpsharpC



#inaccuracy in Delta10 and Delta1E is solved when summing 
#inconsistency in Delta1C and DeltaEsharpC is solved when summing
## WHY ##?????

IGR <- sum(Delta1star, Delta10, Delta1E, Delta1C, Delta1EC)
# = r1bar-r2bar
IGRa <- sum(Delta1star, Delta10, Delta1E, Delta1C, Delta1EsharpC, Delta1ECpar)
IGRb <- sum(Delta1star, Delta10, Delta1E, Delta1C, Delta1EsharpC, Delta1ECbrk, Delta1EpsharpC)

IGR_origFuns <- c(Delta1star, Delta10, Delta1E, Delta1C, Delta1EsharpC, Delta1ECbrk, Delta1EpsharpC)
IGR_modFuns <- c(Delta1star, Delta10, Delta1E, Delta1C, Delta1EsharpC, Delta1ECbrk, Delta1EpsharpC)

#same 
IGR > Delta1ECbrk
Delta1ECbrk>IGR #TRUE
IGR - Delta1ECbrk #IGR would be negative without 

barplot(c(Delta1star, Delta10, Delta1E, Delta1C, Delta1EsharpC, Delta1ECbrk, Delta1EpsharpC, IGR), names.arg=c("Delta*","Delta'","Delta^E","Delta^C", "Delta^(E#C)","Delta^[EC]", "Delta^[E||C]","IGR") )

#plotting to check ATAs and correlation
pv <- function(EC){
	E <- EC[,1]; C <- EC[,2]
	lim <- c(-4,4)
	P <- cor(E,C)
	blue <- rgb(0,0,0.545,0.2)
	plot(E,C, xlim=lim, ylim=lim, xlab='E', ylab='C', col=blue, cex=0.8)
	title(main=bquote(~ rho == .(round(P,4))),line=0.7, cex=0.8)
}

E1norm <- qnorm(pobs(E1))
C1norm <- qnorm(pobs(C1))
#check histograms
EC <- cbind(E1norm, C1norm)
pv(EC)
ECleft <- EC[EC[,1]<0,]
pv(ECleft)
ECright <- EC[EC[,1]>0,]
pv(ECright)

ECupp <- EC[EC[,2]>0,]
pv(ECupp)
EClow <- EC[EC[,2]<0,]
pv(EClow)



pdf("ECnorm.pdf", width=5.5, height=5.5)
par(mar=c(2,2,2,2), oma=c(2,2,3,2))
layout(matrix(c(1,1,2,
				1,1,3,
				4,5,6), nrow=3, byrow=T))

pv(EC)
pv(ECleft); mtext("left", side=1, line=-1.3, at=4, adj=1, cex=.8)
pv(ECright); mtext("right", side=3, line=-1.3, at=4, adj=1, cex=.8)
pv(ECupp); mtext("upper", side=3, line=-1.3, at=4, adj=1, cex=.8)
pv(EClow); mtext("lower", side=3, line=-1.3, at=-4, adj=0, cex=.8)

mtext("Normalized rank distributions", outer=TRUE, side=3, line=0.5, cex.lab=1.25)
mtext("E", outer=TRUE, side=1, line=0.8, cex.lab=1.25)
mtext("C", outer=TRUE, side=2, line=0.5, cex.lab=1.25)
dev.off()

EC <- EC1psharp[,,1]
plot(EC[,1], EC[,2])
Epnorm<-qnorm(pobs(EC[,1]))
Cpnorm<-qnorm(pobs(EC[,2]))
EC <- cbind(Epnorm, Cpnorm)
pv(EC)

normcor(E1, C1)
normcor(EC1psharp[,1,1],EC1psharp[,2,1])























