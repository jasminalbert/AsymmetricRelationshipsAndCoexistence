# a     amplitude
# P     Period (days)
# Tbar  mean temperature (degrees C); theta_0
# time  total time to simulate invasion
# reps  number repititions for simulating symmertric EC (for computing ATA contribution)
# sp: calculate epsilons for 
  # 1: sp 1                    2: sp 2
# invader: species that starts at 0 (other starts at 20)
  # 1: sp 1                    2: sp 2 
# methods: E and r definitions
  # 1: E=temp, r=V(E)/C - D    2: E=V(temp), r=E/C - D
  # use method 1 in publication but use method 2 to check against results in Ellner_2019 SI

getep <- function(a, P, Tbar, time, reps, sp, invader=1, method=1){
  
  parms <- c(Tbar=Tbar, a=a, P=P, D=0.09, S=35)
  #Time <- time 
  times <- seq(0,Time,by=.1)
  
  y0 <- c(R=.1,x1=0,x2=20)
  if (invader==2){ #simulate with 2 as invader, 1 as resident
    y0["x1"]=20; y0["x2"]=0
  }
  
  out <- ode(y0,times,func=forceChemo,parms=parms) 
  cut <- time*(2/3)
  burn <- times > time - round((time - cut)/P,0)*P #integer multiple of P
  out <- out[burn,]
  temp <- out[,5]
  R <- out[,2]
  #use forceChemo to get R to get C 
  
  if (length(temp)%%P!=0){stop("dims are not an integer multiple of P")}
  
  if (sp==2){
    Vfun <- V2modfun; Kfun <- K2flatfun
  } else {Vfun <- V1quad; Kfun <- K1flatfun}
  
  if (method==2){
    r <- function(E, C, parms){E/C - parms["D"]}
    E <- Vfun(temp)
  } else {
    r <- function(E, C, parms){Vfun(E)/C - parms["D"]}
    E <- temp
  }
  
  C <- (Kfun(temp)+R)/R
  
  #special E's and C's
  E0 <- mean(E)
  C0 <- mean(C)
  
  EC <- expand.grid(E, C)
  Esharp <- EC[,1]; Csharp <- EC[,2]
  rm(EC)
  
  ECpsharp <- makePsharp(E, C, reps)
  Epsharp <- ECpsharp[,1,]; Cpsharp <- ECpsharp[,2,]
  
  #special r's
  rsharp <- mean(r(Esharp, Csharp, parms))
  
  rpsharpsims <- apply(ECpsharp, MARGIN=3, function(EC)			{ mean(r(EC[,1], EC[,2], parms)) } );
  rpsharp <- median(rpsharpsims)
  
  rbar <- mean(r(E, C, parms))
  
  #epsilons
  eps0 <- r(E0, C0, parms)
  epsE <- mean(r(E, C0, parms)) - eps0
  epsC <- mean(r(E0, C, parms)) - eps0
  #epsEC <- rbar - epsE - epsC - eps0
  epsEsharpC <- rsharp - epsE - epsC - eps0
  #epsECpar <- rbar - rsharp
  epsECbrk <- rbar - rpsharp
  epsEpsharpC <- rpsharp - rsharp
  
  return(data.frame(eps0=eps0, epsE=epsE, epsC=epsC, epsEsharpC=epsEsharpC, epsECbrk=epsECbrk, epsEpsharpC=epsEpsharpC, time=time))
  
}

getDelt <- function(a, P, Tbar, time, sims, invader=1){
  ep1 <- getep(a, P, Tbar, time, reps=sims, sp=1, invader=invader)
  ep2 <- getep(a, P, Tbar, time, reps=sims, sp=2, invader=invader)
  
  if (invader==2){Delta1 <- ep2-ep1} else {Delta1 <- ep1-ep2}

  Delta1 <- Delta1[,-7]
  Delta1$IGR <- sum(Delta1)
  Delta1$map <- Delta1$epsECbrk/Delta1$IGR
  Delta1$time <- ep1$time
  return(Delta1)
}
m2_1 <- getDelt(6, 60, 18, 3000, 200, invader=1)
m2_2 <- getDelt(6, 60, 18, 3000, 200, invader=2)

m1_1 <- getDelt(6, 60, 18, 3000, 200, invader=1)
m1_2 <- getDelt(6, 60, 18, 3000, 200, invader=2)
print(m1_1); print(m1_2)

dat_loc <- "../results_numeric/fig5dat/"
if(dir.exists(dat_loc)==FALSE){
  dir.create(dat_loc)
}
#amplitude
a <- seq(1,6,length.out=100)
Delta_a2 <- getDelt(a[1], 60, 18, 3000, 200, invader=2)

for (i in 2:length(a)){
  Delta_a2[i,] <- getDelt(a[i], 60, 18, 3000, 200, invader=2)
  print(i)
}
Delta_a2$a <- a
saveRDS(Delta_a2, paste0(dat_loc,'varyAmplitude2.RDS'))
################################################################
#Period
P <- seq(51,199.5,length.out=100)
DeltaP2 <- getDelt(6, P[1], 18, 3000, 200, invader=2)

for (i in 2:length(P)){
  print(paste(i, P[i]))
  DeltaP2[i,] <- getDelt(6, P[i], 18, 3000, 200, invader=2)
}
DeltaP2$P <- P
saveRDS(DeltaP2, paste(dat_loc,'varyPeriod2.RDS',sep=''))
################################################################
#mean temp
Tbar <- seq(16,18,length.out = 100)
DeltaT2 <- getDelt(6, 60, Tbar[1], 3000, 200, invader=2)

for (i in 2:length(Tbar)){
  DeltaT2[i,] <- getDelt(6, 60, Tbar[i], 3000, 200, invader=2)
  print(i)
}
DeltaT2$Tbar <- Tbar
saveRDS(DeltaT2, paste(dat_loc,'varyTbar2.RDS',sep=''))
################################################################
fig_loc <- "../results_figs/"
if(dir.exists(fig_loc)==FALSE){
  dir.create(fig_loc)
}
fig5_2 <- paste(fig_loc,"fig5_2.pdf",sep="")

vAmp <- readRDS(paste(dat_loc,'varyAmplitude2.RDS',sep='')) #from makefig5dat.R
vPer <- readRDS(paste(dat_loc,'varyPeriod2.RDS',sep=''))
vTbr <- readRDS(paste(dat_loc,'varyTbar2.RDS',sep=''))

col <- c(rep("black",4), "blue", "red", "orange")
line <- c(1:4, 1,1,1)

varylist <- list(vAmp, vPer, vTbr)
range <- sapply(varylist, function(X){range(X[,1:7])})
ymin <- min(range)
ymax <- max(range)
orig <- c(a=6, P=60, Tbar=18)

#default is c(5, 4, 4, 2) + 0.1
pdf(fig5_2, height=5, width=15)
par(mfrow=c(1,3), oma=c(0,3,0,0), mar=c(5,1,2,1) )
for (i in 1:3){
  vdat <- varylist[[i]]
  
  #empty box
  plot(0, yaxt='n', xlim=range(vdat[,10]), ylim=c(ymin, ymax), col='white', xlab=colnames(varylist[[i]])[10], ylab='', cex.lab=1.3, cex.axis=1.3)
  
  #axis
  if (i==1){
    axis(2)
  }
  
  #loop for lines
  for (j in 1:7){
    lines(vdat[,10], vdat[,j], col=col[j], lty=line[j])	
  }
  abline(h=0, col="lightgrey") #zero
  abline(v=orig[i], col='magenta', lty=3)
  
  #label
  mtext(paste0("(", letters[i],")"), 3, -1.5, adj=0.985)
  
  #legend
  if (i==1){
    legend("topleft", legend=c(expression(Delta[i]^0), expression(Delta[i]^E), 				expression(Delta[i]^C), expression(Delta[i]^"(E#C)"), expression(Delta[i]^"[E||C]"), expression(Delta[i]^"[EC]") ,expression(IGR)), col = col,lty = line, bty="n", cex=1.2, inset=c(0,0))
  }
}
title(ylab="coexistence", outer=T, line=1.5, cex.lab=1.3)

dev.off()












