#reverse invader resident scenario
# E=temp; r=V(E)/C-D
getep3 <- function(a, P, Tbar, time, reps, sp, invader=1){
  
  parms <- c(Tbar=Tbar, a=a, P=P, D=0.09, S=35)
  Time <- time 
  times <- seq(0,Time,by=.1)
  
  y0 <- c(R=.1,x1=0,x2=20)
  if (invader==2){ #simulate with 2 as invader, 1 as resident
    y0["x1"]=20; y0["x2"]=0
  }
  
  out1i <- ode(y0,times,func=forceChemo,parms=parms) 
  cut <- Time*(2/3)
  burn <- times > Time - round((Time - cut)/P,0)*P #integer multiple of P
  out1i <- out1i[burn,]
  temp <- out1i[,5]
  R <- out1i[,2]
  
  if (length(temp)%%P!=0){stop("dims are not an integer multiple of P")}
  
  if (sp==2){
    Vfun <- V2modfun; Kfun <- K2flatfun
  } else {Vfun <- V1quad; Kfun <- K1flatfun}
  
  
  r <- function(E, C, parms){Vfun(E)/C - parms["D"]}
  
  E <- temp #if after burning, the times is not an integer multiple P, temp is also not
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
  ep1 <- getep3(a, P, Tbar, time, reps=sims, sp=1, invader=invader)
  ep2 <- getep3(a, P, Tbar, time, reps=sims, sp=2, invader=invader)
  
  if (invader==2){Delta1 <- ep2-ep1}
  else {Delta1 <- ep1-ep2}

  Delta1 <- Delta1[,-7]
  Delta1$IGR <- sum(Delta1)
  Delta1$map <- Delta1$epsECbrk/Delta1$IGR
  Delta1$time <- ep1$time
  return(Delta1)
}
#getDelt(6, 60, 18, 3000, 200, invader=1)
#getDelt(6, 60, 18, 3000, 200, invader=2)


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












