#setwd("~/Documents/AsymmetricRelationshipsAndCoexistence/code")
source("diatom/DataAndVKQfuns.R")
source("diatom/partialSharp_fxns.R")
require(deSolve)



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
  times <- seq(0,time,by=.1)
  
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
  
  return(data.frame('0'=eps0, E=epsE, C=epsC, EsharpC=epsEsharpC, ECbrk=epsECbrk, EpsharpC=epsEpsharpC, time=time))
}

getDelt <- function(a, P, Tbar, time, sims, invader=1){
  ep1 <- getep(a, P, Tbar, time, reps=sims, sp=1, invader=invader)
  ep2 <- getep(a, P, Tbar, time, reps=sims, sp=2, invader=invader)
  
  if (invader==2){Delta1 <- ep2-ep1} else {Delta1 <- ep1-ep2}
  
  Delta1 <- Delta1[,-7]
  Delta1$GWR <- sum(Delta1)
  Delta1$map <- Delta1$ECbrk/Delta1$GWR
  Delta1$time <- ep1$time
  return(Delta1)
}

wrapDelt <- function(args){
  Deltas <- getDelt(a=args[1], P=args[2], Tb=args[3], time=args[4], sims=args[5], invader=args[6])
  return(Deltas)
}

