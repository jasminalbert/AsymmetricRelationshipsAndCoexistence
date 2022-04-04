#functions for decomposition: computing epsilons and Deltas
  #getep, getDelt, wrapDelt

##libraries used (invoked with ::): deSolve, stats

### sourcing ###
source("diatom/diatomModel_andVKQ.R")
source("diatom/partialSharp_fxns.R")

### getep ###
# computes epsilons given parameters and definitions
#ARGS:
  #a        amplitude
  #P        Period (days)
  #Tbar     mean temperature (degrees C); theta_0
  #time     total time to simulate invasion
  #reps     number repititions for simulating symmertric EC (for computing ATA contribution)
  #sp       1 or 2; calculate epsilons for species 1 or 2 
  #invader  1 or 2; species that starts at 0 (other starts at 20)
  #methods  1 or 2; E and r definitions:
            #1: E=temp, r=V(E)/C - D;    2: E=V(temp), r=E/C - D
            #paper uses method 1 but use method 2 to check against results in Ellner_2019 SI
#OUT:
  #data.frame of epsilons: 0, E, C, (E#C), [EC], [E||C]; and time 
getep <- function(a, P, Tbar, time, reps, sp, invader=1, method=1, fig6=FALSE){
  
  parms <- c(Tbar=Tbar, a=a, P=P, D=0.09, S=35)
  
  times <- seq(0,time,by=.1)
  
  y0 <- c(R=.1,x1=0,x2=20)
  if (invader==2){ #simulate with sp 2 as rare, sp 1 as common
    y0["x1"]=20; y0["x2"]=0
  }
  
  #simulate competition model (see sourced script)
  #to get R (resouce concentration) 
  out <- deSolve::ode(y0,times,func=forceChemo,parms=parms)
  
  #cut out burn in time
  cut <- time*(2/3)
  burn <- times > time - round((time - cut)/P,0)*P #integer multiple of P
  out <- out[burn,]
  
  temp <- out[,5]
  R <- out[,2]
  
  #make sure simulation length after burning is a multiple of P
  if (length(temp)%%P!=0){stop("dims are not an integer multiple of P")}
  
  if (sp==2){ #species specific responses to temperature
    Vfun <- V2modfun; Kfun <- K2flatfun
  } else {Vfun <- V1quad; Kfun <- K1flatfun}
  
  if (method==2){ #methods of defining r and E
    r <- function(E, C, parms){E/C - parms["D"]}
    E <- Vfun(temp)
  } else {
    r <- function(E, C, parms){Vfun(E)/C - parms["D"]}
    E <- temp
  }
  
  #competition, C
  C <- (Kfun(temp)+R)/R
  
  if (fig6==FALSE){
    ### special E's and C's ###
    #mean; no variation
    E0 <- mean(E)
    C0 <- mean(C)
    
    #variation *per se*; E and C both vary but independently 
    #decouple by getting all possible combinations
    EC <- expand.grid(E, C)
    Esharp <- EC[,1]; Csharp <- EC[,2]
    rm(EC)
    
    ### special r's ###
    #rbar when E and C vary independently 
    rsharp <- mean(r(Esharp, Csharp, parms))
  }
  
  ### more special E and C ###
  #correlation *per se*; E and C covary but perfectly symmetrically
  #see source script; methods defined in SI section 9
  ECpsharp <- makePsharp(E, C, reps)
  Epsharp <- ECpsharp[,1,]; Cpsharp <- ECpsharp[,2,]
  
  ### more special r ###
  #get median of mean or normalized ranking simulations
  #rbar when E and C covary symmetrically 
  rpsharpsims <- apply(ECpsharp, MARGIN=3, function(EC)			{ mean(r(EC[,1], EC[,2], parms)) } );
  rpsharp <- stats::median(rpsharpsims)
  #get SE from 100 bootstrapped samples
  se <- getsdbst(rpsharpsims, q=100)
  
  #true rbar
  rbar <- mean(r(E, C, parms))
  
  #epsilons
  epsECbrk <- rbar - rpsharp
  
  if (fig6==FALSE){
    eps0 <- r(E0, C0, parms)
    epsE <- mean(r(E, C0, parms)) - eps0
    epsC <- mean(r(E0, C, parms)) - eps0
    #epsEC <- rbar - epsE - epsC - eps0
    epsEsharpC <- rsharp - epsE - epsC - eps0
    #epsECpar <- rbar - rsharp; #storage effect
    epsEpsharpC <- rpsharp - rsharp
    
    res <- data.frame('0'=eps0, E=epsE, C=epsC, EsharpC=epsEsharpC, ECbrk=epsECbrk, EpsharpC=epsEpsharpC, rbar=rbar, time=time, SEpsharp=se)
  } else {
    res <- data.frame(ECbrk=epsECbrk, rbar=rbar, time=time, SEpsharp=se)
  }
  
  return(res)
}


### getDelt ###
# get Deltas by subracting species epsilons term by term
#ARGS:
  #a        amplitude, single numeric
  #P        Period, single numeric
  #Tbar     Tbar, single numeric
  #time     time of competition model simulation
  #sims     number of normalized rank simulations (for EC partial sharp; correlation per se)
  #invader  1 or 2; which species is rare?
#OUT:
  #data.frame of Deltas: 0, E, C, (E#C), [EC], [E||C]; and time and mapping ratio (fig 6)
getDelt <- function(a, P, Tbar, time, sims, invader=1, fig6=FALSE){
  
  #compute epsilons for both species
  ep1 <- getep(a, P, Tbar, time, reps=sims, sp=1, invader=invader, fig6=fig6)
  ep2 <- getep(a, P, Tbar, time, reps=sims, sp=2, invader=invader, fig6=fig6)
  
  #subtract one from the other, depending on which one is invader
  if (invader==2){Delta1 <- ep2-ep1} else {Delta1 <- ep1-ep2}
  
  ncol <- ncol(Delta1)
  Delta1 <- Delta1[,-(ncol-1):-ncol]
  
  if (fig6==FALSE){
    Delta1$GWR <- sum(Delta1[,-7]) #GWR
  } else {Delta1$GWR <- Delta1$rbar}
  Delta1$map <- Delta1$ECbrk/Delta1$GWR #ATA/GWR (for fig 6)
  Delta1$time <- ep1$time
  Delta1$SE <- sqrt((ep1$SEpsharp^2)+(ep2$SEpsharp^2))
  return(Delta1)
}

### wrapDelt ###
# wrapper function for getDelt to make possible to use in mclapply 
#ARGS:
  #args   named vector of arguements 
wrapDelt <- function(args){
  fig6 <- ifelse(args[7]==1, TRUE, FALSE)
  Deltas <- getDelt(a=args[1], P=args[2], Tb=args[3], time=args[4], sims=args[5], invader=args[6], fig6=fig6)
  cat(".")
  return(Deltas)
}

