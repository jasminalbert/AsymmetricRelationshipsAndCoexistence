source("DataAndVKQfuns.R")
source("partialSharp_fxns.R")
require(deSolve)

#make into function
#ARGS
#E #C #a #P #reps 
#run model and get Ej and Cj first
## !!only for invader is species 1 (Frag) right now!!

get_epsilons <- function(a, P, Emethod, reps, sp, invader=1){
	
	parms <- c(Tbar=18, a=a, P=P, D=0.09, S=35)
	
	y0 <- c(R=.1,x1=0,x2=20)
	times <- seq(0,3600,by=.1)
	out1i <- ode(y0,times,func=forceChemo,parms=parms) 
	e <- times > max(times-1200) 
	out1i <- out1i[e,] 
	temp <- out1i[,5]
	R <- out1i[,2]
	
	
	if (sp==2){
		Vfun <- V2modfun; Kfun <- K2flatfun
	} else {Vfun <- V1quad; Kfun <- K1flatfun}

	
	if (Emethod=="temp"){
		E <- temp
		r <- function(E, C, parms){Vfun(E)/C - parms["D"]}
	}
	if (Emethod=="Voftemp"){
		E <- Vfun(temp)
		r <- function(E, C, parms){E/C - parms["D"]}
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
	
	epsEC <- rbar - epsE - epsC - eps0
	
	epsEsharpC <- rsharp - epsE - epsC - eps0
	
	epsECpar <- rbar - rsharp
	
	epsECbrk <- rbar - rpsharp
	
	epsEpsharpC <- rpsharp - rsharp

	
	return(data.frame(eps0=eps0, epsE=epsE, epsC=epsC, epsEC=epsEC, epsEsharpC=epsEsharpC, epsECpar=epsECpar, epsECbrk=epsECbrk, epsEpsharpC=epsEpsharpC))
}

#test
#ep1<- get_epsilons(6, 60, "Voftemp", 500, 1)
#ep2<- get_epsilons(6, 60, "Voftemp", 500, 2)
#Delta1 <- ep1-ep2
#Delta1$eps0+Delta1$epsE+Delta1$epsC+Delta1$epsEsharpC+Delta1$epsECbrk+Delta1$epsEpsharpC

#ep1<- get_epsilons(6, 60, "temp", 500, 1)
#ep2<- get_epsilons(6, 60, "temp", 500, 2)
#Delta1t <- ep1-ep2
#Delta1t$eps0+Delta1t$epsE+Delta1t$epsC+Delta1t$epsEsharpC+Delta1t$epsECbrk+Delta1t$epsEpsharpC



getep <- function(a, P, reps, sp, invader=1){
	
	parms <- c(Tbar=18, a=a, P=P, D=0.09, S=35)
	
	y0 <- c(R=.1,x1=0,x2=20)
	times <- seq(0,2000,by=.1)
	out1i <- ode(y0,times,func=forceChemo,parms=parms) 
	e <- times > max(times-660) 
	out1i <- out1i[e,] 
	temp <- out1i[,5]
	R <- out1i[,2]
	
	
	if (sp==2){
		Vfun <- V2modfun; Kfun <- K2flatfun
	} else {Vfun <- V1quad; Kfun <- K1flatfun}

	E <- temp
	r <- function(E, C, parms){Vfun(E)/C - parms["D"]}
	
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

	
	return(data.frame(eps0=eps0, epsE=epsE, epsC=epsC, epsEsharpC=epsEsharpC, epsECbrk=epsECbrk, epsEpsharpC=epsEpsharpC))	
	
}

getep2 <- function(a, P, Tbar, reps, sp, invader=1){
	
	parms <- c(Tbar=Tbar, a=a, P=P, D=0.09, S=35)
	
	y0 <- c(R=.1,x1=0,x2=20)
	times <- seq(0,3000,by=.1)
	out1i <- ode(y0,times,func=forceChemo,parms=parms) 
	e <- times > max(times-1000) 
	out1i <- out1i[e,] 
	temp <- out1i[,5]
	R <- out1i[,2]
	
	
	if (sp==2){
		Vfun <- V2modfun; Kfun <- K2flatfun
	} else {Vfun <- V1quad; Kfun <- K1flatfun}

	E <- temp
	r <- function(E, C, parms){Vfun(E)/C - parms["D"]}
	
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

	
	return(data.frame(eps0=eps0, epsE=epsE, epsC=epsC, epsEsharpC=epsEsharpC, epsECbrk=epsECbrk, epsEpsharpC=epsEpsharpC))
	
}

getep3 <- function(a, P, Tbar, time=3600, reps, sp, invader=1){
  
  parms <- c(Tbar=Tbar, a=a, P=P, D=0.09, S=35)
  
  y0 <- c(R=.1,x1=0,x2=20)
  times <- seq(0,time,by=.1)
  out1i <- ode(y0,times,func=forceChemo,parms=parms) 
  e <- times > max(times-(signif(floor(max(time)/3),2))) 
  out1i <- out1i[e,] 
  temp <- out1i[,5]
  R <- out1i[,2]
  
  
  if (sp==2){
    Vfun <- V2modfun; Kfun <- K2flatfun
  } else {Vfun <- V1quad; Kfun <- K1flatfun}
  
  E <- temp
  r <- function(E, C, parms){Vfun(E)/C - parms["D"]}
  
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
  
  
  return(data.frame(eps0=eps0, epsE=epsE, epsC=epsC, epsEsharpC=epsEsharpC, epsECbrk=epsECbrk, epsEpsharpC=epsEpsharpC))
  
}