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

wrapDelt <- function(args){
  Deltas <- getDelt(a=args[1], P=args[2], Tb=args[3], time=args[4], sims=args[5], invader=args[6])
  return(Deltas)
}

require(parallel)

mapspace <- function(parmlist, sims, time, invader){
	
	#names(parmlist) <- c('a','P','Tbar')
	#len <- unlist(lapply(parmlist, length))
	#varlen <- len[len!=1]
	#rows <- varlen[[1]]*varlen[[2]]
	#res <- matrix(NA, ncol=6, nrow=rows)
	
	a <- parmlist[[1]]
	P <- parmlist[[2]]
	Tb <- parmlist[[3]]
	
	argsList <- list()
	r <- 1
	for (i in 1:length(a)){
		for (j in 1:length(P)){
			for (k in 1:length(Tb)){
				#Delta1 <- getDelt(a[i], P[j], Tb[k], time, sims)
				#res[r,] <- c(a[i], P[j], Tb[k], Delta1$IGR, Delta1$epsECbrk, Delta1$time)
				#print(c(r,i,j,k))
			  argsList[[r]] <- c(a[i], P[j], Tb[k], time, sims, invader)
				r <- r+1
			}	
		}
	}
	#res <- data.frame(res)
	resList <- mclapply(argsList, wrapDelt, mc.cores=7) ###change cores here###
	resdf <- data.frame(t(sapply(resList, function(X){c(X$IGR, X$epsECbrk, X$time, X$map)})))
	parmsdf <- data.frame(matrix(unlist(argsList), ncol=6, byrow=T)[,-4:-6])
	
	res <- cbind(parmsdf, resdf)
	colnames(res) <- c("a", "P", "Tbar", "IGR", "ATA", "time", "map")
	return(res)
}

require(plot3D)

cmMap <- function(dat){
	dat$map[dat$map>1] <- rep(1, length(dat$map[dat$map>1]))
	dat$map[dat$map<(-1)] <- rep(-1, length(dat$map[dat$map<(-1)]))
	
	a <- unique(dat$a); P <- unique(dat$P); Tbar <- unique(dat$Tbar)
	dims <- list(a=round(a,3), P=round(P,3), Tbar=round(Tbar,3))
	
	len <- unlist(lapply(dims, length))
	vars <- dims[len!=1]
	
	map <- matrix(dat$map, nrow=length(vars[[1]]), ncol=length(vars[[2]]), dimnames=vars, byrow=T)
	
	return(map)
}

cmContour <- function(map, ncolor=51, colkey=NULL,...){
	x <- dimnames(map)[2]
	y <- dimnames(map)[1]
	#main <- expression(paste(Delta[i]^"[EC]","/IGR"))
	cm <- cm.colors(ncolor)
	
	if(all(range(map)!=c(-1,1))){
	  by <- 2/(ncolor-1)
	  fullr <- seq(-1,1,by)
	  fullr <- round(fullr,2)
	  
	  roundby <- function(x, by){
	    m <- round(x)
	    if (m<x){
	      i <- seq(m, m+1, by)
	    } else {i <- seq(m, m-1, by)}
	    
	    dd <- tail(i[x>i],1)
	    du <- i[x<i][1]
	    if (x-dd < du-x){
	      x <- dd
	    } else {x <- du}
	    return(x)
	  }
	  range <- range(map)
	  crange <- c(roundby(range[1],by), roundby(range[2],by))
	  cm <- cm[match(crange,fullr)[1]:match(crange,fullr)[2]]
	}
	
	image2D(z=t(map), y=as.numeric(y[[1]]), x=as.numeric(x[[1]]), contour=FALSE, col=cm, colkey=colkey, xlab='', ylab='',...)
	title(xlab=names(x), ylab=ifelse(names(y)=='Tbar', expression(theta[0]), names(y)), line=-1)
	contour(z=t(map), y=as.numeric(y[[1]]), x=as.numeric(x[[1]]),add=TRUE, col='grey50')
}


