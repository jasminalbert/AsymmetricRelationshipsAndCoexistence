setwd("~/Documents/AsymmetricRelationshipsAndCoexistence/code")
source("diatom/DataAndVKQfuns.R")
source("diatom/partialSharp_fxns.R")
require(deSolve)


getep3 <- function(a, P, Tbar, time, reps, sp, invader=1){
	
	parms <- c(Tbar=Tbar, a=a, P=P, D=0.09, S=35)
	time <- round(time/P,0)*P
	
	y0 <- c(R=.1,x1=0,x2=20)
	times <- seq(0,time,by=.1)
	out1i <- ode(y0,times,func=forceChemo,parms=parms) 
	e <- times > max(times-signif(floor(max(times)/2.5),2)) 
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
	
	return(data.frame(eps0=eps0, epsE=epsE, epsC=epsC, epsEsharpC=epsEsharpC, epsECbrk=epsECbrk, epsEpsharpC=epsEpsharpC, time=time))
	
}

getDelt <- function(a, P, Tbar, time, sims){
	ep1 <- getep3(a, P, Tbar, time, reps=sims, sp=1)
	ep2 <- getep3(a, P, Tbar, time, reps=sims, sp=2)
	Delta1 <- ep1-ep2
	Delta1$IGR <- sum(Delta1[,-7])
	Delta1$map <- Delta1$epsECbrk/Delta1$IGR
	Delta1$time <- ep1$time
	return(Delta1)
}

mapspace <- function(parmlist, sims, time){
	
	names(parmlist) <- c('a','P','Tbar')
	len <- unlist(lapply(parmlist, length))
	varlen <- len[len!=1]
	
	rows <- varlen[[1]]*varlen[[2]]
	res <- matrix(NA, ncol=6, nrow=rows)
	r <- 1
	
	a <- parmlist[[1]]
	P <- parmlist[[2]]
	Tb <- parmlist[[3]]
	
	for (i in 1:length(a)){
		for (j in 1:length(P)){
			for (k in 1:length(Tb)){
				Delta1 <- getDelt(a[i], P[j], Tb[k], time, sims)
				
				res[r,] <- c(a[i], P[j], Tb[k], Delta1$IGR, Delta1$epsECbrk, Delta1$time)
				print(c(r,i,j,k))
				r <- r+1
			}	
		}
	}
	res <- data.frame(res)
	colnames(res) <- c("a", "P", "Tbar", "IGR", "ATA", "time")
	res$eff <- res$IGR - res$ATA
	res$map <- res$ATA/res$IGR
	
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
	
	image2D(z=t(map), y=as.numeric(y[[1]]), x=as.numeric(x[[1]]), contour=FALSE, col=cm, colkey=colkey, xlab='', ylab='',...)
	title(xlab=names(x), ylab=names(y), line=-1)
	contour(z=t(map), y=as.numeric(y[[1]]), x=as.numeric(x[[1]]),add=TRUE, col='grey50')
}


