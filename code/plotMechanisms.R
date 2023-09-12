#This script contains functions:
#	makeDeltas
#	makeDeltasB
#	plotco

##libraries used (invoked with ::): graphics, grDevices

source("./decompose_LB.R")
source("./lotteryDecomp_fxns.R")

### makeDeltas function ###
# loops through sigmas to decompose into coexistence mechanisms for logNormal fecundity
#ARGS:
# noise_loc		location of noise
# sigma			vector of sigma values 
# mudif			numeric mu1-mu2
# delta			numeric delta
# qij			T/F, qij!=1?, default is FALSE
#OUT:
# a list of two dataframes: one for Deltas and one for standard errors of Delta
# where rows of dataframes are different sigma values 

#depends on decompose from lotteryDecomp_fxns.R
makeDeltas <- function(noise_loc, sigma, mudif, delta,qij=FALSE){
	load(noise_loc)
	
	message <- "computing coexistence mechanisms..."
	cat(message)
		
	#empty list to store decomposition data.frames
	store <- vector(mode='list')
		
	#iterate across variance parameter to compute decomposition (see sourced scripts)
	for (i in 1:length(sigma)){
    	store[[i]] <- decompose(mudif,sigma[i],delta,b,u)
    	#from lotteryDecomp_fxns.R
    	cat(".")
    }
    cat("done")
    
    if (qij == TRUE){
    	Deltas <- data.frame(t(sapply(store, function(X){X$Dq})))
    	se <- data.frame(t(sapply(store, function(X){X$Dq_se})))
    } else {
    	Deltas <- data.frame(t(sapply(store, function(X){X$D})))
    	se <- data.frame(t(sapply(store, function(X){X$D_se})))
    }
 
    
	colnames(se) <- colnames(Deltas) <- rownames(store[[1]])
	return(list(D=Deltas, se=se))
}
#can probably combine below and above functions

### makeDeltasB function ###
# loops through sigmas to decompose into coexistence mechanisms for Beta fecundity
#ARGS:
# noise_loc		location of noise
# sigma			vector of sigma values 
# mudif			numeric mu1-mu2
# delta			numeric delta
# qij			T/F, qij!=1?, default is FALSE
#OUT:
# a list of two dataframes: one for Deltas and one for standard errors of Delta
# where rows of dataframes are different sigma values 

#depends on decomposeB from decompose_LB.R
makeDeltasB <- function(noise_loc,lb_i,lb_j,ub_i,ub_j,delta,...){
	load(noise_loc)
	
	message <- "\ncomputing coexistence mechanisms..."
	cat(message)
		
	#empty list to store decomposition data.frames
	store <- vector(mode='list')
		
	#iterate across variance parameter to compute decomposition (see sourced scripts)
	for (j in 1:length(ub_j)){
    	store[[j]] <- decomposeB(lb_i,lb_j,ub_i,ub_j[j], delta,b,...)
    	#from decompose_LB.R
    	cat(".")
    }
    cat("done")
    Deltas <- data.frame(t(sapply(store, function(X){X$D})))
	colnames(Deltas) <- rownames(store[[1]])
	return(Deltas)
}

### plotco function ###
# loops through Delta dataframe to plot coexistence mechanisms for lognormal fec
#ARGS:
# Deltas_loc	location of Deltas
# xvar			variable on x axis, vector length nrow Deltas dataframe 
# ylim			vector length 2 or "range" 
#OUT:
# a line graph of coexistence mechanims values over a parameter range (sigma in Fig3)  

plotco <- function(Deltas_loc,xvar,ylim,...){
	
	if (file.exists(Deltas_loc)==FALSE){ #check that file exists
		stop("Deltas do not exist in directory: \ngo make Deltas using appropriate makeDeltas function")
	} else {Deltas <- readRDS(Deltas_loc)}
	
    range <- range(Deltas[,1:7])
    if (ylim[1] == "range"){
    	ylim <- range*1.1	
    }
    
    #ylim <- c(-0.6,0.2)
    
    #ATA effect
	res_vec <- xvar[Deltas$r>0 & Deltas$'[EC]' > Deltas$r]
	ex_vec <- xvar[Deltas$r<0 & Deltas$'[EC]' < Deltas$r]
	
	#plotting set up
	rescol <- grDevices::rgb(227/255, 211/255, 148/255,.5)
	excol <- grDevices::rgb(38/255, 38/255, 38/255,.5)
	wid <- c(2,2,2,2,4,4,5)
	typ <- c(1,2,3,4,1,2,1)
	ncol <- 5
	yor <- grDevices::hcl.colors(ncol, palette = "YlOrRd")
	col <- c("black","black","black","black", yor[2], yor[1], yor[3])
	
	
	##plot##
	#box	
	graphics::plot(xvar, xlab="", ylab="", ylim=ylim, xlim=range(xvar), type="n",xaxt="n", yaxt="n")
	#zero line
	graphics::lines(xvar, rep(0, length(xvar)), col="grey", lwd=0.5)
	#ATA effects
	if (length(res_vec)>0){ #rescue
		res_points <- range(res_vec)
		graphics::polygon(x=c(rep(res_points[1],2), rep(res_points[2],2)), y=c(ylim, rev(ylim))*2, xpd=FALSE, col=rescol, border=NA)
	}
	if (length(ex_vec)>0){ #exclusion
		ex_points <- range(ex_vec)
		graphics::polygon(x=c(rep(ex_points[1],2), rep(ex_points[2],2)), y=c(ylim, rev(ylim))*2, xpd=FALSE, col=excol, border=NA)
	}
	#lines
	for (m in 1:(length(Deltas)-1)){
		graphics::lines(xvar, Deltas[,m], lwd=wid[m], lty=typ[m], col=col[m], xpd=F)
	}
	
}

########
 	
### setting up legend ###
terms <- c(expression(Delta[i]^0), expression(Delta[i]^E),
           expression(Delta[i]^C), expression(Delta[i]^"(E#C)"),
           expression(Delta[i]^"[E||C]"), expression(Delta[i]^"[EC]")
           ,expression(GWR))
ncol <- 5
yor <- grDevices::hcl.colors(ncol, palette = "YlOrRd")
cols <- c("black","black","black","black",yor[1:3])
ltys <- c(1,2,3,4,2,1,1)
lwds <- c(rep(1.5,4),2,2,3)

rescol <- grDevices::rgb(227/255, 211/255, 148/255,.5)
excol <- grDevices::rgb(38/255, 38/255, 38/255,.5)  	
  	

