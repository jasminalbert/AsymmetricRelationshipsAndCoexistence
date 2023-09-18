#This script makes figure 3 in main text and S1 in supplement: 
#plotting LOGNORMAL FECUNDITIES LOTTERY MODEL
#both qij=1 (main) and qij!=1 (supplement)
#decomposed mechanisms contributions to coexistence
#against sigma values
#for different deltas and mu1-mu2

##libraries used (invoked with ::): graphics, grDevices

### source function ###
#lottery model decomposition and plotting
source("./plotMechanisms.R")
	#contains functions: makeDeltas, plotco


### location to save results ###
numRes_loc <- "../results_numeric/"
if(dir.exists(numRes_loc)==FALSE){
  dir.create(numRes_loc)
}
fig_loc <- "../results_figs/"
if(dir.exists(fig_loc)==FALSE){
  dir.create(fig_loc)
}
fig3_loc <- paste(fig_loc,"decompFigLLN.pdf",sep="")
fig3qij_loc <- paste(fig_loc,"decompFigLLNqij.pdf",sep="")

### files to load ###
noise_loc <- paste0(numRes_loc, "noise.RData")
params_loc <- paste0(numRes_loc, "params.RData")
res_loc <- paste0(numRes_loc,"fig3dat/")
if(dir.exists(res_loc)==FALSE){
  dir.create(res_loc)
}
fig3_dat_loc <- paste0(res_loc,"Deltas_d")

load(params_loc)
sigma <- seq(0,7,len=50)
delta <- delta[-4]
### setting dimensions ###
nd <- length(delta) 
nm <- length(mudif)
NP <- nm*nd #number of panels
ph <- 3 #panel height
pw <- 2.5 #panel widht
ht <- c(1,rep(ph, nm)) #vector of panel heights
wd <- c(rep(pw, nd),2) #vector of panel widths
midx <- (pw*(nd/2))/sum(wd) #middle of figure panels; horizontal
midy <- (ph*(nm/2))/sum(ht) #middle of figure panels; vertical 

### setting up legend ###
terms <- c(expression(Delta[i]^0), expression(Delta[i]^E),
           expression(Delta[i]^C), expression(Delta[i]^"(E#C)"),
           expression(Delta[i]^"[E||C]"), expression(Delta[i]^"[EC]")
           ,expression(GWR))
ncol <- 5
yor <- hcl.colors(ncol, palette = "YlOrRd")
cols <- c("black","black","black","black",yor[1:3])
ltys <- c(1,2,3,4,2,1,1)
lwds <- c(rep(1.5,4),2,2,3)

### make fig ###
#qij=1
grDevices::pdf(fig3_loc, height=8, width=8)
graphics::par(mfrow=c(4,3), oma=c(3.2,4,4,10), mar=c(1,1,1,1), mgp=c(3,1,0), xpd=NA)
ylims <- list(c(-0.5,1.25),c(-0.5,1.25),c(-2,1.5), c(-4,2.5))
maxse <- {}

n <- 1
for (m in seq_along(mudif)){ #iterates across mudif values
  for (d in seq_along(delta)){ #iterates across delta values
  	
  	panel <- letters[n]
  	cat(n"/",NP,". panel (",panel,"): mu_1-mu_2=", mudif[m],";delta=",delta[d],"\n")
  	
  	file <- paste0(fig3_dat_loc,delta[d],"_md",mudif[m],".RDS")
  	sefile <- paste0(fig3_dat_loc,delta[d],"_md",mudif[m],"_SE.RDS")
  	#make files?
  	if(file.exists(file)==FALSE){
  		Delts <- makeDeltas(noise_loc, sigma, mudif[m], delta[d], qij=FALSE)
  		#from plotMechanisms.R
  		saveRDS(Delts$D, file=file)
  		saveRDS(Delts$se, file=sefile)
  	}
  	#standard error #need to get rid of NA col
  	maxse[n] <- max(readRDS(sefile),na.rm=T)
  	
  	#plot (from plotMechanisms.R)
  	plotco(file,sigma,ylim=ylims[[m]])
  	#axes
  	if (d==1){
  		graphics::axis(2, cex.axis=1.8, tck=-0.035, lwd.ticks=2)} else {graphics::axis(2, cex.axis=1.8, tck=-0.035, lwd.ticks=2, labels=F)}
  	if (n>9){
  		graphics::axis(1, cex.axis=1.8,tck=-0.028, lwd.ticks=2)} else {graphics::axis(1, cex.axis=1.8, tck=-0.028, lwd.ticks=2, labels=F)}
  	#panel labels
    graphics::mtext(paste0("(", panel,")"), side=3, line=-1.7, at=0, cex=1.3, adj=0)
    
    if (n==2){
    	graphics::mtext("adult death rate, ", line=2.5, font=2, cex=1.5, col="gray40", at=-1, adj=0)
    	graphics::mtext(expression(delta), line=2.5, font=2, cex=2.3, col="gray40", side=3, at=8.25, adj=1)
    }
    if (n<4){
      graphics::mtext(paste(delta[d]), side=3, line=0.2, font=2, cex=1.25, col="gray30")
    }
    if (n%%3==0){
      graphics::text(paste(mudif[m]), srt=-90, x=7.8,y=mean(ylims[[m]]),  font=2, cex=2, col="gray30")
    }
    if (n==6){
    	graphics::text("mean log fecundity difference, ", x=9,y=2.45, font=2, cex=2.5, srt=-90, adj=0, col="gray40")
    	graphics::text(expression(mu[1]-mu[2]), x=9,y=-3.95, font=2, cex=3, srt=-90, adj=1, col="gray40")
    }
    if (n==11){
    	graphics::mtext("log fecundity standard deviation, ",side=1, line=2.8, cex=1.25, at=-3.85, adj=0,col="gray40", font=2)
    	graphics::mtext(expression(sigma),side=1, line=2.8, cex=1.5, at=10.75, adj=1,col="gray40", font=2)
    }
    n <- n+1
  }
}
graphics::mtext("contribution to coexistence", side=2, outer=TRUE, line=2, font=2, cex=1.7, col="gray40")
#legend
graphics::legend("topright", legend=terms, col = cols, lty = ltys, bty="n", cex=2, inset=c(-0.83,-2.75), y.intersp = 1.35, x.intersp = 0.1, seg.len=0.8, lwd=lwds)
graphics::legend("topright",legend=c("ATA \nexcl.", "ATA \nresc."), fill=c(excol,rescol), cex=1.8, bty="n",border=NA, inset=c(-0.76,-0.8), x.intersp = 0.1,y.intersp = 1.7)
grDevices::dev.off() #finish plotting

### get standard error and save ###
fig3maxse <- max(maxse, na.rm=T)
cat("maximum standard error in figure three is", fig3maxse, "\n(M=", M, ")\n")
fig3maxse_loc <- paste(numRes_loc, "fig3maxse.RDS", sep="")
saveRDS(fig3maxse, file=fig3maxse_loc)
###########################################

 ###qij!=1###
res_loc <- paste0(numRes_loc,"fig3qijdat/")
if(dir.exists(res_loc)==FALSE){
  dir.create(res_loc)
}
fig3_dat_loc <- paste0(res_loc,"qDeltas_d")

###make fig###
grDevices::pdf(fig3qij_loc, height=8, width=8)
graphics::par(mfrow=c(4,3), oma=c(3.2,4,4,10), mar=c(1,1,1,1), mgp=c(3,1,0), xpd=NA)
ylims <- list(c(-0.5,1.25),c(-0.5,1.25),c(-2,1.5), c(-4,2.5))
maxse <- {}

n <- 1
for (m in seq_along(mudif)){ #iterates across mudif values
  for (d in seq_along(delta)){ #iterates across delta values
  	
  	file <- paste0(fig3_dat_loc,delta[d],"_md",mudif[m],".RDS")
  	sefile <- paste0(fig3_dat_loc,delta[d],"_md",mudif[m],"_SE.RDS")
  	#make files?
  	if(file.exists(file)==FALSE){
  		Delts <- makeDeltas(noise_loc, sigma, mudif[m], delta[d], qij=TRUE)
  		#from plotMechanims.R
  		saveRDS(Delts$D, file=file)
  		saveRDS(Delts$se, file=sefile)
  	}
  	#standard error #need to get rid of NA col
  	maxse[n] <- max(readRDS(sefile), na.rm=T)
  	
  	#plot (from plotMechanism.R)
  	plotco(file,sigma,ylim=ylims[[m]])
  	#axes
  	if (d==1){
  		graphics::axis(2, cex.axis=1.8, tck=-0.035, lwd.ticks=2)} else {graphics::axis(2, cex.axis=1.8, tck=-0.035, lwd.ticks=2, labels=F)}
  	if (n>9){
  		graphics::axis(1, cex.axis=1.8,tck=-0.028, lwd.ticks=2)} else {graphics::axis(1, cex.axis=1.8, tck=-0.028, lwd.ticks=2, labels=F)}
	#panel labels
    graphics::mtext(paste0("(", letters[n],")"), side=3, line=-1.7, at=0, cex=1.3, adj=0)
    
    if (n==2){
    	graphics::mtext("adult death rate, ", line=2.5, font=2, cex=1.5, col="gray40", at=-1, adj=0)
    	graphics::mtext(expression(delta), line=2.5, font=2, cex=2.3, col="gray40", side=3, at=8.25, adj=1)
    }
    if (n<4){
      graphics::mtext(paste(delta[d]), side=3, line=0.2, font=2, cex=1.25, col="gray30")
    }
    if (n%%3==0){
      graphics::text(paste(mudif[m]), srt=-90, x=7.8,y=mean(ylims[[m]]),  font=2, cex=2, col="gray30")
    }
    if (n==6){
    	graphics::text("mean log fecundity difference, ", x=9,y=2.45, font=2, cex=2.5, srt=-90, adj=0, col="gray40")
    	graphics::text(expression(mu[1]-mu[2]), x=9,y=-3.95, font=2, cex=3, srt=-90, adj=1, col="gray40")
    }
    if (n==11){
    	graphics::mtext("log fecundity standard deviation, ",side=1, line=2.8, cex=1.25, at=-3.25, adj=0,col="gray40")
    	graphics::mtext(expression(sigma),side=1, line=2.8, cex=1.5, at=10.15, adj=1,col="gray40")
    }
    n <- n+1
  }
} 
graphics::mtext("contribution to coexistence", side=2, outer=TRUE, line=2, font=2, cex=1.7, col="gray40")
graphics::legend("topright", legend=terms, col = cols, lty = ltys, bty="n", cex=2, inset=c(-0.83,-2.75), y.intersp = 1.35, x.intersp = 0.1, seg.len=0.8, lwd=lwds)
graphics::legend("topright",legend=c("ATA \nexcl.", "ATA \nresc."), fill=c(excol,rescol), cex=1.8, bty="n",border=NA, inset=c(-0.76,-0.8), x.intersp = 0.1,y.intersp = 1.7)
grDevices::dev.off() #finish plotting

### get standard error and save ###
fig3qijmaxse <- max(maxse, na.rm=T)
cat("maximum standard error in figure three (qij!=1) is", fig3qijmaxse, "\n(M=", M, ")\n")
fig3qijmaxse_loc <- paste(numRes_loc, "fig3qijmaxse.RDS", sep="")
saveRDS(fig3qijmaxse, file=fig3qijmaxse_loc)

