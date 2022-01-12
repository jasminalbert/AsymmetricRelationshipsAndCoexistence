# Figure 5
source("diatom/workflow.R")
dat_loc <- "../results_numeric/fig5dat/"
if (dir.exists(dat_loc)==FALSE){
  source("diatom/makefig5dat.R")
}
fig_loc <- "../results_figs/"
if(dir.exists(fig_loc)==FALSE){
  dir.create(fig_loc)
}
fig5 <- paste(fig_loc,"fig5.pdf",sep="")

vAmp <- readRDS(paste(dat_loc,'varyingAmplitude.RDS',sep='')) #from makefig5dat.R
vPer <- readRDS(paste(dat_loc,'varyingPeriod.RDS',sep=''))
vTbr <- readRDS(paste(dat_loc,'varyTbar.RDS',sep=''))


col <- c(rep("black",4), "red", "blue", "orange")
line <- c(1:4, 1,1,1)

varylist <- list(vAmp, vPer, vTbr)
range <- sapply(varylist, function(X){range(X[,1:7])})
ymin <- min(range)
ymax <- max(range)
orig <- c(a=6, P=60, Tbar=18)

#default is c(5, 4, 4, 2) + 0.1
pdf(fig5, height=5, width=15)
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





