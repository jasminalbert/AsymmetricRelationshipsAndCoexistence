source('./diatom/fig5and6_fxns.R')
# Figure 5
parms <- c('a'=6, 'P'=60, 'Tbar'=18, 'time'=3000, 'reps'=200, 'invader'=1)
dat_loc <- "../results_numeric/fig5dat/"
if (dir.exists(dat_loc)==FALSE){
  #source("diatom/makefig5dat.R")
  a <- seq(1,6,length.out=100)
  Tbar <- seq(16,18,length.out = 100)
  P <- seq(51,199.5,length.out=100)
  dir.create(dat_loc)
  mak5(dat_loc, a_vec=a, P_vec=P, T_vec=Tbar, parms=parms)
}
fig_loc <- "../results_figs/"
if(dir.exists(fig_loc)==FALSE){
  dir.create(fig_loc)
}

Fig5 <- paste(fig_loc,"fig5.pdf",sep="")
fig5(Fig5, dat_loc, invader=1) 


#vAmp <- readRDS(paste(dat_loc,'varyingAmplitude.RDS',sep='')) #from makefig5dat.R
#vPer <- readRDS(paste(dat_loc,'varyingPeriod.RDS',sep=''))
#vTbr <- readRDS(paste(dat_loc,'varyTbar.RDS',sep=''))
#col <- c(rep("black",4), "red", "blue", "orange")
#line <- c(1:4, 1,1,1)
#lwd <- c(rep(2,6),3)
#varylist <- list(vAmp, vPer, vTbr)
#range <- sapply(varylist, function(X){range(X[,1:7])})
#ymin <- min(range)
#ymax <- max(range)
#orig <- c(a=6, P=60, Tbar=18)
#xlab <- c("a", "P", expression(theta[0]))
#default is c(5, 4, 4, 2) + 0.1
#pdf(Fig5, height=5, width=15)
#par(mfrow=c(1,3), oma=c(0,4,0,0), mar=c(5,1,2,1) )
#for (i in 1:3){
#  vdat <- varylist[[i]]
  #empty box
 # plot(0, yaxt='n', xlim=range(vdat[,10]), ylim=c(ymin, ymax), col='white', xlab="", ylab='', cex.lab=1.8, cex.axis=1.8)
#  title(xlab=xlab[i], cex.lab=2.3, line=3.5)
  #axis
 # if (i==1){ axis(2, cex.axis=1.8)}
  #loop for lines
  #for (j in 1:7){
   # lines(vdat[,10], vdat[,j], col=col[j], lty=line[j], lwd=lwd[j])	
  #}
  #abline(h=0, col="lightgrey") #zero
  #abline(v=orig[i], col='magenta', lty=3, lwd=2)
  #label
  #mtext(paste0("(", letters[i],")"), 3, -2.5, adj=0.985, cex=1.8)
  #legend
  #if (i==1){
    #legend("topleft", legend=c(expression(Delta[i]^0), expression(Delta[i]^E), expression(Delta[i]^C), expression(Delta[i]^"(E#C)"), expression(Delta[i]^"[E||C]"), expression(Delta[i]^"[EC]") ,expression(IGR)), 
     #      col = c(rep("black",4), "blue", "red", "orange"),lty = line, bty="n", cex=1.9, inset=c(-0.02,-0.03), 
      #     y.intersp = 1.1, x.intersp = 0.1, seg.len=0.8, lwd=1.8)
  #}
#}
#title(ylab="contribution to coexistence", outer=T, line=2.1, cex.lab=2, font.lab=2)
#dev.off()





