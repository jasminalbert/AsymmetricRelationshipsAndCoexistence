source("./diatom/DataAndVKQfuns.R")

fig_loc <- "../results_figs/"
if(dir.exists(fig_loc)==FALSE){
  dir.create(fig_loc)
}
VKQfig<- paste0(fig_loc,"VKQfig.pdf")

#t <- 1:120
#temp <- 18 + 6*sin(2*pi*t/60)
#plot(t, temp, type='l')

plotVKQ <- function(var, Tvals, dat1, dat2, dat1mod, dat2mod, fun1, fun2, fun1mod, fun2mod){
  ylab <- bquote(.(var)["i"])
  ylim <- range(c(dat1, dat2), na.rm=T)
  t <- seq(min(Tvals),max(Tvals), .1)
  par(mgp=c(2, 0.5, 0),tcl=-0.3)
  
  plot(0, xlab='', ylab=ylab, xlim=range(Tvals), ylim=ylim, xaxt='n', cex.lab=1.5, cex.axis=1.3)
  col <- "#2a2a8c" #sp1
  points(Tvals,dat1, cex=1.3, col=col)
  points(Tvals,dat1mod, pch=19, cex=0.5, col=col)
  lines(t, fun1(t), lty=2, col="#6060bf")
  lines(t, fun1mod(t), col=col)
  
  col <- "#8c2a2a" #sp2
  points(Tvals, dat2, col=col, cex=1.3)
  points(Tvals,dat2mod, pch=19, cex=0.5, col=col)
  lines(t, fun2(t), lty=2, col='#bf6060')
  lines(t, fun2mod(t), col=col)
}

pdf(VKQfig, height=7, width=5)
par(mfrow=c(3,1), mar=c(0.7, 4, 1, 1), oma=c(3, 0,0,0))
plotVKQ("V", Tvals, V1data, V2data, V1data, V2dataMod, V1fun, V2fun, V1quad, V2modfun)
mtext("(a)", 3, -1.5, at=6)
legend("bottomleft", legend=c("sp1","sp2"),  fill=c("#2a2a8c", "#8c2a2a"), inset=c(0,.3),
       bty='n', border='white')
legend("bottomleft", legend=c("original data","modified data","original function","modified function"), 
       lty=c(NA, NA, 2, 1), pch=c(1, 19, NA, NA), bty='n', pt.cex=c(1.3, .5,NA,NA))
plotVKQ("K", Tvals, K1data, K2data, K1flat, K2flat, K1fun, K2fun, K1flatfun, K2flatfun)
mtext("(b)", 3, -1.5, at=6)
plotVKQ("Q", Tvals, Q1data, Q2data, Q1data, Q2mod, Q1fun, Q2fun, Q1fun, Q2fun)
mtext("(c)", 3, -1.5, at=24,)
axis(1, at=Tvals, tick=F, cex.axis=1.3)
axis(1, at=min(Tvals):max(Tvals), labels=F)
title(ylab="K"[2])
title(outer=T, xlab="temperature (°C)", line=1.85, cex.lab=1.5)
dev.off()
##########################################################
# DO THE MODIFIED FUNCTIONS RESULT IN DIFFERENT RESULTS? #
##########################################################
#getDelt(6, 60, 18, 3000, 200)
#getDelt(6, 60, 18, 3000, 200)
