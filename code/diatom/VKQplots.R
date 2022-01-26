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
  points(Tvals,dat1, cex=1.3)
  points(Tvals,dat1mod, pch=19, cex=0.5)
  lines(t, fun1(t), lty=2)
  lines(t, fun1mod(t))

  points(Tvals, dat2, col='red', cex=1.3)
  points(Tvals,dat2mod, pch=19, cex=0.5, col='red')
  lines(t, fun2(t), lty=2, col='red')
  lines(t, fun2mod(t), col='red')
}

pdf(VKQfig, height=7, width=5)
par(mfrow=c(3,1), mar=c(0.7, 4, 1, 1), oma=c(3, 0,0,0))
plotVKQ("V", Tvals, V1data, V2data, V1data, V2dataMod, V1fun, V2fun, V1quad, V2modfun)
plotVKQ("K", Tvals, K1data, K2data, K1flat, K2flat, K1fun, K2fun, K1flatfun, K2flatfun)
plotVKQ("Q", Tvals, Q1data, Q2data, Q1data, Q2data, Q1fun, Q2fun, Q1fun, Q2fun)
axis(1, at=Tvals, tick=F, cex.axis=1.3)
axis(1, at=min(Tvals):max(Tvals), labels=F)
title(ylab="K"[2])
title(outer=T, xlab="temperature (°C)", line=1.85, cex.lab=1.5)
dev.off()