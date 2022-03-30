#This script makes plots of species species responses (V, K, and Q) to
#tempature fluctuations
#plots show original VQK point estimates from Gonzalez(2005) and 
#modified points and functions 

##libraries used (invoked with ::): 

### source VKQ functions ###
source("./diatom/DataAndVKQfuns.R")

### location to save figure ###
fig_loc <- "../results_figs/"
if(dir.exists(fig_loc)==FALSE){
  dir.create(fig_loc)
}
VKQfig<- paste0(fig_loc,"VKQfig.pdf")

### temperature ###
#t <- 1:120
#temp <- 18 + 6*sin(2*pi*t/60)
#plot(t, temp, type='l')


### function ###
### plotVKQ ###
# plot oiriginal estimates and modifications and 
# interpolated functions and modified functions
#ARGS:
	#var			single character, name of varaiable ("V","K","Q"; for axis label)
	#Tvals		vector, temperatures that species responses were evaluated 
	#dat1		vector, point values for sp 1
	#dat2		vector, point values for sp 2
	#dat1mod		vector, modified points for sp 1
	#dat2mod		vector, modified points for sp 2
	#fun1		function of Tvals, based of original data (iterpolation of points) for sp1
	#fun2		same as above for sp2
	#fun1mod		function of Tvals, based of modified data and may not be iterpolation; for sp1
	#fun2mod		same as above for sp2
plotVKQ <- function(var, Tvals, dat1, dat2, dat1mod, dat2mod, fun1, fun2, fun1mod, fun2mod){

	#yaxis
	ylab <- bquote(.(var)["i"])
 	ylim <- range(c(dat1, dat2), na.rm=T)
 	
 	#x values
 	t <- seq(min(Tvals),max(Tvals), .1)
 	
 	#plot
 	graphics::par(mgp=c(2, 0.5, 0),tcl=-0.3)
 	graphics::plot(0, xlab='', ylab=ylab, xlim=range(Tvals), ylim=ylim, xaxt='n', cex.lab=1.5, cex.axis=1.3)
 	
 	col <- "#2a2a8c" #sp1
 	graphics::points(Tvals,dat1, cex=1.3, col=col)
 	graphics::points(Tvals,dat1mod, pch=19, cex=0.5, col=col)
 	graphics::lines(t, fun1(t), lty=2, col="#6060bf")
 	graphics::lines(t, fun1mod(t), col=col)
  
 	col <- "#8c2a2a" #sp2
 	graphics::points(Tvals, dat2, col=col, cex=1.3)
 	graphics::points(Tvals,dat2mod, pch=19, cex=0.5, col=col)
 	graphics::lines(t, fun2(t), lty=2, col='#bf6060')
 	graphics::lines(t, fun2mod(t), col=col)
}


### use function for each V K and Q ###
# makes three panel figure
grDevices::pdf(VKQfig, height=7, width=5)
graphics::par(mfrow=c(3,1), mar=c(0.7, 4, 1, 1), oma=c(3, 0,0,0))

## V 
plotVKQ("V", Tvals, V1data, V2data, V1data, V2dataMod, V1fun, V2fun, V1quad, V2modfun)
graphics::mtext("(a)", 3, -1.5, at=6)

## legend
graphics::legend("bottomleft", legend=c("sp1","sp2"),  fill=c("#2a2a8c", "#8c2a2a"), inset=c(0,.3),
       bty='n', border='white')
graphics::legend("bottomleft", legend=c("original data","modified data","original function","modified function"), 
       lty=c(NA, NA, 2, 1), pch=c(1, 19, NA, NA), bty='n', pt.cex=c(1.3, .5,NA,NA))
       
## K       
plotVKQ("K", Tvals, K1data, K2data, K1flat, K2flat, K1fun, K2fun, K1flatfun, K2flatfun)
graphics::mtext("(b)", 3, -1.5, at=6)

## Q
plotVKQ("Q", Tvals, Q1data, Q2data, Q1data, Q2mod, Q1fun, Q2fun, Q1fun, Q2fun)
graphics::mtext("(c)", 3, -1.5, at=24,)

graphics::axis(1, at=Tvals, tick=F, cex.axis=1.3)
graphics::axis(1, at=min(Tvals):max(Tvals), labels=F)
graphics::title(outer=T, xlab="temperature (°C)", line=1.85, cex.lab=1.5)
grDevices::dev.off() 

