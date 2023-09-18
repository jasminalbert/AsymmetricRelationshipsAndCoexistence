##libraries used (invoked with ::): grDevices, graphics  
  #### theoryfigpanel ####
  # function that makes a panel in figure 2 of paper
  # requires functions from CopulaToolForPedagog.R
## ARGS:
# d			pdf (density function) of marginals; 
#			specified in makTheoryFigs.R
# p			cdf of marginals;
#			specified in makTheoryFigs.R
# q			inverse cdf (quantile function) of marginals;
#			specified in makTheoryFigs.R
# cop		copula object of desired shape;
#			specified in makTheoryFigs.R
# x			x values
# y 		y values
# method	"inv_cdf" (default), "gnrt_var", or "var_perse"
#			for deciding how to make plotted points of random variable
#			1. "inv_cdf" makes vars using the inverse of the desired 
#			distribution (q) of the cdf (oldp) of the 
#			random vars to be transformed (samps)  
#			2. "gnrt_var" makes new random var of desired shape
#			3. "ver_perse" transforms random var (samps) to be decoupled
# samps		a random variable for methods "inv_cdf" and "var_perse"
#			or NULL (default)
# oldp		cdf of samps; for method "inv_cdf", default NULL
#			specified in makTheoryFigs.R			
# col		T/F for coloring points of outputted random var default F
# log_points T/F log random var for plotting? default F
# plot		T/F default T
## OUT:
# new random variable 
# optionally plotted on contour plot of desired pdf  
		
theoryfigpanel <- function(d,p,q,cop,x,y, method="inv_cdf", samps=NULL, oldp=NULL,col=F,log_points=F,plot=T,...){
	
	#makepdf with some copula structure and marginals
	probdisf <- makepdf(cop,d,d,p,p) #from CopulaToolForPedagog.R
	xy <- expand.grid(x,y)
	z <- probdisf(xy)
	z <- matrix(z, length(x), length(y))
	if (plot==T){plot(x,y,type="n", xlab="",ylab="", bty="l")}
	

	if (method == "gnrt_var"){
		#make new random vars for points
		rand_var <- makerandgenrtr(cop,q,q) #CopulaToolForPedagog.R
		samps <- rand_var(1000)
	} else if (method=="inv_cdf"){
		#or use inv cdf to transform from existing rvs
		samps <- tdis(q, oldp, samps) #coded below
	} else if (method=="var_perse"){
		samps <- expand.grid(samps[,1], samps[,2])
		samps <- apply(samps,2,sample,1000)
	}
	if (plot==T){
		if (col==TRUE){
	  		ncol <- 100
	  	if (log_points==T){
	   		z[z==0]<-NA; zvar <- log10(z)
	  	} else { zvar <- z }
	  
	  	z_ran <- range(zvar, na.rm=T)
	  	zbin <- .bincode(zvar, breaks=seq(z_ran[1],z_ran[2],len=ncol+1),TRUE,TRUE)
	  	zbin <- matrix(zbin, length(x), length(y))
	  	cols <- grDevices::hcl.colors(ncol, palette = "YlOrRd",rev=T)
	  	xbin <- .bincode(samps[,1], breaks=seq(min(x),max(x),len=nrow(z)+1),TRUE,TRUE)
	  	ybin <- .bincode(samps[,2], breaks=seq(min(y),max(y),len=ncol(z)+1),TRUE,TRUE)
	  
	  	colindex <- {}
	  	for (i in 1:nrow(samps)) {colindex[i]<-zbin[xbin[i],ybin[i]]}
	  	rgb<-sapply(cols[colindex],FUN=grDevices::col2rgb)
	 	rgba<-rbind(rgb/255,.5)
    	rgbcols <- apply(rgba,2,function(X){rgb(X[1],X[2],X[3],X[4])})
	  #points(samps[,1],samps[,2], col=cols[colindex],pch=20, cex=1)
	  	graphics::points(samps[,1],samps[,2], col=rgbcols,pch=20,cex=1)
		} else {graphics::points(samps[,1],samps[,2],type="p",pch=20, cex=1, col=grDevices::rgb(.74,.74,.74,.4))}
	graphics::contour(x,y,log10(z),nlevels=10, bty='l', add=T)
	pdfhist(samps[,1],x,y,d) #coded below
	pdfhist(samps[,2],x,y,d,horiz=TRUE) #coded below
	graphics::plot.new() 
	}
	return(samps)
}

# transform distribution function
# see arguments of theoryfigpanel
tdis <- function(q, oldp, samps){
	samps <- q(oldp(samps))
}

# makes histograms of marginals
# see arguments of theoryfigpanel
pdfhist <- function(samps,x,y,d, horiz=FALSE,xpd=FALSE){
	if (sum(x==y)!=length(x)){stop("x!=y")}
	
	nbrks <- 15
	breaks <- seq(min(x), max(x), length.out=nbrks)
	
	if (min(samps) < min(x)){samps[samps==min(samps)] <- min(x)}
	h <- graphics::hist(samps, breaks, plot=FALSE)
	
	ylines <- d(x)*length(samps)* diff(h$breaks)[1]
	#xlines <- seq(0, length(h$counts), len=length(x))
	xlines <-x
	
	if(max(x)>1){xpd <- NA}

	if (horiz==TRUE){
		graphics::plot(c(0,max(h$counts)), range(x), type="n",xlab='', ylab='', main='', yaxt='n', xaxt='n', xpd=NA, bty="n")
		graphics::rect(0, h$breaks[-1], h$counts, h$breaks[-nbrks],col = rgb(.74,.74,.74,.4), border="white", lwd=0.5)
		graphics::lines(ylines,xlines,xpd=xpd)
	} else{
		graphics::plot(range(x), c(0,max(h$counts)), type="n",xlab='', ylab='', main='', yaxt='n', xaxt='n', xpd=NA, bty="n")
		graphics::rect(h$breaks[-nbrks], 0, h$breaks[-1], h$counts,col = rgb(.74,.74,.74,.4), border="white", lwd=0.5)
		#barplot(h$counts, space=0, horiz=horiz, col = rgb(.74,.74,.74,.4), border="white", xlab='', ylab='', main='', yaxt='n', xaxt='n', xpd=NA)
		graphics::lines(xlines, ylines,xpd=xpd)
		}
}
#pdfhist(sampsB[,1],dnmarg)
#pdfhist(sampsC[,1],dnmarg)
#pdfhist(sampsD[,1],dbmarg)
#pdfhist(sampsA[,1],dbmarg, xpd=NA)
#cop = ccop; d=dbmarg; p=pbmarg;q=qbmarg

