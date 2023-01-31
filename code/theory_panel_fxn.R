# d		pdf (density function) of marginals
# p		cdf of marginals
# q		inverse cdf (quantile function) of marginals
# cop	copula object of desired shape
# x		x values
# y 	y values
		
theoryfigpanel <- function(d,p,q,cop,x,y, gnrt_var=FALSE, samps=NULL, oldp=NULL,col=F,log_points=F,...){
	
	#makepdf with some copula structure and marginals
	probdisf <- makepdf(cop,d,d,p,p) 
	xy <- expand.grid(x,y)
	z <- probdisf(xy)
	z <- matrix(z, length(x), length(y))
	
	
	plot(x,y,type="n", xlab="",ylab="", bty="l")
	if (gnrt_var == TRUE){
		#make new random vars for points
		rand_var <- makerandgenrtr(cop,q,q)
		samps <- rand_var(1000)
	} else{
		#or use inv cdf to transform from existing rvs
		samps <- tdis(q, oldp, samps)
	}
	if (col==TRUE){
	  ncol <- 100
	  if (log_points==T){
	    z[z==0]<-NA; zvar <- log10(z)
	  } else { zvar <- z }
	  
	  z_ran <- range(zvar, na.rm=T)
	  zbin <- .bincode(zvar, breaks=seq(z_ran[1],z_ran[2],len=ncol+1),TRUE,TRUE)
	  zbin <- matrix(zbin, length(x), length(y))
	  cols <- hcl.colors(ncol, palette = "YlOrRd",rev=T)
	  xbin <- .bincode(samps[,1], breaks=seq(min(x),max(x),len=nrow(z)+1),TRUE,TRUE)
	  ybin <- .bincode(samps[,2], breaks=seq(min(y),max(y),len=ncol(z)+1),TRUE,TRUE)
	  
	  colindex <- {}
	  for (i in 1:nrow(samps)) {colindex[i]<-zbin[xbin[i],ybin[i]]}
	  rgb<-sapply(cols[colindex],FUN=col2rgb)
	  rgba<-rbind(rgb/255,.5)
    rgbcols <- apply(rgba,2,function(X){rgb(X[1],X[2],X[3],X[4])})
	  #points(samps[,1],samps[,2], col=cols[colindex],pch=20, cex=1)
	  points(samps[,1],samps[,2],col=rgbcols,pch=20,cex=1)
	} else {
	  points(samps[,1],samps[,2],type="p",pch=20, cex=1, col=rgb(.74,.74,.74,.4))
	}
	contour(x,y,log10(z),nlevels=10, bty='l', add=T)
	pdfhist(samps[,1],x,y,d)
	pdfhist(samps[,2],x,y,d,horiz=TRUE)
	plot.new() 
	return(samps)
}


tdis <- function(q, oldp, samps){
	samps <- q(oldp(samps))
}

pdfhist <- function(samps,x,y,d, horiz=FALSE,xpd=FALSE){
	min <- min(x)
	max <- max(y)
	h <- hist(samps, seq(min(samps),max(samps),length.out=10), plot=FALSE)
	ylines <- d(x)*length(samps)* diff(h$breaks)[1]
	
	if( max(x)>1){ xent <- 1; xpd <- NA;
		} else{xent <- length(h$counts)}
	barplot(h$counts, space=0, horiz=horiz, col = rgb(.74,.74,.74,.4), border="white", xlab='', ylab='', main='', yaxt='n', xaxt='n', xpd=NA)

	if (horiz==TRUE){
		lines(ylines,(x-min(x))*xent,xpd=xpd)
	} else{
		lines((x-min(x))*xent, ylines,xpd=xpd)
		}
}
#pdfhist(sampsB[,1],dnmarg)
#pdfhist(sampsC[,1],dnmarg)
#pdfhist(sampsD[,1],dbmarg)
#pdfhist(sampsA[,1],dbmarg, xpd=NA)
#cop = ccop; d=dbmarg; p=pbmarg;q=qbmarg

