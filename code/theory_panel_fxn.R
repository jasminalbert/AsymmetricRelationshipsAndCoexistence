
theoryfigpanel <- function(d,p,q,cop,x,y, gnrt_var=FALSE, samps=NULL, oldp=NULL,...){
	probdisf <- makepdf(cop,d,d,p,p)
	xy <- expand.grid(x,y)
	z <- probdisf(xy)
	z <- matrix(z, length(x), length(y))
	
	#layout(matrix(c(2,4,1,3), byrow=T, nrow=2), heights=c(.3,1), widths=c(1,.3))
	#par(mar=c(0.5,0.5,0.5,0.5), oma=c(2,2,1,1), xpd=NA)
	levs<-c(-.5,0,1,1.5)
	contour(x,y,log10(z),xlim=range(x),ylim=range(y),nlevels=10, bty='l')
	if (gnrt_var == TRUE){
		rand_var <- makerandgenrtr(cop,q,q)
		samps <- rand_var(500)
	} else{
		samps <- tdis(q, oldp, samps)
	}
	points(samps[,1],samps[,2],type="p",pch=20, cex=1, col=rgb(.74,.74,.74,.4))
	pdfhist(samps[,1],d,...)
	pdfhist(samps[,2],d,horiz=TRUE,...)
	plot.new() 
	return(samps)
}

tdis <- function(q, oldp, samps){
	samps <- q(oldp(samps))
}

pdfhist <- function(samps, d, horiz=FALSE,xpd=FALSE){
	min <- min(x)
	max <- max(y)
	h <- hist(samps, seq(min(samps),max(samps),length.out=10), plot=FALSE)
	ylines <- d(x)*length(samps)* diff(h$breaks)[1]
	xent <- length(h$counts)
	barplot(h$counts, space=0, horiz=horiz, col = rgb(.74,.74,.74,.4), border="white", xlab='', ylab='', main='', yaxt='n', xaxt='n')

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

