source("./CopulaToolForPedagogFig.R")

#this figure will:
	# have four panels A-D
	#A) E and C with ATA
	#B) E and C transformed to standard normal but still has ATA
	#C) remove ATA by replacing [E,C] with biv normal with same parms
	#D) invert marginals to [E||,C||]
	

#E and C have uniform marginals and clayton copula

#***make an example with beta marginals and normal copula, and then make a contour
shape1<-0.25
shape2<-0.25
dbmarg<-function(x){return(dbeta(x+.5,shape1=shape1,shape2=shape2))}
pbmarg<-function(x){return(pbeta(x+.5,shape1=shape1,shape2=shape2))}
pdf_claycop_betamargs<-makepdf(ccop,dbmarg,dbmarg,pbmarg,pbmarg)
x<-seq(from=-0.51,to=0.49,by=0.01)
y<-seq(from=-0.51,to=0.49,by=0.01)
xy<-expand.grid(x,y)
z<-pdf_claycop_betamargs(xy)
z<-matrix(z,length(x),length(y))


pdf("panelA_0.25.pdf") 
layout(matrix(c(2,4,1,3), byrow=T, nrow=2), heights=c(.3,1), widths=c(1,.3))
par(mar=c(0.5,0.5,0.5,0.5), oma=c(2,2,1,1), xpd=NA)
levs<-c(-.5,0,1,1.5)
#contour(x,y,log10(z),xlim=range(x),ylim=range(y),levels=levs,labels=levs)
contour(x,y,log10(z),xlim=range(x),ylim=range(y),nlevels=10, bty='l')
#marginals
plot(x,dbmarg(x),type="l", bty='n', xaxt='n', xlab="", yaxt="n") #one marginal
plot(dbmarg(x),x,type="l", bty='n', yaxt='n', ylab="", xaxt="n") #the other (which is the same in this case)
plot.new() 
#end
dev.off()




shape1<-0.5
shape2<-0.5

pdf_claycop_betamargs<-makepdf(ccop,dbmarg,dbmarg,pbmarg,pbmarg)
x<-seq(from=-0.51,to=0.49,by=0.01)
y<-seq(from=-0.51,to=0.49,by=0.01)
xy<-expand.grid(x,y)
z<-pdf_claycop_betamargs(xy)
z<-matrix(z,length(x),length(y))

pdf("panelA_0.5.pdf") 
layout(matrix(c(2,4,1,3), byrow=T, nrow=2), heights=c(.3,1), widths=c(1,.3))
par(mar=c(0.5,0.5,0.5,0.5), oma=c(2,2,1,1), xpd=NA)
levs<-c(-.5,0,1,1.5)
contour(x,y,log10(z),xlim=range(x),ylim=range(y),nlevels=10, bty='l')
#marginals
plot(x,dbmarg(x),type="l", bty='n', xaxt='n', xlab="", yaxt="n") #one marginal
plot(dbmarg(x),x,type="l", bty='n', yaxt='n', ylab="", xaxt="n") #the other (which is the same in this case)
plot.new() 
dev.off()

pdf("panelA_0.5_filled.pdf") #issue - filled.con not going into layout
layout(matrix(c(2,4,1,3), byrow=T, nrow=2), heights=c(.3,1), widths=c(1,.3))
par(mar=c(0.5,0.5,0.5,0.5), oma=c(2,2,1,1), xpd=NA)
levs<-c(-.5,0,1,1.5)
filled.contour(x,y,log10(z),xlim=range(x),ylim=range(y),nlevels=40, bty='l')
#marginals
plot(x,dbmarg(x),type="l", bty='n', xaxt='n', xlab="", yaxt="n") #one marginal
plot(dbmarg(x),x,type="l", bty='n', yaxt='n', ylab="", xaxt="n") #the other (which is the same in this case)
plot.new() 
dev.off()

#E and C with normal marginals and clayton copula
dnmarg<-function(x){return(dnorm(x,mean=0,sd=1))}
pnmarg<-function(x){return(pnorm(x,mean=0,sd=1))}
pdf_claycop_normmargs<-makepdf(ccop,dnmarg,dnmarg,pnmarg,pnmarg)
x<-seq(from=-0.51,to=0.49,by=0.01)
y<-seq(from=-0.51,to=0.49,by=0.01)
xy<-expand.grid(x,y)
z<-pdf_claycop_normmargs(xy)
z<-matrix(z,length(x),length(y))

pdf("panelB.pdf")
layout(matrix(c(2,4,1,3), byrow=T, nrow=2), heights=c(.3,1), widths=c(1,.3))
par(mar=c(0.5,0.5,0.5,0.5), oma=c(2,2,1,1), xpd=NA)
levs<-c(-.5,0,1,1.5)
contour(x,y,log10(z),xlim=range(x),ylim=range(y),nlevels=10, bty='l')
#marginals
plot(x,dnmarg(x),type="l", bty='n', xaxt='n', xlab="", yaxt="n") #one marginal
plot(dnmarg(x),x,type="l", bty='n', yaxt='n', ylab="", xaxt="n") #the other (which is the same in this case)
plot.new() 
#end
dev.off()

#normal cop, normal marg
pdf_normcop_normmargs<-makepdf(ncop,dnmarg,dnmarg,pnmarg,pnmarg)
x<-seq(from=-0.51,to=0.49,by=0.01)
y<-seq(from=-0.51,to=0.49,by=0.01)
xy<-expand.grid(x,y)
z<-pdf_normcop_normmargs(xy)
z<-matrix(z,length(x),length(y))

pdf("panelC.pdf")
layout(matrix(c(2,4,1,3), byrow=T, nrow=2), heights=c(.3,1), widths=c(1,.3))
par(mar=c(0.5,0.5,0.5,0.5), oma=c(2,2,1,1), xpd=NA)
levs<-c(-.5,0,1,1.5)
contour(x,y,log10(z),xlim=range(x),ylim=range(y),nlevels=10, bty='l')
#marginals
plot(x,dnmarg(x),type="l", bty='n', xaxt='n', xlab="", yaxt="n") #one marginal
plot(dnmarg(x),x,type="l", bty='n', yaxt='n', ylab="", xaxt="n") #the other (which is the same in this case)
plot.new() 
#end
dev.off()

#normal cop, beta marg
pdf_normcop_betamargs<-makepdf(ncop,dbmarg,dbmarg,pbmarg,pbmarg)
x<-seq(from=-0.51,to=0.49,by=0.01)
y<-seq(from=-0.51,to=0.49,by=0.01)
xy<-expand.grid(x,y)
z<-pdf_normcop_betamargs(xy)
z<-matrix(z,length(x),length(y))

pdf("panelD_0.5.pdf")
layout(matrix(c(2,4,1,3), byrow=T, nrow=2), heights=c(.3,1), widths=c(1,.3))
par(mar=c(0.5,0.5,0.5,0.5), oma=c(2,2,1,1), xpd=NA)
levs<-c(-.5,0,1,1.5)
contour(x,y,log10(z),xlim=range(x),ylim=range(y),nlevels=10, bty='l')
#marginals
plot(x,dbmarg(x),type="l", bty='n', xaxt='n', xlab="", yaxt="n") #one marginal
plot(dbmarg(x),x,type="l", bty='n', yaxt='n', ylab="", xaxt="n") #the other (which is the same in this case)
plot.new() 
#end
dev.off()


