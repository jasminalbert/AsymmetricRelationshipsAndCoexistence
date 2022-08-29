# pedagogical figure
#numeric_results_loc <- "../results_numeric"
#betanoise_loc <- paste0(numeric_results_loc, "/betanoise.RDS")
#b<-readRDS(betanoise_loc)
#removing asymmetries
#start with ATA noise
#blue <- grDevices::rgb(0,0,0.545,0.3) #col2rgb("darkblue")
purp <- grDevices::rgb(.39,.09,.61,0.3)
pdf("pedagogy.pdf",width=7, height=4)
layout(matrix(c(7, 5, 8, 4, 12,12,16,15,15,15,19,
				11,11,11,11,14,14,16,15,15,15,19,
				2, 2, 2, 4, 13, 4,16,18,18,18,19,
				1, 1, 1, 3, 13, 9,16,17,17,17,20,
				1, 1, 1, 3, 13, 6,16,17,17,17,20,
				1, 1, 1, 3, 13,10,16,17,17,17,20), nrow=6, byrow=T), 
				heights=c(0.33,0.33,0.25,0.33, 0.33, 0.33), widths=c(0.33,0.33,0.33,0.25,0.33,0.33,.3,.33,.33,.33,.25))
				
graphics::par(mar=c(0.5,0.5,0,0), oma=c(5,5,2,2), xpd=NA)

plot(b$l[1:800,1], b$l[1:800,2], col=blue,ylab="E", xlab="C", xlim=c(-0.5,0.5), ylim=c(-0.5,0.5), pch=16, cex=1.5, font.lab=2, xaxt="n", yaxt="n")
axis(1, at=c(-0.5,0,0.5))
axis(2, at=c(-0.5,0,0.5))

  d1 <- stats::density(b$l[,1])
  d2 <- stats::density(b$l[,2])
  graphics::plot(d1$x, d1$y, xlab=NA, ylab=NA, sub=NA, bty='n',type='l', xaxt='n', yaxt='n', col="darkblue")
  graphics::plot(d2$y, d2$x, xlab=NA, ylab=NA, sub=NA, bty='n', type='l',xaxt='n', yaxt='n', col="darkblue")
  
plot.new() #4
q <- seq(-0.5,0.5,0.01)
cdf<-pbeta(p, 0.5, 0.5)
plot(q, cdf, type="l", col="darkblue", bty="l", xaxt='n', yaxt='n',xlab='',ylab='',main=expression('F'['E'['i']]) )#5
axis(1, at=c(-0.5,0,0.5),mgp=c(2.5,0.55,0),tck=-0.1)
axis(2, at=c(0,1),mgp=c(2.5,0.55,0),tck=-0.1)
arrows(0.45,0.5,2.6,0.5, col="lightgrey", length=.1)
text(x=1.5, y=0.7, labels=c(expression('F'['E'['i']]("E"["i"]))))

plot(q, cdf, type="l", col="darkblue", bty="l", xlab='',ylab='',xaxt='n',yaxt='n',main=expression('F'['C'['i/i']])) #6
axis(2, at=c(0,1),mgp=c(2.5,0.55,0),tck=-0.1)
axis(1, at=c(-0.5,0,0.5),mgp=c(2.5,0.55,0),tck=-0.1)
arrows(0.0,1.2,0,3, col="lightgrey", length=.1)
#arrows(0.9,0,2.63,0, col="lightgrey", length=.1)
text(x=0.2, y=2, labels=c(expression('F'['C'['i/i']]("C"["i/i"]))), srt=-90)


#7-10
plot.new();plot.new();plot.new();plot.new()
#11
plot.new()
arrows(0.375,-.5,0.375,0.6,length=.15)
arrows(0.375,-.49,0.375,-.5, length=.15)

#12 - cdf of E of E
cdfE <- pbeta(b$l[,1]+0.5,0.5,0.5)
d_cdfE <- stats::density(cdfE)
plot(d_cdfE$x-0.5, d_cdfE$y, xlab=NA, ylab=NA, sub=NA, bty='n',type='l',  col="#27dbbd", main='',xaxt="n",yaxt="n")
axis(2, at=c(0,1),mgp=c(2.5,0.55,0),tck=-0.1)

#13 try to make arrows?
plot.new()
arrows(-.5,.375,0.6,.375, length=.15)
arrows(-.49,.375,-.5,.375, length=.15)
#14
cdfC <- pbeta(b$l[,2]+0.5,0.5,0.5)
d_cdfC <- stats::density(cdfC)
plot(d_cdfC$x-0.5, d_cdfC$y, xlab=NA, ylab=NA, sub=NA, bty='n',type='l',  col="#27dbbd", main='',xaxt="n",yaxt="n")
axis(2, at=c(0,1),mgp=c(2.5,0.55,0),tck=-0.1)
axis(1, at=c(-0.5,0,0.5),mgp=c(2.5,0.55,0),tck=-0.1);


#15 phi: cdf of standard normal 
par(mar=c(2.5,1,0,0))
q <- seq(-3,3,0.05)
ncdf <- pnorm(q)
plot(q, ncdf, type="l", col="darkred", bty="l", xlab='',ylab='', main=expression(phi), xaxt='n', yaxt='n',cex.main=1.5)#5
axis(1, at=c(-3,0,3),mgp=c(2.5,0.55,0),tck=-0.05)
axis(2, at=c(0,1),mgp=c(2.5,0.55,0),tck=-0.05)

plot.new() #16
par(mar=c(0.5,1,0,0))
#17
nE <- qnorm(cdfE[1:800]); nC <- qnorm(cdfC[1:800]) 
plot(nE,nC, xlab='', ylab='', col=purp,pch=16, cex=1.5)

#18
d1 <- stats::density(nE)
graphics::plot(d1$x, d1$y, xlab=NA, ylab=NA, sub=NA, bty='n',type='l', xaxt='n', yaxt='n', col="#59068a")


#19
plot.new()

#20
d2 <- stats::density(nC)
graphics::plot(d2$y, d2$x, xlab=NA, ylab=NA, sub=NA, bty='n', type='l',xaxt='n', yaxt='n', col="#59068a")
dev.off()
