#source("./Figure5.R")

#This script makes figure 5 in main text:
#

##libraries used (invoked with ::): plot3D

### source diatom competition model coexistence decomposition functions ###
source('./diatom/diatomDecomp_fxns.R')
### source functions for computing data for and plotting panels a-c
source("./diatom/plotdiatomdecomp.R")
### source functions for computing data for and plotting panels d-f
source("./mapping_fig6.R")

### location of results ###
numRes_loc <- "../results_numeric/"
dat_loc <- paste0(numRes_loc, "fig5dat/")

# if results are missing, make them
if (dir.exists(dat_loc)==FALSE){
  
  #make folder
  dir.create(dat_loc)
  
  #define original parameters
  parms <- c('a'=6, 'P'=60, 'Tbar'=18, 'time'=3000, 'reps'=200, 'invader'=1)
  
  #define sets of variables
  a <- seq(1,6,length.out=100)
  Tbar <- seq(16,18,length.out = 100)
  P <- seq(51,199.5,length.out=100)
  
  #function that produces results and saves .RDS's to folder
  dat5(dat_loc, a_vec=a, P_vec=P, T_vec=Tbar, parms=parms)
}

### location to save figure ###
fig_loc <- "../results_figs/"
if(dir.exists(fig_loc)==FALSE){
  dir.create(fig_loc)
}
Fig5 <- paste0(fig_loc,"fig5_new.pdf")
xat <- list(a=c(4,5,6), P=seq(50, 120, len=5), T=c(16, 16.5, 17, 17.5, 18))
xlabs <- list(a=c(4,NA,6), P=c(50, 67.5, NA, 102.5, 120), T=c(16, NA, NA, NA, 18))
labs <- c("amplitude (a)", "Period (P)", expression(paste("mean temperature (",theta[0],")" ) ))

pdf(Fig5)
par(mfcol=c(1,1), mar=c(1,3,1.5,0), oma=c(1,1,0,1.5),mgp=c(2,.7,0), cex.axis=1.3, tck=-0.03, lheight=.7)
layout(matrix(c(1,4,7,
				2,5,7,
				3,6,7),ncol=3,byrow=T),widths=c(.7,1,.2))

fig5(dat_loc, invader=1) 
graphics::title(ylab="contribution to coexistence", outer=T, line=-0.6, cex.lab=2, font.lab=2, col.lab="gray40")

ncol <- 100
#col <- hcl.colors(n,"Geyser");col1 <- col[(n/2 + 1):n]
col1 <- hcl.colors(ncol,"YlOrRd", rev=T, alpha=.8)
GWRmat1 <- MAT(aTb, aTb$GWR)#[60:100,]
ATAmat1 <- MAT(aTb, aTb$ATA)#[60:100,]
mat1 <- ATAmat1/abs(GWRmat1)
GWRmat2 <- t(MAT(aP, aP$GWR))#[60:100,1:50] 
ATAmat2 <- t(MAT(aP, aP$ATA))#[60:100,1:50]
mat2 <- ATAmat2/abs(GWRmat2)
GWRmat3 <- (MAT(PTb, PTb$GWR))#[,1:50] 
ATAmat3 <- (MAT(PTb, PTb$ATA))#[,1:50]
mat3 <- ATAmat3/abs(GWRmat3)
matlist <- list(mat1, mat2, mat3)
ranges <- log(rbind(sapply(matlist, min), rep(10,3)))
diffs <- apply(ranges, 2, diff)
standiff <- round(diffs/max(diffs),2)
stanRan <- rbind(ranges[1,]/min(ranges[1,]), ranges[2,]/max(ranges[2,])); stanRan <- round(stanRan,2)

col1 <- col1; col2 <- col1[(ncol*(1-standiff[2])):ncol]; col3 <- col1[(ncol*(1-standiff[3])):ncol]

atmap <- map1(ATAmat1, GWRmat1, col=col1)
graphics::title(ylab=labs[1], cex.lab=1.3, line=.8, xpd=NA)
graphics::title(xlab=labs[3], cex.lab=1.3, line=1.1, xpd=NA)
axis(1, lwd.ticks=2, at=xat$T, labels=xlabs$T); axis(2, lwd.ticks=2,mgp=c(2,.5,0), at=xat$a, labels=xlabs$a)
graphics::mtext(paste0("(", letters[4],")"), 3, -1.5, adj=0.985, cex=1.1)

map1(ATAmat2, GWRmat2, col=col2) #all postive 
graphics::title(ylab=labs[1], cex.lab=1.3, line=.8, xpd=NA)
graphics::title(xlab=labs[2], cex.lab=1.3, line=1.1, xpd=NA)
axis(1, lwd.ticks=2, at=seq(4,6,len=5),labels=xlabs$P); axis(2, lwd.ticks=2,mgp=c(2,.5,0), at=seq(50,120,len=3), labels=c(4,NA,6))
graphics::mtext(paste0("(", letters[5],")"), 3, -1.5, adj=0.985, cex=1.1)

map1(ATAmat3, GWRmat3, col=col3) #all positive
graphics::title(ylab=labs[3], cex.lab=1.3, line=.8, xpd=NA)
graphics::title(xlab=labs[2], cex.lab=1.3, line=1.1, xpd=NA)
axis(1, lwd.ticks=2, at=xat$T, labels=xlabs$P); axis(2, lwd.ticks=2,mgp=c(2,.5,0), at=xat$P, labels=xlabs$T)
graphics::mtext(paste0("(", letters[6],")"), 3, -1.5, adj=0.985, cex=1.1)

par(mar=c(0,0,0,0))
plot.new()
at <- rep(c(5,1),3)*10^(c(-2,-1,-1,0,0,1))
colkey(col1, clog=T, add=T, clim=atmap$lim1, length=0.7, dist=-0.85, shift=0, at=at, labels=c(as.character(signif(at[1:5],1)),NA),cex.axis=2, mgp=c(3,.7,0),, tck=-.3, width=8, side=4)

text(x=par("usr")[1]+c(.65, .9), y=par("usr")[4]-.161, xpd=NA, labels=c(expression("">=""),10), cex=2)
#
#colkey(atmap$col2, clog=T, add=T, clim=atmap$lim2, length=0.065, dist=-0.85, shift=-0.3445, at=c(r2[1],10^(seq(-6,-4,1))), labels=F,cex.axis=1.2, mgp=c(3,.5,0), tck=-.3, width=8, side=4)

#text(x=par("usr")[2]-.6, y=par("usr")[3]+.15, xpd=NA, adj=0, labels=as.character(signif(10^-5,1)*-1), srt=35, cex=1.6)
text(x=par("usr")[1]+.6, y=par("usr")[4]-.07,labels=expression(paste(over(Delta[i]^"[EC]",paste("|",GWR,"|")))), cex=1.8 )

rect(par("usr")[2]-.85,par("usr")[3]+.08,par("usr")[2]-.65,par("usr")[3]+.1, border=rgb(0,0,1,0.7), lwd=2)
text(x=par("usr")[2]-.6, y=par("usr")[3]+.09, xpd=NA, adj=0, labels="coexistence", srt=35, cex=1)
rect(par("usr")[2]-.85,par("usr")[3]+.04,par("usr")[2]-.65,par("usr")[3]+.06, density=18, col=rgb(0,0,0,0.6),border=NA)
text(x=par("usr")[2]-.6, y=par("usr")[3]+.05, xpd=NA, adj=0, labels="ATA rescue", srt=35, cex=1)





dev.off()




### get and save standard error ###
fig5se <- getSE(dat_loc)
fig5maxse <- max(fig5se)
cat("\nmaximum standard error in figure five is", fig5maxse)
fig5maxse_loc <- paste0(numRes_loc, "fig5maxse.RDS")
saveRDS(fig5maxse, file=fig5maxse_loc)





