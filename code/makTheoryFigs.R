#this figure will:
	# have four panels A-D
	#A) E and C with ATA
	#B) E and C transformed to standard normal but still has ATA
	#C) remove ATA by replacing [E,C] with biv normal with same parms
	#D) invert marginals to [E||,C||]
	 

source("./CopulaToolForPedagogFig.R")
source("./theory_panel_fxn.R")

#figure file location and names
fig_loc <- "../results_figs"
if (dir.exists(fig_loc)==F){
  dir.create(fig_loc)
}
theoryfig_col <- paste0(fig_loc,"/theory_fig_col.pdf")
theoryfig_col2 <- paste0(fig_loc,"/theory_fig_col2.pdf")


#copulas
ncop<-normalCopula(.7)
shcop <- normalCopula(0.001)
ccop0<-claytonCopula(2)
ccop<-claytonCopula(iRho(ccop0,rho(ncop)))
shcop <-claytonCopula(iRho(ccop0,rho(shcop)))

#distribution functions
##beta
dbmarg<-function(x){return(dbeta(x,shape1=shape1,shape2=shape2))}
pbmarg<-function(x){return(pbeta(x,shape1=shape1,shape2=shape2))}
qbmarg<-function(x){return(qbeta(x,shape1=shape1,shape2=shape2))}
##normal
dnmarg<-function(x){return(dnorm(x,mean=0,sd=1))}
pnmarg<-function(x){return(pnorm(x,mean=0,sd=1))}
qnmarg<-function(x){return(qnorm(x,mean=0,sd=1))}


shape1<-0.55;shape2<-0.55

### set up ###
hts <- c(.3,1)
wds <- c(rep(c(1,.3,.3),3), c(1,.3))
lblc <- matrix(c(-2.1,-4.8,-4.3,-4.4,-1.7,0.2), nrow=3,dimnames=list(c("E","C","lab"),c("x","y")))

#single row layout
ly_mat<- matrix(c (	2,4,17,6,8,17,10,12,17,14,16,
				1,3,17,5,7,17,9,11,17,13,15), byrow=T, nrow=2)

#colored points
#### start pdf ####
#pdf(theoryfig_col, height=4, width=15)
#layout(ly_mat, heights=hts, widths=wds)#;layout.show(n=16)
#par(mar=c(0,0.25,0,0.25), oma=c(6,5,2,3), xpd=NA)
#source("./theoryFig_script.R")
#dev.off() 
########

#new fig for reviewers
#### start pdf ####
new <- TRUE
ly_mat <- cbind(ly_mat, c(17,17), c(19,18), c(21,20))
lblc[2:3,"x"] <- lblc[2:3,"x"]+.17
lblc["E","y"] <- lblc["E","y"]+.1  
lblc["C","y"] <- lblc["C","y"]-.18  
wds[(1:3)*3] <- 0.2
pdf(theoryfig_col2, height=4, width=19)
layout(ly_mat, heights=hts, widths=c(wds,wds[3:5]))#;layout.show(n=16)
par(mar=c(0,0.25,0,0.25), oma=c(6,5,2,3), xpd=NA)
source("./theoryFig_script.R")
dev.off() 
########

#quad layout
#theoryfig_qd <- paste0(fig_loc,"/theory_fig_qd.pdf")
#lblc[2,1] <- -4.5; lblc[1,2] <- -4.15; lblc[3,] <- c(-4.1,0.15)
#ly_mat<- matrix(c(2,4,17,6,8,
#				1,3,17,5,7,
#				17,17,17,17,17,
#				14,16,17,10,12,
#				13,15,17,9,11), byrow=T, nrow=5)
#hts <- c(.3,1,.2,.3,1); wds <- c(1,.3,.2,1,.3)
#### start pdf ####
#pdf(theoryfig_qd, height=10, width=10)
#layout(ly_mat, heights=hts, widths=wds);#layout.show(n=17)
#par(mar=c(0,0.25,0,0.25), oma=c(6,5,2,3), xpd=NA)
#source("./theoryfig_meat.R")
#dev.off() 
########

#no color
#theoryfig <- paste0(fig_loc,"/theory_fig_rw.pdf")
#### start pdf ####
#pdf(theoryfig, height=4, width=15)
#layout(ly_mat, heights=hts, widths=wds)#;layout.show(n=16)
#par(mar=c(0,0.25,0,0.25), oma=c(6,5,2,3), xpd=NA)
#source("./theoryfig_meat.R")
#dev.off() 
########












#pdf("panelA_0.5_filled.pdf") #issue - filled.con not going into layout
#layout(matrix(c(2,4,1,3), byrow=T, nrow=2), heights=c(.3,1), widths=c(1,.3))
#par(mar=c(0.5,0.5,0.5,0.5), oma=c(2,2,1,1), xpd=NA)
#levs<-c(-.5,0,1,1.5)
#filled.contour(x,y,log10(z),xlim=range(x),ylim=range(y),nlevels=40, bty='l')
#marginals
#plot(x,dbmarg(x),type="l", bty='n', xaxt='n', xlab="", yaxt="n") #one marginal
#plot(dbmarg(x),x,type="l", bty='n', yaxt='n', ylab="", xaxt="n") #the other (which is the same in this case)
#plot.new() 
#dev.off()






