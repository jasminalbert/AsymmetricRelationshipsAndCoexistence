source('./diatom/fig5and6_fxns.R')
# Figure 6
parms <- c('a'=6, 'P'=60, 'Tbar'=18, 'time'=3000, 'reps'=200, 'invader'=1)
dat_loc <- "../results_numeric/fig6dat/"
if (dir.exists(dat_loc)==FALSE){
  #source("diatom/makefig6dat.R")
  a <- seq(1,6,length.out=100)
  Tbar <- seq(16,18,length.out = 100)
  P <- seq(51,199.5,length.out=100)
  dir.create(dat_loc)
  mak6(dat_loc, a_vec=a, P_vec=P, T_vec=Tbar, parms)
}
fig_loc <- "../results_figs/"
if(dir.exists(fig_loc)==FALSE){
  dir.create(fig_loc)
}
Fig6 <- paste0(fig_loc,"fig6.pdf")
fig6(Fig6, dat_loc, invader=1) 

#aP <- readRDS(paste0(dat_loc,'aP.RDS')) #from makefig6dat.R
#aTbar <- readRDS(paste0(dat_loc,'aTb.RDS'))
#PTbar <- readRDS(paste0(dat_loc,'pTb.RDS'))

#maps <- list(aP, aTbar, PTbar)

#pdf(fig6, height=15, width=5)
#par(mfrow=c(1,1), oma=c(3,0,1,0), mar=c(2,3,1,1), bty='n', xpd=T)
#layout(matrix(c(1,2,3,4), byrow=T), heights=c(1,1,1,0.1))
#for (i in 1:length(maps)){
#	cmContour(maps[[i]], colkey=F, xaxt='n', yaxt='n')
#	axis(side=1, mgp=c(3,0.5,0.2), col='gray50')
#	axis(side=2, mgp=c(3,0.5,0.2), col='gray50')
#	mtext(paste0("(", letters[i],")"), side=3, line=-1.5, adj=0.985)
#}
#colkey(col=cm.colors(51), clim=c(-1,1),side=1, width=10,)

#title(main=expression(paste(Delta[i]^"[EC]","/IGR")), outer=T, line=-111.5, cex.main=1.5, adj=0.53)
#dev.off()

