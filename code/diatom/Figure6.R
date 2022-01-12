# Figure 6
source("diatom/workflow.R")
dat_loc <- "../results_numeric/fig6dat/"
if (dir.exists(dat_loc)==FALSE){
  source("diatom/makefig6dat.R")
}
fig_loc <- "../results_figs/"
if(dir.exists(fig_loc)==FALSE){
  dir.create(fig_loc)
}
fig6 <- paste0(fig_loc,"fig6.pdf")

aP <- readRDS(paste0(dat_loc,'aP.RDS')) #from makefig6dat.R
aTbar <- readRDS(paste0(dat_loc,'aTb.RDS'))
PTbar <- readRDS(paste0(dat_loc,'pTb.RDS'))

maps <- list(aP, aTbar, PTbar)

pdf(fig6, height=15, width=5)
par(mfrow=c(1,1), oma=c(3,0,1,0), mar=c(2,3,1,1), bty='n', xpd=T)
layout(matrix(c(1,2,3,4), byrow=T), heights=c(1,1,1,0.1))
for (i in 1:length(maps)){
	cmContour(maps[[i]], colkey=F, xaxt='n', yaxt='n')
	axis(side=1, mgp=c(3,0.5,0.2), col='gray50')
	axis(side=2, mgp=c(3,0.5,0.2), col='gray50')
	mtext(paste0("(", letters[i],")"), side=3, line=-1.5, adj=0.985)
}
colkey(col=cm.colors(51), clim=c(-1,1),side=1, width=10,)

title(main=expression(paste(Delta[i]^"[EC]","/IGR")), outer=T, line=-111.5, cex.main=1.5, adj=0.53)
dev.off()

