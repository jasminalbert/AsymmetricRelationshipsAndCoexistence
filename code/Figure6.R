# Figure 6
source("diatom/mapping_fxns.R")
dat_loc <- "../results_numeric/"
fig_loc <- "../results_figs/"
if(dir.exists(fig_loc)==FALSE){
  dir.create(fig_loc)
}
fig6 <- paste(fig_loc,"fig6.pdf",sep="")

aP <- readRDS(paste(dat_loc,'.RDS',sep='')) #from .R
aTbar <- readRDS(paste(dat_loc,'.RDS',sep=''))
PTbar <- readRDS(paste(dat_loc,'.RDS',sep=''))

datlist <- list(aP, aTbar, PTbar)

#transform dataframes into matrices filled with ratio value
maps <- lapply(datlist, FUN=cmMap)


pdf(fig6, height=15, width=5)
par(mfrow=c(1,1), oma=c(3,0,3,0), mar=c(2,3,1,3), bty='n', xpd=T)
layout(matrix(c(1,2,3,4), byrow=T), heights=c(1,1,1,0.1))
for (i in 1:length(maps)){
	cmContour(maps[[i]], colkey=F, xaxt='n', yaxt='n')
	axis(side=1, mgp=c(3,0.5,0.2), col='gray50')
	axis(side=2, mgp=c(3,0.5,0.2), col='gray50')
}
colkey(col=cm.colors(51), clim=c(-1,1),side=1, width=10)

title(main=expression(paste(Delta[i]^"[EC]","/IGR")), outer=T, line=1, cex.main=1.5)
dev.off()

