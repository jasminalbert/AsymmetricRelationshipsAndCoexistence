#qij !=1 version for beta distributed fecundities lottery model
source("./decompose_LB.R")
numRes_loc <- "../results_numeric/"
betanoise_loc <- paste0(numRes_loc, "betanoise.RDS")
B<-readRDS(betanoise_loc)


fig_loc <- "../results_figs/"
res_loc <- paste0(numRes_loc,"LBqij_figdat/")
if(dir.exists(fig_loc)==FALSE){
  dir.create(fig_loc)
}
if(dir.exists(res_loc)==FALSE){
  dir.create(res_loc)
}

fig3_LB_loc <- paste0(fig_loc,"decompFigLB_qij.pdf")
fig3_LB_dat_loc <- paste0(res_loc,"Deltas_dr")

### setting up legend ###
terms <- c(expression(Delta[i]^0), expression(Delta[i]^E),
           expression(Delta[i]^C), expression(Delta[i]^"(E#C)"),
           expression(Delta[i]^"[E||C]"), expression(Delta[i]^"[EC]")
           ,expression(GWR))
ncol <- 5
yor <- hcl.colors(ncol, palette = "YlOrRd")
cols <- c("black","black","black","black",yor[1:3])
ltys <- c(1,2,3,4,2,1,1)
lwds <- c(rep(1.5,4),2,2,3)

#params
etaj <- seq(1,5,0.1)
etai<-1; 
delta<-c(0.2,0.4,0.6)

#start fig
pdf(fig3_LB_loc, height=8, width=10)
par(mfrow=c(2,3), oma=c(3,4,4,7), mar=c(1,1,1,1), mgp=c(3,1,0), xpd=NA)

#panel a-c (left tail)
maxseLT <- {}
for (d in seq_along(delta)){
	
	file <- paste0(fig3_LB_dat_loc,delta[d],"_LT.RDS")
	sefile <- paste0(fig3_LB_dat_loc,delta[d],"SE_LT.RDS")
	DeltasLT <- plotco(etai,etaj,delta[d], B,Deltas_loc=file, qij=T)
	maxseLT[d] <- max(DeltasLT$SE,na.rm=T)
	mtext(paste(delta[d]), side=3, line=0.2, font=2, cex=1.5, col="gray30")
	#axis(1, cex.axis=2,tck=-0.028, lwd.ticks=2)
	if (d==1){
		axis(2, cex.axis=1.8, tck=-0.035, lwd.ticks=2)
	}
	if (d==2){
		graphics::mtext("adult death rate, ", line=2.5, font=2, cex=1.5, col="gray40", at=1.25, adj=0)
    	graphics::mtext(expression(delta), line=2.5, font=2, cex=2.3, col="gray40", side=3, at=5, adj=1)
	} 
	if (d==3){
		text("LEFT-TAILED", x=5.4, y=0.23, font=2, srt=-90, cex=1.8,adj=0, col="grey")
	}
	
	#save
	if (file.exists(file)==FALSE){saveRDS(DeltasLT$D,file=file)}
	if (file.exists(sefile)==FALSE){saveRDS(DeltasLT$SE,file=sefile)}
}


#d-f (right-tail)
maxseRT <- {}
for (d in 1:length(delta)){
	file <- paste0(fig3_LB_dat_loc,delta[d],"_RT.RDS")
	sefile <- paste0(fig3_LB_dat_loc,delta[d],"SE_RT.RDS")	
	DeltasRT <- plotco(etai,etaj,delta[d],B, Deltas_loc=file, dir="RIGHT", qij=T)
	maxseRT[d] <- max(DeltasRT$SE,na.rm=T)
	axis(1, cex.axis=1.8,tck=-0.028, lwd.ticks=2)
	if (d==1){
		axis(2, cex.axis=1.8, tck=-0.035, lwd.ticks=2)
	}
	if (d==2){
		title(xlab="upper bound ratio", line=2.75, font.lab=2, cex.lab=2, col.lab="gray40")
	}	
	if (d==3){
		text("RIGHT-TAILED", x=5.4, y=-0.62, font=2, srt=-90, cex=1.8,adj=1, col="grey")
	}
	#save
	if (file.exists(file)==FALSE){saveRDS(DeltasRT$D,file=file)}
	if (file.exists(sefile)==FALSE){saveRDS(DeltasRT$SE,file=sefile)}
}


#7
rescol <- rgb(227/255, 211/255, 148/255,.5)
excol <- rgb(38/255, 38/255, 38/255,.5)
#graphics::plot.new() #1: legend panel
graphics::legend("topright", legend=terms, col = cols, lty = ltys, bty="n", cex=2, inset=c(-.36,-.72), y.intersp = 1.35, x.intersp = 0.1, seg.len=0.8, lwd=lwds)
graphics::legend("topright",legend=c("ATA \nexlusion", "ATA \nrescue"), fill=c(excol,rescol), cex=1.8, bty="n",border=NA, inset=c(-.427,0.15), x.intersp = 0.1,y.intersp = 1.7)

mtext("contribution to coexistence", side=2, outer=TRUE, line=2, font=2, cex=1.7, col="gray40")

dev.off()

maxSE_LBdecomp <- max(c(maxseLT, maxseRT))
saveRDS(maxSE_LBdecomp,file=paste0(numRes_loc,"fig4maxse.RDS"))






