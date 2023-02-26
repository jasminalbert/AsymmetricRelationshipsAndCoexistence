source("./decompose_LB.R")
numRes_loc <- "../results_numeric/"
betanoise_loc <- paste0(numRes_loc, "betanoise.RDS")
B<-readRDS(betanoise_loc)

plotco <- function(etai,etaj,delta, noise,Deltas_loc=0,...){
	
	if (Deltas_loc!=0){ #check that file exists
		
		if (file.exists(Deltas_loc)==FALSE){ 
			Deltas_loc<-0
			} else {Deltas <- readRDS(Deltas_loc)}
		}
		
	if (Deltas_loc==0){
		message <- "computing coexistence mechanisms..."
		cat(message)
		store <- vector(mode='list')

		for (j in 1:length(etaj)){
    		store[[j]] <- decomp(etai,etaj[j], delta,B,...)
    		cat(".")
    	}
    	cat("done")
    	Deltas <- data.frame(t(sapply(store, function(X){X$Delta_i})))
		colnames(Deltas) <- rownames(store[[1]])
		SE <- data.frame(t(sapply(store, function(X){X$stanErr_D_i})))
		colnames(SE) <- colnames(Deltas)

	}
	
    range <- range(Deltas)
    #ylim <- range*1.1
    ylim <- c(-0.6,0.2)
    
    #ATA effect
	res_vec <- etaj[Deltas$r>0 & Deltas$ATA > Deltas$r]
	ex_vec <- etaj[Deltas$r<0 & Deltas$ATA < Deltas$r]
	
	#plotting set up
	rescol <- rgb(227/255, 211/255, 148/255,.5)
	excol <- rgb(38/255, 38/255, 38/255,.5)
	wid <- c(2,2,2,2,4,4,5)
	typ <- c(1,2,3,4,2,1,1)
	ncol <- 5
	yor <- hcl.colors(ncol, palette = "YlOrRd")
	col <- c("black","black","black","black", yor[1], yor[2], yor[3])
	
	
	##plot##
	#box	
	plot(etaj, xlab="", ylab="", ylim=ylim, 			xlim=range(etaj), type="n",xaxt="n", yaxt="n")
	#zero line
	lines(etaj, rep(0, length(etaj)), col="grey", lwd=0.5)
	#ATA effects
	if (length(res_vec)>0){ #rescue
		res_points <- range(res_vec)
		polygon(x=c(rep(res_points[1],2), rep(res_points[2],2)), y=c(ylim, rev(ylim))*2, xpd=FALSE, col=rescol, border=NA)
	}
	if (length(ex_vec)>0){ #exclusion
		ex_points <- range(ex_vec)
		polygon(x=c(rep(ex_points[1],2), rep(ex_points[2],2)), y=c(ylim, rev(ylim))*2, xpd=FALSE, col=excol, border=NA)
	}
	#lines
	for (m in 1:length(Deltas)){
		lines(etaj, Deltas[,m], lwd=wid[m], lty=typ[m], col=col[m], xpd=F)
	}
	
	return(list(D=Deltas, SE=SE))
	#title(main=paste("mu_j=",params["mu_j"], "sigma_i=",params["sigma_i"], "sigma_j=",params["sigma_j"], "delta=",params["delta"]))
}
fig_loc <- "../results_figs/"
res_loc <- paste0(numRes_loc,"fig3LB/")
if(dir.exists(fig_loc)==FALSE){
  dir.create(fig_loc)
}
if(dir.exists(res_loc)==FALSE){
  dir.create(res_loc)
}

fig3_LB_loc <- paste0(fig_loc,"decompFigLB.pdf")
fig3_LB_dat_loc <- paste0(res_loc,"fig3LB_Deltas_dr")

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
	DeltasLT <- plotco(etai,etaj,delta[d], B,Deltas_loc=file)
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
#text("LEFT-TAILED ASYMMETRY", x=5.4, y=-0.2, font=2, srt=-90, cex=1.5)

#d-f (right-tail)
maxseRT <- {}
for (d in 1:length(delta)){
	file <- paste0(fig3_LB_dat_loc,delta[d],"_RT.RDS")
	sefile <- paste0(fig3_LB_dat_loc,delta[d],"SE_RT.RDS")	
	DeltasRT <- plotco(etai,etaj,delta[d],B, Deltas_loc=file, dir="RIGHT")
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
#text("RIGHT-TAILED ASYMMETRY", x=5.4, y=-0.2, font=2, srt=-90, cex=1.5)


#7
rescol <- rgb(227/255, 211/255, 148/255,.5)
excol <- rgb(38/255, 38/255, 38/255,.5)
#graphics::plot.new() #1: legend panel
graphics::legend("topright", legend=terms, col = cols, lty = ltys, bty="n", cex=2, inset=c(-.36,-.72), y.intersp = 1.35, x.intersp = 0.1, seg.len=0.8, lwd=lwds)
graphics::legend("topright",legend=c("ATA \nexlusion", "ATA \nrescue"), fill=c(excol,rescol), cex=1.8, bty="n",border=NA, inset=c(-.427,0.15), x.intersp = 0.1,y.intersp = 1.7)

#mtext("adult death rate", side=3, outer=T, line=0.2, font=2, cex=1.5)
mtext("contribution to coexistence", side=2, outer=TRUE, line=2, font=2, cex=1.7, col="gray40")
#mtext("variance metric", side=1, outer=T, line=1, font=2, cex=1.5)
dev.off()

maxSE_LBdecomp <- max(c(maxseLT, maxseRT))
saveRDS(maxSE_LBdecomp,file=paste0(numRes_loc,"fig4maxse.RDS"))






