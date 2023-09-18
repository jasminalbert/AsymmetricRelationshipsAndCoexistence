#This script generates figure 1 from the paper
#combines ATA examples and lottery pop sims using that noise
#and plankton ATA examples
#Also computes plankton partial correlation stats (see Ghosh2020)


##libraries used (invoked with ::): stats, graphics, grDevices

### locations of inputs and outputs needed for this script
numRes_loc <- "../results_numeric/"
noise_loc <- paste0(numRes_loc,"noise.RData") #from makenoise_LB.R
popsim_loc <- paste0(numRes_loc,"betapopsim.RData") #from beta_popsims.R
plankton_loc <- paste0(numRes_loc,"ceratium1x2.RData") #from plankton.R
plankStats_loc <- paste0(numRes_loc, "plankStats.RDS")

#locations where result figs are stored
fig_loc <- "../results_figs/"
if(dir.exists(fig_loc)==FALSE){
  dir.create(fig_loc)
}
fig1_loc <- paste0(fig_loc,"fig1_a-h.pdf")

 ##load the noise and data we'll need
load(noise_loc) #loads list called b which has noise in it
load(popsim_loc) #list, popsim, w 3 kinds of popsim from beta noise
load(plankton_loc) #list, ceratium1x2, from 2 locations
####
####
#### Plankton Stats ####
#plankton ATA stats
source("partCor.R")

cerfus23 <- ceratium1x2$loc23$cerfus
cerfur23 <- ceratium1x2$loc23$cerfur
cerfus40 <- ceratium1x2$loc40$cerfus
cerfur40 <- ceratium1x2$loc40$cerfur
#ceratium partial correlation for panel D
cPCD <- partCor(cerfus23, cerfur23, 1/2, "average")
#ceratium partial correlation for panel E
cPCE <- partCor(cerfus40, cerfur40, 1/2, "average")
plankStats <- list(d=cPCD, e=cPCE)
#plankStatsATA <- list(d=cPCD$diff, e=cPCE$diff)
saveRDS(plankStats ,plankStats_loc)
#partCor(cerfus23, cerfur23, 1/3, "average")
#partCor(cerfus40, cerfur40, 1/3, 'average')
#cor(cerfus40, cerfur40, method="spearman")


#### Make the figure ###

nt <- c("l", "s", "r") #nt = noise type
blue <- grDevices::rgb(0,0,0.545,0.3) #col2rgb("darkblue")

grDevices::pdf(fig1_loc, width=7.5, height=5)

graphics::par(mar=c(0.5,0.5,0,0), oma=c(4,3,2,1),mgp=c(2.2,0.6,0), xpd=NA, cex.lab=1.7)

#graphics::layout(matrix(c(
#				2 ,13,13, 4,14,15,
 #               1 ,3 ,13, 4,14,15,
  #              6 ,13,13, 8,14,15,
   #             5 ,7, 13, 8,14,15,
    #            5 ,7, 13, 8,14,13,
     #           5 ,7, 13, 8,14,16,
      #          10,13,13,12,14,16,
       #         9,11,13,12,14,16),nrow=8,byrow=T),
       #heights=c(0.22,1,0.22,0.36,0.03,0.61,0.22,1), widths=c(1,0.25,0.4,2.5,0.45,1.85))
graphics::layout(matrix(c(
				2 ,13,13, 15,14,4,
                1 ,3 ,13, 15,14,4,
                6 ,13,13, 15,14,8,
                5 ,7, 13, 15,14,8,
                5 ,7, 13, 13,14,8,
                5 ,7, 13, 16,14,8,
                10,13,13,16,14,12,
                9,11,13,16,14,12),nrow=8,byrow=T),
       heights=c(0.22,1,0.22,0.36,0.03,0.61,0.22,1), widths=c(1,0.25,0.4,1.85,0.45,2.5))       

for (panel in 1:3){
  
  noise <- b[[nt[panel]]]
  P <- stats::cor(noise, method='pearson')[1,2]
  
  #noise
  graphics::par(mar=c(0.5,0.25,0,0))
  graphics::plot(noise[1:800,1], noise[1:800,2], col=blue, 
  ylab=NA, xlab=ifelse(panel==3,"variable 2",NA), xlim=c(-4,4), ylim=c(-4,4), pch=16, cex=1.5, col.lab="darkblue", xaxt=ifelse(panel==3,"s","n"))
  graphics::title(ylab=ifelse(panel==2,"variable 1",NA), col.lab="darkblue", line=1.8)
  graphics::text(x=4,y=-3.8,labels=bquote(~ rho == .(round(P,5))), adj=1)
  graphics::mtext(paste0("(",letters[panel],")"), side=3, line=-1.3, at=-3)
  
  #marginals
  d1 <- stats::density(noise[,1])
  d2 <- stats::density(noise[,2])
  graphics::par(mar=c(0.25,0.25,0,0))
  graphics::plot(d1$x, d1$y, xlab=NA, ylab=NA, sub=NA, main=NA, bty='n',type='l', xaxt='n', yaxt='n')
  graphics::par(mar=c(0.25,0.25,0,0))
  graphics::plot(d2$y, d2$x, xlab=NA, ylab=NA, sub=NA, main=NA, bty='n', type='l',xaxt='n', yaxt='n')
  
  #popsims
  graphics::par(mar=c(0.5,0.25,1,0))
  graphics::plot(popsim[[panel]]$N1[1:250], type="l", ylim=c(0,50), xlab=ifelse(panel==3,"time",NA), ylab=NA, xaxt=ifelse (panel==3,"s","n"), bty="l", lwd=1.5)
  graphics::title(ylab=ifelse(panel==2,"species 1 population",NA), line=1.7)
  graphics::mtext(paste0("(",letters[panel+5],")"), side=3, line=-0.8, at=7)
}

graphics::par(mar=c(0.5,0.25,0.3,0))
graphics::plot.new();graphics::plot.new()
# plankton example - ranked density
graphics::title(ylab="C. furca density", line=17)
for (panel in 1:2){
	#cerfus <- ceratium1x2[[panel]]$cerfus
	#cerfur <- ceratium1x2[[panel]]$cerfur
	#Tot <- dim(ceratium1x2[[panel]])[1]
	pStats <- plankStats[[panel]]
	plankton <- pStats$uv
	graphics::plot(plankton, type="n", bty="n", ylim=c(0,1),xlim=c(0,1), xlab=ifelse(panel==2,"C. fusus density",NA), ylab=NA, xaxt=ifelse(panel==2,"s","n"))
	#abline(v=0.5, col="lightgrey", lwd=0.5, xpd=F)
	#abline(h=0.5, col="lightgrey", lwd=0.5, xpd=F)
	graphics::lines(0:1,0:1, lwd=0.5, col="darkgrey")
	graphics::points(plankton[pStats$bounds[,"left"],], col="darkgrey")
	graphics::points(plankton[pStats$bounds[,"right"],], col="darkgrey", pch=19)
	#graphics::plot(rank(cerfus)/Tot, rank(cerfur)/Tot, bty="n",pch=20, col="darkgrey", cex=2, ylim=c(0,1),xlim=c(0,1), xlab=ifelse(panel==2,"C. fusus density",NA), ylab=NA, xaxt=ifelse(panel==2,"s","n"))
	graphics::mtext(paste0("(",letters[panel+3],")"), side=3, line=-1.5, at=0.05)
}
#mtext("C. furca density", side=4,outer=T, line=-14, cex=1)


grDevices::dev.off()



