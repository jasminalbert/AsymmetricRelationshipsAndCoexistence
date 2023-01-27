#This script generates figure 1 from the paper
#combines ATA examples and lottery pop sims using that noise
#and plankton ATA examples

### Preamble ###

##libraries used (invoked with ::): stats, graphics, grDevices

#locations of inputs needed for this script
numRes_loc <- "../results_numeric/"
noise_loc <- paste0(numRes_loc,"noise.RData") #from makenoise_LB.R
popsim_loc <- paste0(numRes_loc,"betapopsim.RData") #from beta_popsims.R
plankton_loc <- paste0(numRes_loc,"ceratium1x2.RData") #from plankton.R

#locations where result figs are stored
fig_loc <- "../results_figs/"
if(dir.exists(fig_loc)==FALSE){
  dir.create(fig_loc)
}
fig1_loc <- paste(fig_loc,"fig1_a-h.pdf",sep="")

#load the noise we'll need
load(noise_loc) #this loads an obtect called b which has noise in it
load(popsim_loc) #list called popsim with three kinds of popsim from beta noise
load(plankton_loc) #list called ceratium1x2 from two locations

####

#### Make the figure ###

nt <- c("l", "s", "r") #nt = noise type
blue <- grDevices::rgb(0,0,0.545,0.3) #col2rgb("darkblue")

grDevices::pdf(fig1_loc, width=7.5, height=5)

graphics::par(mar=c(0.5,0.5,0,0), oma=c(4,3,2,1),mgp=c(2.2,0.6,0), xpd=NA, cex.lab=1.7)

graphics::layout(matrix(c(
				2 ,13,13, 4,14,15,
                1 ,3 ,13, 4,14,15,
                6 ,13,13, 8,14,15,
                5 ,7, 13, 8,14,15,
                5 ,7, 13, 8,14,13,
                5 ,7, 13, 8,14,16,
                10,13,13,12,14,16,
                9,11,13,12,14,16),nrow=8,byrow=T),
       heights=c(0.22,1,0.22,0.36,0.03,0.61,0.22,1), widths=c(1,0.25,0.4,2.5,0.45,1.85))

for (panel in 1:3){
  
  noise <- b[[nt[panel]]]
  P <- stats::cor(noise, method='pearson')[1,2]
  
  #noise
  graphics::par(mar=c(0.5,0.25,0,0))
  graphics::plot(noise[1:800,1], noise[1:800,2], col=blue, 
  ylab=NA, xlab=ifelse(panel==3,"variable 2",NA), xlim=c(-4,4), ylim=c(-4,4), pch=16, cex=1.5, col.lab="darkblue", xaxt=ifelse(panel==3,"s","n"))
  title(ylab=ifelse(panel==2,"variable 1",NA), col.lab="darkblue", line=1.8)
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
  graphics::mtext(paste0("(",letters[panel+3],")"), side=3, line=-0.8, at=7)
}

graphics::par(mar=c(0.5,0.25,0.3,0))
plot.new();plot.new()
# plankton example - ranked density
title(ylab="C. furca density", line=-1.6)
for (panel in 1:2){
	cerfus <- ceratium1x2[[panel]]$cerfus
	cerfur <- ceratium1x2[[panel]]$cerfur
	Tot <- dim(ceratium1x2[[panel]])[1]
	graphics::plot(rank(cerfus)/Tot, rank(cerfur)/Tot, bty="n",pch=20, col="darkgrey", cex=2, ylim=c(0,1),xlim=c(0,1), xlab=ifelse(panel==2,"C. fusus density",NA), ylab=NA, xaxt=ifelse(panel==2,"s","n"))
	graphics::mtext(paste0("(",letters[panel+6],")"), side=3, line=-1.5, at=0.05)
}
#mtext("C. furca density", side=4,outer=T, line=-14, cex=1)


grDevices::dev.off()




