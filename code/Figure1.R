#This script generates figure 1 from the paper
#combines ATA examples and lottery pop sims using that noise
#and plankton ATA examples

### Preamble ###

##libraries used (invoked with ::): stats, graphics, grDevices

#locations of inputs needed for this script
numRes_loc <- "../results_numeric/"
noise_loc <- paste0(numRes_loc,"noiseB.RData") #from makenoise_LB.R
popsim_loc <- paste0(numRes_loc,"betapopsim.RData") #from betapopsim.R
plankton_loc <- paste0(numRes_loc,"ceratium1x2.RData") #from plankton.R

#locations where result figs are stored
fig_loc <- "../results_figs/"
if(dir.exists(fig_loc)==FALSE){
  dir.create(fig_loc)
}
fig1 <- paste(fig_loc,"fig1_a-h.pdf",sep="")

#load the noise we'll need
load(noise_loc) #this loads an obtect called v which has beta noise in it
load(popsim_loc) #list called popsim with three kinds of popsim from beta noise
load(plankton_loc) #list called ceratium1x2 from two locations

####

#### Make the figure ###

nt <- c("l", "s", "r") #nt = noise type
blue <- grDevices::rgb(0,0,0.545,0.3) #col2rgb("darkblue")
labels <- c("(a)","(b)", "(c)")

grDevices::pdf(fig1_vertical_loc, width=2.8, height=7)

graphics::par(mar=c(0.5,0.5,0,0), oma=c(5,5,2,2))

graphics::layout(matrix(c(2,10,
                1,3 ,
                13,13,
                5,11,
                4,6,
                14,14,
                8,12,
                7,9),nrow=8,byrow=T),
       heights=c(0.25,1,0.2,0.25,1,0.2,0.25,1), widths=c(1,0.25))

for (panel in 1:3){
  
  noise <- b_tilde[[nt[panel]]]
  P <- stats::cor(noise, method='pearson')[1,2]
  
  #noise
  graphics::plot(noise[1:800,1], noise[1:800,2], col=blue, 
       ylab=NA, xlab=NA, xlim=c(-4,4), ylim=c(-4,4), pch=16, cex=1.5)
  graphics::text(x=4,y=-3.8,labels=bquote(~ rho == .(round(P,5))), adj=1)
  graphics::mtext(labels[panel], side=3, line=-1.6, at=-3.5)
  
  #marginals
  d1 <- stats::density(noise[,1])
  d2 <- stats::density(noise[,2])
  graphics::plot(d1$x, d1$y, xlab=NA, ylab=NA, sub=NA, main=NA, bty='n',type='l', xaxt='n', yaxt='n')
  graphics::plot(d2$y, d2$x, xlab=NA, ylab=NA, sub=NA, main=NA, bty='n', type='l',xaxt='n', yaxt='n')
}
graphics::mtext("variable 1", outer=TRUE, side=1, line=2.5, cex.lab=1.25, at=0.43)
graphics::mtext("variable 2", outer=TRUE, side=2, line=2.5, cex.lab=1.25, at=0.48)
grDevices::dev.off()

### Make the figure again, same but with panels arranged horizontally, as a possible alternative

grDevices::pdf(fig1_horiz_loc, width=7, height=3)

graphics::par(mar=c(0.5,0.5,0,0), oma=c(5,5,2,2))
graphics::layout(matrix(c(2,10,10,5,11,11,8,12,
                1,3 ,13,4,6 ,14,7,9 ),nrow=2,byrow=T),
       heights=c(0.25,1), widths=c(1,0.25,0.25,1,0.25,0.25,1,0.25))

for (panel in 1:3){
  
  noise <- b_tilde[[nt[panel]]]
  P <- stats::cor(noise, method='pearson')[1,2]
  
  #noise
  graphics::plot(noise[1:800,1], noise[1:800,2], col=blue, 
       ylab=NA, xlab=NA, xlim=c(-4,4), ylim=c(-4,4), pch=16, cex=1.5)
  graphics::text(x=4,y=-3.8,labels=bquote(~ rho == .(round(P,5))), adj=1)
  graphics::mtext(labels[panel], side=3, line=-1.6, at=-3.5)
  
  #marginals
  d1 <- stats::density(noise[,1])
  d2 <- stats::density(noise[,2])
  graphics::plot(d1$x, d1$y, xlab=NA, ylab=NA, sub=NA, main=NA, bty='n',type='l', xaxt='n', yaxt='n')
  graphics::plot(d2$y, d2$x, xlab=NA, ylab=NA, sub=NA, main=NA, bty='n', type='l',xaxt='n', yaxt='n')
}

graphics::title(xlab="variable 1", outer=TRUE, line=3, cex.lab=1.5)
graphics::title(ylab="variable 2", outer=TRUE, line=3, cex.lab=1.5)

grDevices::dev.off()














