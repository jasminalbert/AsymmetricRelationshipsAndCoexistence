#This script contains functions for making data for figure 5 and plotting

##libraries used (invoked with ::): graphics, grDevices

### dat5 ###
# computes Deltas for plotting 
#ARGS:
  #dat_loc  #folder to put data
  #a_vec    vector containg set of amplitudes
  #P_vec    vector containing set of Periods
  #T_vec    vector containg set of mean temperatures
  #parms    orignal parameters of diatom coexistence to stay constant
            #parms = c('a'=6, 'P'=60, 'Tbar'=18, 'time'=3000, 'reps'=200, 'invader'=1)
#OUT:
  #folder containing three .RDS files for each panel of figure 5
dat5 <- function(dat_loc, a_vec, P_vec, T_vec, parms){
  
  #a,P, and T vectors should be same length
  if((length(a_vec)==length(P_vec) & length(P_vec)==length(T_vec)) == FALSE){
    paste("vectors are not equal length")
  }
  
  len <- length(a_vec)
  names <- c("amplitude","Period","Tbar")
  plots <- 1:3
  aPT <- cbind(a_vec, P_vec, T_vec)
  
  for (p in plots){ #for each panel of figure
    
    #fill a matrix with original paramters
    Mat <- matrix(c(parms, FALSE), nrow=len, ncol=length(parms)+1, byrow=T)
    
    #change appropriate column to one of the vectors (allowing one parm to vary)
    Mat[,p] <- aPT[,p]
    
    #function from diatomDecomp_fxns.R; computes Deltas
    cat("making data for panel",p)
    Delta <- wrapDelt(args=Mat[1,])
    
    #iterate over rows of Mat to compute Deltas
    for (i in 2:len){
      Delta[i,] <- wrapDelt(Mat[i,])
    }
    Delta$var <- aPT[,p]
    
    #save
    saveRDS(Delta, paste0(dat_loc, "vary", names[p],parms["invader"],".RDS"))
  }
} 

#### fig5 ###
# function that plots fig 5 data and produces figure as PDF
#ARGS:
  #filename   name of PDF
  #dat_loc    name of folder containing the .RDS data files
  #invader    1 or 2, which species is the rare species?
#OUT:
  #PDF file of figure 5
fig5 <- function(dat_loc, invader){
  names <- c("amplitude","Period","Tbar")
  plots <- 1:3
  xat <- list(a=c(1,2,3,4,5,6), P=c(50, 80, 110, 140, 170, 200), T=c(16, 16.5, 17, 17.5, 18))
  xlabs <- list(a=c(1,2,NA,NA,5,6), P=c(50, 80, NA, NA, 170, 200), T=c(16, NA, NA, NA, 18))
  #empty list
  varylist <- list()
  
  for (p in plots){ #for each panel of figure
    
    #load data from folder in directory
    varylist[[p]] <- readRDS(paste0(dat_loc, "vary", names[p],invader,".RDS"))
  }
  
  #colors and line styles
  col <- c(rep("black",4), "red", "blue", "orange")
  line <- c(1,2,4,3,1,1,1); lwd <- c(rep(2,6),3)
  #plotting set up
	rescol <- grDevices::rgb(227/255, 211/255, 148/255,.5)
	excol <- grDevices::rgb(38/255, 38/255, 38/255,.5)
	lwd <- c(2,2,2,2,4,4,5)
	line <- c(1,2,3,4,1,2,1)
	ncol <- 5
	yor <- grDevices::hcl.colors(ncol, palette = "YlOrRd")
	col <- c("black","black","black","black", yor[2], yor[1], yor[3])
  
  #get range for y axis
  range <- sapply(varylist, function(X){range(X[,1:7])})
  ymin <- min(range); ymax <- max(range) + (max(range))*0
  
  #set original diatom coexistence paramters (Gonzalez2005) to denote in panels
  orig <- c(a=6, P=60, Tbar=18)
  xlab <- c("amplitude (a)", "Period (P)", expression(paste("mean temperature (",theta[0],")" ) ))
  
  #start fig
  #grDevices::pdf(filename, height=5.5, width=14)
  #graphics::par(mfrow=c(1,3), oma=c(1,6,1,0), mar=c(6,1,2,1))
  
  for (i in 1:3){
    vdat <- varylist[[i]]
    #empty box
    graphics::plot(0, yaxt='n', xaxt='n', xlim=range(vdat$var), ylim=c(ymin, ymax), col='white', ylab='', xlab='')
    graphics::title(xlab=xlab[i], cex.lab=1.3, line=1.1, xpd=NA)
    #axis
    graphics::axis(1, lwd.ticks=2, at=xat[[i]], labels=xlabs[[i]])
    #if (i==1){
      graphics::axis(2, lwd.ticks=2,mgp=c(2,.6,0))
    #}
    #loop for lines
    for (j in 1:7){
      graphics::lines(vdat$var, vdat[,j], col=col[j], lty=line[j], lwd=lwd[j])	
    }
    graphics::abline(h=0, col="lightgrey") #zero
    graphics::abline(v=orig[i], col='darkgrey', lty=3, lwd=2)
    
    #ATA effect blocking
    resc <- vdat[vdat$ECbrk>vdat$GWR & vdat$GWR>0,]
    graphics::rect(resc[1,10], ymin*1.2, tail(resc,1)[1,10], ymax*1.2, col=rescol, border=NA)
    
    excl <- vdat[vdat$ECbrk<vdat$GWR & vdat$GWR<0,]
    graphics::rect(excl[1,10], ymin, tail(excl,1)[1,10], ymax, col=excol, border=NA)
    
    #label
    graphics::mtext(paste0("(", letters[i],")"), 3, -1.5, adj=0.985, cex=1.1)
    #legend
    if (i==1){
      graphics::legend("topleft", legend=c(expression(Delta[i]^0), expression(Delta[i]^E), expression(Delta[i]^C), expression(Delta[i]^"(E#C)")), 
             col ="black",lty = line[1:4], bty="n", cex=1.5, inset=c(-0.01,-0.055), 
             y.intersp = 1.1, x.intersp = 0.1, seg.len=0.8, lwd=1.5)
      graphics::legend("topleft", legend=c(expression(Delta[i]^"[E||C]"), expression(Delta[i]^"[EC]") ,expression(GWR)), col = yor[1:3], lty=c(2,1,1), bty='n', cex=1.5, inset=c(0.2, -0.055), x.intersp=0.2, seg.len=0.8, lwd=c(2, 2, 3), y.intersp=1.1)
      graphics::legend("topleft", legend=paste("ATA\nrescue"), fill=rescol, bty='n', border=NA,cex=1.4, inset=c(0.47, -0.055), x.intersp=0.2, y.intersp=1.1)	
    }
  }
  #grDevices::dev.off()
}

