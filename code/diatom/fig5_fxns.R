#This script contains functions for making data for figure 5 and plotting

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
    Mat <- matrix(parms, nrow=len, ncol=length(parms), byrow=T)
    
    #change appropriate column to one of the vectors (allowing one parm to vary)
    Mat[,p] <- aPT[,p]
    
    #function from diatomDecomp_fxns.R; computes Deltas
    Delta <- wrapDelt(args=Mat[1,])
    
    #iterate over rows of Mat to compute Deltas
    for (i in 2:len){
      Delta[i,] <- wrapDelt(Mat[i,])
      print(i)
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
fig5 <- function(filename, dat_loc, invader){
  names <- c("amplitude","Period","Tbar")
  plots <- 1:3
  
  #empty list
  varylist <- list()
  
  for (p in plots){ #for each panel of figure
    
    #load data from folder in directory
    varylist[[p]] <- readRDS(paste0(dat_loc, "vary", names[p],invader,".RDS"))
  }
  
  #colors and line styles
  col <- c(rep("black",4), "red", "blue", "orange")
  line <- c(1:4, 1,1,1); lwd <- c(rep(2,6),3)
  
  #get range for y axis
  range <- sapply(varylist, function(X){range(X[,1:7])})
  ymin <- min(range); ymax <- max(range)
  
  #set original diatom coexistence paramters (Gonzalez2005) to denote in panels
  orig <- c(a=6, P=60, Tbar=18)
  xlab <- c("a", "P", expression(theta[0]))
  
  #start fig
  pdf(filename, height=5, width=15)
  par(mfrow=c(1,3), oma=c(0,4,0,0), mar=c(5,1,2,1) )
  
  for (i in 1:3){
    vdat <- varylist[[i]]
    #empty box
    plot(0, yaxt='n', xlim=range(vdat[,10]), ylim=c(ymin, ymax), col='white', xlab="", ylab='', cex.lab=1.8, cex.axis=1.8)
    title(xlab=xlab[i], cex.lab=2.3, line=3.5)
    #axis
    if (i==1){
      axis(2, cex.axis=1.8)
    }
    #loop for lines
    for (j in 1:7){
      lines(vdat[,10], vdat[,j], col=col[j], lty=line[j], lwd=lwd[j])	
    }
    abline(h=0, col="lightgrey") #zero
    abline(v=orig[i], col='magenta', lty=3, lwd=2)
    
    #ATA effect blocking
    resc <- varylist[[i]][varylist[[i]]$epsECbrk>varylist[[i]]$IGR & varylist[[i]]$IGR>0,]
    rect(resc[1,10], ymin*1.2, tail(resc,1)[1,10], ymax*1.2, col="darkgoldenrod2", border=NA, density=40, lty=3)
    
    excl <- varylist[[i]][varylist[[i]]$epsECbrk<varylist[[i]]$IGR & varylist[[i]]$IGR<0,]
    rect(excl[1,10], ymin, tail(excl,1)[1,10], ymax, col="hotpink2", border=NA, density=40, lty=3)
    
    #label
    mtext(paste0("(", letters[i],")"), 3, -2.5, adj=0.985, cex=1.8)
    #legend
    if (i==1){
      legend("topleft", legend=c(expression(Delta[i]^0), expression(Delta[i]^E), expression(Delta[i]^C), expression(Delta[i]^"(E#C)"), expression(Delta[i]^"[E||C]"), expression(Delta[i]^"[EC]") ,expression(IGR)), 
             col = c(rep("black",4), "blue", "red", "orange"),lty = line, bty="n", cex=1.9, inset=c(-0.02,-0.03), 
             y.intersp = 1.1, x.intersp = 0.1, seg.len=0.8, lwd=1.8)
    }
  }
  title(ylab="contribution to coexistence", outer=T, line=2.1, cex.lab=2, font.lab=2)
  dev.off()
}