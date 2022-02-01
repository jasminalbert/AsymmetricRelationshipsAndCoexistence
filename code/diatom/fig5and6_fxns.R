source('./diatom/workflow.R')
#***FUNCTIONS FOR MAKING DATA AND PLOTTING OF FIGURE 5 AND FIGURE 6***#
#######################################################################

# make dat functions #
# fig5
#parms = c('a'=6, 'P'=60, 'Tbar'=18, 'time'=3000, 'reps'=200, 'invader'=2)
mak5 <- function(dat_loc, a_vec, P_vec, T_vec, parms){
  if((length(a_vec)==length(P_vec) & length(P_vec)==length(T_vec)) == FALSE){
    paste("vectors are not equal length")
  }
  len <- length(a_vec)
  names <- c("amplitude","Period","Tbar")
  plots <- 1:3
  aPT <- cbind(a_vec, P_vec, T_vec)
  
  for (p in plots){
    Mat <- matrix(parms, nrow=len, ncol=length(parms), byrow=T)
    Mat[,p] <- aPT[,p]
    Delta <- wrapDelt(args=Mat[1,])
    
    for (i in 2:len){
      Delta[i,] <- wrapDelt(Mat[i,])
      print(i)
    }
    Delta$var <- aPT[,p]
    saveRDS(Delta, paste0(dat_loc, "vary", names[p],parms["invader"],".RDS"))
  }
} ###############################################################################
# fig 6
#argsList = list(a_vec, P_vec, T_vec, sims, time)
mak6_wrap <- function(argsList){
  mapdf <- mapspace(parmlist=argsList[1:3], sims=argsList[[5]], time=argsList[[4]], invader=argsList[[6]])
  mapmat <- cmMap(mapdf)
}
mak6 <- function(dat_loc, a_vec, P_vec, T_vec, parms){
  if((length(a_vec)==length(P_vec) & length(P_vec)==length(T_vec)) == FALSE){
    paste("vectors are not equal length")
  }
  plots <- 1:3
  aPT <- list(a=a_vec, P=P_vec, Tb=T_vec)
  varcom <- list(c(1,2), c(1,3), c(2,3))
  names <- c('a','P','Tb')
  
  for (p in plots){
    List <- as.list(parms)
    List[varcom[[p]]] <- aPT[varcom[[p]]]
    
    map <- mak6_wrap(argsList=List)
    saveRDS(map,paste0(dat_loc,paste0(names[varcom[[p]]], collapse=''),parms['invader'],'.RDS'))
  }
}#################test############
a_vec=seq(3.5,6,length.out=2); P_vec=seq(51,53,length.out=2); T_vec=seq(16,19,length.out=2)
parms = c('a'=6, 'P'=60, 'Tbar'=18, 'time'=3000, 'reps'=200, 'invader'=1)
dat_loc <- "../results_numeric/fig6dat2/"
######################################################
######################################################

# plotting functions #
# fig 5
fig5 <- function(filename, dat_loc, invader){
  names <- c("amplitude","Period","Tbar")
  plots <- 1:3
  varylist <- list()
  
  for (p in plots){
    varylist[[p]] <- readRDS(paste0(dat_loc, "vary", names[p],parms["invader"],".RDS"))
  }
  col <- c(rep("black",4), "red", "blue", "orange")
  line <- c(1:4, 1,1,1); lwd <- c(rep(2,6),3)
  
  range <- sapply(varylist, function(X){range(X[,1:7])})
  ymin <- min(range); ymax <- max(range)
  orig <- c(a=6, P=60, Tbar=18)
  xlab <- c("a", "P", expression(theta[0]))
  
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
}#########################################################################
# fig 6
fig6 <- function(filename, dat_loc, invader){
  names <- c('a','P','Tb')
  varcom <- list(c(1,2), c(1,3), c(2,3))
  plots <- 1:3
  maps <- list()
  
  for (p in plots){
    maps[[p]] <- readRDS(paste0(dat_loc,paste0(names[varcom[[p]]], collapse=''),invader,'.RDS'))
  }
  
  pdf(filename, height=15, width=5)
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
}