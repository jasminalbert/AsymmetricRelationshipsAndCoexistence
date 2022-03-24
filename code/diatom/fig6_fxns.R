#This script contains functions for making data for figure 5 and plotting

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