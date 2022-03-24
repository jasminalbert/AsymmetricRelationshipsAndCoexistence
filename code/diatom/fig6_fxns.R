#This script contains functions for making data for figure 6 and plotting


require(parallel)

mapspace <- function(parmlist, sims, time, invader){
  
  #names(parmlist) <- c('a','P','Tbar')
  #len <- unlist(lapply(parmlist, length))
  #varlen <- len[len!=1]
  #rows <- varlen[[1]]*varlen[[2]]
  #res <- matrix(NA, ncol=6, nrow=rows)
  
  a <- parmlist[[1]]
  P <- parmlist[[2]]
  Tb <- parmlist[[3]]
  
  argsList <- list()
  r <- 1
  for (i in 1:length(a)){
    for (j in 1:length(P)){
      for (k in 1:length(Tb)){
        #Delta1 <- getDelt(a[i], P[j], Tb[k], time, sims)
        #res[r,] <- c(a[i], P[j], Tb[k], Delta1$IGR, Delta1$epsECbrk, Delta1$time)
        #print(c(r,i,j,k))
        argsList[[r]] <- c(a[i], P[j], Tb[k], time, sims, invader)
        r <- r+1
      }	
    }
  }
  #res <- data.frame(res)
  resList <- mclapply(argsList, wrapDelt, mc.cores=7) ###change cores here###
  resdf <- data.frame(t(sapply(resList, function(X){c(X$IGR, X$epsECbrk, X$time, X$map)})))
  parmsdf <- data.frame(matrix(unlist(argsList), ncol=6, byrow=T)[,-4:-6])
  
  res <- cbind(parmsdf, resdf)
  colnames(res) <- c("a", "P", "Tbar", "IGR", "ATA", "time", "map")
  return(res)
}

require(plot3D)

cmMap <- function(dat){
  dat$map[dat$map>1] <- rep(1, length(dat$map[dat$map>1]))
  dat$map[dat$map<(-1)] <- rep(-1, length(dat$map[dat$map<(-1)]))
  
  a <- unique(dat$a); P <- unique(dat$P); Tbar <- unique(dat$Tbar)
  dims <- list(a=round(a,3), P=round(P,3), Tbar=round(Tbar,3))
  
  len <- unlist(lapply(dims, length))
  vars <- dims[len!=1]
  
  map <- matrix(dat$map, nrow=length(vars[[1]]), ncol=length(vars[[2]]), dimnames=vars, byrow=T)
  
  return(map)
}

cmContour <- function(map, ncolor=51, colkey=NULL,...){
  x <- dimnames(map)[2]
  y <- dimnames(map)[1]
  #main <- expression(paste(Delta[i]^"[EC]","/IGR"))
  cm <- cm.colors(ncolor)
  
  if(all(range(map)!=c(-1,1))){ #need to match actual range to full range to adj color index
    by <- 2/(ncolor-1) 
    ran0 <- matrix(seq(-1,1,by)) #full range
    ran <- range(map) #actual range
    rMatch <- c(round(ran[1]/by)*by, round(ran[2]/by)*by) #round to match 
    
    min <- apply(ran0, 1, function(X){all.equal(X, rMatch[1])}) #why does R have to lie about floats
    max <- apply(ran0, 1, function(X){all.equal(X, rMatch[2])}) 
    cm <- cm[which(min==TRUE):which(max==TRUE)] #cropped color palette..was there an easier way to do this?
    
    
    #silly function I wrote before I realized the round(a/b)*b trick 
    #also match() fucntion doesnt work because R lies about floats
    
    #roundby <- function(x, by){
    #m <- round(x)
    #if (m<x){i <- seq(m, m+1, by)
    #} else {i <- seq(m, m-1, -by)}
    #dd <- tail(i[x>i],1);du <- i[x<i][1]
    #if (x-dd < du-x){x <- dd
    #} else {x <- du}
    #return(x)
    #}
    #ran <- range(map)
    #rMatch <- c(roundby(ran[1],by), roundby(ran[2],by))
    #cm <- cm[match(rMatch,ran0)[1]:match(rMatch,ran0)[2]]
  }
  
  image2D(z=t(map), y=as.numeric(y[[1]]), x=as.numeric(x[[1]]), contour=FALSE, col=cm, colkey=colkey, xlab='', ylab='',...)
  title(xlab=names(x), ylab=ifelse(names(y)=='Tbar', expression(theta[0]), names(y)), line=-1)
  contour(z=t(map), y=as.numeric(y[[1]]), x=as.numeric(x[[1]]),add=TRUE, col='grey50')
}



### dat6_wrap ###
# a wrapper function *to use* in dat6
#ARGS:
  #argsList   a list of parameters *made inside dat6 function
              #argsList = list(a_vec, P_vec, T_vec, sims, time)
dat6_wrap <- function(argsList){
  #calls functions from diatomDecomp
  mapdf <- mapspace(parmlist=argsList[1:3], sims=argsList[[5]], time=argsList[[4]], invader=argsList[[6]])
  mapmat <- cmMap(mapdf)
}
dat6 <- function(dat_loc, a_vec, P_vec, T_vec, parms){
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
    
    map <- dat6_wrap(argsList=List)
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