#This script contains functions for making data for figure 6 and plotting
#6 functions:
  #mapDat, cmMap, cmContour, dat6_wrap, dat6, fig6

##libraries used (invoked with ::): parallel, plot3D, grDevices, graphics

### mapDat ###
# description
#ARGS:
  #parmlist   list length three containing a, P, and Tbar, intended that one element is a single numeric
  #sims       number of simulations for gerating partial sharp noise and ATA contribution
  #time       length of competition simulation
  #invader    1 or 2; which species is the rare species?
#OUT:
  #results data.frame with aPT combinations, GWR, ATA contribution, time, ATA/GWR
mapDat <- function(parmlist, sims, time, invader){
  
  #separate list elements
  a <- parmlist[[1]]
  P <- parmlist[[2]]
  Tb <- parmlist[[3]]
  
  #empty list
  argsList <- list()
  r <- 1
  for (i in 1:length(a)){
    for (j in 1:length(P)){
      for (k in 1:length(Tb)){
        #put aPT in string with other parameters
        #each unique set is a list element
        argsList[[r]] <- c(a[i], P[j], Tb[k], time, sims, invader)
        r <- r+1
      }	
    }
  }
  
  #use mclapply to quicken computation of Deltas
  resList <- parallel::mclapply(argsList, wrapDelt, mc.cores=7) ###change cores here###
  
  #extract results of interest into one data.frame
  #map is ATA/GWR
  resdf <- data.frame(t(sapply(resList, function(X){c(X$GWR, X$ECbrk, X$time, X$map)})))
  
  #extract aPT combinations and 
  parmsdf <- data.frame(matrix(unlist(argsList), ncol=6, byrow=T)[,-4:-6])
  
  #ammend to results data.frame
  res <- cbind(parmsdf, resdf)
  colnames(res) <- c("a", "P", "Tbar", "GWR", "ATA", "time", "map")
  return(res)
}

### cmMap ###
# puts mapping values (ATA/GWR) in matrix for contour color map plotting 
#ARGS:
  #dat  data.frame obtained from mapDat
#OUT:
  #matrix of ratios with rows and columns and the two varying paramter sets
cmMap <- function(dat){
  # $map is ATA/GWR
  #constraining to [-1,1]
  dat$map[dat$map>1] <- rep(1, length(dat$map[dat$map>1]))
  dat$map[dat$map<(-1)] <- rep(-1, length(dat$map[dat$map<(-1)]))
  
  #get param sets back
  a <- unique(dat$a); P <- unique(dat$P); Tbar <- unique(dat$Tbar)
  
  #get dimensions for matrix
  dims <- list(a=round(a,3), P=round(P,3), Tbar=round(Tbar,3))
  vars <- dims[len!=1]
  
  map <- matrix(dat$map, nrow=length(vars[[1]]), ncol=length(vars[[2]]), dimnames=vars, byrow=T)
  
  return(map)
}

### cmContour ###
# description
#ARGS:
  #map      matrix obtained from cmMap
  #ncolor   number of colors, should be odd
  #colkey   color bar key, default to NULL
  #...      other arguments to image2D
cmContour <- function(map, ncolor=51, colkey=NULL,...){
  x <- dimnames(map)[2]
  y <- dimnames(map)[1]
  
  cm <- grDevices::cm.colors(ncolor)
  
  #aligning color index if range of values is not [-1,1] so that not all colors are used
  if(all(range(map)!=c(-1,1))){ 
    by <- 2/(ncolor-1) 
    ran0 <- matrix(seq(-1,1,by)) #full range
    ran <- range(map) #actual range
    rMatch <- c(round(ran[1]/by)*by, round(ran[2]/by)*by) #round to match 
    
    #working around method of subsetting because R lie about floats
    min <- apply(ran0, 1, function(X){all.equal(X, rMatch[1])}) 
    max <- apply(ran0, 1, function(X){all.equal(X, rMatch[2])}) 
    cm <- cm[which(min==TRUE):which(max==TRUE)] #cropped color palette..was there an easier way to do this?
  }
  
  #color map plotting, plotting values of matrix in color
  plot3D::image2D(z=t(map), y=as.numeric(y[[1]]), x=as.numeric(x[[1]]), contour=FALSE, col=cm, colkey=colkey, xlab='', ylab='',...)
  graphics::title(xlab=names(x), ylab=ifelse(names(y)=='Tbar', expression(theta[0]), names(y)), line=-1)
  
  #contour lines
  graphics::contour(z=t(map), y=as.numeric(y[[1]]), x=as.numeric(x[[1]]),add=TRUE, col='grey50')
}

### dat6_wrap ###
# a wrapper function *to use* in dat6
# combines mapDat and cmMap, get mapping matrix from list of parameters
#ARGS:
  #argsList   a list of parameters *made inside dat6 function
              #argsList = list(a_vec, P_vec, T_vec, sims, time)
#OUT:
  #mapping matrix 
dat6_wrap <- function(argsList){
  mapdf <- mapDat(parmlist=argsList[1:3], sims=argsList[[5]], time=argsList[[4]], invader=argsList[[6]])
  mapmat <- cmMap(mapdf)
}

### dat6 ###
# uses above functions to generate and save results for making figure 6
#ARGS:
  #dat_loc    folder to put data
  #a_vec      vector containing set of amplitudes
  #P_vec      vector containing set of Periods
  #T_vec      vector containing set of mean temperatures
  #parms      named vector of original diatom coexistence parameters and such
              #parms = c('a'=6, 'P'=60, 'Tbar'=18, 'time'=3000, 'reps'=200, 'invader'=1)
#OUT:
  #folder containg .RDS files each containg data for fig 6
dat6 <- function(dat_loc, a_vec, P_vec, T_vec, parms){
  
  #a,P, and T vectors should be same length
  if((length(a_vec)==length(P_vec) & length(P_vec)==length(T_vec)) == FALSE){
    paste("vectors are not equal length")
  }
  plots <- 1:3
  aPT <- list(a=a_vec, P=P_vec, Tb=T_vec)
  varcom <- list(c(1,2), c(1,3), c(2,3))
  names <- c('a','P','Tb')
  
  for (p in plots){ #for each panel of figure
    List <- as.list(parms)
    
    #replace two of the constant variables with two varying sets
    List[varcom[[p]]] <- aPT[varcom[[p]]]
    
    #get mapping matrix from those parametrs
    map <- dat6_wrap(argsList=List)
    
    #save
    saveRDS(map,paste0(dat_loc,paste0(names[varcom[[p]]], collapse=''),parms['invader'],'.RDS'))
  }
}

### fig6 ###
# produce figure 6 as PDF from .RDS files
#ARGS:
  #filename   name of PDF
  #dat_loc    folder location of .RDS files
  #invader    1 or 2; which species is rare?
fig6 <- function(filename, dat_loc, invader){
  names <- c('a','P','Tb')
  varcom <- list(c(1,2), c(1,3), c(2,3))
  plots <- 1:3
  maps <- list()
  
  #load .RDS files
  for (p in plots){
    maps[[p]] <- readRDS(paste0(dat_loc,paste0(names[varcom[[p]]], collapse=''),invader,'.RDS'))
  }
  
  #start figure
  grDevices::pdf(filename, height=15, width=5)
  graphics::par(mfrow=c(1,1), oma=c(3,0,1,0), mar=c(2,3,1,1), bty='n', xpd=T)
  graphics::layout(matrix(c(1,2,3,4), byrow=T), heights=c(1,1,1,0.1))
  
  #make countor color plots
  for (i in 1:length(maps)){
    cmContour(maps[[i]], colkey=F, xaxt='n', yaxt='n')
    graphics::axis(side=1, mgp=c(3,0.5,0.2), col='gray50')
    graphics::axis(side=2, mgp=c(3,0.5,0.2), col='gray50')
    graphics::mtext(paste0("(", letters[i],")"), side=3, line=-1.5, adj=0.985)
  }
  #color bar key
  plot3D::colkey(col=cm.colors(51), clim=c(-1,1),side=1, width=10,)
  
  graphics::title(main=expression(paste(Delta[i]^"[EC]","/IGR")), outer=T, line=-111.5, cex.main=1.5, adj=0.53)
  grDevices::dev.off()
}