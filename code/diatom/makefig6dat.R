#make files for figure 6
source("./diatom/workflow.R")

dat_loc <- "../results_numeric/fig6dat/"
if(dir.exists(dat_loc)==FALSE){
  dir.create(dat_loc)
}
#################################################################################
# make into function#
#################################################################################
# argsList = list(a_vec, P_vec, T_vec, sims, time)
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
    saveRDS(map,paste0(dat_loc,paste0(names[varcom[[p]]], collapse=''),'.RDS'))
  }
}
#################test############
a_vec=seq(3.5,6,length.out=2); P_vec=seq(51,53,length.out=2); T_vec=seq(16,19,length.out=2)
parms = c('a'=6, 'P'=60, 'Tbar'=18, 'time'=3000, 'reps'=200, 'invader'=1)
dat_loc <- "../results_numeric/fig6dat2/"

##################################
# amplitude x period
a <- seq(3.5,6,length.out=100)
P <- seq(51,199.5,length.out=100)
Tb <- 18
parmlist <- list(a,P, Tb)

aP <- mapspace(parmlist, sims=200, time=3000) #makes dataframe 
aP <- cmMap(aP)
saveRDS(aP, paste0(dat_loc,'aP.RDS'))

# amplitude x Tbar
a <- seq(3.5,6,length.out=100)
P<- 60
Tb <- seq(12, 24, length.out=100)
parmlist <- list(a ,P, Tb)

aTb <- mapspace(parmlist, sims=200, time=3000)
aTb <- cmMap(aTb)
saveRDS(aTb, paste0(dat_loc,'aTb.RDS'))

# Tbar x period
a <- 6
P <- seq(51,199.5,length.out=100)
Tb <- seq(12, 24, length.out=100)
parmlist <- list(a,P, Tb)

pTb <- mapspace(parmlist, sims=200, time=3000)
pTb <- t(cmMap(pTb))
saveRDS(pTb, paste0(dat_loc,'pTb.RDS'))

