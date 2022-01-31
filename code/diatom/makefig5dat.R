source("diatom/workflow.R")
dat_loc <- "../results_numeric/fig5dat/"
if(dir.exists(dat_loc)==FALSE){
  dir.create(dat_loc)
}
#amplitude
a <- seq(1,6,length.out=100)
Delta_a <- getDelt(a[1], 60, 18, 3000, 200)

for (i in 2:length(a)){
	Delta_a[i,] <- getDelt(a[i], 60, 18, 3000, 200)
	print(i)
}

Delta_a$a <- a
saveRDS(Delta_a, paste(dat_loc,'varyingAmplitude.RDS',sep=''))


#######################################################################

#vary P at a=6
P <- seq(51,199.5,length.out=100)
DeltaP <- getDelt(6, P[1], 18, 3000, 200)

for (i in 2:length(P)){
  print(paste(i, P[i]))
  DeltaP[i,] <- getDelt(6, P[i], 18, 3000, 200)
}

DeltaP$P <- P

saveRDS(DeltaP, paste(dat_loc,'varyingPeriod.RDS',sep=''))

#######################################################################

Tbar <- seq(16,19,length.out = 100)

DeltaT <- getDelt(6, 60, Tbar[1], 3000, 200)

for (i in 2:length(Tbar)){
	DeltaT[i,] <- getDelt(6, 60, Tbar[i], 3000, 200)
	print(i)
}


DeltaT$Tbar <- Tbar

saveRDS(DeltaT, paste(dat_loc,'varyTbar.RDS',sep=''))

#######################################################################
# make into function#
#make fig 5 dat function

# args
# dat_loc
# a_vec
# P_vec
# T_vec
# len

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
    saveRDS(Delta, paste0(dat_loc, "vary", names[p], ".RDS"))
  }
}
################################################################################
dat_loc <- "../results_numeric/fig5dat2/"
if(dir.exists(dat_loc)==FALSE){
  dir.create(dat_loc)
}
parms = c('a'=6, 'P'=60, 'Tbar'=18, 'time'=3000, 'reps'=200, 'invader'=1)
mak5(dat_loc, a_vec=a, P_vec=P, T_vec=Tbar, parms=parms)

