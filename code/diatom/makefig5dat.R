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

Tbar <- seq(16,18,length.out = 100)

DeltaT <- getDelt(6, 60, Tbar[1], 3000, 200)

for (i in 2:length(Tbar)){
	DeltaT[i,] <- getDelt(6, 60, Tbar[i], 3000, 200)
	print(i)
}


DeltaT$Tbar <- Tbar

saveRDS(DeltaT, paste(dat_loc,'varyTbar.RDS',sep=''))

#######################################################################


