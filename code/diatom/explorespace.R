source("getEpsilons_fxn.R")
dat_loc <- "../results_numeric/"
#loop through amplitude? 
ep1<- get_epsilons(3, 60, "Voftemp", 500, 1)
ep2<- get_epsilons(3, 60, "Voftemp", 500, 2)
Delta1 <- ep1-ep2
sum(Delta1)-(Delta1$epsEC+Delta1$epsECpar)

a <- seq(0.5,7,0.5)
ep1 <- get_epsilons(a[1], 60, "temp", 500, 1)
ep2 <- get_epsilons(a[1], 60, "temp", 500, 2)
Delta1 <- ep1-ep2

for (i in 2:length(a)){
	ep1 <- get_epsilons(a[i], 60, "temp", 500, 1)
	ep2 <- get_epsilons(a[i], 60, "temp", 500, 2)
	Delta1[i,] <- ep1-ep2
	print(i)
}
Delta1 <- Delta1[,-c(4,6)]
Delta1$IGR <- rowSums(Delta1)
Delta1$amplitude <- a

saveRDS(Delta1, paste(dat_loc,'varyingAmplitude.RDS',sep=''))
round(Delta1[,c(5,7,8)],5)
#######################################################################

ep1 <- get_epsilons(4.415, 60, "temp", 500, 1)
ep2 <- get_epsilons(4.415, 60, "temp", 500, 2)
D <- ep1-ep2
D <- D[-c(4,6)]
D[7] <- sum(D)
D[7]-D[5];D[7]
#4.29-4.41 have ATA rescue
#######################################################################


#vary P at a=6
P <- seq(20,120,10)
ep1 <- get_epsilons(6, P[1], "temp", 500, 1)
ep2 <- get_epsilons(6, P[1], "temp", 500, 2)
DeltaP <- ep1-ep2

for (i in 2:length(P)){
	ep1 <- get_epsilons(6, P[i], "temp", 500, 1)
	ep2 <- get_epsilons(6, P[i], "temp", 500, 2)
	DeltaP[i,] <- ep1-ep2
	print(i)
}
DeltaP <- DeltaP[,-c(4,6)]
DeltaP$IGR <- rowSums(DeltaP)

round(DeltaP[,c(5,7,8)],5)
DeltaP$Period <- P
saveRDS(DeltaP, paste(dat_loc,'varyingPerioda-6.RDS',sep=''))

#######################################################################

Tbar <- seq(16,19,0.25)
ep1 <- getep2(6, 60, Tbar[1], 200, 1)
ep2 <- getep2(6, 60, Tbar[1], 200, 2)
Delta1 <- ep1-ep2

for (i in 2:length(Tbar)){
	ep1 <- getep2(6, 60, Tbar[i], 200, 1)
	ep2 <- getep2(6, 60, Tbar[i], 200, 2)
	Delta1[i,] <- ep1-ep2
	print(i)
}
Delta1$IGR <- rowSums(Delta1)
Delta1$Tbar <- Tbar
saveRDS(Delta1, paste(dat_loc,'varyTbar.RDS',sep=''))



