#make files for figure 6
source("./diatom/workflow.R")

dat_loc <- "../results_numeric/fig6dat/"
if(dir.exists(dat_loc)==FALSE){
  dir.create(dat_loc)
}

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
Tb <- seq(50, 220, length.out=100)

parmlist <- list(a,P, Tb)
aTb <- mapspace(parmlist, sims=200, time=3000)
aTb <- cmMap(aTb)
saveRDS(aTb, paste0(dat_loc,'aTb.RDS'))

# Tbar x period
a <- 6
P <- seq(51,199.5,length.out=100)
Tb <- seq(50, 220, length.out=100)

parmlist <- list(a,P, Tb)
pTb <- mapspace(parmlist, sims=200, time=3000)
pTb <- t(cmMap(pTb))
saveRDS(pTb, paste0(dat_loc,'pTb.RDS'))

