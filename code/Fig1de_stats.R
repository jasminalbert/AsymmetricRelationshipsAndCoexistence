#This scripts computes the statistics associated with Fig 1 d and e, which are the
#panels showing ATAs visually for plankton population data. We are showing those 
#panels just as a demo of ATAs early in the paper.
#
#Reuman
#Begun 2023 07 24

#***
#Setup
#***

rm(list=ls())

source("partCor.R")

nsims<-10000
resloc<-"../results_numeric/"

set.seed(101)

#***
#load data and take a look
#***

dat <- readRDS("../results_numeric/plankStats.RDS")

# class(dat)
# names(dat)
# class(dat$e)
# names(dat$e)
# class(dat$e$uv)
# dim(dat$e$uv)
# plot(dat$e$uv[,1],dat$e$uv[,2],type="p")
# #this looks like the right tail associated plot
# 
# class(dat$d)
# names(dat$d)
# class(dat$d$uv)
# dim(dat$d$uv)
# plot(dat$d$uv[,1],dat$d$uv[,2],type="p")
# #this looks like the left tail associated plot

panelD<-dat$d$uv
panelE<-dat$e$uv

#***
#Get the difference of left minus right partial correlation for the real data. So
#positive values correspond to left tail association, negative to right
#***

resD<-partCor(panelD[,1],panelD[,2],.5,"average")$diff
#resD
resE<-partCor(panelE[,1],panelE[,2],.5,"average")$diff
#resE

#***
#Now compare with an appropriate symmetric null hypothesis
#***

#***start with the d panel

#get the copula to be simulated
ncop1<-copula::normalCopula(.5)
parmD<-copula::iRho(ncop1,cor(panelD[,1],panelD[,2],method="spearman"))
ncop2<-copula::normalCopula(parmD)

#simulate it a bunch of times and for each time compute the difference of partial Spearmans
resDsurr<-NA*numeric(nsims)
for (counter in 1:nsims)
{
  datsurr<-copula::rCopula(dim(panelD)[1],ncop2)
  resDsurr[counter]<-partCor(datsurr[,1],datsurr[,2],.5,"average")$diff
}

#now get a p-value
#hist(resDsurr)
#points(resD,0,col="red")
pD<-sum(resDsurr>resD)/nsims

#***now do the e panel

#get the copula to be simulated
ncop1<-copula::normalCopula(.5)
parmE<-copula::iRho(ncop1,cor(panelE[,1],panelE[,2],method="spearman"))
ncop2<-copula::normalCopula(parmE)

#simulate it a bunch of times and for each time compute the difference of partial Spearmans
resEsurr<-NA*numeric(nsims)
for (counter in 1:nsims)
{
  datsurr<-copula::rCopula(dim(panelE)[1],ncop2)
  resEsurr[counter]<-partCor(datsurr[,1],datsurr[,2],.5,"average")$diff
}

#now get a p-value
#hist(resEsurr)
#points(resE,0,col="red")
pE<-sum(resEsurr<resE)/nsims
#Note: we are using 1-tailed tests here

#***
#Save the various results into RDS files for later uptake into the latex
#***

saveRDS(resD,paste0(resloc,"DiffPartialSpearmanPanelD.Rds"))
saveRDS(resE,paste0(resloc,"DiffPartialSpearmanPanelE.Rds"))
saveRDS(pD,paste0(resloc,"DiffPartialSpearman_pval_PanelD.Rds"))
saveRDS(pE,paste0(resloc,"DiffPartialSpearman_pval_PanelE.Rds"))
saveRDS(parmD,paste0(resloc,"DiffPartialSpearman_parm_PanelD.Rds"))
saveRDS(parmE,paste0(resloc,"DiffPartialSpearman_parm_PanelE.Rds"))
saveRDS(nsims,paste0(resloc,"DiffPartialSpearmannsims.Rds"))
