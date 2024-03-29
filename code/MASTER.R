#Running this script should reproduce all analyses for the paper. Assumes the R working directory is the 
#directory this file is in, i.e., the code directory for the project.

rm(list=ls())

### libraries used (invoked with ::):
  # MASS, stats, graphics, grDevices, 
  # copula, parallel, plot3D, deSolve

#Make standard environmental noise (transformed to use in different simulations), change M=length of 
#simulations here, if desired 
source("./makenoise.R") 

#transform b to beta distributed B
#includes plotting checks
source("./beta_transform.R")

#beta pop sim
#can change parameters for pop sim here
source("./beta_popsims.R")

#plankton data for fig 1
source("./plankton.R")

#Set parameters sigma, mu1, mu2, delta for the lottery model
source("./parameters.R") 

#FIGURE 1
source("./Figure1.R")

#FIGURE 2
source("./makTheoryFigs.R")

#FIGURE 3
source("./Figure3.R")

#FIGURE 4
source("./LBdecompFig.R")
source("./LBdecompfig_qij.R")

#FIGURE 5
#source("./diatom/Figure5.R")
#source("./diatom/Figure5def.R")
source("./diatom/fig5NEW.R")

#SI versions of 5&6 - not used in paper but created in support of statements in text
#source("./diatom/SuppFigs5and6.R")

#VKQ figure (SI), shows assumptions of how V, K and Q depend on temperature 
source("./diatom/VKQplots.R")

#Stats the referee asked us to add regarding Fig 1d, e
source("Fig1de_stats.R")

#Figurs supporting new explanatory material referee asked for
source("ATAEffectIntuitionDemo_betaFecund.R")
source("ATAEffectIntuitionDemo_lognormalFecund.R")