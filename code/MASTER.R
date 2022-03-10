#Running this script should reproduce all analyses for the paper. Assumes the R working directory is the 
#directory this file is in, i.e., the code directory for the project.

rm(list=ls())

#Make standard environmental noise (transformed to use in different simulations), change M=length of 
#simulations here, if desired 
source("./makenoise.R") 

#Set parameters sigma, mu1, mu2, delta for the lottery model
source("./parameters.R") 

#FIGURE 1
source("./Figure1.R")

#FIGURE 2
source("./Figure2.R")

#FIGURE 3
source("./Figure3.R")

#FIGURE 4
source("./Figure4.R")

#FIGURE 5
source("./diatom/Figure5.R")

#FIGURE 6
source("./diatom/Figure6.R")

#SI versions of 5&6 - not used in paper but created in support of statements in text
source("./diatom/SuppFigs5and6.R")

#VKQ figure (SI), shows assumptions of how V, K and Q depend on temperature 
source("./diatom/VKQplots.R")



