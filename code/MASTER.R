#MASTER - assumes the R working directory is the directory this file is in

rm(list=ls())

source("./makenoise.R") #make standard noise, change M here 

source("./parameters.R") #set parameters sigma, mu1, mu2, delta

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




