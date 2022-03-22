#This script makes figure 2 in main text: population simulations with and without ATAs 

### source funtions ###
#script containing functions to make plot
source("./popsim_fxns.R")

### location to save results ###
noise_loc <- "../results_numeric/noise.RData"
fig_loc <- "../results_figs/"
if(dir.exists(fig_loc)==FALSE){
  dir.create(fig_loc)
}
fig2_loc <- paste(fig_loc,"fig2.pdf",sep="")
fig2data_loc <- "../results_numeric/fig2data.RDS"

### make fig ###
pdf(fig2_loc) 
par(mfrow=c(2,1), mar=c(3,3,1,1), oma=c(2,1,1,1))

fig2 <- simsPlot(noise_loc, mudif=-0.8, delta=0.8, sigma=5, end=1000) # r_woATA = 0.098, r = -0.0373, [E||C] = -0.1354

dev.off()

### save data ###
saveRDS(fig2, file=fig2data_loc)


