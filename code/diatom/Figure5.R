#This script makes figure 5 in main text:
#

##libraries used (invoked with ::): none

### source diatom competition model coexistence decomposition functions ###
source('./diatom/diatomDecomp_fxns.R')
### source functions for computing data for and plotting figure 5
source('./diatom/fig5_fxns.R')

#define original parameters
parms <- c('a'=6, 'P'=60, 'Tbar'=18, 'time'=3000, 'reps'=200, 'invader'=1)

### location of results ###
numeric_results_loc <- "../results_numeric/"
dat_loc <- paste0(numeric_results_loc, "fig5dat/")

# if results are missing, make them
if (dir.exists(dat_loc)==FALSE){
  
  #make folder
  dir.create(dat_loc)
  
  #define sets of variables
  a <- seq(1,6,length.out=100)
  Tbar <- seq(16,18,length.out = 100)
  P <- seq(51,199.5,length.out=100)
  
  #function that produces results and saves .RDS's to folder
  dat5(dat_loc, a_vec=a, P_vec=P, T_vec=Tbar, parms=parms)
}

### location to save figure ###
fig_loc <- "../results_figs/"
if(dir.exists(fig_loc)==FALSE){
  dir.create(fig_loc)
}
Fig5 <- paste(fig_loc,"fig5.pdf",sep="")

### function that makes figure ###
fig5(Fig5, dat_loc, invader=1) 

### get and save standard error ###
fig5se <- getSE(dat_loc)
fig5maxse <- max(fig5se)
cat("\nmaximum standard error in figure five is", fig5maxse)
fig5maxse_loc <- paste0(numeric_results_loc, "fig5maxse.RDS")
saveRDS(fig5maxse, file=fig5maxse_loc)










