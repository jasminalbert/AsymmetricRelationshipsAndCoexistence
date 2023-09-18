#This script computes results for 
#panels d-f of figure 5 in main text:

##libraries used (invoked with ::): none

### source diatom model coexistence decomposition functions ###
source('./diatom/diatomDecomp_fxns.R')
### source functions for computing data and plotting of figure 6 ###
source('./diatom/fig5def_fxns.R')

#define original parameters
parms <- c('a'=6, 'P'=60, 'Tbar'=18, 'time'=3000, 'reps'=200, 'invader'=1)

### location to save results ##
numRes_loc <- "../results_numeric/"
dat_loc <- paste0(numRes_loc,"fig6dat/")

#if folder is missing, make it
if (dir.exists(dat_loc)==FALSE){
  #make folder
  dir.create(dat_loc)
}

#define sets of variables
a <- seq(4,6,length.out=101)
Tbar <- seq(16,18,length.out = 101)
P <- seq(50,120,length.out=101)

#function that makes and saves results
dat6(dat_loc, a_vec=a, P_vec=P, T_vec=Tbar, parms)

#### location to store figure ###
#fig_loc <- "../results_figs/"
#if(dir.exists(fig_loc)==FALSE){
#  dir.create(fig_loc)
#}
#Fig6 <- paste0(fig_loc,"fig6.pdf")

### function that makes figure ###
#fig6(Fig6, dat_loc, invader=1) 

### get and save standard error ###
#fig6se <- getSE(dat_loc)
#fig6maxse <- max(fig6se)
#cat("\nmaximum standard error in figure six is", fig6maxse)
#fig6maxse_loc <- paste0(numeric_results_loc, "fig6maxse.RDS")
#saveRDS(fig6maxse, file=fig6maxse_loc)


