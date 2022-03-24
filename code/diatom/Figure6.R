#This script makes figure 6 in main text:
#

### source diatom model coexistence decomposition functions ###
source('./diatomDecomp_fxns.R')
### source functions for computing data and plotting of figure 6 ###
source('./diatom/fig6_fxns.R')

#define original parameters
parms <- c('a'=6, 'P'=60, 'Tbar'=18, 'time'=3000, 'reps'=200, 'invader'=1)

#location of data files
dat_loc <- "../results_numeric/fig6dat/"

#if data is missing, make them
if (dir.exists(dat_loc)==FALSE){
  
  #make folder
  dir.create(dat_loc)
  
  #define sets of variables
  a <- seq(1,6,length.out=100)
  Tbar <- seq(16,18,length.out = 100)
  P <- seq(51,199.5,length.out=100)
  
  #function that makes and saves results
  dat6(dat_loc, a_vec=a, P_vec=P, T_vec=Tbar, parms)
}

#### location to store figure ###
fig_loc <- "../results_figs/"
if(dir.exists(fig_loc)==FALSE){
  dir.create(fig_loc)
}
Fig6 <- paste0(fig_loc,"fig6.pdf")

### function that makes figure ###
fig6(Fig6, dat_loc, invader=1) 



