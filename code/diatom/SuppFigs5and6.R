#If desired, this script can be used to produce versions of 
#fig 5 and fig 6 where Cyclotella psuedosteligaria 
#is the rare species instead of Fragilaria crotonensis

##libraries used (invoked with ::): none

### sourcing ###
source("./diatom/diatomDecomp_fxns.R")
source("./diatom/fig5_fxns.R")
source("./diatom/fig6_fxns.R")

### location to save figs ###
fig_loc <- "../results_figs/"
if(dir.exists(fig_loc)==FALSE){
  dir.create(fig_loc)
}

### define parameters (invader=2 this time) ###
parms <- c('a'=6, 'P'=60, 'Tbar'=18, 'time'=3000, 'reps'=200, 'invader'=2)


### supplemental version of fig 5 ###

## location of numeric results ## 
dat_loc <- "../results_numeric/fig5dat2/"

#if data is missing, make them
if (dir.exists(dat_loc)==FALSE){
  dir.create(dat_loc)
  
  #define sets of variables
  a <- seq(1,6,length.out=100)
  Tbar <- seq(16,18,length.out = 100)
  P <- seq(51,199.5,length.out=100)
  
  #makes results as folder of .RDS files
  dat5(dat_loc, a_vec=a, P_vec=P, T_vec=Tbar, parms=parms)
}
## name of figure file ##
fig5_2 <- paste(fig_loc,"fig5_2.pdf",sep="")

## makes figure ##
fig5(fig5_2, dat_loc, invader=2) 



### supplemental version of fig 6 ###

## location of numeric results ## 
dat_loc <- "../results_numeric/fig6dat2/"

#if data is missing, make them
if (dir.exists(dat_loc)==FALSE){
  dir.create(dat_loc)
  
  #define sets of variables
  a <- seq(3.5,6,length.out=100)
  P <- seq(51,199.5,length.out=100)
  Tbar <- seq(16,19,length.out=100)
  
   #makes results as folder of .RDS files
  dat6(dat_loc, a_vec=a, P_vec=P, T_vec=Tbar, parms)
}
## name of figure file ##
fig6_2 <- paste0(fig_loc,"fig6_2.pdf")

## makes figure ##
fig6(fig6_2, dat_loc, invader=2) 









