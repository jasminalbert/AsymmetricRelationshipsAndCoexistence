source("./diatom/fig5and6_fxns.R")
#       SI figs       #
#######################
fig_loc <- "../results_figs/"
if(dir.exists(fig_loc)==FALSE){
  dir.create(fig_loc)
}
parms <- c('a'=6, 'P'=60, 'Tbar'=18, 'time'=3000, 'reps'=200, 'invader'=2)

#######################
# SI version of fig 5 #
#######################
# make data here 
dat_loc <- "../results_numeric/fig5dat2/"
if (dir.exists(dat_loc)==FALSE){
  a <- seq(1,6,length.out=100)
  Tbar <- seq(16,18,length.out = 100)
  P <- seq(51,199.5,length.out=100)
  mak5(dat_loc, a_vec=a, P_vec=P, T_vec=Tbar, parms=parms)
}
fig5_2 <- paste(fig_loc,"fig5_2.pdf",sep="")

fig5(fig5_2, dat_loc, invader=2) #makes figure


#######################
# SI version of fig 6 #
#######################
# make data here
dat_loc <- "../results_numeric/fig6dat2/"
if (dir.exists(dat_loc)==FALSE){
  a <- seq(3.5,6,length.out=100)
  P <- seq(51,199.5,length.out=100)
  Tbar <- seq(16,19,length.out=100)
  mak6(dat_loc, a_vec=a, P_vec=P, T_vec=Tbar, parms)
}
fig6_2 <- paste0(fig_loc,"fig6_2.pdf")

fig6(fig6_2, dat_loc, invader=2) #makes figure








