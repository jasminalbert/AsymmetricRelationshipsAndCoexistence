#This is a new script as of 2023 07 06, hoping to demo to a referee
#and other readers some intuition behind the influence of ATAs.
#Reuman
#begun 2023 07 06 

#***
#Setup
#***

rm(list=ls())

#data locations, copied from the scripts that generated the noise
numRes_loc <- "../results_numeric/"

noise_loc <- paste0(numRes_loc, "noise.RData")
M_loc <- paste0(numRes_loc, "M.RDS")
rho_loc <- paste0(numRes_loc, "rho.RDS")

#locations to save resulting plots in
resloc<-"../results_figs/"

#***
#Load noise previously generated (standard normal marginals)
#***

#load and rename
load(noise_loc)
rm(u)
nl<-b[[1]]
nr<-b[[2]]
ns<-b[[3]]
rm(b)

# #check
# plot(nl[1:1000,1],nl[1:1000,2],type="p")
# plot(nr[1:1000,1],nr[1:1000,2],type="p")
# plot(ns[1:1000,1],ns[1:1000,2],type="p")
# cor(nl)
# cor(nr)
# cor(ns)

#cut from 1000000 to numpoints points
numpoints<-1000
nl<-nl[1:1000,]
nr<-nr[1:1000,]
ns<-ns[1:1000,]

#***
#Parameters for a particular beta-fecundities model
#***

eta2<-1.2
eta1<-1
bf_delta<-0.6

#***
#Calculate E1 and C_{2\1} for the cases of right-tail, left-tail, and symmetric noise, beta-fecundities
#***

bf_E1l<-eta1*qbeta(pnorm(nl[,1]),0.5,0.5)
bf_C1bs1l<-eta2*qbeta(pnorm(nl[,2]),0.5,0.5)/bf_delta

bf_E1s<-eta1*qbeta(pnorm(ns[,1]),0.5,0.5)
bf_C1bs1s<-eta2*qbeta(pnorm(ns[,2]),0.5,0.5)/bf_delta

bf_E1r<-eta1*qbeta(pnorm(nr[,1]),0.5,0.5)
bf_C1bs1r<-eta2*qbeta(pnorm(nr[,2]),0.5,0.5)/bf_delta

#ATA contributions to coexistence, for checking against Fig. 3 
mean(log(1-bf_delta+bf_E1l/bf_C1bs1l))-mean(log(1-bf_delta+bf_E1s/bf_C1bs1s)) #ATA contribution to coexistence for the left-tail noise
mean(log(1-bf_delta+bf_E1r/bf_C1bs1r))-mean(log(1-bf_delta+bf_E1s/bf_C1bs1s)) #ATA contribution to coexistence for the right-tail noise

#***
#Parameters for a particular log-normal-fecundities model
#***

sig<-1
lnf_delta<-0.6
mu1<-0
mu2<-0.5

#***
#Calculate E1 and C_{2\1} for the cases of right-tail, left-tail, and symmetric noise, log-normal-fecundities
#***

lnf_E1l<-exp(sig*nl[,1]+mu1)
lnf_C1bs1l<-exp(sig*nl[,2]+mu2)/lnf_delta

lnf_E1r<-exp(sig*nr[,1]+mu1)
lnf_C1bs1r<-exp(sig*nr[,2]+mu2)/lnf_delta

lnf_E1s<-exp(sig*ns[,1]+mu1)
lnf_C1bs1s<-exp(sig*ns[,2]+mu2)/lnf_delta

#ATA contributions to coexistence, for checking against Fig. 2 
mean(log(1-lnf_delta+lnf_E1l/lnf_C1bs1l))-mean(log(1-lnf_delta+lnf_E1s/lnf_C1bs1s)) #ATA contribution to coexistence 
mean(log(1-lnf_delta+lnf_E1r/lnf_C1bs1r))-mean(log(1-lnf_delta+lnf_E1s/lnf_C1bs1s)) #ATA contribution to coexistence, should be almost the same 

#***
#Set up the overall plot
#***

#figure dimensions, inches
xaxht<-0.5
yaxwd<-0.5
gap<-0.3
totwd<-6.5
panwd<-(totwd-4*gap-yaxwd)/4
panht<-panwd
totht<-2*xaxht+2*panht+gap
cex_lab<-1.5
cex_lab_cont<-0.5

jpeg(file=paste0(resloc,"PostHocExplanatoryFig.jpg"),width=totwd,height=totht,units="in",quality=90,res=400)

#***
#Now generate the distributions of r_{i\i} for left tail noise and symmetric noise,
#and plot on top of a color image showing the function ln(1-delta+E1/Cjbs1),
#for lognormal-fecundities
#***

E1ta<-lnf_E1l
C1bs1ta<-lnf_C1bs1l
E1s<-lnf_E1s
C1bs1s<-lnf_C1bs1s

#***linear scale first

#set up the plotting panel
i<-1
j<-1
par(fig=c((yaxwd+(j-1)*(panwd+gap))/totwd,
          (yaxwd+(j-1)*(panwd+gap)+panwd)/totwd,
          (xaxht+(2-i)*(panht+xaxht))/totht,
          (xaxht+(2-i)*(panht+xaxht)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25)

#first make the color image
x<-seq(from=0.01,to=120.01,by=0.1)
y<-x
z<-outer(x,y,function(x,y){log(1-bf_delta+x/y)})
breaks<-pretty(z,10)
image(x,y,z,xlab="E1",ylab="C1bs1",col = hcl.colors(length(breaks)-1, "YlOrRd", rev = TRUE),
      breaks=breaks) 
contour(x,y,z,nlevels=10,add=T,labcex=cex_lab_cont)
mtext(expression(E[1]),1,1.2)
mtext(expression(C[1/1]),2,1.2)
text(max(x),max(y),"A",cex=cex_lab,adj=c(1,1))

#now add the points for the left-tail noise
points(E1ta,C1bs1ta,type="p",pch=20,cex=.1,col="green")

#now add the points for the symmetric noise, using a different color
points(E1s,C1bs1s,type="p",pch=20,cex=.1)

#***then log scale

#set up the plotting panel
i<-2
j<-1
par(fig=c((yaxwd+(j-1)*(panwd+gap))/totwd,
          (yaxwd+(j-1)*(panwd+gap)+panwd)/totwd,
          (xaxht+(2-i)*(panht+xaxht))/totht,
          (xaxht+(2-i)*(panht+xaxht)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)

#first make the color image
log10x<-seq(from=-2,to=log10(120),by=0.01)
log10y<-log10x
z<-outer(10^log10x,10^log10y,function(x,y){log(1-bf_delta+x/y)})
breaks<-pretty(z,10)
image(log10x,log10y,z,xlab="log10E1",ylab="log10C1bs1",col = hcl.colors(length(breaks)-1, "YlOrRd", rev = TRUE),
      breaks=breaks) 
contour(log10x,log10y,z,nlevels=10,add=T,labcex=cex_lab_cont)
mtext(expression(log[10](E[1])),1,1.2)
mtext(expression(log[10](C[1/1])),2,1.2)
text(max(log10x),max(log10y),"B",cex=cex_lab,adj=c(1,1))

#now add the points for the symmetric noise, using a different color
points(log10(E1s),log10(C1bs1s),type="p",pch=20,cex=.1)

#now add the points for the left-tail noise
points(log10(E1ta),log10(C1bs1ta),type="p",pch=20,cex=.1,col="green")

#***
#Now generate the distributions of r_{i\i} for right tail noise and symmetric noise,
#and plot on top of a color image showing the function ln(1-delta+E1/Cjbs1),
#for lognormal-fecundities
#***

E1ta<-lnf_E1r
C1bs1ta<-lnf_C1bs1r
E1s<-lnf_E1s
C1bs1s<-lnf_C1bs1s

#***linear scale first

#set up the plotting panel
i<-1
j<-2
par(fig=c((yaxwd+(j-1)*(panwd+gap))/totwd,
          (yaxwd+(j-1)*(panwd+gap)+panwd)/totwd,
          (xaxht+(2-i)*(panht+xaxht))/totht,
          (xaxht+(2-i)*(panht+xaxht)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)

#first make the color image
x<-seq(from=0.01,to=120.01,by=0.1)
y<-x
z<-outer(x,y,function(x,y){log(1-bf_delta+x/y)})
breaks<-pretty(z,10)
image(x,y,z,xlab="E1",ylab="C1bs1",col = hcl.colors(length(breaks)-1, "YlOrRd", rev = TRUE),
      breaks=breaks) 
contour(x,y,z,nlevels=10,add=T,labcex=cex_lab_cont)
mtext(expression(E[1]),1,1.2)
text(max(x),max(y),"C",cex=cex_lab,adj=c(1,1))

#now add the points for the symmetric noise, using a different color
points(E1s,C1bs1s,type="p",pch=20,cex=.1)

#now add the points for the left-tail noise
points(E1ta,C1bs1ta,type="p",pch=20,cex=.1,col="green")

#***then log scale

#set up the plotting panel
i<-2
j<-2
par(fig=c((yaxwd+(j-1)*(panwd+gap))/totwd,
          (yaxwd+(j-1)*(panwd+gap)+panwd)/totwd,
          (xaxht+(2-i)*(panht+xaxht))/totht,
          (xaxht+(2-i)*(panht+xaxht)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)

#first make the color image
log10x<-seq(from=-2,to=log10(120),by=0.01)
log10y<-log10x
z<-outer(10^log10x,10^log10y,function(x,y){log(1-bf_delta+x/y)})
breaks<-pretty(z,10)
image(log10x,log10y,z,xlab="log10E1",ylab="log10C1bs1",col = hcl.colors(length(breaks)-1, "YlOrRd", rev = TRUE),
      breaks=breaks) 
contour(log10x,log10y,z,nlevels=10,add=T,labcex=cex_lab_cont)
mtext(expression(log[10](E[1])),1,1.2)
text(max(log10x),max(log10y),"D",cex=cex_lab,adj=c(1,1))

#now add the points for the symmetric noise, using a different color
points(log10(E1s),log10(C1bs1s),type="p",pch=20,cex=.1)

#now add the points for the left-tail noise
points(log10(E1ta),log10(C1bs1ta),type="p",pch=20,cex=.1,col="green")

#***
#Now generate the distributions of r_{i\i} for left tail noise and symmetric noise,
#and plot on top of a color image showing the function ln(1-delta+E1/Cjbs1),
#for beta-fecundities
#***

E1ta<-bf_E1l
C1bs1ta<-bf_C1bs1l
E1s<-bf_E1s
C1bs1s<-bf_C1bs1s

#***linear scale first

#set up the plotting panel
i<-1
j<-3
par(fig=c((yaxwd+(j-1)*(panwd+gap))/totwd,
          (yaxwd+(j-1)*(panwd+gap)+panwd)/totwd,
          (xaxht+(2-i)*(panht+xaxht))/totht,
          (xaxht+(2-i)*(panht+xaxht)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)

#first make the color image
mx<-max(E1ta,C1bs1ta)
mx<-1.05*mx
x<-seq(from=0.01,to=mx,by=0.01)
y<-x
z<-outer(x,y,function(x,y){log(1-bf_delta+x/y)})
breaks<-pretty(z,10)
image(x,y,z,xlab="E1",ylab="C1bs1",col = hcl.colors(length(breaks)-1, "YlOrRd", rev = TRUE),
      breaks=breaks) 
contour(x,y,z,nlevels=10,add=T,labcex=cex_lab_cont)
mtext(expression(E[1]),1,1.2)
text(max(x),max(y),"E",cex=cex_lab,adj=c(1,1))

#now add the points for the left-tail noise
points(E1ta,C1bs1ta,type="p",pch=20,cex=.1,col="green")

#now add the points for the symmetric noise, using a different color
points(E1s,C1bs1s,type="p",pch=20,cex=.1)

#***then log scale

#set up the plotting panel
i<-2
j<-3
par(fig=c((yaxwd+(j-1)*(panwd+gap))/totwd,
          (yaxwd+(j-1)*(panwd+gap)+panwd)/totwd,
          (xaxht+(2-i)*(panht+xaxht))/totht,
          (xaxht+(2-i)*(panht+xaxht)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)

#first make the color image
mx<-max(E1ta,C1bs1ta)
mx<-1.5*mx
log10x<-seq(from=-2,to=log10(mx),by=0.01)
log10y<-log10x
z<-outer(10^log10x,10^log10y,function(x,y){log(1-bf_delta+x/y)})
breaks<-pretty(z,10)
image(log10x,log10y,z,xlab="log10E1",ylab="log10C1bs1",col = hcl.colors(length(breaks)-1, "YlOrRd", rev = TRUE),
      breaks=breaks) 
contour(log10x,log10y,z,nlevels=10,add=T,labcex=cex_lab_cont)
mtext(expression(log[10](E[1])),1,1.2)
text(max(log10x),max(log10y),"F",cex=cex_lab,adj=c(1,1))

#now add the points for the symmetric noise, using a different color
points(log10(E1s),log10(C1bs1s),type="p",pch=20,cex=.1)

#now add the points for the left-tail noise
points(log10(E1ta),log10(C1bs1ta),type="p",pch=20,cex=.1,col="green")

#***
#Now generate the distributions of r_{i\i} for right tail noise and symmetric noise,
#and plot on top of a color image showing the function ln(1-delta+E1/Cjbs1),
#for beta fecundities
#***

E1ta<-bf_E1r
C1bs1ta<-bf_C1bs1r
E1s<-bf_E1s
C1bs1s<-bf_C1bs1s

#***linear scale first

#set up the plotting panel
i<-1
j<-4
par(fig=c((yaxwd+(j-1)*(panwd+gap))/totwd,
          (yaxwd+(j-1)*(panwd+gap)+panwd)/totwd,
          (xaxht+(2-i)*(panht+xaxht))/totht,
          (xaxht+(2-i)*(panht+xaxht)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)

#first make the color image
mx<-max(E1ta,C1bs1ta)
mx<-1.05*mx
x<-seq(from=0.01,to=mx,by=0.01)
y<-x
z<-outer(x,y,function(x,y){log(1-bf_delta+x/y)})
breaks<-pretty(z,10)
image(x,y,z,xlab="E1",ylab="C1bs1",col = hcl.colors(length(breaks)-1, "YlOrRd", rev = TRUE),
      breaks=breaks) 
contour(x,y,z,nlevels=10,add=T,labcex=cex_lab_cont)
mtext(expression(E[1]),1,1.2)
text(max(x),max(y),"G",cex=cex_lab,adj=c(1,1))

#now add the points for the left-tail noise
points(E1ta,C1bs1ta,type="p",pch=20,cex=.1,col="green")

#now add the points for the symmetric noise, using a different color
points(E1s,C1bs1s,type="p",pch=20,cex=.1)

#***then log scale

#set up the plotting panel
i<-2
j<-4
par(fig=c((yaxwd+(j-1)*(panwd+gap))/totwd,
          (yaxwd+(j-1)*(panwd+gap)+panwd)/totwd,
          (xaxht+(2-i)*(panht+xaxht))/totht,
          (xaxht+(2-i)*(panht+xaxht)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)

#first make the color image
mx<-max(E1ta,C1bs1ta)
mx<-1.5*mx
log10x<-seq(from=-2,to=log10(mx),by=0.01)
log10y<-log10x
z<-outer(10^log10x,10^log10y,function(x,y){log(1-bf_delta+x/y)})
breaks<-pretty(z,10)
image(log10x,log10y,z,xlab="log10E1",ylab="log10C1bs1",col = hcl.colors(length(breaks)-1, "YlOrRd", rev = TRUE),
      breaks=breaks) 
contour(log10x,log10y,z,nlevels=10,add=T,labcex=cex_lab_cont)
mtext(expression(log[10](E[1])),1,1.2)
text(max(log10x),max(log10y),"H",cex=cex_lab,adj=c(1,1))

#now add the points for the symmetric noise, using a different color
points(log10(E1s),log10(C1bs1s),type="p",pch=20,cex=.1)

#now add the points for the left-tail noise
points(log10(E1ta),log10(C1bs1ta),type="p",pch=20,cex=.1,col="green")

dev.off()