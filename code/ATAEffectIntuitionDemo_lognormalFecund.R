#This script to help explain ATA intuitively. Uses the lognormal fecundities lottery 
#model case to make a figure supporting text added to the Discussion and SI
#following a referee request.
#
#Reuman

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

#graphical params
numhistbins<-20
numpoints<-1000

#***
#Functions
#***

#A function that makes it relatively easier for me to make the barplots
#I want to make. Adds to a plot assumed to already be open.
#
#Args
#m      A vector of counts
#bks    Breaks. Length is length(m)+1
#col    A colors
#base   bottoms of rectangles
#top    top of highest rectangle
#axis   Use 1 for plotting the bars vertically, using the x axis for breaks.
#         Use 2 for plotting the bars horizontally, using the y axis for breaks.
#reddot Put a red dot at this location on the base axis
#
#Output - none, but adds to a plot
#
#Notes. 
#If the arguments bks, base, top do not comport with the axes already
#open then the barplot will appear (invisibly) off the plotting axes.
#
mybarplot<-function(m,bk,col,base,top,axis,reddot)
{
  Mx<-max(m)
  
  if (axis==1)
  {
    for (rc in 1:length(m))
    {
      rect(bk[rc],base,bk[rc+1],base+(top-base)*m[rc]/Mx,col=col,border=NA)
    }
    points(reddot,base,col="red",pch=20,cex=0.5)
  }
  
  if (axis==2)
  {
    for (rc in 1:length(m))
    {
      rect(base,bk[rc],base+(top-base)*m[rc]/Mx,bk[rc+1],col=col,border=NA)
    }
    points(base,reddot,col="red",pch=20,cex=0.5)
  }
}

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
nl<-nl[1:1000,]
nr<-nr[1:1000,]
ns<-ns[1:1000,]

#***
#Parameters for a particular lognormal fecundities model
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
yaxwd<-0.6
gap<-0.3
totwd<-6.5
panwd<-(totwd-3*gap-3*yaxwd)/3
panht<-panwd
totht<-2*xaxht+2*panht+gap
cex_lab<-1.25
cex_lab_cont<-0.5

jpeg(file=paste0(resloc,"PostHocExplanatoryFig_lognormalFecund.jpg"),width=totwd,height=totht,units="in",quality=90,res=400)
#pdf(file=paste0(resloc,"PostHocExplanatoryFig.pdf"),width=totwd,height=totht)

#***
#Now generate the distributions for left tail noise, and plot on top of a 
#color image showing the function ln(1-delta+E1/C1bs1), for lognormal fecundities
#***

E1ta<-lnf_E1l
C1bs1ta<-lnf_C1bs1l

#set up the plotting panel
i<-1
j<-1
par(fig=c((yaxwd+(j-1)*(panwd+yaxwd))/totwd,
          (yaxwd+(j-1)*(panwd+yaxwd)+panwd)/totwd,
          (xaxht+(2-i)*(panht+xaxht))/totht,
          (xaxht+(2-i)*(panht+xaxht)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25)

#first make the color image
lims<-c(-1.9,2.3)
dlims<-diff(lims)
log10x<-seq(from=-1.9,to=2.3,by=0.01)
log10y<-log10x
z<-outer(10^log10x,10^log10y,function(x,y){log(1-lnf_delta+x/y)})
breaks<-pretty(z,10)
image(log10x,log10y,z,xlab="log10E1",ylab="log10C1bs1",col = hcl.colors(length(breaks)-1, "YlOrRd", rev = TRUE),
      breaks=breaks,xlim=lims,ylim=lims) 
contour(log10x,log10y,z,nlevels=10,add=T,labcex=cex_lab_cont)

#now add a reference line
lines(c(min(log10x),max(log10x)),c(min(log10y),max(log10y))-log10(lnf_delta)+log10(exp(mu2))-log10(exp(mu1)),type="l",lty="dashed",lwd=1.5)

#a divider line between the "tails"
#lnf_E1l<-exp(sig*nl[,1]+mu1)
#lnf_C1bs1l<-exp(sig*nl[,2]+mu2)/lnf_delta
#line thru log10(exp(mu1)),log10(exp(mu2))-log10(lnf_delta), with slope -1
#(y-log10(exp(mu2))+log10(lnf_delta))=-(x-log0(exp(mu1)))
#y=-x+log10(exp(mu1))+log10(exp(mu2))-log10(lnf_delta)
abline(log10(exp(mu1))+log10(exp(mu2))-log10(lnf_delta),-1,lty="dashed",lwd=1.5)

#label the "tails"
text(-1.35,1.8,"right",adj=c(0,0),cex=0.65)
text(-1.35,1.8,"left",adj=c(1,1),cex=0.65)

#label the dashed lines
text(-0.55,0.9,expression(L[2]),cex=0.65,adj=c(0,0))
text(-1.4,-1,expression(L[1]),cex=0.65,adj=c(1,0))

#now add the points for the left-tail noise
points(log10(E1ta),log10(C1bs1ta),type="p",pch=20,cex=.1,col="green")

#panel label
rect(.88*dlims+lims[1],.865*dlims+lims[1],lims[2],lims[2],col="white")
text(lims[2],lims[2],"A",cex=cex_lab,adj=c(1.05,1.05))

#axes
mtext(expression(log[10](E[1])),1,1.4)
mtext(expression(log[10](C["1\\1"])),2,1)
mtext("left-tail assoc.",3,0.2)

#***
#Now make the corresponding plot that shows storage effects
#***

#set up the plotting panel
i<-2
j<-1
par(fig=c((yaxwd+(j-1)*(panwd+yaxwd))/totwd,
          (yaxwd+(j-1)*(panwd+yaxwd)+panwd)/totwd,
          (xaxht+(2-i)*(panht+xaxht))/totht,
          (xaxht+(2-i)*(panht+xaxht)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)

#Get the x axis values for the line plot for this and the two analogous panels,
#as well as x axis limits, and same for y
xDEF_rg<-range(log10(lnf_E1l/lnf_C1bs1l),log10(lnf_E1s/lnf_C1bs1s),log10(lnf_E1r/lnf_C1bs1r))
xDEF<-seq(from=xDEF_rg[1],to=xDEF_rg[2],length.out=100)
yDEF<-log(1-lnf_delta+10^xDEF)
yDEF_rg<-range(yDEF)

#expand the axis limits to accommodate the two histograms on each axis, also get common breaks
bk_x<-seq(from=xDEF_rg[1],to=xDEF_rg[2],length.out=numhistbins+1)
bk_y<-seq(from=yDEF_rg[1],to=yDEF_rg[2],length.out=numhistbins+1)

yhistrg_l<-c(xDEF_rg[1]-.4*diff(xDEF_rg),xDEF_rg[1]-.25*diff(xDEF_rg))
yhistrg_r<-c(xDEF_rg[1]-.15*diff(xDEF_rg),xDEF_rg[1])
xhistrg_l<-c(yDEF_rg[1]-.4*diff(yDEF_rg),yDEF_rg[1]-.25*diff(yDEF_rg))
xhistrg_r<-c(yDEF_rg[1]-.15*diff(yDEF_rg),yDEF_rg[1])
xDEF_rg[1]<-yhistrg_l[1]
yDEF_rg[1]<-xhistrg_l[1]

#plot the line
plot(xDEF,yDEF,type="l",xlim=xDEF_rg,ylim=yDEF_rg)
mtext(expression(log[10](E[1]/C["1\\1"])),1,1.4)
mtext(expression(ln(1-delta+E[1]/C["1\\1"])),2,1)

#now assemble the data for the x axis histograms and make them
x<-log10(E1ta/C1bs1ta)
inds_l<-(-log10(E1ta)+log10(exp(mu1))+log10(exp(mu2))-log10(lnf_delta)>log10(C1bs1ta))
inds_r<-(-log10(E1ta)+log10(exp(mu1))+log10(exp(mu2))-log10(lnf_delta)<=log10(C1bs1ta))
x_l<-x[inds_l]
x_r<-x[inds_r]
cts_l<-hist(x_l,breaks=bk_x,plot=FALSE)$counts
cts_r<-hist(x_r,breaks=bk_x,plot=FALSE)$counts
mybarplot(m=cts_l,bk=bk_x,col="green",base=xhistrg_l[1],top=xhistrg_l[2],axis=1,reddot=mean(x_l))
mybarplot(m=cts_r,bk=bk_x,col="green",base=xhistrg_r[1],top=xhistrg_r[2],axis=1,reddot=mean(x_r))

#now assemble the data for the y axis histograms and make them
y<-log(1-lnf_delta+10^x)
y_l<-y[inds_l]
y_r<-y[inds_r]
cts_l<-hist(y_l,breaks=bk_y,plot=FALSE)$counts
cts_r<-hist(y_r,breaks=bk_y,plot=FALSE)$counts
mybarplot(m=cts_l,bk=bk_y,col="green",base=yhistrg_l[1],top=yhistrg_l[2],axis=2,reddot=mean(y_l))
mybarplot(m=cts_r,bk=bk_y,col="green",base=yhistrg_r[1],top=yhistrg_r[2],axis=2,reddot=mean(y_r))

#panel label
xlim<-c(xDEF_rg[1]-.04*diff(xDEF_rg),xDEF_rg[2]+.04*diff(xDEF_rg))
ylim<-c(yDEF_rg[1]-.04*diff(yDEF_rg),yDEF_rg[2]+.04*diff(yDEF_rg))
rect(.88*diff(xlim)+xlim[1],.865*diff(ylim)+ylim[1],
     xlim[2],ylim[2],col="white")
text(xlim[2],ylim[2],"D",cex=cex_lab,adj=c(1.05,1.05))

#growth
text(xlim[1],ylim[2],paste0("growth=",round(mean(y),4)),cex=0.75,adj=c(-0.05,1.2))

#label the "tails"
text(xlim[2],xhistrg_l[1],"left",adj=c(1.25,0),cex=0.65)
text(xlim[2],xhistrg_r[1],"right",adj=c(1.2,0),cex=0.65)
text(yhistrg_l[1],min(yDEF),"left",adj=c(0,0),srt=270,cex=0.65)
text(yhistrg_r[1],min(yDEF),"right",adj=c(0,0),srt=270,cex=0.65)

#***
#Now generate the distributions for symmetric noise, and plot on top of a 
#color image showing the function ln(1-delta+E1/C1bs1), for lognormal fecundities
#***

E1s<-lnf_E1s
C1bs1s<-lnf_C1bs1s

#set up the plotting panel
i<-1
j<-2
par(fig=c((yaxwd+(j-1)*(panwd+yaxwd))/totwd,
          (yaxwd+(j-1)*(panwd+yaxwd)+panwd)/totwd,
          (xaxht+(2-i)*(panht+xaxht))/totht,
          (xaxht+(2-i)*(panht+xaxht)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)

#first make the color image
lims<-c(-1.9,2.3)
dlims<-diff(lims)
log10x<-seq(from=-1.9,to=2.3,by=0.01)
log10y<-log10x
z<-outer(10^log10x,10^log10y,function(x,y){log(1-lnf_delta+x/y)})
breaks<-pretty(z,10)
image(log10x,log10y,z,xlab="log10E1",ylab="log10C1bs1",col = hcl.colors(length(breaks)-1, "YlOrRd", rev = TRUE),
      breaks=breaks,xlim=lims,ylim=lims) 
contour(log10x,log10y,z,nlevels=10,add=T,labcex=cex_lab_cont)

#now add a reference line
lines(c(min(log10x),max(log10x)),c(min(log10y),max(log10y))-log10(lnf_delta)+log10(exp(mu2))-log10(exp(mu1)),type="l",lty="dashed",lwd=1.5)

#a divider line between the "tails"
#lnf_E1l<-exp(sig*nl[,1]+mu1)
#lnf_C1bs1l<-exp(sig*nl[,2]+mu2)/lnf_delta
#line thru log10(exp(mu1)),log10(exp(mu2))-log10(lnf_delta), with slope -1
#(y-log10(exp(mu2))+log10(lnf_delta))=-(x-log0(exp(mu1)))
#y=-x+log10(exp(mu1))+log10(exp(mu2))-log10(lnf_delta)
abline(log10(exp(mu1))+log10(exp(mu2))-log10(lnf_delta),-1,lty="dashed",lwd=1.5)

#label the "tails"
text(-1.35,1.8,"right",adj=c(0,0),cex=0.65)
text(-1.35,1.8,"left",adj=c(1,1),cex=0.65)

#label the dashed lines
text(-0.55,0.9,expression(L[2]),cex=0.65,adj=c(0,0))
text(-1.4,-1,expression(L[1]),cex=0.65,adj=c(1,0))

#now add the points for the left-tail noise
points(log10(E1s),log10(C1bs1s),type="p",pch=20,cex=.1,col="black")

#panel label
rect(.88*dlims+lims[1],.865*dlims+lims[1],lims[2],lims[2],col="white")
text(lims[2],lims[2],"B",cex=cex_lab,adj=c(1.05,1.05))

#axes
mtext(expression(log[10](E[1]^{"||"})),1,1.4)
mtext(expression(log[10](C["1\\1"]^{"||"})),2,1)

#***
#Now make the corresponding plot that shows storage effects
#***

#set up the plotting panel
i<-2
j<-2
par(fig=c((yaxwd+(j-1)*(panwd+yaxwd))/totwd,
          (yaxwd+(j-1)*(panwd+yaxwd)+panwd)/totwd,
          (xaxht+(2-i)*(panht+xaxht))/totht,
          (xaxht+(2-i)*(panht+xaxht)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)

#plot the line
plot(xDEF,yDEF,type="l",xlim=xDEF_rg,ylim=yDEF_rg)
mtext(expression(log[10](E[1]/C["1\\1"])),1,1.4)
mtext(expression(ln(1-delta+E[1]/C["1\\1"])),2,1)

#now assemble the data for the x axis histograms and make them
x<-log10(E1s/C1bs1s)
inds_l<-(-log10(E1s)+log10(exp(mu1))+log10(exp(mu2))-log10(lnf_delta)>log10(C1bs1s))
inds_r<-(-log10(E1s)+log10(exp(mu1))+log10(exp(mu2))-log10(lnf_delta)<=log10(C1bs1s))
x_l<-x[inds_l]
x_r<-x[inds_r]
cts_l<-hist(x_l,breaks=bk_x,plot=FALSE)$counts
cts_r<-hist(x_r,breaks=bk_x,plot=FALSE)$counts
mybarplot(m=cts_l,bk=bk_x,col="black",base=xhistrg_l[1],top=xhistrg_l[2],axis=1,reddot=mean(x_l))
mybarplot(m=cts_r,bk=bk_x,col="black",base=xhistrg_r[1],top=xhistrg_r[2],axis=1,reddot=mean(x_r))

#now assemble the data for the y axis histograms and make them
y<-log(1-lnf_delta+10^x)
y_l<-y[inds_l]
y_r<-y[inds_r]
cts_l<-hist(y_l,breaks=bk_y,plot=FALSE)$counts
cts_r<-hist(y_r,breaks=bk_y,plot=FALSE)$counts
mybarplot(m=cts_l,bk=bk_y,col="black",base=yhistrg_l[1],top=yhistrg_l[2],axis=2,reddot=mean(y_l))
mybarplot(m=cts_r,bk=bk_y,col="black",base=yhistrg_r[1],top=yhistrg_r[2],axis=2,reddot=mean(y_r))

#panel label
xlim<-c(xDEF_rg[1]-.04*diff(xDEF_rg),xDEF_rg[2]+.04*diff(xDEF_rg))
ylim<-c(yDEF_rg[1]-.04*diff(yDEF_rg),yDEF_rg[2]+.04*diff(yDEF_rg))
rect(.88*diff(xlim)+xlim[1],.865*diff(ylim)+ylim[1],
     xlim[2],ylim[2],col="white")
text(xlim[2],ylim[2],"E",cex=cex_lab,adj=c(1.05,1.05))

#growth
text(xlim[1],ylim[2],paste0("growth=",round(mean(y),4)),cex=0.75,adj=c(-0.05,1.2))

#label the "tails"
text(xlim[2],xhistrg_l[1],"left",adj=c(1.25,0),cex=0.65)
text(xlim[2],xhistrg_r[1],"right",adj=c(1.2,0),cex=0.65)
text(yhistrg_l[1],min(yDEF),"left",adj=c(0,0),srt=270,cex=0.65)
text(yhistrg_r[1],min(yDEF),"right",adj=c(0,0),srt=270,cex=0.65)

#***
#Now generate the distributions for right tail noise, and plot on top of a 
#color image showing the function ln(1-delta+E1/C1bs1), for lognormal fecundities
#***

E1ta<-lnf_E1r
C1bs1ta<-lnf_C1bs1r

#set up the plotting panel
i<-1
j<-3
par(fig=c((yaxwd+(j-1)*(panwd+yaxwd))/totwd,
          (yaxwd+(j-1)*(panwd+yaxwd)+panwd)/totwd,
          (xaxht+(2-i)*(panht+xaxht))/totht,
          (xaxht+(2-i)*(panht+xaxht)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)

#first make the color image
lims<-c(-1.9,2.3)
dlims<-diff(lims)
log10x<-seq(from=-1.9,to=2.3,by=0.01)
log10y<-log10x
z<-outer(10^log10x,10^log10y,function(x,y){log(1-lnf_delta+x/y)})
breaks<-pretty(z,10)
image(log10x,log10y,z,xlab="log10E1",ylab="log10C1bs1",col = hcl.colors(length(breaks)-1, "YlOrRd", rev = TRUE),
      breaks=breaks,xlim=lims,ylim=lims) 
contour(log10x,log10y,z,nlevels=10,add=T,labcex=cex_lab_cont)

#now add a reference line
lines(c(min(log10x),max(log10x)),c(min(log10y),max(log10y))-log10(lnf_delta)+log10(exp(mu2))-log10(exp(mu1)),type="l",lty="dashed",lwd=1.5)

#a divider line between the "tails"
#lnf_E1l<-exp(sig*nl[,1]+mu1)
#lnf_C1bs1l<-exp(sig*nl[,2]+mu2)/lnf_delta
#line thru log10(exp(mu1)),log10(exp(mu2))-log10(lnf_delta), with slope -1
#(y-log10(exp(mu2))+log10(lnf_delta))=-(x-log0(exp(mu1)))
#y=-x+log10(exp(mu1))+log10(exp(mu2))-log10(lnf_delta)
abline(log10(exp(mu1))+log10(exp(mu2))-log10(lnf_delta),-1,lty="dashed",lwd=1.5)

#label the "tails"
text(-1.35,1.8,"right",adj=c(0,0),cex=0.65)
text(-1.35,1.8,"left",adj=c(1,1),cex=0.65)

#label the dashed lines
text(-0.55,0.9,expression(L[2]),cex=0.65,adj=c(0,0))
text(-1.4,-1,expression(L[1]),cex=0.65,adj=c(1,0))

#now add the points for the left-tail noise
points(log10(E1ta),log10(C1bs1ta),type="p",pch=20,cex=.1,col="green")

#panel label
rect(.88*dlims+lims[1],.865*dlims+lims[1],lims[2],lims[2],col="white")
text(lims[2],lims[2],"C",cex=cex_lab,adj=c(1.05,1.05))

#axes
mtext(expression(log[10](E[1])),1,1.4)
mtext(expression(log[10](C["1\\1"])),2,1)
mtext("right-tail assoc.",3,0.2)

#***
#Now make the corresponding plot that shows storage effects
#***

#set up the plotting panel
i<-2
j<-3
par(fig=c((yaxwd+(j-1)*(panwd+yaxwd))/totwd,
          (yaxwd+(j-1)*(panwd+yaxwd)+panwd)/totwd,
          (xaxht+(2-i)*(panht+xaxht))/totht,
          (xaxht+(2-i)*(panht+xaxht)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)

#plot the line
plot(xDEF,yDEF,type="l",xlim=xDEF_rg,ylim=yDEF_rg)
mtext(expression(log[10](E[1]/C["1\\1"])),1,1.4)
mtext(expression(ln(1-delta+E[1]/C["1\\1"])),2,1)

#now assemble the data for the x axis histograms and make them
x<-log10(E1ta/C1bs1ta)
inds_l<-(-log10(E1ta)+log10(exp(mu1))+log10(exp(mu2))-log10(lnf_delta)>log10(C1bs1ta))
inds_r<-(-log10(E1ta)+log10(exp(mu1))+log10(exp(mu2))-log10(lnf_delta)<=log10(C1bs1ta))
x_l<-x[inds_l]
x_r<-x[inds_r]
cts_l<-hist(x_l,breaks=bk_x,plot=FALSE)$counts
cts_r<-hist(x_r,breaks=bk_x,plot=FALSE)$counts
mybarplot(m=cts_l,bk=bk_x,col="green",base=xhistrg_l[1],top=xhistrg_l[2],axis=1,reddot=mean(x_l))
mybarplot(m=cts_r,bk=bk_x,col="green",base=xhistrg_r[1],top=xhistrg_r[2],axis=1,reddot=mean(x_r))

#now assemble the data for the y axis histograms and make them
y<-log(1-lnf_delta+10^x)
y_l<-y[inds_l]
y_r<-y[inds_r]
cts_l<-hist(y_l,breaks=bk_y,plot=FALSE)$counts
cts_r<-hist(y_r,breaks=bk_y,plot=FALSE)$counts
mybarplot(m=cts_l,bk=bk_y,col="green",base=yhistrg_l[1],top=yhistrg_l[2],axis=2,reddot=mean(y_l))
mybarplot(m=cts_r,bk=bk_y,col="green",base=yhistrg_r[1],top=yhistrg_r[2],axis=2,reddot=mean(y_r))

#panel label
xlim<-c(xDEF_rg[1]-.04*diff(xDEF_rg),xDEF_rg[2]+.04*diff(xDEF_rg))
ylim<-c(yDEF_rg[1]-.04*diff(yDEF_rg),yDEF_rg[2]+.04*diff(yDEF_rg))
rect(.88*diff(xlim)+xlim[1],.865*diff(ylim)+ylim[1],
     xlim[2],ylim[2],col="white")
text(xlim[2],ylim[2],"F",cex=cex_lab,adj=c(1.05,1.05))

#growth
text(xlim[1],ylim[2],paste0("growth=",round(mean(y),4)),cex=0.75,adj=c(-0.05,1.2))

#label the "tails"
text(xlim[2],xhistrg_l[1],"left",adj=c(1.25,0),cex=0.65)
text(xlim[2],xhistrg_r[1],"right",adj=c(1.2,0),cex=0.65)
text(yhistrg_l[1],min(yDEF),"left",adj=c(0,0),srt=270,cex=0.65)
text(yhistrg_r[1],min(yDEF),"right",adj=c(0,0),srt=270,cex=0.65)

dev.off()
