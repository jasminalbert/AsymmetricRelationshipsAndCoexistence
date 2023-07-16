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

#graphical params
numhistbins<-12

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
yaxwd<-0.6
gap<-0.3
totwd<-6.5
panwd<-(totwd-4*gap-yaxwd)/4
panht<-panwd
totht<-3*xaxht+3*panht+gap
cex_lab<-1.5
cex_lab_cont<-0.5

jpeg(file=paste0(resloc,"PostHocExplanatoryFig.jpg"),width=totwd,height=totht,units="in",quality=90,res=400)
#pdf(file=paste0(resloc,"PostHocExplanatoryFig.pdf"),width=totwd,height=totht)

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
          (xaxht+(3-i)*(panht+xaxht))/totht,
          (xaxht+(3-i)*(panht+xaxht)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25)

#first make the color image
lims<-c(0,120)
dlims<-diff(lims)
x<-seq(from=0.01,to=120.01,by=0.1)
y<-x
z<-outer(x,y,function(x,y){log(1-bf_delta+x/y)})
breaks<-pretty(z,10)
image(x,y,z,xlab="E1",ylab="C1bs1",col = hcl.colors(length(breaks)-1, "YlOrRd", rev = TRUE),
      breaks=breaks,xlim=lims,ylim=lims) 
contour(x,y,z,nlevels=10,add=T,labcex=cex_lab_cont)
mtext(expression(E[1]),1,1.2)
mtext(expression(C[1/1]),2,1.2)

#now add the points for the left-tail noise
points(E1ta,C1bs1ta,type="p",pch=20,cex=.1,col="green")

#now add the points for the symmetric noise, using a different color
points(E1s,C1bs1s,type="p",pch=20,cex=.1)

#now add the panel label
rect(.84*dlims+lims[1],.82*dlims+lims[1],lims[2],lims[2],col="white")
text(lims[2],lims[2],"A",cex=cex_lab,adj=c(1.05,1.05))

#***then log scale

#set up the plotting panel
i<-2
j<-1
par(fig=c((yaxwd+(j-1)*(panwd+gap))/totwd,
          (yaxwd+(j-1)*(panwd+gap)+panwd)/totwd,
          (xaxht+(3-i)*(panht+xaxht))/totht,
          (xaxht+(3-i)*(panht+xaxht)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)

#first make the color image
lims<-c(-1.9,2.3)
dlims<-diff(lims)
log10x<-seq(from=lims[1],to=lims[2],by=0.01)
log10y<-log10x
z<-outer(10^log10x,10^log10y,function(x,y){log(1-bf_delta+x/y)})
breaks<-pretty(z,10)
image(log10x,log10y,z,xlab="log10E1",ylab="log10C1bs1",col = hcl.colors(length(breaks)-1, "YlOrRd", rev = TRUE),
      breaks=breaks,xlim=lims,ylim=lims) 
contour(log10x,log10y,z,nlevels=10,add=T,labcex=cex_lab_cont)
mtext(expression(log[10](E[1])),1,1.2)
mtext(expression(log[10](C[1/1])),2,1.2)

#now add a reference line
lines(c(min(log10x),max(log10x)),c(min(log10y),max(log10y))-log10(lnf_delta)+log10(exp(mu2))-log10(exp(mu1)),type="l",lty="dashed",lwd=1.5)

#a divider line between the "tails"
#lnf_E1l<-exp(sig*nl[,1]+mu1)
#lnf_C1bs1l<-exp(sig*nl[,2]+mu2)/lnf_delta
#line thru log10(exp(mu1)),log10(exp(mu2))-log10(lnf_delta), with slope -1
#(y-log10(exp(mu2))+log10(lnf_delta))=-(x-log0(exp(mu1)))
#y=-x+log10(exp(mu1))+log10(exp(mu2))-log10(lnf_delta)
abline(log10(exp(mu1))+log10(exp(mu2))-log10(lnf_delta),-1,lty="dashed",lwd=1.5)

#now add the points for the symmetric noise, using a different color
points(log10(E1s),log10(C1bs1s),type="p",pch=20,cex=.1)

#now add the points for the left-tail noise
points(log10(E1ta),log10(C1bs1ta),type="p",pch=20,cex=.1,col="green")

#now the panel label
rect(.84*dlims+lims[1],.82*dlims+lims[1],lims[2],lims[2],col="white")
text(lims[2],lims[2],"B",cex=cex_lab,adj=c(1.05,1.05))

#***now deviations from the central reference line

#x axis quantities for points on the new panel - signed distance from the central reference line
xs<-((log10(E1s)-log10(C1bs1s)-log10(lnf_delta)+log10(exp(mu2))-log10(exp(mu1)))/sqrt(2))
xta<-((log10(E1ta)-log10(C1bs1ta)-log10(lnf_delta)+log10(exp(mu2))-log10(exp(mu1)))/sqrt(2))

#y axis quantities for the same
ys<-log(1-bf_delta+E1s/C1bs1s)
yta<-log(1-bf_delta+E1ta/C1bs1ta)

#now separate these into those in the left v right tails
#y=-x+log10(exp(mu1))+log10(exp(mu2))-log10(lnf_delta)
inds_l<-(-log10(E1s)+log10(exp(mu1))+log10(exp(mu2))-log10(lnf_delta)>log10(C1bs1s))
inds_r<-(-log10(E1s)+log10(exp(mu1))+log10(exp(mu2))-log10(lnf_delta)<=log10(C1bs1s))
xs_l<-xs[inds_l]
xs_r<-xs[inds_r]
ys_l<-ys[inds_l]
ys_r<-ys[inds_r]
inds_l<-(-log10(E1ta)+log10(exp(mu1))+log10(exp(mu2))-log10(lnf_delta)>log10(C1bs1ta))
inds_r<-(-log10(E1ta)+log10(exp(mu1))+log10(exp(mu2))-log10(lnf_delta)<=log10(C1bs1ta))
xta_l<-xta[inds_l]
xta_r<-xta[inds_r]
yta_l<-yta[inds_l]
yta_r<-yta[inds_r]

#coordinates for the plotted line which represents the value of 
#log(1-bf_delta+x/y) as a function of distance from the central 
#reference line
xol<-log10x
yol<-rev(log10y)
yl<-log(1-bf_delta+(10^xol)/(10^yol))
xl<-((xol-yol-log10(lnf_delta)+log10(exp(mu2))-log10(exp(mu1))))/sqrt(2)  
yl<-yl[xl>=min(xs,xta) & xl<=max(xs,xta)]
xl<-xl[xl>=min(xs,xta) & xl<=max(xs,xta)]

#set up the plotting panel
i<-3
j<-1
par(fig=c((yaxwd+(j-1)*(panwd+gap))/totwd,
          (yaxwd+(j-1)*(panwd+gap)+panwd)/totwd,
          (xaxht+(3-i)*(panht+xaxht))/totht,
          (xaxht+(3-i)*(panht+xaxht)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)

#plot the line
xlim<-range(xl)
ylim<-range(yl)
xhistrg_l<-c(xlim[1]-.4*diff(xlim),xlim[1]-.25*diff(xlim))
xhistrg_r<-c(xlim[1]-.15*diff(xlim),xlim[1])
yhistrg_l<-c(ylim[1]-.4*diff(ylim),ylim[1]-.25*diff(ylim))
yhistrg_r<-c(ylim[1]-.15*diff(ylim),ylim[1])
xlim[1]<-xhistrg_l[1]
ylim[1]<-yhistrg_l[1]
plot(xl,yl,type="l",xlim=xlim,ylim=ylim)
mtext("Perp. distance",1,1.2)
mtext(expression(log(1-delta+E[1]/C[1/1])),2,1.2)

#plot the mini histograms on the x axis
bks<-seq(from=min(xs,xta),to=max(xs,xta),length.out=numhistbins+1)
cts_l<-hist(xs_l,breaks=bks,plot=F)$counts
cts_r<-hist(xs_r,breaks=bks,plot=F)$counts
ctta_l<-hist(xta_l,breaks=bks,plot=F)$counts
ctta_r<-hist(xta_r,breaks=bks,plot=F)$counts

#A function that makes it relatively easier for me to make the barplots
#I want to make. For plotting two interdigitated bar plots for two variables.
#Adds to a plot assumed to already be open.
#
#Args
#m      A 2 by N matrix of counts
#bks    Breaks. Length is dim(m)[2]+1
#cols   Length-2 vector of colors
#base   bottoms of rectangles
#top    top of highest rectangle
#axis   Use 1 for plotting the bars vertically, using the x axis for breaks.
#         Use 2 for plotting the bars horizontally, using the y axis for breaks.
#
#Output - none, but adds to a plot
#
#Notes. 
#If the arguments bks, base, top do not comport with the axes already
#open then the barplot will appear (invisibly) off the plotting axes.
#
mybarplot<-function(m,bk,cols,base,top,axis)
{
  Mx<-max(m)
  
  if (axis==1)
  {
    for (rc in 1:(dim(m)[2]))
    {
      rect(bk[rc],base,mean(bk[rc:(rc+1)]),base+(top-base)*m[1,rc]/Mx,col=cols[1],border=NA)
      rect(mean(bk[rc:(rc+1)]),base,bk[rc+1],base+(top-base)*m[2,rc]/Mx,col=cols[2],border=NA)
    }
  }
  
  if (axis==2)
  {
    for (rc in 1:(dim(m)[2]))
    {
      rect(base,bk[rc],base+(top-base)*m[1,rc]/Mx,mean(bk[rc:(rc+1)]),col=cols[1],border=NA)
      rect(base,mean(bk[rc:(rc+1)]),base+(top-base)*m[2,rc]/Mx,bk[rc+1],col=cols[2],border=NA)
    }
  }
}

mybarplot(m=matrix(c(cts_l,ctta_l),2,length(cts_l),byrow=TRUE),bk=bks,cols=c("black","green"),
          base=yhistrg_l[1],top=yhistrg_l[2],axis=1)
mybarplot(m=matrix(c(cts_r,ctta_r),2,length(cts_r),byrow=TRUE),bk=bks,cols=c("black","green"),
          base=yhistrg_r[1],top=yhistrg_r[2],axis=1)

#plot the mini histogram on the y axis
bks<-seq(from=min(ys,yta),to=max(ys,yta),length.out=numhistbins+1)
cts_l<-hist(ys_l,breaks=bks,plot=F)$counts
cts_r<-hist(ys_r,breaks=bks,plot=F)$counts
ctta_l<-hist(yta_l,breaks=bks,plot=F)$counts
ctta_r<-hist(yta_r,breaks=bks,plot=F)$counts

mybarplot(m=matrix(c(cts_l,ctta_l),2,length(cts_l),byrow=TRUE),bk=bks,cols=c("black","green"),
          base=xhistrg_l[1],top=xhistrg_l[2],axis=2)
mybarplot(m=matrix(c(cts_r,ctta_r),2,length(cts_r),byrow=TRUE),bk=bks,cols=c("black","green"),
          base=xhistrg_r[1],top=xhistrg_r[2],axis=2)

#panel label
xlim<-c(xlim[1]-.04*diff(xlim),xlim[2]+.04*diff(xlim))
ylim<-c(ylim[1]-.04*diff(ylim),ylim[2]+.04*diff(ylim))
rect(.84*diff(xlim)+xlim[1],.82*diff(ylim)+ylim[1],
     xlim[2],ylim[2],col="white")
text(xlim[2],ylim[2],"C",cex=cex_lab,adj=c(1.05,1.05))

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
          (xaxht+(3-i)*(panht+xaxht))/totht,
          (xaxht+(3-i)*(panht+xaxht)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)

#first make the color image
lims<-c(0,120)
dlims<-diff(lims)
x<-seq(from=0.01,to=120.01,by=0.1)
y<-x
z<-outer(x,y,function(x,y){log(1-bf_delta+x/y)})
breaks<-pretty(z,10)
image(x,y,z,xlab="E1",ylab="C1bs1",col = hcl.colors(length(breaks)-1, "YlOrRd", rev = TRUE),
      breaks=breaks,xlim=lims,ylim=lims) 
contour(x,y,z,nlevels=10,add=T,labcex=cex_lab_cont)
mtext(expression(E[1]),1,1.2)

#now add the points for the symmetric noise, using a different color
points(E1s,C1bs1s,type="p",pch=20,cex=.1)

#now add the points for the left-tail noise
points(E1ta,C1bs1ta,type="p",pch=20,cex=.1,col="green")

#panel label
rect(.84*dlims+lims[1],.82*dlims+lims[1],lims[2],lims[2],col="white")
text(lims[2],lims[2],"D",cex=cex_lab,adj=c(1.05,1.05))

#***then log scale

#set up the plotting panel
i<-2
j<-2
par(fig=c((yaxwd+(j-1)*(panwd+gap))/totwd,
          (yaxwd+(j-1)*(panwd+gap)+panwd)/totwd,
          (xaxht+(3-i)*(panht+xaxht))/totht,
          (xaxht+(3-i)*(panht+xaxht)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)

#first make the color image
lims<-c(-1.9,2.3)
dlims<-diff(lims)
log10x<-seq(from=lims[1],to=lims[2],by=0.01)
log10y<-log10x
z<-outer(10^log10x,10^log10y,function(x,y){log(1-bf_delta+x/y)})
breaks<-pretty(z,10)
image(log10x,log10y,z,xlab="log10E1",ylab="log10C1bs1",col = hcl.colors(length(breaks)-1, "YlOrRd", rev = TRUE),
      breaks=breaks,xlim=lims,ylim=lims) 
contour(log10x,log10y,z,nlevels=10,add=T,labcex=cex_lab_cont)
mtext(expression(log[10](E[1])),1,1.2)

#now add a reference line
lines(c(min(log10x),max(log10x)),c(min(log10y),max(log10y))-log10(lnf_delta)+log10(exp(mu2))-log10(exp(mu1)),type="l",lty="dashed",lwd=1.5)

#a divider line between the "tails"
abline(log10(exp(mu1))+log10(exp(mu2))-log10(lnf_delta),-1,lty="dashed",lwd=1.5)

#now add the points for the symmetric noise, using a different color
points(log10(E1s),log10(C1bs1s),type="p",pch=20,cex=.1)

#now add the points for the left-tail noise
points(log10(E1ta),log10(C1bs1ta),type="p",pch=20,cex=.1,col="green")

#panel label
rect(.84*dlims+lims[1],.82*dlims+lims[1],lims[2],lims[2],col="white")
text(lims[2],lims[2],"E",cex=cex_lab,adj=c(1.05,1.05))

#***now deviations from the central reference line

#x axis quantities for points on the new panel - signed distance from the central reference line
xs<-((log10(E1s)-log10(C1bs1s)-log10(lnf_delta)+log10(exp(mu2))-log10(exp(mu1)))/sqrt(2))
xta<-((log10(E1ta)-log10(C1bs1ta)-log10(lnf_delta)+log10(exp(mu2))-log10(exp(mu1)))/sqrt(2))

#y axis quantities for the same
ys<-log(1-bf_delta+E1s/C1bs1s)
yta<-log(1-bf_delta+E1ta/C1bs1ta)

#now separate these into those in the left v right tails
#y=-x+log10(exp(mu1))+log10(exp(mu2))-log10(lnf_delta)
inds_l<-(-log10(E1s)+log10(exp(mu1))+log10(exp(mu2))-log10(lnf_delta)>log10(C1bs1s))
inds_r<-(-log10(E1s)+log10(exp(mu1))+log10(exp(mu2))-log10(lnf_delta)<=log10(C1bs1s))
xs_l<-xs[inds_l]
xs_r<-xs[inds_r]
ys_l<-ys[inds_l]
ys_r<-ys[inds_r]
inds_l<-(-log10(E1ta)+log10(exp(mu1))+log10(exp(mu2))-log10(lnf_delta)>log10(C1bs1ta))
inds_r<-(-log10(E1ta)+log10(exp(mu1))+log10(exp(mu2))-log10(lnf_delta)<=log10(C1bs1ta))
xta_l<-xta[inds_l]
xta_r<-xta[inds_r]
yta_l<-yta[inds_l]
yta_r<-yta[inds_r]

#coordinates for the plotted line which represents the value of 
#log(1-bf_delta+x/y) as a function of distance from the central 
#reference line
xol<-log10x
yol<-rev(log10y)
yl<-log(1-bf_delta+(10^xol)/(10^yol))
xl<-((xol-yol-log10(lnf_delta)+log10(exp(mu2))-log10(exp(mu1))))/sqrt(2)  
yl<-yl[xl>=min(xs,xta) & xl<=max(xs,xta)]
xl<-xl[xl>=min(xs,xta) & xl<=max(xs,xta)]

#set up the plotting panel
i<-3
j<-2
par(fig=c((yaxwd+(j-1)*(panwd+gap))/totwd,
          (yaxwd+(j-1)*(panwd+gap)+panwd)/totwd,
          (xaxht+(3-i)*(panht+xaxht))/totht,
          (xaxht+(3-i)*(panht+xaxht)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)

#plot the line
xlim<-range(xl)
ylim<-range(yl)
xhistrg_l<-c(xlim[1]-.4*diff(xlim),xlim[1]-.25*diff(xlim))
xhistrg_r<-c(xlim[1]-.15*diff(xlim),xlim[1])
yhistrg_l<-c(ylim[1]-.4*diff(ylim),ylim[1]-.25*diff(ylim))
yhistrg_r<-c(ylim[1]-.15*diff(ylim),ylim[1])
xlim[1]<-xhistrg_l[1]
ylim[1]<-yhistrg_l[1]
plot(xl,yl,type="l",xlim=xlim,ylim=ylim)
mtext("Perp. distance",1,1.2)

#plot the mini histograms on the x axis
bks<-seq(from=min(xs,xta),to=max(xs,xta),length.out=numhistbins+1)
cts_l<-hist(xs_l,breaks=bks,plot=F)$counts
cts_r<-hist(xs_r,breaks=bks,plot=F)$counts
ctta_l<-hist(xta_l,breaks=bks,plot=F)$counts
ctta_r<-hist(xta_r,breaks=bks,plot=F)$counts

mybarplot(m=matrix(c(cts_l,ctta_l),2,length(cts_l),byrow=TRUE),bk=bks,cols=c("black","green"),
          base=yhistrg_l[1],top=yhistrg_l[2],axis=1)
mybarplot(m=matrix(c(cts_r,ctta_r),2,length(cts_r),byrow=TRUE),bk=bks,cols=c("black","green"),
          base=yhistrg_r[1],top=yhistrg_r[2],axis=1)

#plot the mini histogram on the y axis
bks<-seq(from=min(ys,yta),to=max(ys,yta),length.out=numhistbins+1)
cts_l<-hist(ys_l,breaks=bks,plot=F)$counts
cts_r<-hist(ys_r,breaks=bks,plot=F)$counts
ctta_l<-hist(yta_l,breaks=bks,plot=F)$counts
ctta_r<-hist(yta_r,breaks=bks,plot=F)$counts

mybarplot(m=matrix(c(cts_l,ctta_l),2,length(cts_l),byrow=TRUE),bk=bks,cols=c("black","green"),
          base=xhistrg_l[1],top=xhistrg_l[2],axis=2)
mybarplot(m=matrix(c(cts_r,ctta_r),2,length(cts_r),byrow=TRUE),bk=bks,cols=c("black","green"),
          base=xhistrg_r[1],top=xhistrg_r[2],axis=2)

#panel label
xlim<-c(xlim[1]-.04*diff(xlim),xlim[2]+.04*diff(xlim))
ylim<-c(ylim[1]-.04*diff(ylim),ylim[2]+.04*diff(ylim))
rect(.84*diff(xlim)+xlim[1],.82*diff(ylim)+ylim[1],
     xlim[2],ylim[2],col="white")
text(xlim[2],ylim[2],"F",cex=cex_lab,adj=c(1.05,1.05))

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
          (xaxht+(3-i)*(panht+xaxht))/totht,
          (xaxht+(3-i)*(panht+xaxht)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)

#first make the color image
lims<-c(0,2.1)
dlims<-diff(lims)
x<-seq(from=0.01,to=2.11,by=0.01)
y<-x
z<-outer(x,y,function(x,y){log(1-bf_delta+x/y)})
breaks<-pretty(z,10)
image(x,y,z,xlab="E1",ylab="C1bs1",col = hcl.colors(length(breaks)-1, "YlOrRd", rev = TRUE),
      breaks=breaks,xlim=lims,ylim=lims) 
contour(x,y,z,nlevels=10,add=T,labcex=cex_lab_cont)
mtext(expression(E[1]),1,1.2)

#now add the points for the left-tail noise
points(E1ta,C1bs1ta,type="p",pch=20,cex=.1,col="green")

#now add the points for the symmetric noise, using a different color
points(E1s,C1bs1s,type="p",pch=20,cex=.1)

#panel label
rect(.83*dlims+lims[1],.82*dlims+lims[1],lims[2],lims[2],col="white")
text(lims[2],lims[2],"G",cex=cex_lab,adj=c(1.01,1.07))

#***then log scale

#set up the plotting panel
i<-2
j<-3
par(fig=c((yaxwd+(j-1)*(panwd+gap))/totwd,
          (yaxwd+(j-1)*(panwd+gap)+panwd)/totwd,
          (xaxht+(3-i)*(panht+xaxht))/totht,
          (xaxht+(3-i)*(panht+xaxht)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)

#first make the color image
lims<-c(-5,1.6)
dlims<-diff(lims)
log10x<-seq(from=-5,to=1.6,by=0.01)
log10y<-log10x
z<-outer(10^log10x,10^log10y,function(x,y){log(1-bf_delta+x/y)})
breaks<-pretty(z,10)
image(log10x,log10y,z,xlab="log10E1",ylab="log10C1bs1",col = hcl.colors(length(breaks)-1, "YlOrRd", rev = TRUE),
      breaks=breaks,xlim=lims,ylim=lims) 
contour(log10x,log10y,z,nlevels=10,add=T,labcex=cex_lab_cont)
mtext(expression(log[10](E[1])),1,1.2)

#now add a reference line
lines(c(min(log10x),max(log10x)),c(min(log10y),max(log10y))-log10(bf_delta)+log10(eta2/eta1),type="l",lty="dashed",lwd=1.5)

#now add a divider between tails
#bf_E1l<-eta1*qbeta(pnorm(nl[,1]),0.5,0.5)
#bf_C1bs1l<-eta2*qbeta(pnorm(nl[,2]),0.5,0.5)/bf_delta
#line thru log10(eta1*qbeta(pnorm(0),0.5,0.5)), log10(eta2*qbeta(pnorm(0),0.5,0.5))-log10(bf_delta), 
#of slope -1
#(y-log10(eta2*qbeta(pnorm(0),0.5,0.5))+log10(bf_delta))=-1(x-log10(eta1*qbeta(pnorm(0),0.5,0.5)))
#y=-x+log10(eta1*qbeta(pnorm(0),0.5,0.5))+log10(eta2*qbeta(pnorm(0),0.5,0.5))-log10(bf_delta)
abline(log10(eta1*qbeta(pnorm(0),0.5,0.5))+log10(eta2*qbeta(pnorm(0),0.5,0.5))-log10(bf_delta),-1,
      lty="dashed",lwd=1.5)

#now add the points for the symmetric noise, using a different color
points(log10(E1s),log10(C1bs1s),type="p",pch=20,cex=.1)

#now add the points for the left-tail noise
points(log10(E1ta),log10(C1bs1ta),type="p",pch=20,cex=.1,col="green")

#panel label
rect(.84*dlims+lims[1],.82*dlims+lims[1],lims[2],lims[2],col="white")
text(lims[2],lims[2],"H",cex=cex_lab,adj=c(1.05,1.05))

#***now deviations from the central reference line

#x axis quantities for points on the new panel - signed distance from the central reference line
xs<-((log10(E1s)-log10(C1bs1s)-log10(bf_delta)+log10(eta2)-log10(eta1))/sqrt(2))
xta<-((log10(E1ta)-log10(C1bs1ta)-log10(bf_delta)+log10(eta2)-log10(eta1))/sqrt(2))

#y axis quantities for the same
ys<-log(1-bf_delta+E1s/C1bs1s)
yta<-log(1-bf_delta+E1ta/C1bs1ta)

#now separate these into those in the left v right tails
#y=-x+log10(eta1*qbeta(pnorm(0),0.5,0.5))+log10(eta2*qbeta(pnorm(0),0.5,0.5))-log10(bf_delta)
inds_l<-(-log10(E1s)+log10(eta1*0.5)+log10(eta2*0.5)-log10(bf_delta)>log10(C1bs1s))
inds_r<-(-log10(E1s)+log10(eta1*0.5)+log10(eta2*0.5)-log10(bf_delta)<=log10(C1bs1s))
xs_l<-xs[inds_l]
xs_r<-xs[inds_r]
ys_l<-ys[inds_l]
ys_r<-ys[inds_r]
inds_l<-(-log10(E1ta)+log10(eta1*0.5)+log10(eta2*0.5)-log10(bf_delta)>log10(C1bs1ta))
inds_r<-(-log10(E1ta)+log10(eta1*0.5)+log10(eta2*0.5)-log10(bf_delta)<=log10(C1bs1ta))
xta_l<-xta[inds_l]
xta_r<-xta[inds_r]
yta_l<-yta[inds_l]
yta_r<-yta[inds_r]

#coordinates for the plotted line which represents the value of 
#log(1-bf_delta+x/y) as a function of distance from the central 
#reference line
xol<-log10x
yol<-rev(log10y)
yl<-log(1-bf_delta+(10^xol)/(10^yol))
xl<-((xol-yol-log10(bf_delta)+log10(eta2)-log10(eta1)))/sqrt(2)  
yl<-yl[xl>=min(xs,xta) & xl<=max(xs,xta)]
xl<-xl[xl>=min(xs,xta) & xl<=max(xs,xta)]

#set up the plotting panel
i<-3
j<-3
par(fig=c((yaxwd+(j-1)*(panwd+gap))/totwd,
          (yaxwd+(j-1)*(panwd+gap)+panwd)/totwd,
          (xaxht+(3-i)*(panht+xaxht))/totht,
          (xaxht+(3-i)*(panht+xaxht)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)

#plot the line
xlim<-range(xl)
ylim<-range(yl)
xhistrg_l<-c(xlim[1]-.4*diff(xlim),xlim[1]-.25*diff(xlim))
xhistrg_r<-c(xlim[1]-.15*diff(xlim),xlim[1])
yhistrg_l<-c(ylim[1]-.4*diff(ylim),ylim[1]-.25*diff(ylim))
yhistrg_r<-c(ylim[1]-.15*diff(ylim),ylim[1])
xlim[1]<-xhistrg_l[1]
ylim[1]<-yhistrg_l[1]
plot(xl,yl,type="l",xlim=xlim,ylim=ylim)
mtext("Perp. distance",1,1.2)

#plot the mini histograms on the x axis
bks<-seq(from=min(xs,xta),to=max(xs,xta),length.out=numhistbins+1)
cts_l<-hist(xs_l,breaks=bks,plot=F)$counts
cts_r<-hist(xs_r,breaks=bks,plot=F)$counts
ctta_l<-hist(xta_l,breaks=bks,plot=F)$counts
ctta_r<-hist(xta_r,breaks=bks,plot=F)$counts

mybarplot(m=matrix(c(cts_l,ctta_l),2,length(cts_l),byrow=TRUE),bk=bks,cols=c("black","green"),
          base=yhistrg_l[1],top=yhistrg_l[2],axis=1)
mybarplot(m=matrix(c(cts_r,ctta_r),2,length(cts_r),byrow=TRUE),bk=bks,cols=c("black","green"),
          base=yhistrg_r[1],top=yhistrg_r[2],axis=1)

#plot the mini histogram on the y axis
bks<-seq(from=min(ys,yta),to=max(ys,yta),length.out=numhistbins+1)
cts_l<-hist(ys_l,breaks=bks,plot=F)$counts
cts_r<-hist(ys_r,breaks=bks,plot=F)$counts
ctta_l<-hist(yta_l,breaks=bks,plot=F)$counts
ctta_r<-hist(yta_r,breaks=bks,plot=F)$counts

mybarplot(m=matrix(c(cts_l,ctta_l),2,length(cts_l),byrow=TRUE),bk=bks,cols=c("black","green"),
          base=xhistrg_l[1],top=xhistrg_l[2],axis=2)
mybarplot(m=matrix(c(cts_r,ctta_r),2,length(cts_r),byrow=TRUE),bk=bks,cols=c("black","green"),
          base=xhistrg_r[1],top=xhistrg_r[2],axis=2)

#panel label
xlim<-c(xlim[1]-.04*diff(xlim),xlim[2]+.04*diff(xlim))
ylim<-c(ylim[1]-.04*diff(ylim),ylim[2]+.04*diff(ylim))
rect(.84*diff(xlim)+xlim[1],.82*diff(ylim)+ylim[1],
     xlim[2],ylim[2],col="white")
text(xlim[2],ylim[2],"I",cex=cex_lab,adj=c(1.7,1.05))

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
          (xaxht+(3-i)*(panht+xaxht))/totht,
          (xaxht+(3-i)*(panht+xaxht)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)

#first make the color image
lims<-c(0,2.1)
dlims<-diff(lims)
x<-seq(from=0.01,to=2.11,by=0.01)
y<-x
z<-outer(x,y,function(x,y){log(1-bf_delta+x/y)})
breaks<-pretty(z,10)
image(x,y,z,xlab="E1",ylab="C1bs1",col = hcl.colors(length(breaks)-1, "YlOrRd", rev = TRUE),
      breaks=breaks,xlim=lims,ylim=lims) 
contour(x,y,z,nlevels=10,add=T,labcex=cex_lab_cont)
mtext(expression(E[1]),1,1.2)

#now add the points for the left-tail noise
points(E1ta,C1bs1ta,type="p",pch=20,cex=.1,col="green")

#now add the points for the symmetric noise, using a different color
points(E1s,C1bs1s,type="p",pch=20,cex=.1)

#panel label
rect(.84*dlims+lims[1],.82*dlims+lims[1],lims[2],lims[2],col="white")
text(lims[2],lims[2],"J",cex=cex_lab,adj=c(1.05,1.05))

#***then log scale

#set up the plotting panel
i<-2
j<-4
par(fig=c((yaxwd+(j-1)*(panwd+gap))/totwd,
          (yaxwd+(j-1)*(panwd+gap)+panwd)/totwd,
          (xaxht+(3-i)*(panht+xaxht))/totht,
          (xaxht+(3-i)*(panht+xaxht)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)

#first make the color image
lims<-c(-5,1.6)
dlims<-diff(lims)
log10x<-seq(from=-5,to=1.6,by=0.01)
log10y<-log10x
z<-outer(10^log10x,10^log10y,function(x,y){log(1-bf_delta+x/y)})
breaks<-pretty(z,10)
image(log10x,log10y,z,xlab="log10E1",ylab="log10C1bs1",col = hcl.colors(length(breaks)-1, "YlOrRd", rev = TRUE),
      breaks=breaks,xlim=lims,ylim=lims) 
contour(log10x,log10y,z,nlevels=10,add=T,labcex=cex_lab_cont)
mtext(expression(log[10](E[1])),1,1.2)

#now add a reference line
lines(c(min(log10x),max(log10x)),c(min(log10y),max(log10y))-log10(bf_delta)+log10(eta2/eta1),type="l",lty="dashed",lwd=1.5)

#now add a divider between tails
abline(log10(eta1*qbeta(pnorm(0),0.5,0.5))+log10(eta2*qbeta(pnorm(0),0.5,0.5))-log10(bf_delta),-1,
       lty="dashed",lwd=1.5)

#now add the points for the symmetric noise, using a different color
points(log10(E1s),log10(C1bs1s),type="p",pch=20,cex=.1)

#now add the points for the left-tail noise
points(log10(E1ta),log10(C1bs1ta),type="p",pch=20,cex=.1,col="green")

#panel label
rect(.84*dlims+lims[1],.82*dlims+lims[1],lims[2],lims[2],col="white")
text(lims[2],lims[2],"K",cex=cex_lab,adj=c(1.05,1.05))

#***now deviations from the central reference line

#x axis quantities for points on the new panel - signed distance from the central reference line
xs<-((log10(E1s)-log10(C1bs1s)-log10(bf_delta)+log10(eta2)-log10(eta1))/sqrt(2))
xta<-((log10(E1ta)-log10(C1bs1ta)-log10(bf_delta)+log10(eta2)-log10(eta1))/sqrt(2))

#y axis quantities for the same
ys<-log(1-bf_delta+E1s/C1bs1s)
yta<-log(1-bf_delta+E1ta/C1bs1ta)

#now separate these into those in the left v right tails
#y=-x+log10(eta1*qbeta(pnorm(0),0.5,0.5))+log10(eta2*qbeta(pnorm(0),0.5,0.5))-log10(bf_delta)
inds_l<-(-log10(E1s)+log10(eta1*0.5)+log10(eta2*0.5)-log10(bf_delta)>log10(C1bs1s))
inds_r<-(-log10(E1s)+log10(eta1*0.5)+log10(eta2*0.5)-log10(bf_delta)<=log10(C1bs1s))
xs_l<-xs[inds_l]
xs_r<-xs[inds_r]
ys_l<-ys[inds_l]
ys_r<-ys[inds_r]
inds_l<-(-log10(E1ta)+log10(eta1*0.5)+log10(eta2*0.5)-log10(bf_delta)>log10(C1bs1ta))
inds_r<-(-log10(E1ta)+log10(eta1*0.5)+log10(eta2*0.5)-log10(bf_delta)<=log10(C1bs1ta))
xta_l<-xta[inds_l]
xta_r<-xta[inds_r]
yta_l<-yta[inds_l]
yta_r<-yta[inds_r]

#coordinates for the plotted line which represents the value of 
#log(1-bf_delta+x/y) as a function of distance from the central 
#reference line
xol<-log10x
yol<-rev(log10y)
yl<-log(1-bf_delta+(10^xol)/(10^yol))
xl<-((xol-yol-log10(bf_delta)+log10(eta2)-log10(eta1)))/sqrt(2)  
yl<-yl[xl>=min(xs,xta) & xl<=max(xs,xta)]
xl<-xl[xl>=min(xs,xta) & xl<=max(xs,xta)]

#set up the plotting panel
i<-3
j<-4
par(fig=c((yaxwd+(j-1)*(panwd+gap))/totwd,
          (yaxwd+(j-1)*(panwd+gap)+panwd)/totwd,
          (xaxht+(3-i)*(panht+xaxht))/totht,
          (xaxht+(3-i)*(panht+xaxht)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)

#plot the line
xlim<-range(xl)
ylim<-range(yl)
xhistrg_l<-c(xlim[1]-.4*diff(xlim),xlim[1]-.25*diff(xlim))
xhistrg_r<-c(xlim[1]-.15*diff(xlim),xlim[1])
yhistrg_l<-c(ylim[1]-.4*diff(ylim),ylim[1]-.25*diff(ylim))
yhistrg_r<-c(ylim[1]-.15*diff(ylim),ylim[1])
xlim[1]<-xhistrg_l[1]
ylim[1]<-yhistrg_l[1]
plot(xl,yl,type="l",xlim=xlim,ylim=ylim)
mtext("Perp. distance",1,1.2)

#plot the mini histograms on the x axis
bks<-seq(from=min(xs,xta),to=max(xs,xta),length.out=numhistbins+1)
cts_l<-hist(xs_l,breaks=bks,plot=F)$counts
cts_r<-hist(xs_r,breaks=bks,plot=F)$counts
ctta_l<-hist(xta_l,breaks=bks,plot=F)$counts
ctta_r<-hist(xta_r,breaks=bks,plot=F)$counts

mybarplot(m=matrix(c(cts_l,ctta_l),2,length(cts_l),byrow=TRUE),bk=bks,cols=c("black","green"),
          base=yhistrg_l[1],top=yhistrg_l[2],axis=1)
mybarplot(m=matrix(c(cts_r,ctta_r),2,length(cts_r),byrow=TRUE),bk=bks,cols=c("black","green"),
          base=yhistrg_r[1],top=yhistrg_r[2],axis=1)

#plot the mini histogram on the y axis
bks<-seq(from=min(ys,yta),to=max(ys,yta),length.out=numhistbins+1)
cts_l<-hist(ys_l,breaks=bks,plot=F)$counts
cts_r<-hist(ys_r,breaks=bks,plot=F)$counts
ctta_l<-hist(yta_l,breaks=bks,plot=F)$counts
ctta_r<-hist(yta_r,breaks=bks,plot=F)$counts

mybarplot(m=matrix(c(cts_l,ctta_l),2,length(cts_l),byrow=TRUE),bk=bks,cols=c("black","green"),
          base=xhistrg_l[1],top=xhistrg_l[2],axis=2)
mybarplot(m=matrix(c(cts_r,ctta_r),2,length(cts_r),byrow=TRUE),bk=bks,cols=c("black","green"),
          base=xhistrg_r[1],top=xhistrg_r[2],axis=2)

#panel label
xlim<-c(xlim[1]-.04*diff(xlim),xlim[2]+.04*diff(xlim))
ylim<-c(ylim[1]-.04*diff(ylim),ylim[2]+.04*diff(ylim))
rect(.84*diff(xlim)+xlim[1],.82*diff(ylim)+ylim[1],
     xlim[2],ylim[2],col="white")
text(xlim[2],ylim[2],"L",cex=cex_lab,adj=c(1.05,1.05))

dev.off()