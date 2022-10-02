source("./CopulaToolForPedagogFig.R")

#########examples using only makepdf()##############


#***make an example with uniform marginals and normal copula, and then make a contour 
#plot of the pdf - looks like Fig 3a of Ghosh et al, 2020, Adv Ecol Res
pdf_normcop_unifmargs<-makepdf(ncop,dumarg,dumarg,pumarg,pumarg)
x<-seq(from=0.01,to=0.99,by=0.01)
y<-seq(from=0.01,to=0.99,by=0.01)
xy<-expand.grid(x,y)
z<-pdf_normcop_unifmargs(xy)
z<-matrix(z,length(x),length(y))
levs<-c(-2,-1,0)
contour(x,y,log10(z),xlim=range(x),ylim=range(y),levels=levs,labels="")

#then the marginals could be plotted, as in the side panels of Fig. 3a:
plot(x,dumarg(x),type="l") #one marginal
plot(x,dumarg(x),type="l") #the other (which is the same in this case)

#***make an example with normal marginals and a normal copula - looks like Fig 3c 
#from Ghosh et al, 2020, Adv Ecol Res
pdf_normcop_normmargs<-makepdf(ncop,dnmarg,dnmarg,pnmarg,pnmarg)
x<-seq(from=-3,to=3,by=0.05)
y<-x
xy<-expand.grid(x,y)
z<-pdf_normcop_normmargs(xy)
z<-matrix(z,length(x),length(y))
levs<-c(-8,-6,-4,-2,-1)
contour(x,y,log10(z),xlim=range(x),ylim=range(y),levels=levs,labels="")

#then the marginals could be plotted, as in the side panels 
plot(x,dnmarg(x),type="l") #one marginal
plot(x,dnmarg(x),type="l") #the other (which is the same in this case)

#***make an example with gamma marginals and a normal copula - looks like Fig. 3e
pdf_normcop_gammargs<-makepdf(ncop,dgmarg,dgmarg,pgmarg,pgmarg)
x<-seq(from=0.1,to=16,by=0.1)
y<-x
xy<-expand.grid(x,y)
z<-pdf_normcop_gammargs(xy)
z<-matrix(z,length(x),length(y))
levs<-c(-10,-5,-4,-3,-2)
contour(x,y,log10(z),xlim=range(x),ylim=range(y),levels=levs,labels="")

#then the marginals could be plotted, as in the side panels 
plot(x,dgmarg(x),type="l") #one marginal
plot(x,dgmarg(x),type="l") #the other (which is the same in this case)

#***make an example with gamma marginals and clayton copula - looks like Fig. 3f
pdf_claycop_gammargs<-makepdf(ccop,dgmarg,dgmarg,pgmarg,pgmarg)
x<-seq(from=0.1,to=16,by=0.1)
y<-x
xy<-expand.grid(x,y)
z<-pdf_claycop_gammargs(xy)
z<-matrix(z,length(x),length(y))
levs<-c(-10,-5,-4,-3,-2)
contour(x,y,log10(z),xlim=range(x),ylim=range(y),levels=levs,labels="")

#then the marginals could be plotted, as in the side panels, and same as above
#since the marginals have not changed this time
plot(x,dgmarg(x),type="l") #one marginal
plot(x,dgmarg(x),type="l") #the other (which is the same in this case)


############examples using makepdf() and makerandgenrtr()#################### 

#***Do an example with uniform marginals and normal copula. Plot it on the pdf
#for comparison
pdf_normcop_unifmargs<-makepdf(ncop,dumarg,dumarg,pumarg,pumarg)
x<-seq(from=0.01,to=0.99,by=0.01)
y<-seq(from=0.01,to=0.99,by=0.01)
xy<-expand.grid(x,y)
z<-pdf_normcop_unifmargs(xy)
z<-matrix(z,length(x),length(y))
levs<-c(-2,-1,0)
contour(x,y,log10(z),xlim=range(x),ylim=range(y),levels=levs,labels="")

rand_normcop_unifmargs<-makerandgenrtr(ncop,qumarg,qumarg)
samps<-rand_normcop_unifmargs(10000)
points(samps[,1],samps[,2],type="p",pch=20,cex=.1)

#***Do an example with gamma marginals and clayton copula. Plot it on the pdf
#for comparison.
pdf_claycop_gammargs<-makepdf(ccop,dgmarg,dgmarg,pgmarg,pgmarg)
x<-seq(from=0.1,to=16,by=0.1)
y<-x
xy<-expand.grid(x,y)
z<-pdf_claycop_gammargs(xy)
z<-matrix(z,length(x),length(y))
levs<-c(-10,-5,-4,-3,-2)
contour(x,y,log10(z),xlim=range(x),ylim=range(y),levels=levs,labels="")

rand_claycop_gammargs<-makerandgenrtr(ccop,qgmarg,qgmarg)
samps<-rand_claycop_gammargs(10000)
points(samps[,1],samps[,2],type="p",pch=20,cex=.1)
