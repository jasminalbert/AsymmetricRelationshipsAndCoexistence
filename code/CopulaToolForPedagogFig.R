#The tools in this script can be used to get the pdf of a bivariate random variable 
#with various choices for copula and marginals made separately. This is intended
#as a guide/help to Jasmin in making a new pedagogical figure for our tail association
#and competition theory paper. Perhaps more tools will be needed but this is a first 
#step. Reuman, 2022 09 09.

library(copula)

#***Some copulas, to use as building blocks of an overall dsitribution, see section below.
#Can also be used on their own. 

#Gaussian copula - by default dimension 2, the parameters controls the strength
#of relationship between the two components.
ncop<-normalCopula(.7)

#Clayton copula with same Spearman as the above normal copula
ccop<-claytonCopula(2)
ccop<-claytonCopula(iRho(ccop,rho(ncop))) 
#iRho converts rho to appropriate parameter of certain copula

#There are many copulas encoded in the copula and VineCopula packages in R. Jasmin,
#if you want one with different properties, let me know.

#***Some marginals, to use as building blocks of an overall dsitribution, see section below.
#Can also be used on their own.

#normal marginals, convenience functions, codifies parameters now so no one has to
#think about them below - pdf, then cdf, then inverse cdf  
dnmarg<-function(x){return(dnorm(x,mean=0,sd=1))}
pnmarg<-function(x){return(pnorm(x,mean=0,sd=1))}
qnmarg<-function(x){return(qnorm(x,mean=0,sd=1))}

#similar but for a gamma distribution
shape<-5
scale<-1
dgmarg<-function(x){return(dgamma(x,shape=shape,scale=scale))}
pgmarg<-function(x){return(pgamma(x,shape=shape,scale=scale))}
qgmarg<-function(x){return(qgamma(x,shape=shape,scale=scale))}

#similar but for a beta distribution
shape1<-0.25
shape2<-0.25
dbmarg<-function(x){return(dbeta(x+.5,shape1=shape1,shape2=shape2))}
pbmarg<-function(x){return(pbeta(x+.5,shape1=shape1,shape2=shape2))}
qbmarg<-function(x){return(qbeta(x,shape1=shape1,shape2=shape2))}

#make sure this beta looks like the one we want
#hist(rbeta(10000,shape1=shape1,shape2=shape2))

#One can play around with shape1 and shape2, keeping them equal, to 
#get different degrees of bimodality. Smaller values mean less 
#probability mass in the middle.

#Similar but for uniform distribution - not strictly necessary but may help 
dumarg<-function(x){return(dunif(x))}
pumarg<-function(x){return(punif(x))}
qumarg<-function(x){return(qunif(x))}

#***Here is how you combine a copula with two marginals and get the pdf of the 
#result

#***JASMIN: You will make your figure partly using the tool makepdf below, and 
#taking inspiration from the examples below the function.

#This function takes a copula and a couple marginals and returns a function which
#is the pdf of the bivariate random variable determined by that copula and those 
#marginals. One can then call that function at different values and use the results
#to make a contour plot of the pdf, or whatever.
#
#Args
#cop          The copula. A copula object from the copula package. Examples given
#               in a section above.
#dmarg1       The pdf of first marginal. Examples given above.
#dmarg2       The pdf of the second marginal.
#pmarg1       The cdf of the first marginal. Examples above.
#pmarg2       The cdf of the second marginal. 
#
#Output - A function which takes a single argument which is an n by 2 matrix, each 
#row of which has a pair of x and y coordinates at which to evaluate the pdf.
#
#Notes
#If dmargi and pmargi are not from the same marginal distribution, then who knows
#what will happen but it likely won't be good. Applies for i=1,2. So make sure they
#match!
#
makepdf<-function(cop,dmarg1,dmarg2,pmarg1,pmarg2)
{
  #the argument m is an n by 2 matrix
  res<-function(m)
  {
    return(dCopula(cbind(pmarg1(m[,1]),pmarg2(m[,2])),cop,log=FALSE)*dmarg1(m[,1])*dmarg2(m[,2]))
  }
  
  return(res)
}



#***Now make function that allows samples to be drawn from a distribution specified
#in the same way.

#***JASMIN: You will make your figure partly using the tool makerandgenrtr below, and 
#taking inspiration from the examples below the function.

#This function takes a copula and a couple marginals and returns a function which
#can generate random bivariate numbers from the bivariate random variable determined
#by that copula and those marginals. 
#
#Args
#cop          The copula. A copula object from the copula package. Examples given
#               in a section above.
#qmarg1       The inverse cdf of the first marginal. Examples given above.
#qmarg2       The inverse cdf of the second marginal.
#
#Output - A function which takes a single argument which is the number of samples
#desired, and returns that many samples from the random variable.
#
makerandgenrtr<-function(cop,qmarg1,qmarg2)
{
  res<-function(numsamps)
  {
    samps_cop<-rCopula(numsamps,cop)
    samps<-cbind(qmarg1(samps_cop[,1]),qmarg2(samps_cop[,2]))
    return(samps)
  }
  
  return(res)
}

